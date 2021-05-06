#' Perform automated comparisons of continuous variables between two levels of a categorical variable
#'
#' @description get_2c_comp() performs automated comparisons between the two levels of a categorical variable on a set of continuous variables.
#'
#' @details For independent samples: normality is tested independently on the two levels of the categorical variable with Shapiro–Wilk test (sample size must be >= 3). Homoscedasticity is tested with Bartlett's test in case of normality, or modified Levene's test in case of non-normality. Comparisons between the two levels of the categorical variable are performed with independent Student's t-test for normal/homoscedastic data, independent Welch's t-test for normal/heteroscedastic data, or Mann-Whitney U-test for non-normal/homoscedastic data. Comparison is not performed in case of non-normal/heteroscedastic data.
#'
#' For paired samples: normality is test on the difference between the two levels of the categorical variable with Shapiro–Wilk test (sample size must be >= 3). Comparisons between the two levels of the categorical variable are performed with paired Student's t-test for normal data, or Wilcoxon signed-rank test for non-normal data.
#'
#' Statistical significance is adjusted by false discovery rate (FDR), using the Benjamini-Hochberg procedure.
#'
#' @param tibble a tibble.
#' @param grp a string indicating the column which contains the names of the continuous variables.
#' @param dep_var a string indicating the column which contains the values of the continuous variables.
#' @param comp_nm a string indicating the column which contains the values of the categorical variable.
#' @param comp_lvl1 a string providing the first level of the categorical variable.
#' @param comp_lvl2 a string providing the second level of the categorical variable.
#' @param paired a logical indicating whether the samples are independent or paired, default is FALSE.
#' @param pairing_key if paired = TRUE, a string indicating the column to use as pairing key. Default is NULL.
#' @param FDR a numeric indicating the q-value threshold to use for FDR.
#' @param grp_as_label a logical indicating if the names of the continuous variables should be used to label the y axis on the plots. If FALSE, the value of the dep_var argument is systematically used as label. Default is FALSE.
#' @param base_size a numeric provided to the base_size argument of theme_pubr() for plotting.
#' @param multi_diff a numeric provided to the nudge_y argument of geom_text(). The higher the multi_diff value, the greater the distance between the p value label and the data points on the graph. Default is 1 (no adjustment).
#'
#' @return A tibble containing, for each continuous variable, a statistical comparison between the two levels of the categorical variable and the corresponding plot.
#'
#' @export
#'
#' @importFrom car leveneTest
#' @importFrom dplyr arrange
#' @importFrom dplyr bind_rows
#' @importFrom dplyr case_when
#' @importFrom dplyr desc
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr group_modify
#' @importFrom dplyr if_else
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr ungroup
#' @importFrom forcats as_factor
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 position_jitterdodge
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggpubr theme_pubr
#' @importFrom magrittr %>%
#' @importFrom mcp.project fdr
#' @importFrom purrr map
#' @importFrom purrr map2_chr
#' @importFrom purrr map2_dbl
#' @importFrom purrr map2_lgl
#' @importFrom purrr map_chr
#' @importFrom purrr map_dbl
#' @importFrom purrr map_lgl
#' @importFrom purrr pmap
#' @importFrom purrr pmap_chr
#' @importFrom purrr pmap_dbl
#' @importFrom purrr pmap_lgl
#' @importFrom rlang !!
#' @importFrom rlang :=
#' @importFrom rlang enexpr
#' @importFrom rlang enexprs
#' @importFrom rlang sym
#' @importFrom rlang syms
#' @importFrom stats bartlett.test
#' @importFrom stats shapiro.test
#' @importFrom stats t.test
#' @importFrom stats wilcox.test
#' @importFrom stringr str_c
#' @importFrom tidyr nest
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr unnest
#' @importFrom tidyselect all_of
#' @importFrom tidyselect everything
#' @importFrom tidyselect matches
#' @importFrom tidyselect starts_with
#' @importFrom generics tidy
#'
#' @examples
#'
get_2c_comp <- function(tibble, grp = NULL, dep_var, comp_nm, comp_lvl1, comp_lvl2, paired = FALSE, pairing_key = NULL, FDR = 0.2, grp_as_label = FALSE, base_size = 12, multi_diff = 1) {

  if (isTRUE(paired) & is.null(pairing_key)) {
    stop("A pairing key is required for comparison of paired data")
  }

  if (is.null(grp) & isTRUE(grp_as_label)) {
    stop("Grp cannot be used as y_label if NULL")
  }

  dep_var <- enexpr(dep_var)
  comp <- enexprs(comp_nm, comp_lvl1, comp_lvl2)
  pairing_key <- enexpr(pairing_key)
  y_label <- syms(grp)

  if (paired == FALSE) {

    output <- tibble %>%
      group_by(.dots = grp) %>%
      group_modify(~ nest(.data = .x, data = everything())) %>%
      ungroup() %>%
      mutate(!!str_c("mean_", comp[[2]]) := map_dbl(data,
                                                    ~ mean(filter(.x, !!sym(comp[[1]]) == !!comp[[2]]) %>% select(!!sym(dep_var)) %>% {.[[1]]}, na.rm = TRUE)),
             !!str_c("mean_", comp[[3]]) := map_dbl(data,
                                                    ~ mean(filter(.x, !!sym(comp[[1]]) == !!comp[[3]]) %>% select(!!sym(dep_var)) %>% {.[[1]]}, na.rm = TRUE)),
             !!str_c("difference_mean_", comp[[2]], "-", comp[[3]]) := map2_dbl(!!sym(str_c("mean_", comp[[2]])), !!sym(str_c("mean_", comp[[3]])), ~ .x - .y),
             !!str_c("shapiro_", comp[[2]]) := map(data,
                                                   ~ if (length(.x[[comp[[1]]]][.x[[comp[[1]]]] == comp[[2]]]) >= 3) {
                                                     shapiro.test(filter(.x, !!sym(comp[[1]]) == !!comp[[2]]) %>% select(!!sym(dep_var)) %>% {.[[1]]})
                                                   } else {
                                                     NULL
                                                   }),
             !!str_c("statistic_shapiro_", comp[[2]]) := map(!!sym(str_c("shapiro_", comp[[2]])), ~ .x$statistic[[1]]),
             !!str_c("p.value_shapiro_", comp[[2]]) := map(!!sym(str_c("shapiro_", comp[[2]])), ~ .x$p.value),
             !!str_c("shapiro_", comp[[3]]) := map(data,
                                                   ~ if (length(.x[[comp[[1]]]][.x[[comp[[1]]]] == comp[[3]]]) >= 3) {
                                                     shapiro.test(filter(.x, !!sym(comp[[1]]) == !!comp[[3]]) %>% select(!!sym(dep_var)) %>% {.[[1]]})
                                                   } else {
                                                     NULL
                                                   }),
             !!str_c("statistic_shapiro_", comp[[3]]) := map(!!sym(str_c("shapiro_", comp[[3]])), ~ .x$statistic[[1]]),
             !!str_c("p.value_shapiro_", comp[[3]]) := map(!!sym(str_c("shapiro_", comp[[3]])), ~ .x$p.value),
             normal = map2_lgl(!!sym(str_c("p.value_shapiro_", comp[[2]])), !!sym(str_c("p.value_shapiro_", comp[[3]])),
                               ~ case_when(.x > 0.05 & .y > 0.05 ~ TRUE,
                                           is.na(.x) | is.na(.y) ~ NA,
                                           TRUE ~ FALSE)),
             bartlett = map2(normal, data,
                             ~ if (isTRUE(.x)) {
                               bartlett.test(formula = !!sym(dep_var) ~ !!sym(comp[[1]]), data = .y)
                             } else {
                               NULL
                             }),
             statistic_bartlett = map_dbl(bartlett, ~ if_else(is.null(.x), NA_real_, .x$statistic[[1]])),
             p.value_bartlett = map_dbl(bartlett, ~ if_else(is.null(.x), NA_real_, .x$p.value)),
             levene = map2(normal, data,
                           ~ if (isFALSE(.x)) {
                             leveneTest(y = !!sym(dep_var) ~ !!sym(comp[[1]]),
                                        data = .y %>% mutate(!!sym(comp[[1]]) := as_factor(!!sym(comp[[1]]))),
                                        center = median)
                           } else {
                             NULL
                           }),
             statistic_levene = map_dbl(levene, ~ if_else(is.null(.x), NA_real_, .x[[2]][[1]])),
             p.value_levene = map_dbl(levene, ~ if_else(is.null(.x), NA_real_, .x[[3]][[1]])),
             var.equal = pmap_lgl(list(normal,
                                       p.value_bartlett,
                                       p.value_levene),
                                  ~ case_when((isTRUE(..1) & ..2 > 0.05) | (isFALSE(..1) & ..3 > 0.05) ~ TRUE,
                                              is.na(..1) ~ NA,
                                              TRUE ~ FALSE)),
             indep.student = pmap(list(normal,
                                       var.equal,
                                       data),
                                  ~ if (isTRUE(..1) & isTRUE(..2)) {
                                    t.test(formula = !!sym(dep_var) ~ !!sym(comp[[1]]),
                                           data = ..3,
                                           alternative = "two.sided",
                                           paired = FALSE,
                                           var.equal = TRUE)
                                  } else {
                                    NULL
                                  }),
             statistic_indep.student = map_dbl(indep.student, ~ if_else(is.null(.x), NA_real_, .x$statistic[[1]])),
             p.value_indep.student = map_dbl(indep.student, ~ if_else(is.null(.x), NA_real_, .x$p.value)),
             indep.welch = pmap(list(normal,
                                     var.equal,
                                     data),
                                ~ if (isTRUE(..1) & isFALSE(..2)) {
                                  t.test(formula = !!sym(dep_var) ~ !!sym(comp[[1]]),
                                         data = ..3,
                                         alternative = "two.sided",
                                         paired = FALSE,
                                         var.equal = FALSE)
                                } else {
                                  NULL
                                }),
             statistic_indep.welch = map_dbl(indep.welch, ~ if_else(is.null(.x), NA_real_, .x$statistic[[1]])),
             p.value_indep.welch = map_dbl(indep.welch, ~ if_else(is.null(.x), NA_real_, .x$p.value)),
             mann.whitney = pmap(list(normal,
                                      var.equal,
                                      data),
                                 ~ if (isFALSE(..1) & isTRUE(..2)) {
                                   wilcox.test(formula = !!sym(dep_var) ~ !!sym(comp[[1]]),
                                               data = ..3,
                                               alternative = "two.sided",
                                               paired = FALSE)
                                 } else {
                                   NULL
                                 }),
             statistic_mann.whitney = map_dbl(mann.whitney, ~ if_else(is.null(.x), NA_real_, .x$statistic[[1]])),
             p.value_mann.whitney = map_dbl(mann.whitney, ~ if_else(is.null(.x), NA_real_, .x$p.value)),
             type_pairwise.comp = map2_chr(normal, var.equal,
                                           ~ case_when(isTRUE(.x) & isTRUE(.y) ~ "indep.student",
                                                       isTRUE(.x) & isFALSE(.y) ~ "indep.welch",
                                                       isFALSE(.x) & isTRUE(.y) ~ "mann.whitney",
                                                       TRUE ~ NA_character_)),
             p.value_pairwise.comp = pmap_dbl(list(normal,
                                                   var.equal,
                                                   p.value_indep.student,
                                                   p.value_indep.welch,
                                                   p.value_mann.whitney),
                                              ~ case_when(isTRUE(..1) & isTRUE(..2) ~ ..3,
                                                          isTRUE(..1) & isFALSE(..2) ~ ..4,
                                                          isFALSE(..1) & isTRUE(..2) ~ ..5,
                                                          TRUE ~ NA_real_)))

    out_comp <- output %>%
      filter(!is.na(p.value_pairwise.comp)) %>%
      nest(data = everything()) %>%
      mutate(fdr = map(data, ~ fdr(.x$p.value_pairwise.comp, q = FDR, method = "BH")$Pvals)) %>%
      unnest(cols = c(data, fdr)) %>%
      rename(pass_fdr = rejected) %>%
      mutate(fdr_result = map_chr(pass_fdr,
                                  ~ if_else(isTRUE(.x),
                                            str_c("sign (FDR = ", FDR, ")"),
                                            str_c("ns (FDR = ", FDR, ")"))))

    out_no_comp <- output %>%
      filter(is.na(p.value_pairwise.comp))

    output <- bind_rows(out_comp, out_no_comp)

    if (isTRUE(grp_as_label)) {

      output <- output %>%
        mutate(y_label = pmap_chr(list(!!!y_label), ~ str_c(..., sep = " - ")),
               plot = pmap(list(data, y_label, p.value_pairwise.comp, fdr_result),
                           ~ ggplot(data = ..1, mapping = aes(x = .data[[comp[[1]]]], y = .data[[dep_var]])) +
                             geom_boxplot(aes(fill = .data[[comp[[1]]]]), alpha = 0.4, outlier.shape = NA) +
                             geom_point(aes(color = .data[[comp[[1]]]]), position = position_jitterdodge(jitter.width = 0.2), alpha = 0.7) +
                             scale_colour_manual(values = c("#5d76cb", "#ca3767", "#3bb08f", "#ffcf48"), aesthetics = c("fill", "colour")) +
                             labs(x = NULL, y = ..2) +
                             theme_pubr(legend = "none", base_size = base_size) +
                             if (is.na(..3)) {
                               geom_text(data = ..1 %>%
                                           filter(!!sym(dep_var) %in% c(max(.[.[[comp[[1]]]] == comp[[2]],][[!!dep_var]]), max(.[.[[comp[[1]]]] == comp[[3]],][[!!dep_var]]))) %>%
                                           arrange(!!sym(dep_var)) %>%
                                           {.[1,]} %>%
                                           mutate(lab = "p not computed") %>%
                                           select(!!sym(comp[[1]]), !!sym(dep_var), lab),
                                         aes(x = !!sym(comp[[1]]), y = !!sym(dep_var), label = lab),
                                         nudge_y = ..1 %>%
                                           select(!!sym(comp[[1]]), !!sym(dep_var)) %>%
                                           group_by(!!sym(comp[[1]])) %>%
                                           nest() %>%
                                           mutate(data = map(data, ~ .x %>% arrange(desc(!!sym(dep_var))) %>% {.[1,]})) %>%
                                           unnest(cols = data) %>%
                                           pivot_wider(names_from = comp[[1]], values_from = dep_var) %>%
                                           mutate(abs_diff = abs(!!sym(comp[[2]]) - !!sym(comp[[3]]))) %>%
                                           select(abs_diff) %>% {.[[1]] * multi_diff})
                             } else {
                               geom_text(data = ..1 %>%
                                           filter(!!sym(dep_var) %in% c(max(.[.[[comp[[1]]]] == comp[[2]],][[!!dep_var]]), max(.[.[[comp[[1]]]] == comp[[3]],][[!!dep_var]]))) %>%
                                           arrange(!!sym(dep_var)) %>%
                                           {.[1,]} %>%
                                           mutate(lab = str_c("p = ", round(..3, 4), "\n", ..4)) %>%
                                           select(!!sym(comp[[1]]), !!sym(dep_var), lab),
                                         aes(x = !!sym(comp[[1]]), y = !!sym(dep_var), label = lab),
                                         nudge_y = ..1 %>%
                                           select(!!sym(comp[[1]]), !!sym(dep_var)) %>%
                                           group_by(!!sym(comp[[1]])) %>%
                                           nest() %>%
                                           mutate(data = map(data, ~ .x %>% arrange(desc(!!sym(dep_var))) %>% {.[1,]})) %>%
                                           unnest(cols = data) %>%
                                           pivot_wider(names_from = comp[[1]], values_from = dep_var) %>%
                                           mutate(abs_diff = abs(!!sym(comp[[2]]) - !!sym(comp[[3]]))) %>%
                                           select(abs_diff) %>% {.[[1]] * multi_diff})
                             }))

    } else {

      output <- output %>%
        mutate(plot = pmap(list(data, p.value_pairwise.comp, fdr_result),
                           ~ ggplot(data = ..1, mapping = aes(x = .data[[comp[[1]]]], y = .data[[dep_var]])) +
                             geom_boxplot(aes(fill = .data[[comp[[1]]]]), alpha = 0.4, outlier.shape = NA) +
                             geom_point(aes(color = .data[[comp[[1]]]]), position = position_jitterdodge(jitter.width = 0.2), alpha = 0.7) +
                             scale_colour_manual(values = c("#5d76cb", "#ca3767", "#3bb08f", "#ffcf48"), aesthetics = c("fill", "colour")) +
                             labs(x = NULL) +
                             theme_pubr(legend = "none", base_size = base_size) +
                             if (is.na(..2)) {
                               geom_text(data = ..1 %>%
                                           filter(!!sym(dep_var) %in% c(max(.[.[[comp[[1]]]] == comp[[2]],][[!!dep_var]]), max(.[.[[comp[[1]]]] == comp[[3]],][[!!dep_var]]))) %>%
                                           arrange(!!sym(dep_var)) %>%
                                           {.[1,]} %>%
                                           mutate(lab = "p not computed") %>%
                                           select(!!sym(comp[[1]]), !!sym(dep_var), lab),
                                         aes(x = !!sym(comp[[1]]), y = !!sym(dep_var), label = lab),
                                         nudge_y = ..1 %>%
                                           select(!!sym(comp[[1]]), !!sym(dep_var)) %>%
                                           group_by(!!sym(comp[[1]])) %>%
                                           nest() %>%
                                           mutate(data = map(data, ~ .x %>% arrange(desc(!!sym(dep_var))) %>% {.[1,]})) %>%
                                           unnest(cols = data) %>%
                                           pivot_wider(names_from = comp[[1]], values_from = dep_var) %>%
                                           mutate(abs_diff = abs(!!sym(comp[[2]]) - !!sym(comp[[3]]))) %>%
                                           select(abs_diff) %>% {.[[1]] * multi_diff})
                             } else {
                               geom_text(data = ..1 %>%
                                           filter(!!sym(dep_var) %in% c(max(.[.[[comp[[1]]]] == comp[[2]],][[!!dep_var]]), max(.[.[[comp[[1]]]] == comp[[3]],][[!!dep_var]]))) %>%
                                           arrange(!!sym(dep_var)) %>%
                                           {.[1,]} %>%
                                           mutate(lab = str_c("p = ", round(..2, 4), "\n", ..3)) %>%
                                           select(!!sym(comp[[1]]), !!sym(dep_var), lab),
                                         aes(x = !!sym(comp[[1]]), y = !!sym(dep_var), label = lab),
                                         nudge_y = ..1 %>%
                                           select(!!sym(comp[[1]]), !!sym(dep_var)) %>%
                                           group_by(!!sym(comp[[1]])) %>%
                                           nest() %>%
                                           mutate(data = map(data, ~ .x %>% arrange(desc(!!sym(dep_var))) %>% {.[1,]})) %>%
                                           unnest(cols = data) %>%
                                           pivot_wider(names_from = comp[[1]], values_from = dep_var) %>%
                                           mutate(abs_diff = abs(!!sym(comp[[2]]) - !!sym(comp[[3]]))) %>%
                                           select(abs_diff) %>% {.[[1]] * multi_diff})

                             }))
    }

    output <- output %>%
      select(all_of(grp), starts_with("mean"), matches("^difference.*$"), type_pairwise.comp, p.value_pairwise.comp, fdr_result, everything())

  } else {

    output <- tibble %>%
      mutate(!!pairing_key := as_factor(.[[pairing_key]])) %>%
      arrange(.[[pairing_key]], .[[comp[[1]]]]) %>%
      group_by(.dots = grp) %>%
      group_modify(~ nest(.data = .x, data = everything())) %>%
      ungroup() %>%
      mutate(!!str_c("mean_", comp[[2]]) := map_dbl(data,
                                                    ~ mean(filter(.x, !!sym(comp[[1]]) == !!comp[[2]]) %>% select(!!sym(dep_var)) %>% {.[[1]]}, na.rm = TRUE)),
             !!str_c("mean_", comp[[3]]) := map_dbl(data,
                                                    ~ mean(filter(.x, !!sym(comp[[1]]) == !!comp[[3]]) %>% select(!!sym(dep_var)) %>% {.[[1]]}, na.rm = TRUE)),
             data_difference = map(data,
                                   ~ select(.x, !!pairing_key, !!dep_var, !!comp[[1]]) %>%
                                     pivot_wider(names_from = !!comp[[1]], values_from = !!dep_var) %>%
                                     mutate(diff = .[[comp[[2]]]] - .[[comp[[3]]]])),
             !!str_c("difference_mean_", comp[[2]], "-", comp[[3]]) := map_dbl(data_difference, ~ mean(.x$diff, na.rm = TRUE)),
             shapiro_difference = map(data_difference, ~ shapiro.test(.x$diff)),
             statistic_shapiro_difference = map_dbl(shapiro_difference, ~ .x$statistic),
             p.value_shapiro_difference = map_dbl(shapiro_difference, ~ .x$p.value),
             normal_difference = map_lgl(p.value_shapiro_difference, ~ if_else(.x > 0.05, TRUE, FALSE)),
             paired.student = map2(normal_difference, data,
                                   ~ if (isTRUE(.x)) {
                                     t.test(formula = !!sym(dep_var) ~ !!sym(comp[[1]]),
                                            data = .y,
                                            alternative = "two.sided",
                                            paired = TRUE)
                                   } else {
                                     NULL
                                   }),
             statistic_paired.student = map2_dbl(normal_difference, paired.student,
                                                 ~ if (isTRUE(.x)) {
                                                   .y$statistic[[1]]
                                                 } else {
                                                   NA_real_
                                                 }),
             p.value_paired.student = map2_dbl(normal_difference, paired.student,
                                               ~ if (isTRUE(.x)) {
                                                 .y$p.value
                                               } else {
                                                 NA_real_
                                               }),
             wilcoxon.signed.rank = map2(normal_difference, data,
                                         ~ if (isFALSE(.x)) {
                                           wilcox.test(formula = !!sym(dep_var) ~ !!sym(comp[[1]]),
                                                       data = .y,
                                                       alternative = "two.sided",
                                                       paired = TRUE)
                                         } else {
                                           NULL
                                         }),
             statistic_wilcoxon.signed.rank = map2_dbl(normal_difference, wilcoxon.signed.rank,
                                                       ~ if (isFALSE(.x)) {
                                                         .y$statistic[[1]]
                                                       } else {
                                                         NA_real_
                                                       }),
             p.value_wilcoxon.signed.rank = map2_dbl(normal_difference, wilcoxon.signed.rank,
                                                     ~ if (isFALSE(.x)) {
                                                       .y$p.value
                                                     } else {
                                                       NA_real_
                                                     }),
             p.value_pairwise.comp = pmap_dbl(list(normal_difference,
                                                   p.value_paired.student,
                                                   p.value_wilcoxon.signed.rank),
                                              ~ case_when(isTRUE(..1) ~ ..2,
                                                          isFALSE(..1) ~ ..3,
                                                          TRUE ~ NA_real_)),
             type_pairwise.comp = map_chr(normal_difference,
                                          ~ case_when(isTRUE(.x) ~ "paired.student",
                                                      isFALSE(.x) ~ "wilcoxon.signed.rank",
                                                      TRUE ~ NA_character_)))

    out_comp <- output %>%
      filter(!is.na(p.value_pairwise.comp)) %>%
      nest(data = everything()) %>%
      mutate(fdr = map(data, ~ fdr(.x$p.value_pairwise.comp, q = FDR, method = "BH")$Pvals)) %>%
      unnest(cols = c(data, fdr)) %>%
      rename(pass_fdr = rejected) %>%
      mutate(fdr_result = map_chr(pass_fdr,
                                  ~ if_else(isTRUE(.x),
                                            str_c("sign (FDR = ", FDR, ")"),
                                            str_c("ns (FDR = ", FDR, ")"))))

    out_no_comp <- output %>%
      filter(is.na(p.value_pairwise.comp))

    output <- bind_rows(out_comp, out_no_comp)

    if (isTRUE(grp_as_label)) {

      output <- output %>%
        mutate(y_label = pmap_chr(list(!!!y_label), ~ str_c(..., sep = " - ")),
               plot = pmap(list(data, y_label, p.value_pairwise.comp, fdr_result),
                           ~ ggplot(data = ..1, mapping = aes(x = .data[[comp[[1]]]], y = .data[[dep_var]])) +
                             geom_boxplot(aes(fill = .data[[comp[[1]]]]), alpha = 0.4, outlier.shape = NA) +
                             geom_line(aes(group = .data[[pairing_key]]), color = "grey", alpha = 0.5) +
                             geom_point(aes(color = .data[[comp[[1]]]]), alpha = 0.7) +
                             scale_colour_manual(values = c("#5d76cb", "#ca3767", "#3bb08f", "#ffcf48"), aesthetics = c("fill", "colour")) +
                             labs(x = NULL, y = ..2) +
                             theme_pubr(legend = "none", base_size = base_size) +
                             if (is.na(..3)) {
                               geom_text(data = ..1 %>%
                                           filter(!!sym(dep_var) %in% c(max(.[.[[comp[[1]]]] == comp[[2]],][[!!dep_var]]), max(.[.[[comp[[1]]]] == comp[[3]],][[!!dep_var]]))) %>%
                                           arrange(!!sym(dep_var)) %>%
                                           {.[1,]} %>%
                                           mutate(lab = "p not computed") %>%
                                           select(!!sym(comp[[1]]), !!sym(dep_var), lab),
                                         aes(x = !!sym(comp[[1]]), y = !!sym(dep_var), label = lab),
                                         nudge_y = ..1 %>%
                                           select(!!sym(comp[[1]]), !!sym(dep_var)) %>%
                                           group_by(!!sym(comp[[1]])) %>%
                                           nest() %>%
                                           mutate(data = map(data, ~ .x %>% arrange(desc(!!sym(dep_var))) %>% {.[1,]})) %>%
                                           unnest(cols = data) %>%
                                           pivot_wider(names_from = comp[[1]], values_from = dep_var) %>%
                                           mutate(abs_diff = abs(!!sym(comp[[2]]) - !!sym(comp[[3]]))) %>%
                                           select(abs_diff) %>% {.[[1]] * multi_diff})
                             } else {
                               geom_text(data = ..1 %>%
                                           filter(!!sym(dep_var) %in% c(max(.[.[[comp[[1]]]] == comp[[2]],][[!!dep_var]]), max(.[.[[comp[[1]]]] == comp[[3]],][[!!dep_var]]))) %>%
                                           arrange(!!sym(dep_var)) %>%
                                           {.[1,]} %>%
                                           mutate(lab = str_c("p = ", round(..3, 4), "\n", ..4)) %>%
                                           select(!!sym(comp[[1]]), !!sym(dep_var), lab),
                                         aes(x = !!sym(comp[[1]]), y = !!sym(dep_var), label = lab),
                                         nudge_y = ..1 %>%
                                           select(!!sym(comp[[1]]), !!sym(dep_var)) %>%
                                           group_by(!!sym(comp[[1]])) %>%
                                           nest() %>%
                                           mutate(data = map(data, ~ .x %>% arrange(desc(!!sym(dep_var))) %>% {.[1,]})) %>%
                                           unnest(cols = data) %>%
                                           pivot_wider(names_from = comp[[1]], values_from = dep_var) %>%
                                           mutate(abs_diff = abs(!!sym(comp[[2]]) - !!sym(comp[[3]]))) %>%
                                           select(abs_diff) %>% {.[[1]] * multi_diff})
                             }))

    } else {

      output <- output %>%
        mutate(plot = pmap(list(data, p.value_pairwise.comp, fdr_result),
                           ~ ggplot(data = ..1, mapping = aes(x = .data[[comp[[1]]]], y = .data[[dep_var]])) +
                             geom_boxplot(aes(fill = .data[[comp[[1]]]]), alpha = 0.4, outlier.shape = NA) +
                             geom_line(aes(group = .data[[pairing_key]]), color = "grey", alpha = 0.5) +
                             geom_point(aes(color = .data[[comp[[1]]]]), alpha = 0.7) +
                             scale_colour_manual(values = c("#5d76cb", "#ca3767", "#3bb08f", "#ffcf48"), aesthetics = c("fill", "colour")) +
                             labs(x = NULL) +
                             theme_pubr(legend = "none", base_size = base_size) +
                             if (is.na(..2)) {
                               geom_text(data = ..1 %>%
                                           filter(!!sym(dep_var) %in% c(max(.[.[[comp[[1]]]] == comp[[2]],][[!!dep_var]]), max(.[.[[comp[[1]]]] == comp[[3]],][[!!dep_var]]))) %>%
                                           arrange(!!sym(dep_var)) %>%
                                           {.[1,]} %>%
                                           mutate(lab = "p not computed") %>%
                                           select(!!sym(comp[[1]]), !!sym(dep_var), lab),
                                         aes(x = !!sym(comp[[1]]), y = !!sym(dep_var), label = lab),
                                         nudge_y = ..1 %>%
                                           select(!!sym(comp[[1]]), !!sym(dep_var)) %>%
                                           group_by(!!sym(comp[[1]])) %>%
                                           nest() %>%
                                           mutate(data = map(data, ~ .x %>% arrange(desc(!!sym(dep_var))) %>% {.[1,]})) %>%
                                           unnest(cols = data) %>%
                                           pivot_wider(names_from = comp[[1]], values_from = dep_var) %>%
                                           mutate(abs_diff = abs(!!sym(comp[[2]]) - !!sym(comp[[3]]))) %>%
                                           select(abs_diff) %>% {.[[1]] * multi_diff})
                             } else {
                               geom_text(data = ..1 %>%
                                           filter(!!sym(dep_var) %in% c(max(.[.[[comp[[1]]]] == comp[[2]],][[!!dep_var]]), max(.[.[[comp[[1]]]] == comp[[3]],][[!!dep_var]]))) %>%
                                           arrange(!!sym(dep_var)) %>%
                                           {.[1,]} %>%
                                           mutate(lab = str_c("p = ", round(..2, 4), "\n", ..3)) %>%
                                           select(!!sym(comp[[1]]), !!sym(dep_var), lab),
                                         aes(x = !!sym(comp[[1]]), y = !!sym(dep_var), label = lab),
                                         nudge_y = ..1 %>%
                                           select(!!sym(comp[[1]]), !!sym(dep_var)) %>%
                                           group_by(!!sym(comp[[1]])) %>%
                                           nest() %>%
                                           mutate(data = map(data, ~ .x %>% arrange(desc(!!sym(dep_var))) %>% {.[1,]})) %>%
                                           unnest(cols = data) %>%
                                           pivot_wider(names_from = comp[[1]], values_from = dep_var) %>%
                                           mutate(abs_diff = abs(!!sym(comp[[2]]) - !!sym(comp[[3]]))) %>%
                                           select(abs_diff) %>% {.[[1]] * multi_diff})
                             }))
    }

    output <- output %>%
      select(all_of(grp), starts_with("mean"), matches("^difference.*$"), type_pairwise.comp, p.value_pairwise.comp, fdr_result, everything())

  }
}
