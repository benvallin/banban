#' Perform automated ANOVA and posthoc test
#'
#' @description get_complex_anova() performs automated ANOVA and subsequent posthoc test on a set of dependent variables.
#'
#' @details Normality of residuals is tested with Shapiro–Wilk test. In case of normality, the model is used for all subsequent analyses.
#' In case of non-normality, the dependent variable is subjected to Box-Cox transformation and normality of residuals is tested again with Shapiro–Wilk test. In case of normality, the Box-Cox-transformed model is used for all subsequent analyses.
#' In case of non-normality after Box-Cox transformation, the dependent variable is subjected to log transformation and normality of residuals is tested once again with Shapiro–Wilk test. In case of normality, the log-transformed model is used for all subsequent analyses.
#' If the distribution of residuals remains non-normal after the two types of transformation, the dependent variable is not further analysed.
#'
#' Homoscedasticity is tested with Breusch–Pagan test.
#'
#' In case of normal/homoscedastic data, classic type III ANOVA is performed. In case of normal/heteroscedastic data, type III ANOVA with White adjustment is performed.
#'
#' Factors with a significant effect in ANOVA are further subjected to posthoc test. In case of significant interaction effect(s), only the significant interaction term(s) are subjected to posthoc analysis. For normal/homoscedastic data, posthoc test is performed by comparison of estimated marginal means and one of the following adjustment methods: "tukey", "scheffe", "sidak", "bonferroni", "dunnettx", "mvt", and "none". For normal/heteroscedastic data, posthoc test is performed with Games-Howell test.
#'
#' @param tibble a tibble.
#' @param grp a string indicating the column which contains the names of the dependent variables.
#' @param lhs a string indicating the column which contains the values of the dependent variables.
#' @param rhs a string indicating the name(s) of the independent variable(s). Must follow the classic formula rules (use "+", "*" and ":" as separators).
#' @param base.size a numeric provided to the base_size argument of theme_pubr() for plotting.
#' @param legend a string provided to the legend argument of theme_pubr() for plotting. Must be one of "top", "bottom", "left", "right", "none".
#' @param spe.contr a logical indicating if estimated marginal means are desired over a specific independent variable in case of significant interaction. Default is FALSE.
#' @param specs if spe.contr = TRUE, a string indicating the name of the independent variable over which estimated marginal means are desired. Will be used in case of significant interaction between the independent variables passed to the specs and by arguments. Default is NA_character_.
#' @param by if spe.contr = TRUE, a string indicating the name of the independent variable to condition on. Will be used in case of significant interaction between the independent variables passed to the specs and by arguments. Default is NA_character_.
#' @param adjust.emm.comp a string indicating which adjustment method should be used after computation of estimated marginal means. Must be one of "tukey", "scheffe", "sidak", "bonferroni", "dunnettx", "mvt", and "none". Default is "tukey".
#' @param block a character vector indicating the name(s) of the independent variable(s) for which the ANOVA-associated p values should not be displayed on the plot. Default is NA_character_.
#' @param jumps a numeric indicating the number of new lines that should be added at the end of the caption which provides ANOVA-associated p values on the plot.
#'
#' @return A tibble containing, for each dependent variable, an ANOVA and associated posthoc test as well as the corresponding plot.
#'
#' @export
#'
#' @importFrom MASS boxcox
#' @importFrom car Anova
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
#' @importFrom dplyr select_if
#' @importFrom dplyr ungroup
#' @importFrom emmeans contrast
#' @importFrom emmeans emmeans
#' @importFrom forcats as_factor
#' @importFrom generics tidy
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 position_jitterdodge
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggpubr theme_pubr
#' @importFrom graphics pairs
#' @importFrom lmtest bptest
#' @importFrom magrittr %>%
#' @importFrom purrr keep
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom purrr map2_chr
#' @importFrom purrr map2_lgl
#' @importFrom purrr map_chr
#' @importFrom purrr map_dbl
#' @importFrom purrr map_lgl
#' @importFrom purrr modify_at
#' @importFrom purrr pmap
#' @importFrom purrr pmap_chr
#' @importFrom purrr some
#' @importFrom rlang !!
#' @importFrom rlang enexpr
#' @importFrom rlang is_empty
#' @importFrom rlang parse_expr
#' @importFrom rlang syms
#' @importFrom rstatix games_howell_test
#' @importFrom stats as.formula
#' @importFrom stats lm
#' @importFrom stats residuals
#' @importFrom stats shapiro.test
#' @importFrom stringr str_c
#' @importFrom stringr str_detect
#' @importFrom stringr str_remove_all
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_split
#' @importFrom tibble as_tibble
#' @importFrom tidyr nest
#' @importFrom tidyselect all_of
#' @importFrom tidyselect everything
#' @importFrom tidyselect matches
#'
#' @examples
#'
get_complex_anova <- function(tibble, grp = NULL, lhs, rhs, base.size = 14, legend = "none", spe.contr = FALSE, specs = NA_character_, by = NA_character_, adjust.emm.comp = "tukey", block = NA_character_, jumps = 0) {

  options(contrasts = c("contr.sum", "contr.poly"))

  after.emmeans <- function(input) {
    output <- vector("list", length(input))
    for (i in seq_along(input)) {
      output[[i]] <- as_tibble(summary(input[[i]]))
    }
    output <- bind_rows(output)
  }

  # Process user's arguments to construct formula components
  lhs <- enexpr(lhs)
  rhs <- enexpr(rhs)
  formula <- parse_expr(str_c(lhs, " ~ ", rhs))

  indep.vars <- str_split(string = as.character(rhs), pattern = "\\+|\\*|\\:")[[1]] %>% str_remove_all(pattern = "\\s") %>% unique

  y.label <- syms(grp)

  # Transform input tibble according to user's groups and dependent / independent variables
  out <- tibble %>%
    modify_at(.at = indep.vars, .f = function(x) if(is.character(x)) as_factor(x) else x) %>%
    group_by(.dots = grp) %>%
    group_modify(~ nest(.data = .x, data = everything())) %>%
    ungroup() %>%
    mutate(dependent.variable = lhs,
           independent.variables = rhs,
           # mod
           mod = map(data, ~ lm(!!formula, data = .x)),
           shapiro.mod = map(mod, ~ shapiro.test(residuals(.x))),
           statistic.shapiro.mod = map_dbl(shapiro.mod, ~ .x$statistic[[1]]),
           p.value.shapiro.mod = map_dbl(shapiro.mod, ~ .x$p.value),
           normality.mod = map_lgl(p.value.shapiro.mod, ~ if_else(.x > 0.05, TRUE, FALSE)),
           breusch.pagan.mod = map(mod, ~ bptest(.x)),
           statistic.breusch.pagan.mod = map_dbl(breusch.pagan.mod, ~ .x$statistic[[1]]),
           p.value.breusch.pagan.mod = map_dbl(breusch.pagan.mod, ~ .x$p.value[[1]]),
           homoscedasticity.mod = map_lgl(p.value.breusch.pagan.mod, ~ if_else(.x > 0.05, TRUE, FALSE)),
           norm.homo.mod = map2_lgl(normality.mod,
                                    homoscedasticity.mod,
                                    ~ if_else(isTRUE(.x) & isTRUE(.y), TRUE, FALSE)),
           anova.mod = map(mod, ~ Anova(.x, type = 3)),
           white.anova.mod = map(mod, ~ Anova(.x, type = 3, white.adjust = TRUE)),
           # mod.boxcox
           lambda.boxcox = map_dbl(data, ~ boxcox(!!formula, data = .x) %>% as_tibble() %>% arrange(desc(y)) %>% {.[[1, 1]]}),
           mod.boxcox = map2(lambda.boxcox, data, ~ lm(as.formula(str_c("((", lhs, "^", .x, "-1) / ", .x, ") ~", rhs)), data = .y)),
           shapiro.mod.boxcox = map(mod.boxcox, ~ shapiro.test(residuals(.x))),
           statistic.shapiro.mod.boxcox = map_dbl(shapiro.mod.boxcox, ~ .x$statistic[[1]]),
           p.value.shapiro.mod.boxcox = map_dbl(shapiro.mod.boxcox, ~ .x$p.value),
           normality.mod.boxcox = map_lgl(p.value.shapiro.mod.boxcox, ~ if_else(.x > 0.05, TRUE, FALSE)),
           breusch.pagan.mod.boxcox = map(mod.boxcox, ~ bptest(.x)),
           statistic.breusch.pagan.mod.boxcox = map_dbl(breusch.pagan.mod.boxcox, ~ .x$statistic[[1]]),
           p.value.breusch.pagan.mod.boxcox = map_dbl(breusch.pagan.mod.boxcox, ~ .x$p.value[[1]]),
           homoscedasticity.mod.boxcox = map_lgl(p.value.breusch.pagan.mod.boxcox, ~ if_else(.x > 0.05, TRUE, FALSE)),
           norm.homo.mod.boxcox = map2_lgl(normality.mod.boxcox,
                                           homoscedasticity.mod.boxcox,
                                           ~ if_else(isTRUE(.x) & isTRUE(.y), TRUE, FALSE)),
           anova.mod.boxcox = map(mod.boxcox, ~ Anova(.x, type = 3)),
           white.anova.mod.boxcox = map(mod.boxcox, ~ Anova(.x, type = 3, white.adjust = TRUE)),
           # mod.log
           mod.log = map(data, ~ lm(as.formula(str_c("`log`(", lhs, ")", " ~ ", rhs)), data = .x)),
           shapiro.mod.log = map(mod.log, ~ shapiro.test(residuals(.x))),
           statistic.shapiro.mod.log = map_dbl(shapiro.mod.log, ~ .x$statistic[[1]]),
           p.value.shapiro.mod.log = map_dbl(shapiro.mod.log, ~ .x$p.value),
           normality.mod.log = map_lgl(p.value.shapiro.mod.log, ~ if_else(.x > 0.05, TRUE, FALSE)),
           breusch.pagan.mod.log = map(mod.log, ~ bptest(.x)),
           statistic.breusch.pagan.mod.log = map_dbl(breusch.pagan.mod.log, ~ .x$statistic[[1]]),
           p.value.breusch.pagan.mod.log = map_dbl(breusch.pagan.mod.log, ~ .x$p.value[[1]]),
           homoscedasticity.mod.log = map_lgl(p.value.breusch.pagan.mod.log, ~ if_else(.x > 0.05, TRUE, FALSE)),
           norm.homo.mod.log = map2_lgl(normality.mod.log,
                                        homoscedasticity.mod.log,
                                        ~ if_else(isTRUE(.x) & isTRUE(.y), TRUE, FALSE)),
           anova.mod.log = map(mod.log, ~ Anova(.x, type = 3)),
           white.anova.mod.log = map(mod.log, ~ Anova(.x, type = 3, white.adjust = TRUE)),
           # summary
           best.mod.nm = pmap_chr(list(normality.mod, #1
                                       lambda.boxcox, #2
                                       normality.mod.boxcox, #3
                                       normality.mod.log), #4
                                  ~ case_when(isTRUE(..1) ~ "mod",
                                              (..2 < -0.01 | ..2 > 0.01) & isTRUE(..3) ~ "mod.boxcox",
                                              isTRUE(..4) ~ "mod.log",
                                              TRUE ~ "none")),
           formula.analysis = pmap_chr(list(best.mod.nm, #1
                                            dependent.variable, #2
                                            lambda.boxcox), #3
                                       ~ if (..1 == "none") {
                                         NA_character_
                                       } else if (..1 == "mod") {
                                         str_c(..2, " ~ ", rhs)
                                       } else if (..1 == "mod.boxcox") {
                                         str_c("((", ..2, "^", round(..3, 4), "-1) / ", round(..3, 4), ") ~ ", rhs)
                                       } else {
                                         str_c("log(", ..2, ") ~ ", rhs)
                                       }),
           best.mod.val = pmap(list(best.mod.nm, #1
                                    mod, #2
                                    mod.boxcox, #3
                                    mod.log), #4
                               ~ case_when(..1 == "mod" ~ list(..2),
                                           ..1 == "mod.boxcox" ~ list(..3),
                                           ..1 == "mod.log" ~ list(..4),
                                           TRUE ~ list(NULL))[[1]]),
           best.anova.nm = pmap_chr(list(best.mod.nm, #1
                                         norm.homo.mod, #2
                                         norm.homo.mod.boxcox, #3
                                         norm.homo.mod.log), #4
                                    ~ case_when(..1 == "mod" & isTRUE(..2) ~ "anova.mod",
                                                ..1 == "mod" & isFALSE(..2) ~ "white.anova.mod",
                                                ..1 == "mod.boxcox" & isTRUE(..3) ~ "anova.mod.boxcox",
                                                ..1 == "mod.boxcox" & isFALSE(..3) ~ "white.anova.mod.boxcox",
                                                ..1 == "mod.log" & isTRUE(..4) ~ "anova.mod.log",
                                                ..1 == "mod.log" & isFALSE(..4) ~ "white.anova.mod.log",
                                                ..1 == "none" ~ "none")),
           best.anova.val = pmap(list(best.anova.nm, #1
                                      anova.mod, #2
                                      white.anova.mod, #3
                                      anova.mod.boxcox, #4
                                      white.anova.mod.boxcox, #5
                                      anova.mod.log, #6
                                      white.anova.mod.log), #7
                                 ~ case_when(..1 == "anova.mod" ~ list(tidy(..2)),
                                             ..1 == "white.anova.mod" ~ list(tidy(..3)),
                                             ..1 == "anova.mod.boxcox" ~ list(tidy(..4)),
                                             ..1 == "white.anova.mod.boxcox" ~ list(tidy(..5)),
                                             ..1 == "anova.mod.log" ~ list(tidy(..6)),
                                             ..1 == "white.anova.mod.log" ~ list(tidy(..7)),
                                             TRUE ~ list(NULL))[[1]]),
           sign.terms.best.anova = map_chr(best.anova.val,
                                           ~ if(is.null(.x)) {
                                             NA_character_
                                           } else {
                                             filter(.x, !term %in% c("(Intercept)", "Residuals"), p.value <= 0.05) %>%
                                               select(term) %>% {.[[1]] %>% str_c(collapse = " / ")} %>%
                                               {if_else(. == "", "none", .)}
                                           }),
           sign.interaction.term.best.anova = map_chr(sign.terms.best.anova,
                                                      ~ if_else(is.na(.x), "NA", .x) %>%
                                                        str_split(" / ") %>%
                                                        {.[[1]]} %>%
                                                        as.list() %>%
                                                        keep(.p = function(x) str_detect(x, "^.*\\:.*$")) %>%
                                                        {.[1][[1]]} %>%
                                                        {if_else(is.null(.), NA_character_, .)}),
           target.terms.posthoc = pmap(list(sign.terms.best.anova,
                                            sign.interaction.term.best.anova,
                                            data),
                                       ~ if (is.na(..1) | ..1 == "none") {
                                         NA_character_
                                       } else if (is.na(..2)) {
                                         ..3[, ..1 %>% str_split(" / ") %>% {.[[1]]}] %>%
                                           select_if(is.factor) %>%
                                           colnames() %>%
                                           {if (is_empty(.)) { NA_character_ } else { . }}
                                       } else if (!is.na(..2)) {
                                         ..3[, str_split(..2, ":") %>% {.[[1]]}] %>%
                                           select_if(is.factor) %>%
                                           colnames() %>%
                                           str_c(collapse = ":")
                                       }),
           posthoc.nm = map2_chr(target.terms.posthoc,
                                 best.anova.nm,
                                 ~ if (length(.x) == 1L && is.na(.x)) {
                                   NA_character_
                                 } else {
                                   case_when(.y %in% c("anova.mod", "anova.mod.boxcox", "anova.mod.log") ~ str_c("EMM comparison with ", adjust.emm.comp, " adjustment"),
                                             .y %in% c("white.anova.mod", "white.anova.mod.boxcox", "white.anova.mod.log") ~ "Games-Howell test",
                                             TRUE ~ NA_character_)
                                 }),
           posthoc.val = pmap(list(target.terms.posthoc,
                                   posthoc.nm,
                                   sign.interaction.term.best.anova,
                                   best.mod.val),
                              ~ if(length(..1) == 1L && is.na(..1)) {
                                NULL
                              } else if (str_detect(..2, "EMM") & is.na(..3)) {
                                emmeans(..4, as.list(..1)) %>% pairs(adjust = adjust.emm.comp) %>% after.emmeans
                              } else if (str_detect(..2, "EMM") & !is.na(..3)) {
                                if (isFALSE(spe.contr)) {
                                  emmeans(..4, as.formula(str_c("pairwise ~ ", ..3)), adjust = adjust.emm.comp)$contrasts %>% as_tibble
                                } else {
                                  emmeans(..4, specs = specs, by = by) %>% contrast(method = "pairwise", simple = list(specs), combine = TRUE, adjust = adjust.emm.comp) %>% as_tibble
                                }
                              } else if (is.na(..3)) {
                                pool <- vector("list", length(..1))
                                for (i in seq_along(..1)) {
                                  pool[[i]] <- games_howell_test(data = ..4$model, formula = as.formula(str_c(colnames(..4$model)[[1]], " ~ ", ..1[[i]])))
                                }
                                bind_rows(pool)
                              } else if (!is.na(..3)) {
                                mutate(..4$model,
                                       interaction.term = interaction(..4$model[[..3 %>% str_split(":") %>% {.[[1]][[1]]}]], ..4$model[[..3 %>% str_split(":") %>% {.[[1]][[2]]}]])) %>%
                                  games_howell_test(formula = as.formula(str_c(colnames(.)[[1]], " ~ ", "interaction.term")))
                              }),
           posthoc.summary = pmap(list(target.terms.posthoc, #1
                                       posthoc.nm, #2
                                       posthoc.val, #3
                                       sign.interaction.term.best.anova), #4
                                  ~ if (length(..1) == 1L && is.na(..1)) {
                                    NULL
                                  } else if (str_detect(..2, "EMM")) {
                                    if (isFALSE(spe.contr) | is.na(..4)) {
                                      ..3 %>% select(comparison = contrast, estimate, adj.p.value = p.value)
                                    } else {
                                      ..3 %>% mutate(comparison = str_c(.[[by]], contrast, sep = ": ")) %>% select(comparison, estimate, adj.p.value = p.value)
                                    }
                                  } else if (is.na(..4)) {
                                    ..3 %>% mutate(comparison = str_c(group1, group2, sep = " - ")) %>%
                                      select(comparison, estimate, adj.p.value = p.adj)
                                  } else {
                                    ..3 %>% mutate(comparison = str_c(group1, group2, sep = " - ") %>% str_replace_all("\\.", " ")) %>%
                                      select(comparison, estimate, adj.p.value = p.adj)
                                  }),
           is.sign.posthoc = map_lgl(posthoc.summary,
                                     ~ if (is.null(.x)) {
                                       NA
                                     } else {
                                       some(.x$adj.p.value, ~ .x <= 0.05)
                                     }),
           sign.terms.posthoc = map_chr(posthoc.summary,
                                        ~ if (is.null(.x)) {
                                          NA_character_
                                        } else {
                                          .x %>% filter(adj.p.value <= 0.05) %>%
                                            mutate(adj.p.value = formatC(adj.p.value, format = "e", digits = 2),
                                                   sign.terms = str_c(comparison, ": ", adj.p.value)) %>% select(sign.terms) %>% {.[[1]] %>% str_c(collapse = "; ")} %>% {if_else(. == "", NA_character_, .)}
                                        }),
           y.label = pmap_chr(list(!!!y.label), ~ str_c(..., sep = " - ")),

           p.values.best.anova = map2_chr(sign.terms.best.anova, best.anova.val,
                                          ~ if (is.na(.x) | .x == "none") {
                                            NA_character_
                                          } else {
                                            .y %>%
                                              select(term, p.value) %>%
                                              filter(!term %in% c("(Intercept)", "Residuals", block)) %>%
                                              mutate(sign = case_when(p.value > 0.05 ~ "ns",
                                                                      p.value <= 0.0001 ~ "****",
                                                                      p.value <= 0.001 ~ "***",
                                                                      p.value <= 0.01 ~ "**",
                                                                      p.value <= 0.05 ~ "*"),
                                                     p.val = str_c(term, ": ", formatC(p.value, format = "e", digits = 2), " ", sign)) %>%
                                              select(p.val) %>%
                                              {str_c(.[[1]], collapse = "\n")}
                                          }),

           plot = pmap(list(sign.terms.best.anova, #1
                            mod, #2
                            y.label, #3
                            p.values.best.anova), #4
                       ~ if (is.na(..1) | ..1 == "none") {
                         ggplot(data = ..2$mod, mapping = aes(x = .data[[indep.vars[[1]]]], y = .data[[lhs]])) +
                           geom_boxplot(aes(fill = .data[[indep.vars[[2]]]]), alpha = 0.4, outlier.shape = NA) +
                           geom_point(aes(color = .data[[indep.vars[[2]]]]), position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0), alpha = 0.7) +
                           scale_colour_manual(values = c("#5d76cb", "#ca3767", "#3bb08f", "#ffcf48"), aesthetics = c("fill", "colour")) +
                           labs(x = NULL, y = ..3, caption = str_c("ANOVA not performed", str_c(rep(x = "\n", times = jumps), collapse = ""))) +
                           theme_pubr(legend = legend, base_size = base.size)
                       } else {
                         ggplot(data = ..2$mod, mapping = aes(x = .data[[indep.vars[[1]]]], y = .data[[lhs]])) +
                           geom_boxplot(aes(fill = .data[[indep.vars[[2]]]]), alpha = 0.4, outlier.shape = NA) +
                           geom_point(aes(color = .data[[indep.vars[[2]]]]), position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0), alpha = 0.7) +
                           scale_colour_manual(values = c("#5d76cb", "#ca3767", "#3bb08f", "#ffcf48"), aesthetics = c("fill", "colour")) +
                           labs(x = NULL, y = ..3, caption = str_c("p values:", "\n", ..4)) + # "\n" to add/suppress
                           theme_pubr(legend = legend, base_size = base.size)
                       })) %>%
    select(all_of(grp),
           matches("^.*variable.*$"),
           formula.analysis,
           best.mod.nm,
           best.anova.nm,
           sign.terms.best.anova,
           sign.interaction.term.best.anova,
           sign.terms.posthoc,
           everything())

  options(contrasts = c("contr.treatment", "contr.poly"))

  out

}
