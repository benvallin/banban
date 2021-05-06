#' Blabla
#'
#' @param tibble blabla
#' @param grp blabla
#' @param lhs blabla
#' @param rhs blabla
#' @param base.size blabla
#' @param legend blabla
#' @param spe.contr blabla
#' @param specs blabla
#' @param by blabla
#' @param adjust.emm.comp blabla
#' @param block blabla
#' @param jumps blabla
#'
#' @return blabla
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
get_complex_anova_beta <- function(tibble, grp = NULL, lhs, rhs, base.size = 14, legend = "none", spe.contr = FALSE, specs = NA_character_, by = NA_character_, adjust.emm.comp = "tukey", block = NA_character_, jumps = 0) {

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
           target.terms.posthoc = map2(sign.interaction.term.best.anova,
                                       sign.terms.best.anova,
                                       ~ case_when(is.na(.x) & .y != "none" ~ .y,
                                                   is.na(.x) & .y == "none" ~ NA_character_,
                                                   !is.na(.x) ~ .x) %>%
                                         str_split(" / ") %>% {.[[1]]}),
           posthoc.nm = map_chr(best.anova.nm,
                                ~ case_when(.x %in% c("anova.mod", "anova.mod.boxcox", "anova.mod.log") ~ str_c("EMM comparison with ", adjust.emm.comp, " adjustment"),
                                            .x %in% c("white.anova.mod", "white.anova.mod.boxcox", "white.anova.mod.log") ~ "Games-Howell test",
                                            TRUE ~ NA_character_)),
           posthoc.val = pmap(list(target.terms.posthoc, #1
                                   posthoc.nm, #2
                                   sign.interaction.term.best.anova, #3
                                   best.mod.val), #4
                              ~ if(is.na(..1)) {
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
                                  pool[[i]] <- games_howell_test(data = ..4$mod, formula = as.formula(str_c(colnames(..4$mod)[[1]], " ~ ", .x[[i]])))
                                }
                                bind_rows(pool)
                              } else if (!is.na(..3)) {
                                mutate(..4$mod,
                                       interaction.term = interaction(..4$mod[[..3 %>% str_split(":") %>% {.[[1]][[1]]}]], ..4$mod[[..3 %>% str_split(":") %>% {.[[1]][[2]]}]])) %>%
                                  games_howell_test(formula = as.formula(str_c(colnames(.)[[1]], " ~ ", "interaction.term")))
                              }),
           posthoc.summary = pmap(list(target.terms.posthoc, #1
                                       posthoc.nm, #2
                                       posthoc.val, #3
                                       sign.interaction.term.best.anova), #4
                                  ~ if (is.na(..1)) {
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
