#' Perform automated ANOVA with random effect and posthoc test
#'
#' @description get_mixed_mod_anova() performs automated ANOVA including random effect(s) and subsequent posthoc test on a set of dependent variables.
#'
#' @details Normality of residuals is tested with Shapiro–Wilk test. In case of normality, the model is used for all subsequent analyses.
#' In case of non-normality, the dependent variable is subjected to log transformation and normality of residuals is tested again with Shapiro–Wilk test.
#' Currently, the log-transformed model is used for all subsequent analyses, irrespective of the result of the Shapiro-Wilk test. This should be modified in the future, so that the dependent variable is further analysed only if the distribution of residuals becomes normal after log transformation.
#'
#' Additionally, normality and homoscedasticity can be graphically inspected through 3 types of plot: residuals vs observed values, fitted values vs residuals, and qq plot.
#'
#' Type I ANOVA is performed with the Kenward-Roger's approximation to degrees of freedom.
#'
#' Factors with a significant effect in ANOVA are further subjected to posthoc test.
#' In case of significant interaction effect(s), only the significant interaction term(s) are subjected to posthoc analysis.
#' Posthoc test is performed by comparison of estimated marginal means and one of the following adjustment methods: "tukey", "scheffe", "sidak", "bonferroni", "dunnettx", "mvt", and "none".
#'
#' @param tibble a tibble.
#' @param grp a string indicating the column which contains the names of the dependent variables.
#' @param lhs a string indicating the column which contains the values of the dependent variables.
#' @param rhs a string indicating the name(s) of the independent variable(s) and random factor(s). Must follow the classic formula rules (use "+", "*" and ":" as separators, and see ?lmerTest::lmer for notations of random effects).
#' @param base.size a numeric provided to the base_size argument of theme_pubr() for plotting.
#' @param legend a string provided to the legend argument of theme_pubr() for plotting. Must be one of "top", "bottom", "left", "right", "none".
#' @param spe.contr a logical indicating if estimated marginal means are desired over a specific independent variable in case of significant interaction. Default is FALSE.
#' @param specs if spe.contr = TRUE, a string indicating the name of the independent variable over which estimated marginal means are desired. Will be used in case of significant interaction between the independent variables passed to the specs and by arguments. Default is NA_character_.
#' @param by if spe.contr = TRUE, a string indicating the name of the independent variable to condition on. Will be used in case of significant interaction between the independent variables passed to the specs and by arguments. Default is NA_character_.
#' @param adjust.emm.comp a string indicating which adjustment method should be used after computation of estimated marginal means. Must be one of "tukey", "scheffe", "sidak", "bonferroni", "dunnettx", "mvt", and "none". Default is "tukey".
#' @param block a character vector indicating the name(s) of the independent variable(s) for which the ANOVA-associated p values should not be displayed on the plot. Default is NA_character_.
#'
#' @return A tibble containing, for each dependent variable, an ANOVA (possibly with random factor) and associated posthoc test as well as the corresponding plot.
#' Also included are residuals vs observed values plot, fitted values vs residuals plot, and qq plot.
#'
#' @export
#'
#' @importFrom dplyr bind_rows
#' @importFrom dplyr case_when
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
#' @importFrom ggplot2 stat_qq
#' @importFrom ggplot2 stat_qq_line
#' @importFrom ggpubr theme_pubr
#' @importFrom graphics pairs
#' @importFrom lmerTest lmer
#' @importFrom magrittr %>%
#' @importFrom purrr keep
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom purrr map2_chr
#' @importFrom purrr map_chr
#' @importFrom purrr map_dbl
#' @importFrom purrr map_lgl
#' @importFrom purrr modify_at
#' @importFrom purrr pmap
#' @importFrom purrr pmap_chr
#' @importFrom purrr some
#' @importFrom rlang !!
#' @importFrom rlang :=
#' @importFrom rlang enexpr
#' @importFrom rlang parse_expr
#' @importFrom rlang syms
#' @importFrom stats anova
#' @importFrom stats as.formula
#' @importFrom stats fitted
#' @importFrom stats residuals
#' @importFrom stats shapiro.test
#' @importFrom stringr str_c
#' @importFrom stringr str_detect
#' @importFrom stringr str_remove_all
#' @importFrom stringr str_split
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' @importFrom tidyr nest
#' @importFrom tidyselect all_of
#' @importFrom tidyselect everything
#' @importFrom tidyselect matches
#'
#' @examples
#'
get_mixed_mod_anova <- function(tibble, grp = NULL, lhs, rhs, base.size = 14, legend = "none", spe.contr = FALSE, specs = NA_character_, by = NA_character_, adjust.emm.comp = "tukey", block = NA_character_) {

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
  lhs.log <- str_c("log.", lhs)
  formula.log <- parse_expr(str_c(lhs.log, " ~ ", rhs))

  indep.vars <- str_split(string = as.character(rhs), pattern = "\\+|\\*|\\:")[[1]] %>% str_remove_all(pattern = "\\s") %>% unique

  y.label <- syms(grp)

  out <- tibble %>%
    modify_at(.at = indep.vars, .f = function(x) if(is.character(x)) as_factor(x) else x) %>%
    mutate(!!str_c("log.", lhs) := log(.[[lhs]])) %>%
    group_by(.dots = grp) %>%
    group_modify(~ nest(.data = .x, data = everything())) %>%
    ungroup() %>%
    mutate(dependent.variable = lhs,
           independent.variables = rhs,
           # mod
           mod = map(data, ~ lmer(!!formula, data = .x)),
           shapiro.mod = map(mod, ~ shapiro.test(residuals(.x))),
           p.value.shapiro.mod = map_dbl(shapiro.mod, ~ .x$p.value),
           normality.mod = map_lgl(p.value.shapiro.mod, ~ if_else(.x > 0.05, TRUE, FALSE)),
           p.mod.res.vs.obs = map2(mod, data, ~ ggplot(data = tibble(res = residuals(.x), obs = .y[[lhs]]),
                                                       aes(x = res, y = obs)) + geom_point()),
           p.mod.fit.vs.res = map(mod, ~ ggplot(data = tibble(fit = fitted(.x), res = residuals(.x)),
                                                aes(x = fit, y = res)) + geom_point()),
           p.mod.qq = map(mod, ~ ggplot(data = tibble(res = residuals(.x)),
                                        aes(sample = res))
                          + stat_qq() + stat_qq_line() + labs(x = "Standard normal quantiles", y = "Standardized residuals")),
           anova.mod = map(mod, ~ anova(.x, ddf = "Kenward-Roger")),
           sign.terms.anova.mod = map_chr(anova.mod,
                                          ~ tidy(.x) %>%
                                            filter(p.value <= 0.05) %>% select(term) %>% {.[[1]] %>% str_c(collapse = " / ")} %>%
                                            {if_else(. == "", "none", .)}),
           sign.interaction.term.anova.mod = map_chr(sign.terms.anova.mod,
                                                     ~ if_else(is.na(.x), "NA", .x) %>%
                                                       str_split(" / ") %>%
                                                       {.[[1]]} %>%
                                                       as.list() %>%
                                                       keep(.p = function(x) str_detect(x, "^.*\\:.*$")) %>%
                                                       {.[1][[1]]} %>%
                                                       {if_else(is.null(.), NA_character_, .)}),
           target.terms.posthoc.anova.mod = map2(sign.interaction.term.anova.mod,
                                                 sign.terms.anova.mod,
                                                 ~ case_when(is.na(.x) & .y != "none" ~ .y,
                                                             is.na(.x) & .y == "none" ~ NA_character_,
                                                             !is.na(.x) ~ .x) %>%
                                                   str_split(" / ") %>% {.[[1]]}),
           posthoc.anova.mod = pmap(list(target.terms.posthoc.anova.mod, #1
                                         sign.interaction.term.anova.mod, #2
                                         mod), #3
                                    ~ if(is.na(..1)) {
                                      NULL
                                    } else if (is.na(..2)) {
                                      emmeans(..3, as.list(..1)) %>% pairs(adjust = adjust.emm.comp) %>% after.emmeans
                                    } else {
                                      if (isFALSE(spe.contr)) {
                                        emmeans(..3, as.formula(str_c("pairwise ~ ", ..2)), adjust = adjust.emm.comp)$contrasts %>% as_tibble
                                      } else {
                                        emmeans(..3, specs = specs, by = by) %>% contrast(method = "pairwise", simple = list(specs), combine = TRUE, adjust = adjust.emm.comp) %>% as_tibble
                                      }
                                    }),
           posthoc.anova.mod.summary = pmap(list(target.terms.posthoc.anova.mod, #1
                                                 posthoc.anova.mod, #2
                                                 sign.interaction.term.anova.mod), #3
                                            ~ if (is.na(..1)) {
                                              NULL
                                            } else if (isFALSE(spe.contr) | is.na(..3)) {
                                              ..2 %>% select(comparison = contrast, estimate, adj.p.value = p.value)
                                            } else {
                                              ..2 %>% mutate(comparison = str_c(.[[by]], contrast, sep = ": ")) %>% select(comparison, estimate, adj.p.value = p.value)
                                            }),
           is.sign.posthoc.anova.mod = map_lgl(posthoc.anova.mod.summary,
                                               ~ if (is.null(.x)) {
                                                 NA
                                               } else {
                                                 some(.x$adj.p.value, ~ .x <= 0.05)
                                               }),
           sign.terms.posthoc.anova.mod = map_chr(posthoc.anova.mod.summary,
                                                  ~ if (is.null(.x)) {
                                                    NA_character_
                                                  } else {
                                                    .x %>% filter(adj.p.value <= 0.05) %>%
                                                      mutate(adj.p.value = formatC(adj.p.value, format = "e", digits = 2),
                                                             sign.terms = str_c(comparison, ": ", adj.p.value)) %>% select(sign.terms) %>% {.[[1]] %>% str_c(collapse = "; ")} %>% {if_else(. == "", NA_character_, .)}
                                                  }),
           # log
           mod.log = map(data, ~ lmer(!!formula.log, data = .x)),
           shapiro.mod.log = map(mod.log, ~ shapiro.test(residuals(.x))),
           p.value.shapiro.mod.log = map_dbl(shapiro.mod.log, ~ .x$p.value),
           normality.mod.log = map_lgl(p.value.shapiro.mod.log, ~ if_else(.x > 0.05, TRUE, FALSE)),
           p.mod.res.vs.obs = map2(mod.log, data, ~ ggplot(data = tibble(res = residuals(.x), obs = .y[[lhs.log]]),
                                                           aes(x = res, y = obs)) + geom_point()),
           p.mod.log.fit.vs.res = map(mod.log, ~ ggplot(data = tibble(fit = fitted(.x), res = residuals(.x)),
                                                        aes(x = fit, y = res)) + geom_point()),
           p.mod.log.qq = map(mod.log, ~ ggplot(data = tibble(res = residuals(.x)),
                                                aes(sample = res))
                              + stat_qq() + stat_qq_line() + labs(x = "Standard normal quantiles", y = "Standardized residuals")),
           anova.mod.log = map(mod.log, ~ anova(.x, ddf = "Kenward-Roger")),
           sign.terms.anova.mod.log = map_chr(anova.mod.log,
                                              ~ tidy(.x) %>%
                                                filter(p.value <= 0.05) %>% select(term) %>% {.[[1]] %>% str_c(collapse = " / ")} %>%
                                                {if_else(. == "", "none", .)}),
           sign.interaction.term.anova.mod.log = map_chr(sign.terms.anova.mod.log,
                                                         ~ if_else(is.na(.x), "NA", .x) %>%
                                                           str_split(" / ") %>%
                                                           {.[[1]]} %>%
                                                           as.list() %>%
                                                           keep(.p = function(x) str_detect(x, "^.*\\:.*$")) %>%
                                                           {.[1][[1]]} %>%
                                                           {if_else(is.null(.), NA_character_, .)}),
           target.terms.posthoc.anova.mod.log = map2(sign.interaction.term.anova.mod.log,
                                                     sign.terms.anova.mod.log,
                                                     ~ case_when(is.na(.x) & .y != "none" ~ .y,
                                                                 is.na(.x) & .y == "none" ~ NA_character_,
                                                                 !is.na(.x) ~ .x) %>%
                                                       str_split(" / ") %>% {.[[1]]}),
           posthoc.anova.mod.log = pmap(list(target.terms.posthoc.anova.mod.log, #1
                                             sign.interaction.term.anova.mod.log, #2
                                             mod.log), #3
                                        ~ if(is.na(..1)) {
                                          NULL
                                        } else if (is.na(..2)) {
                                          emmeans(..3, as.list(..1)) %>% pairs(adjust = adjust.emm.comp) %>% after.emmeans
                                        } else {
                                          if (isFALSE(spe.contr)) {
                                            emmeans(..3, as.formula(str_c("pairwise ~ ", ..2)), adjust = adjust.emm.comp)$contrasts %>% as_tibble
                                          } else {
                                            emmeans(..3, specs = specs, by = by) %>% contrast(method = "pairwise", simple = list(specs), combine = TRUE, adjust = adjust.emm.comp) %>% as_tibble
                                          }
                                        }),
           posthoc.anova.mod.log.summary = pmap(list(target.terms.posthoc.anova.mod.log, #1
                                                     posthoc.anova.mod.log, #2
                                                     sign.interaction.term.anova.mod.log), #3
                                                ~ if (is.na(..1)) {
                                                  NULL
                                                } else if (isFALSE(spe.contr) | is.na(..3)) {
                                                  ..2 %>% select(comparison = contrast, estimate, adj.p.value = p.value)
                                                } else {
                                                  ..2 %>% mutate(comparison = str_c(.[[by]], contrast, sep = ": ")) %>% select(comparison, estimate, adj.p.value = p.value)
                                                }),
           is.sign.posthoc.anova.mod.log = map_lgl(posthoc.anova.mod.log.summary,
                                                   ~ if (is.null(.x)) {
                                                     NA
                                                   } else {
                                                     some(.x$adj.p.value, ~ .x <= 0.05)
                                                   }),
           sign.terms.posthoc.anova.mod.log = map_chr(posthoc.anova.mod.log.summary,
                                                      ~ if (is.null(.x)) {
                                                        NA_character_
                                                      } else {
                                                        .x %>% filter(adj.p.value <= 0.05) %>%
                                                          mutate(adj.p.value = formatC(adj.p.value, format = "e", digits = 2),
                                                                 sign.terms = str_c(comparison, ": ", adj.p.value)) %>% select(sign.terms) %>% {.[[1]] %>% str_c(collapse = "; ")} %>% {if_else(. == "", NA_character_, .)}
                                                      }),
           # summary
           best.mod.nm = case_when(normality.mod == TRUE ~ "mod", TRUE ~ "mod.log"),
           formula.analysis = map2_chr(best.mod.nm, #1
                                       dependent.variable, #2
                                       ~ if_else(.x == "mod", str_c(..2, " ~ ", rhs), str_c("log(", ..2, ") ~ ", rhs))),
           best.mod.val = case_when(normality.mod == TRUE ~ mod, TRUE ~ mod.log),
           best.anova.val = case_when(normality.mod == TRUE ~ anova.mod, TRUE ~ anova.mod.log),
           best.anova.sign.terms = case_when(normality.mod == TRUE ~ sign.terms.anova.mod, TRUE ~ sign.terms.anova.mod.log),
           best.posthoc.target.terms = case_when(normality.mod == TRUE ~ target.terms.posthoc.anova.mod, TRUE ~ target.terms.posthoc.anova.mod.log),
           best.posthoc.val = case_when(normality.mod == TRUE ~ posthoc.anova.mod, TRUE ~ posthoc.anova.mod.log),
           best.posthoc.summary = case_when(normality.mod == TRUE ~ posthoc.anova.mod.summary, TRUE ~ posthoc.anova.mod.log.summary),
           best.posthoc.sign.terms = case_when(normality.mod == TRUE ~ sign.terms.posthoc.anova.mod, TRUE ~ sign.terms.posthoc.anova.mod.log),
           y.label = pmap_chr(list(!!!y.label), ~ str_c(..., sep = " - ")),
           p.values.best.anova = map_chr(best.anova.val,
                                         ~ tidy(.x) %>%
                                           select(term, p.value) %>%
                                           filter(!term %in% c(block)) %>%
                                           mutate(sign = case_when(p.value > 0.05 ~ "ns",
                                                                   p.value <= 0.0001 ~ "****",
                                                                   p.value <= 0.001 ~ "***",
                                                                   p.value <= 0.01 ~ "**",
                                                                   p.value <= 0.05 ~ "*"),
                                                  p.val = str_c(term, ": ", formatC(p.value, format = "e", digits = 2), " ", sign)) %>%
                                           select(p.val) %>%
                                           {str_c(.[[1]], collapse = "\n")}),
           plot = pmap(list(data, #1
                            y.label, #2
                            p.values.best.anova), #4
                       ~ ggplot(data = ..1, mapping = aes(x = .data[[indep.vars[[1]]]], y = .data[[lhs]])) +
                         geom_boxplot(aes(fill = .data[[indep.vars[[2]]]]), alpha = 0.4, outlier.shape = NA) +
                         geom_point(aes(color = .data[[indep.vars[[2]]]]), position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0), alpha = 0.7) +
                         scale_colour_manual(values = c("#5d76cb", "#ca3767", "#3bb08f", "#ffcf48"), aesthetics = c("fill", "colour")) +
                         labs(x = NULL, y = ..2, caption = str_c("p values:", "\n", ..3)) + # "\n" to add/suppress
                         theme_pubr(legend = legend, base_size = base.size))) %>%
    select(all_of(grp),
           matches("^.*variable.*$"),
           formula.analysis,
           matches("^.*best.*$"),
           everything())

  options(contrasts = c("contr.treatment", "contr.poly"))

  out
}
