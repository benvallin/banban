#' Adjust choose_on_aic() output for multiple comparisons
#'
#' @description get_sign_mod() adjusts the p value related to the F statistic of each "best model" provided by choose_on_aic().
#'
#' The function first filters out every best model containing only the intercept (1) as independent variable (dependent variable ~ mean(dependent variable)).
#'
#' Adjustment for multiple comparisons then relies either on false discovery rate (FDR) using the Benjamini-Hochberg procedure, or on one of the following adjustment methods: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#'
#' Best models receive a label "best_is_sign" = TRUE if the p value related to their F statistic reaches statistical significance after adjustment (passed FDR correction, or adjusted value <= 0.05) AND the p value of the t statistic related to the independent variable of interest <= 0.05.
#'
#' @param tibble a tibble produced by choose_on_aic().
#' @param method a string indicating the adjustment method to be used. Must be one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and "none".
#' @param FDR a numeric indicating the q-value threshold to use for FDR.
#'
#' @return A tibble containing, for each combination of dependent variable, independent variable of interest and covariates, i) the corresponding best model, ii) its associated statistics (F p value, t p value) and iii) the label "best_is_sign" = TRUE/FALSE.
#'
#' Combinations of dependent variable, independent variable of interest and covariates for which the best model contains only the intercept as independent variable are discarded, so the output may be shorter than the input tibble.
#'
#' @export
#'
#' @importFrom dplyr case_when
#' @importFrom dplyr filter
#' @importFrom dplyr if_else
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom magrittr %>%
#' @importFrom mcp.project fdr
#' @importFrom purrr map
#' @importFrom purrr map_chr
#' @importFrom rlang !!
#' @importFrom rlang enexpr
#' @importFrom stats p.adjust
#' @importFrom stringr str_c
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#' @importFrom tidyselect everything
#'
#' @examples
#'
get_sign_mod <- function(tibble, method, FDR = 0.2) {

  method <- enexpr(method)

  output <- tibble %>% filter(!is.na(best_mod_f_p_value))

  if (method == "BH") {
    output <- output %>%
      nest(data = everything()) %>%
      mutate(fdr = map(data, ~ fdr(.x$best_mod_f_p_value, q = FDR, method = "BH")$Pvals)) %>%
      unnest(cols = c(data, fdr)) %>%
      rename(pass_fdr = rejected) %>%
      mutate(fdr_result = map_chr(pass_fdr, ~ if_else(isTRUE(.x),
                                                      str_c("sign (FDR = ", FDR, ")"),
                                                      str_c("ns (FDR = ", FDR, ")"))),
             best_is_sign = case_when(best_mod_t_p_value <= 0.05 & pass_fdr == TRUE ~ TRUE,
                                      TRUE ~ FALSE))
  } else {
    output <- output %>%
      mutate(adj_best_mod_f_p_value = p.adjust(best_mod_f_p_value, method = !!method),
             best_is_sign = case_when(best_mod_t_p_value <= 0.05 & adj_best_mod_f_p_value <= 0.05 ~ TRUE,
                                      TRUE ~ FALSE))
  }
}
