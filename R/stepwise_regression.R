#' Perform automated stepwise regression
#'
#' @description stepwise_regression() performs automated stepwise regression on a set of models differing by an independent variable of interest.
#' All original full models share the same dependent variable and covariates, but differ by their independent variable of interest.
#'
#' From each original full model, a minimal model (dependent variable ~ independent variable of interest), a forward selection model, a stepwise selection model and a backward elimination model are generated. For each independent variable of interest, 5 distinct models are therefore produced (original full, minimal, forward selection, stepwise selection, backward elimination).
#'
#' For all models, the adjusted R square, the F statistic, and the t statistic related to the independent variable of interest are computed.
#' For stepwise regression models only (forward selection, stepwise selection, backward elimination), the AIC value is also computed.
#' Note that the AIC and F values are computed only if df > 0, which implies n.obs > n.independent.variables +1 (intercept).
#'
#' Models which display an F p value <= 0.05 AND a t p value (related to the independent variable of interest) <= 0.05 are labelled sign = TRUE.
#'
#' @param tibble a tibble containing ONLY i) the column with the values of the dependent variable, ii) the column with the names of the independent variables of interest, iii) the column with the values of the independent variables of interest, iv) the column(s) with the categorical variable(s) on which the independent variables of interest should be conditioned, and v) the column(s) with the values of the covariate(s).
#' @param dep_var a string indicating the column with the values of the dependent variable.
#' @param indep_oi a string indicating the column with the values of the independent variables of interest.
#' @param grps a character vector indicating i) the column with the names of the independent variables of interest, and optionally ii) the column(s) with the categorical variable(s) on which the independent variables of interest should be conditioned.
#'
#' @return A tibble containing, for each (possibly conditioned) independent variable of interest, all five models (original full, minimal, forward selection, stepwise selection, backward elimination) and their respective statistics (AIC, adjusted R square, F statistic, t statistic related to the independent variable of interest, "sign" label).
#'
#' @export
#'
#' @importFrom MASS stepAIC
#' @importFrom dplyr case_when
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr na_if
#' @importFrom dplyr select
#' @importFrom dplyr ungroup
#' @importFrom generics setdiff
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom purrr map2_dbl
#' @importFrom purrr map2_lgl
#' @importFrom purrr map_chr
#' @importFrom purrr map_dbl
#' @importFrom purrr map_int
#' @importFrom rlang enexpr
#' @importFrom stats as.formula
#' @importFrom stats extractAIC
#' @importFrom stats lm
#' @importFrom stats na.omit
#' @importFrom stats pf
#' @importFrom stringr str_c
#' @importFrom stringr str_detect
#' @importFrom stringr str_split
#' @importFrom stringr str_subset
#' @importFrom tidyr nest
#' @importFrom tidyselect all_of
#' @importFrom tidyselect everything
#'
#' @examples
#'
stepwise_regression <- function(tibble, dep_var, indep_oi, grps, hackfull = FALSE) {

  # record user's dependant variable Y
  # dep_var <- enexpr(dep_var) %>% as.character()
  dep_var <- as.character(dep_var)

  # record user's independant variable of interest Xoi
  indep_oi <- enexpr(indep_oi) %>% as.character()

  # create a formula for the minimal model of interest: Y ~ Xoi
  min_oi <- str_c(dep_var, "~", indep_oi) %>% as.formula()

  # create a string expressing the formula for the minimal possible model: Y ~ 1 (Y ~ mean(Y))
  min <- str_c(dep_var, "~", "1")

  # create a formula for the full model: Y ~ Xoi + Xn
  # the covariates Xn correspond to all the variables in the input table which are not defined as dep_var, indep_oi or grps
  # max <- str_c(dep_var,
  #              "~",
  #              str_c(setdiff(setdiff(colnames(tibble), dep_var), grps), collapse = "+")) %>% as.formula()

  # create a string pattern to look for the presence of Xoi in the RHS of the returned stepwise regression models
  look <- str_c("^.*", indep_oi, ".*$")

  # compute all models & related variables in user's table
  output <- tibble %>%
    group_by(.dots = grps) %>%
    nest() %>%
    mutate(dep_var = dep_var,
           indep_var_oi = indep_oi,
           # START NEW
           data = map(data, ~ select_if(.x, .predicate = function(x) if_else(every(x, is.na), FALSE, TRUE))),
           max = map2(dep_var, data, ~ str_c(.x, "~", str_c(setdiff(colnames(.y), dep_var), collapse = "+")) %>% as.formula),
           n_obs = map_int(data, ~ length(na.omit(.)[[1]])),
           n_terms_full = map_int(data, ~ length(setdiff(colnames(.), dep_var))),
           categorical_terms = map(data, ~ select_if(.x, .predicate = function(x) if_else((is.factor(x) | is.character(x)), TRUE, FALSE))),
           n_categorical_terms = map_int(categorical_terms, ~ if_else(is_empty(.x), 0L, length(.x))),
           n_levels_categorical_terms = map2_int(n_categorical_terms,
                                                 categorical_terms,
                                                 ~ if (.x > 0L) { .y %>% mutate_all(n_distinct) %>% unique %>% reduce(`+`) } else { 0L }),
           df_took_by_cat = n_levels_categorical_terms - n_categorical_terms,
           df_took_by_all_terms = (n_terms_full - n_categorical_terms) + df_took_by_cat,
           # resid_df_sup_0 = map2_lgl(n_obs, n_terms_full, ~ case_when(.x > (.y + 1) ~ TRUE, TRUE ~ FALSE)),
           resid_df_sup_0 = map2_lgl(n_obs, df_took_by_all_terms, ~ if_else(.x > (.y + 1), TRUE, FALSE)),
           # END NEW
           minoi_mod = map(data, ~ lm(min_oi, data = .)),
           minoi_terms = map_chr(minoi_mod, ~ as.character(.$terms[[3]])),
           n_terms_minoi = 1,
           indepoi_in_minoi = TRUE,
           minoi_t_p_value = map_dbl(minoi_mod, ~ summary(.)$coefficients[indep_oi, 4]),
           minoi_f_p_value = map_dbl(minoi_mod, ~ summary(.)$fstatistic %>% {pf(.[[1]], .[[2]], .[[3]], lower.tail = FALSE)}),
           minoi_r_squared = map_dbl(minoi_mod, ~ summary(.)$r.squared),
           minoi_sign = case_when(minoi_t_p_value <= 0.05 ~ TRUE,
                                  TRUE ~ FALSE),
           # START NEW
           full_mod = map2(max, data, ~ lm(.x, data = .y)),
           # END NEW
           full_terms = map(full_mod, ~ .$terms[[3]]) %>% as.character(),
           indepoi_in_full = TRUE,
           full_t_p_value = map_dbl(full_mod, ~ summary(.)$coefficients[indep_oi, 4]),
           full_f_p_value = map_dbl(full_mod, ~ summary(.)$fstatistic %>% {pf(.[[1]], .[[2]], .[[3]], lower.tail = FALSE)}),
           full_r_squared = map_dbl(full_mod, ~ summary(.)$r.squared),
           full_adj_r_squared = map_dbl(full_mod, ~ summary(.)$adj.r.squared),
           full_t_sign = case_when(full_t_p_value <= 0.05 ~ TRUE,
                                   TRUE ~ FALSE),
           full_f_sign = case_when(full_f_p_value <= 0.05 ~ TRUE,
                                   TRUE ~ FALSE),
           full_sign = case_when(full_t_sign == TRUE & full_f_sign == TRUE ~ TRUE,
                                 TRUE ~ FALSE),
           # START NEW
           forw_mod = map2(data, max, ~ stepAIC(lm(as.formula(min), data = .x), scope = .y, direction = "forward", trace = FALSE)),
           # END NEW
           forw_terms = map(forw_mod, ~ .$terms[[3]]) %>% as.character(),
           n_terms_forw = map_int(forw_terms, ~ str_split(., " ") %>% unlist() %>% str_subset("[^\\+]") %>% length()),
           indepoi_in_forw = case_when(str_detect(forw_terms, look) ~ TRUE,
                                       TRUE ~ FALSE),
           forw_t_p_value = map2_dbl(forw_mod, indepoi_in_forw,
                                     ~ summary(.x) %>% {ifelse(isTRUE(.y),
                                                               .$coefficients[indep_oi, 4],
                                                               NA_real_)}),
           forw_f_p_value = map2_dbl(forw_mod, forw_terms,
                                     ~ summary(.x) %>% {ifelse(.y == "1",
                                                               NA_real_,
                                                               pf(.$fstatistic[[1]],
                                                                  .$fstatistic[[2]],
                                                                  .$fstatistic[[3]],
                                                                  lower.tail = FALSE))}),
           forw_r_squared = map_dbl(forw_mod, ~ summary(.)$r.squared),
           forw_adj_r_squared = map_dbl(forw_mod, ~ summary(.)$adj.r.squared),
           forw_AIC = map_dbl(forw_mod, ~ extractAIC(.)[[2]]),
           forw_t_sign = case_when(forw_t_p_value <= 0.05 ~ TRUE,
                                   TRUE ~ FALSE),
           forw_f_sign = case_when(forw_f_p_value <= 0.05 ~ TRUE,
                                   TRUE ~ FALSE),
           forw_sign = case_when(forw_t_sign == TRUE & forw_f_sign == TRUE ~ TRUE,
                                 TRUE ~ FALSE),
           # START NEW
           step_mod = map2(data, max, ~ stepAIC(lm(as.formula(min), data = .x), scope = .y, direction = "both", trace = FALSE)),
           # END NEW
           step_terms = map(step_mod, ~ .$terms[[3]]) %>% as.character(),
           n_terms_step = map_int(step_terms, ~ str_split(., " ") %>% unlist() %>% str_subset("[^\\+]") %>% length()),
           indepoi_in_step = case_when(str_detect(step_terms, look) ~ TRUE,
                                       TRUE ~ FALSE),
           step_t_p_value = map2_dbl(step_mod, indepoi_in_step,
                                     ~ summary(.x) %>% {ifelse(isTRUE(.y),
                                                               .$coefficients[indep_oi, 4],
                                                               NA_real_)}),
           step_f_p_value = map2_dbl(step_mod, step_terms,
                                     ~ summary(.x) %>% {ifelse(.y == "1",
                                                               NA_real_,
                                                               pf(.$fstatistic[[1]],
                                                                  .$fstatistic[[2]],
                                                                  .$fstatistic[[3]],
                                                                  lower.tail = FALSE))}),
           step_r_squared = map_dbl(step_mod, ~ summary(.)$r.squared),
           step_adj_r_squared = map_dbl(step_mod, ~ summary(.)$adj.r.squared),
           step_AIC = map_dbl(step_mod, ~ extractAIC(.)[[2]]),
           step_t_sign = case_when(step_t_p_value <= 0.05 ~ TRUE,
                                   TRUE ~ FALSE),
           step_f_sign = case_when(step_f_p_value <= 0.05 ~ TRUE,
                                   TRUE ~ FALSE),
           step_sign = case_when(step_t_sign == TRUE & step_f_sign == TRUE ~ TRUE,
                                 TRUE ~ FALSE),
           # START NEW
           # back_mod = map2(resid_df_sup_0, data,
           #                 ~ ifelse(isTRUE(.x),
           #                          list(stepAIC(lm(formula = max, data = .y),
           #                                             direction = "backward",
           #                                             trace = FALSE)),
           #                          list(NULL))[[1]]),
           back_mod = pmap(list(resid_df_sup_0, max, data),
                           ~ if (isTRUE(..1)) {
                             stepAIC(lm(formula = ..2, data = ..3),
                                     direction = "backward",
                                     trace = FALSE)
                           } else {
                             NULL
                           }),
           # END NEW
           back_terms = map2(resid_df_sup_0, back_mod,
                             ~ ifelse(isTRUE(.x),
                                      list(.y$terms[[3]]),
                                      list(NA))[[1]]) %>% as.character() %>% na_if("NA"),
           n_terms_back = map_int(back_terms, ~ str_split(., " ") %>% unlist() %>% str_subset("[^\\+]") %>% length()),
           indepoi_in_back = case_when(str_detect(back_terms, look) ~ TRUE,
                                       TRUE ~ FALSE),
           back_t_p_value = map2_dbl(back_mod, indepoi_in_back,
                                     ~ summary(.x) %>% {ifelse(!is.null(.x) & isTRUE(.y),
                                                               .$coefficients[indep_oi, 4],
                                                               NA_real_)}),
           back_f_p_value = map2_dbl(back_mod, back_terms,
                                     ~ summary(.x) %>% {ifelse(.y %in% c("1", NA),
                                                               NA_real_,
                                                               pf(.$fstatistic[[1]],
                                                                  .$fstatistic[[2]],
                                                                  .$fstatistic[[3]],
                                                                  lower.tail = FALSE))}),
           back_r_squared = map2_dbl(back_mod, resid_df_sup_0,
                                     ~ summary(.x) %>% {ifelse(isTRUE(.y),
                                                               .$r.squared,
                                                               NA_real_)}),
           back_adj_r_squared = map2_dbl(back_mod, resid_df_sup_0,
                                         ~ summary(.x) %>% {ifelse(isTRUE(.y),
                                                                   .$adj.r.squared,
                                                                   NA_real_)}),
           back_AIC = map_dbl(back_mod, ~ ifelse(!is.null(.x),
                                                 extractAIC(.x)[[2]],
                                                 NA_real_)),
           back_t_sign = case_when(back_t_p_value <= 0.05 ~ TRUE,
                                   TRUE ~ FALSE),
           back_f_sign = case_when(back_f_p_value <= 0.05 ~ TRUE,
                                   TRUE ~ FALSE),
           back_sign = case_when(back_t_sign == TRUE & back_f_sign == TRUE ~ TRUE,
                                 TRUE ~ FALSE)) %>%
    select(dep_var, 1:data, everything()) %>%
    ungroup()

  if (isTRUE(hackfull)) {
    output <- output %>%
      select(dep_var:n_obs, matches("full")) %>%
      mutate(full_mod_coeff_indepoi = map2_dbl(full_mod, indep_var_oi,
                                               ~ .x$coefficients[[.y]]),
             full_mod_corr_indepoi = if_else(full_mod_coeff_indepoi < 0, "negative", "positive")) %>%
      select(dep_var, all_of(grps), model = max, n_obs, n_terms_full,
             full_mod_coeff_indepoi, full_mod_corr_indepoi,
             full_t_p_value:full_sign,
             data, full_mod)
  } else {
    output
  }
}
