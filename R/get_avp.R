#' Add added-variable plots
#'
#' get_avp() draws the added-variable plots for the "best models" determined by choose_on_aic().
#'
#' @param tibble a tibble produced by choose_on_aic() or get_sign_mod().
#'
#' @return A tibble in which the added-variable plots for each best model are stored in a variable called avp.
#'
#' @export
#'
#' @importFrom car avPlots
#' @importFrom dplyr mutate
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggpubr theme_pubr
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom stats formula
#' @importFrom stats lm
#' @importFrom stats setNames
#' @importFrom stringr str_c
#' @importFrom stringr str_replace
#' @importFrom tibble as_tibble
#'
#' @examples
#'
get_avp <- function(tibble) {

  output <- tibble %>%
    mutate(avp = map2(best_mod_val, rq_nm,
                      ~ avPlots(lm(formula(paste(.x$terms %>% deparse %>% str_replace("rq_val", .y), collapse = " ")),
                                   data = .x$model %>% setNames(str_replace(colnames(.), "rq_val", .y))), ask = FALSE) %>%
                        map(~ as_tibble(.x) %>% ggplot(aes(.[[1]], .[[2]])) +
                              geom_smooth(method = "lm", se = FALSE, color = "orange") +
                              geom_point() +
                              theme_pubr() +
                              labs(x = str_c(colnames(.)[[1]], " | others"),
                                   y = str_c(colnames(.)[[2]], " | others")))))
}
