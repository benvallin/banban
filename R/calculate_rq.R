#' Calculate RQ
#'
#' calculate_rq() computes quantification of target amplicons relative to reference amplicons.
#'
#' @param data a tibble provided by calculate_n0().
#' @param tar a character vector indicating the names of the target amplicons.
#' @param ref a character vector indicating the names of the reference amplicons.
#'
#' @return A tibble containing the RQ values for each sample. RQ values are computed as the ratio target mean_n0 / reference mean_n0 and stored in the rq_val variable.
#'
#' @export
#'
#' @importFrom dplyr filter
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom purrr pmap
#' @importFrom stringr str_c
#' @importFrom stringr str_detect
#' @importFrom stringr str_remove
#' @importFrom stringr str_replace
#' @importFrom tibble tibble
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#'
#' @examples
#'
calculate_rq <- function(data, tar, ref) {

  tar_ref <- union(tar, ref)

  out <- data %>%
    filter(amplicon %in% tar_ref) %>%
    nest(data = -sample) %>%
    mutate(rq = list(tibble(rq_nm = map(.x = tar, .f = ~ str_c(.x, ref, sep = "BANBAN")) %>% unlist %>% {.[!str_detect(., "(.*)BANBAN\\1$")]},
                            tar_nm = str_remove(rq_nm, "BANBAN.*$"),
                            ref_nm = str_remove(rq_nm, "^.*BANBAN")) %>%
                       mutate(rq_nm = str_replace(rq_nm, "BANBAN", "_"))) %>%
             {map2(.x = .,
                   .y = data,
                   ~ left_join(.x, .y[, c("amplicon", "mean_n0")] %>% rename(tar_nm = amplicon, tar_mean_n0 = mean_n0), by = "tar_nm")%>%
                     left_join(.y[, c("amplicon", "mean_n0")] %>% rename(ref_nm = amplicon, ref_mean_n0 = mean_n0), by = "ref_nm") %>%
                     mutate(rq_val = tar_mean_n0 / ref_mean_n0))}) %>%
    unnest(cols = c(rq)) %>%
    mutate(data = pmap(.l = list(data, tar_nm, ref_nm),
                       .f = ~ ..1 %>% filter(amplicon %in% c(..2, ..3)))) %>%
    select(sample, rq_nm, rq_val, tar_nm, tar_mean_n0, ref_nm, ref_mean_n0, data)
}
