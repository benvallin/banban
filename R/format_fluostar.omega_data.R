#' Format fluorescence data from FLUOstar Omega plate reader
#'
#' blabla
#'
#' @param data_file blabla
#' @param infos_samples_file blabla
#' @param extend blabla
#'
#' @return blabla
#'
#' @export
#'
#' @importFrom dplyr arrange
#' @importFrom dplyr full_join
#' @importFrom dplyr group_by
#' @importFrom dplyr if_else
#' @importFrom dplyr mutate
#' @importFrom dplyr rename_all
#' @importFrom dplyr select
#' @importFrom dplyr ungroup
#' @importFrom forcats as_factor
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom purrr map_dbl
#' @importFrom readr read_csv
#' @importFrom stats sd
#' @importFrom stringr str_c
#' @importFrom stringr str_remove
#' @importFrom stringr str_split
#' @importFrom tidyr nest
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr unnest
#' @importFrom tidyselect matches
#'
#' @examples
#'
format_fluostar.omega_data <- function(data_file, infos_samples_file, extend = FALSE) {

  output <- read_csv(data_file, col_names = TRUE, skip = 4) %>%
    rename_all(.funs = function(x) str_remove(x, "Raw.*- ")) %>%
    mutate(well_id = str_c(`Well Row`, `Well Col`)) %>%
    pivot_longer(cols = matches("\\d*\\sh.*"), names_to = "hour", values_to = "fluorescence") %>%
    mutate(hour =
             str_remove(hour, "\\smin") %>%
             str_split("\\sh") %>%
             map(.f = function(x) as.integer(x) %>% {if_else(is.na(.), 0L, .)}) %>%
             map_dbl(.f = function(x) (x[[1]]*60 + x[[2]])/60)) %>%
    select(well_id, hour, fluorescence) %>%
    full_join(read_csv(infos_samples_file) %>% mutate(sample_id = as_factor(sample_id))) %>%
    nest(replicate_data = -c(sample_id, hour)) %>%
    mutate(mean_fluorescence = map_dbl(replicate_data, ~ mean(.x$fluorescence)),
           sd_fluorescence = map_dbl(replicate_data, ~ sd(.x$fluorescence))) %>%
    select(sample_id, hour, mean_fluorescence, sd_fluorescence, replicate_data) %>%
    arrange(hour)
  if(extend == TRUE) {
    output <- output %>%
      unnest(cols = replicate_data) %>%
      group_by(sample_id, hour) %>%
      mutate(well_id = seq_along(well_id) %>% {str_c("fluo_rep", .)}) %>%
      pivot_wider(names_from = "well_id", values_from = "fluorescence") %>%
      ungroup()
  } else {
    output
  }
}
