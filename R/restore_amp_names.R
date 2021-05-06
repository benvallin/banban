#' Restore amplicons names in LinRegPCR output
#'
#' @description restore_amp_names() is useful when LinRegPCR output files contain a single amplicon, in which case the amplicon name is replaced by "not_named".
#'
#' The function takes xlsx files produced by LinRegPCR as inputs and returns csv files in which the amplicon names have been restored based on files names.
#'
#' @param path a string specifying the path to the directory containing LinRegPCR output files.
#' @param filenames a character vector containing the names of xlsx files produced by LinRegPCR, with the xlsx extension included. The files must be of xlsx type and their name must match the following regex pattern: "^.*_amplicon.name.xlsx$".
#'
#' @return For each xlsx file provided as input, a csv file in which the amplicon name has been restored. The csv files are written in the directory indicated by the path argument.
#'
#' @export
#'
#' @importFrom dplyr if_else
#' @importFrom dplyr mutate
#' @importFrom purrr walk
#' @importFrom readr write_csv
#' @importFrom readxl read_excel
#' @importFrom stringr str_c
#' @importFrom stringr str_remove
#'
#' @examples
#'
restore_amp_names <- function(path, filenames) {

  walk(.x = filenames,
       .f = ~ read_excel(str_c(path, .x)) %>%
         mutate(`...12` = if_else(`...12` == "not_named",
                                  .x %>% str_remove("^.*_") %>% str_remove(".xlsx"),
                                  `...12`),
                `...14` = if_else(`...14` == "sample",
                                  "Sample",
                                  `...14`)) %>%
         write_csv(str_c(path, .x %>% str_remove(".xlsx"), ".csv")))

}
