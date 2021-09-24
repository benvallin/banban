#' Gather identical amplicons from distinct CFX output files
#'
#' @description gather_amplicons() aims at creating virtual CFX output as if the amplification of a given cDNA, performed on distinct qPCR plates, had been performed on a single plate.
#'
#' It creates a single set of virtual CFX output for each amplicon spread across distinct original CFX output files. Virtual infos_plate and fluo_data variables are created for each amplicon and can then be exported as virtual CFX output files to feed LinRegPCR software.
#'
#' @param path a string specifying the path to the directory containing LinRegPCR output directories.
#' @param datafiles a character vector containing the names of LinRegPCR output directories. Each LinRegPCR output directory must be named after its corresponding CFX data file and must contain the "Quantification Cq Results" and "Quantification Amplification Results" files provided by CFX software. These files must be exported as csv from CFX software and their names contain, as a prefix, the name of the corresponding CFX data file.
#' @param max_out a logical, default is FALSE. If TRUE, all amplicons are gathered on a single set of virtual CFX output. So, LinRegPCR software can be fed with a single pair of virtual infos_plate and fluo_data files.
#'
#' @return A tibble containing virtual infos_plate and fluo_data variables, either for each amplicon (max_out = FALSE) or for all amplicons gathered together (max_out = TRUE).
#'
#' @export
#'
#' @importFrom dplyr filter
#' @importFrom dplyr if_else
#' @importFrom dplyr mutate
#' @importFrom dplyr rename_all
#' @importFrom dplyr rename_at
#' @importFrom dplyr select
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom purrr pmap
#' @importFrom purrr reduce
#' @importFrom readr read_csv
#' @importFrom stringr str_c
#' @importFrom stringr str_detect
#' @importFrom stringr str_remove
#' @importFrom tibble tibble
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#' @importFrom tidyselect everything
#'
#' @examples
#'
gather_amplicons <- function(path, datafiles, max_out = FALSE) {

  pool <- tibble(datafiles = datafiles) %>%
    mutate(infos_plate = map(.x = datafiles,
                             .f = function(dirname) {
                               read_csv(str_c(path, "/", dirname, "/", dirname, " -  Quantification Cq Results_0.csv")) %>%
                                 rename_all(.funs = tolower) %>%
                                 select(well, target, sample) %>%
                                 mutate(well = if_else(str_detect(well, "^[A-H]0.*$"), str_remove(well, "0"), well))
                             }),
           fluo_plate = map(.x = datafiles,
                            .f = function (dirname) {
                              read_csv(str_c(path, "/", dirname, "/", dirname, " -  Quantification Amplification Results_SYBR.csv"))
                            }))

  out <- pool["infos_plate"] %>%
    unnest(cols = c(infos_plate)) %>%
    select(target) %>%
    unique %>%
    mutate(raw_data = map(target,
                          ~ pool %>%
                            mutate(target = .x,
                                   infos_plate = map2(infos_plate,
                                                      target,
                                                      ~ filter(.x, target == .y)),
                                   fluo_plate = pmap(list(fluo_plate,
                                                          infos_plate,
                                                          datafiles),
                                                     ~ select(..1, X1 = `...1`, Cycle, ..2$well) %>%
                                                       rename_at(.vars = -c(1, 2), .funs = function(., ..3) str_c(., ..3, sep = "_"))),
                                   infos_plate = map2(infos_plate, datafiles,
                                                      ~ mutate(.x, well = str_c(well, .y, sep = "_")))) %>%
                            select(target, everything())),
           infos_plate = map(raw_data,
                             ~ select(.x, infos_plate) %>% unnest(cols = c(infos_plate))),
           fluo_data = map(raw_data,
                           ~ select(.x, fluo_plate) %>% {.[[1]]} %>% reduce(full_join)))

  if (max_out == FALSE) {
    out
  } else {
    max_out <- out %>%
      select(infos_plate, fluo_data) %>%
      nest(infos_plate = c(infos_plate), fluo_data = c(fluo_data)) %>%
      mutate(infos_plate = map(infos_plate, ~ unnest(.x, cols = c(infos_plate))),
             fluo_data = map(fluo_data, ~ .x[[1]] %>% reduce(full_join)))
  }
}
