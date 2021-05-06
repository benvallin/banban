#' Calculate mean n0 from LinRegPCR analysis
#'
#'calculate_n0() computes the mean n0 values of qPCR replicates from xls/xlsx/csv files containing LinRegPCR output.
#'
#' @param path a string specifying the path to the directory containing LinRegPCR output files.
#' @param linregpcr_files a character vector containing the names of LinRegPCR output files, with the extension included (.xls, .xlsx or .csv).
#' @param type a string specifying the extension of LinRegPCR output files (either .xls, .xlsx or .csv).
#' @param include_all a logical, default is FALSE. If TRUE, will also include the samples that did not reach the plateau phase or displayed deviating PCR efficiency.
#' @param hkg optional. A named list of character vectors specifying the sets of amplicons to be used for computing reference mean n0 values. For each vector of the list, reference mean n0 values are named after the vector and computed as the geometric mean of the vector elements.
#'
#' @return A tibble containing the mean n0 value of each amplicon (including newly created references) for each RT sample.
#'
#' @export
#'
#' @importFrom EnvStats geoMean
#' @importFrom dplyr arrange
#' @importFrom dplyr filter
#' @importFrom dplyr full_join
#' @importFrom dplyr group_by
#' @importFrom dplyr if_else
#' @importFrom dplyr mutate
#' @importFrom dplyr n_distinct
#' @importFrom dplyr rename_all
#' @importFrom dplyr select
#' @importFrom dplyr ungroup
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom purrr map_dbl
#' @importFrom purrr map_int
#' @importFrom purrr map_lgl
#' @importFrom purrr pmap
#' @importFrom readr read_csv
#' @importFrom readxl read_excel
#' @importFrom rlang enquo
#' @importFrom rlang !!
#' @importFrom rlang :=
#' @importFrom stats sd
#' @importFrom stringr str_c
#' @importFrom stringr str_detect
#' @importFrom stringr str_remove
#' @importFrom tibble tibble
#' @importFrom tidyr nest
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr unnest
#' @importFrom tidyselect everything
#'
#' @examples
#'
calculate_n0 <- function(path,
                         linregpcr_files,
                         type,
                         include_all = FALSE,
                         hkg = NULL) {

  linregpcr_files <- enquo(linregpcr_files)
  path <- enquo(path)
  include_all <- enquo(include_all)

  if (type == "csv") {
      reader <- read_csv
    } else if (type == "excel") {
      reader <- read_excel
    }

  infos <- tibble(file_nm = !!linregpcr_files) %>%
    mutate(data = map(.x = file_nm,
                      .f = function (var) reader(str_c(!!path, var), col_names = FALSE)),
           infos = map(.x = data,
                       .f = function (var) tibble(software_version = var[[1, 2]] %>% str_remove("version:"),
                                                  input_sheet = var[[3, 1]] %>% str_remove("Input Sheet: "),
                                                  analysis_date = var[[2, 1]] %>% str_remove("analysis date:"),
                                                  wol_setting = var[[1, 3]] %>% str_remove("WoL: "),
                                                  points_in_wol = var[[2, 3]] %>% str_remove("points in WoL: "),
                                                  threshold_setting = var[[3, 3]] %>% str_remove("Threshold: "),
                                                  chemistry = var[[1, 10]] %>% str_remove("Chemistry: "),
                                                  input_dna = var[[2, 10]] %>%  str_remove("Input: "),
                                                  mean_eff_criteria = var[c(26, 27), 18][[1]] %>% str_c(collapse = " AND "),
                                                  mean_eff_deviation_criterion = var[[28, 18]],
                                                  log_linear_phase_criterion = var[[30, 18]] %>% str_remove("log-linear phase criterion: ")))) %>%
    select(-data) %>%
    unnest(cols = c(infos))

  out <- tibble(file_nm = !!linregpcr_files) %>%
    mutate(data = map(.x = file_nm,
                      .f = function (var) reader(str_c(!!path, var), skip = 3) %>%
                        select(-c(1, length(.))) %>%
                        filter(!is.na(name)) %>%
                        rename_all(.funs = tolower) %>%
                        mutate(n0_mean_eff =
                                 if (isFALSE(!!include_all)) {
                                   if_else(str_detect(sample_use, "3"), `n0_(mean eff)`, NA_real_)
                                 } else {
                                   if_else(str_detect(quality_checks, "[12]"), NA_real_, threshold / mean_pcr_eff ^ cq)
                                 },
                               use_for_wol_setting = if_else(str_detect(sample_use, "1"), TRUE, FALSE),
                               contributes_to_mean_pcr_efficiency = if_else(str_detect(sample_use, "2"), TRUE, FALSE),
                               n0_value_calculated = if_else(str_detect(sample_use, "3"), TRUE, FALSE)) %>%
                        {if (isTRUE(!!include_all)) {
                          mutate(., n0_value_REcalculated = if_else(!is.na(n0_mean_eff), TRUE, FALSE))
                        } else {
                          .
                        }
                        } %>%
                        mutate(passed_all_checks = if_else(str_detect(quality_checks, "0"), TRUE, FALSE),
                               no_amplification = if_else(str_detect(quality_checks, "1"), TRUE, FALSE),
                               baseline_error = if_else(str_detect(quality_checks, "2"), TRUE, FALSE),
                               no_plateau = if_else(str_detect(quality_checks, "3"), TRUE, FALSE),
                               noisy_sample = if_else(str_detect(quality_checks, "4"), TRUE, FALSE),
                               pcr_efficiency_outside_5_percent = if_else(str_detect(quality_checks, "5"), TRUE, FALSE),
                               excluded_from_mean_eff = if_else(str_detect(quality_checks, "6"), TRUE, FALSE),
                               excluded_by_user = if_else(str_detect(quality_checks, "7"), TRUE, FALSE),
                               included_by_user = if_else(str_detect(quality_checks, "8"), TRUE, FALSE),
                               manual_baseline = if_else(str_detect(quality_checks, "9"), TRUE, FALSE)))) %>%
    unnest(cols = c(data)) %>%
    full_join(infos, by = "file_nm") %>%
    nest(linregpcr_data = -c(sample, amplicon)) %>%
    mutate(replicate_n = map_int(.x = linregpcr_data, .f = function (var) var[!is.na(var$n0_mean_eff), "n0_mean_eff"] %>% {length(.[[1]])}),
           mean_n0 = map_dbl(.x = linregpcr_data, .f = function (var) mean(var$n0_mean_eff, na.rm = TRUE)) %>% {if_else(is.nan(.), NA_real_, .)},
           sd_n0 = map_dbl(.x = linregpcr_data, .f = function (var) sd(var$n0_mean_eff)),
           replicates_on_same_plate = map_lgl(.x = linregpcr_data, .f = function (var) if_else(length(var$file_nm) > 1, n_distinct(var$file_nm) %>% {if_else(. == 1, TRUE, FALSE)}, NA))) %>%
    select(sample, everything()) %>%
    arrange(sample, amplicon)

  if (is.null(hkg)) {

    out

  } else {

    hkg <- enquo(hkg)
    infos_hkg <- tibble(nm = names(!!hkg),
                        val = map(.x = !!hkg,
                                  .f = function (var) str2expression(str_c("c(", str_c(var, collapse = ", "), ")"))) %>% unlist())

    out_hkg <- out %>%
      select(sample, amplicon, mean_n0) %>%
      pivot_wider(names_from = "amplicon", values_from = "mean_n0")

    for (i in seq_along(infos_hkg$nm)) {
      out_hkg <- out_hkg %>%
        group_by(sample) %>%
        mutate(!!str_c("geomean.", infos_hkg$nm[[i]]) := geoMean(!!infos_hkg$val[[i]])) %>%
        ungroup
    }

    out_hkg %>% pivot_longer(cols = -sample, names_to = "amplicon", values_to = "mean_n0") %>%
      full_join(out, by = c("sample", "amplicon", "mean_n0")) %>%
      select(colnames(out)) %>%
      mutate(linregpcr_data = pmap(.l = list(linregpcr_data,
                                             sample,
                                             amplicon),
                                   .f = ~ if (!is.null(..1)) {
                                     ..1
                                   } else {
                                     filter(out,
                                            sample == ..2,
                                            amplicon %in% (!!hkg %>% {.[[str_remove(..3, "geomean.")]]})) %>%
                                       {.[, c("amplicon", "linregpcr_data")] %>%
                                           unnest(cols = c(linregpcr_data))}
                                   }))
  }
}
