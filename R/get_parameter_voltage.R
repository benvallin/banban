#' Get in range flow cytometry data
#'
#' get_parameter_voltage() extracts the voltage value of parameters in each flowFrame of a flowSet. If requested, it also checks that all voltage values of a parameter are identical across the flowFrames.
#'
#' @param fs a flowSet to tidy.
#' @param check_identical a logical indicating if the equality of voltage values should be checked.
#'
#' @return a data frame summarizing the voltage values of parameters included in the flowSet. If check_identical = TRUE, get_parameter_voltage() returns a list composed of the voltage_table and the voltage_check.
#'
#' @export
#'
#' @importFrom flowCore fsApply
#' @importFrom flowCore keyword
#' @importFrom tidyr nest
#'
#'
#' @examples
#'
get_parameter_voltage <- function(fs, check_identical = FALSE) {

  output <- fsApply(x = fs,
                    FUN = function(ff) {
                      param.table <- data.frame(par_number = as.character(seq(as.integer(keyword(ff, c("$PAR"))[[1]]))))
                      param.table$par_name <- vapply(X = param.table$par_number,
                                                     FUN = function(x) { keyword(ff, paste0("$P", x, "N"))[[1]] },
                                                     FUN.VALUE = vector(mode = "character", length = 1))
                      param.table$par_voltage <- vapply(X = param.table$par_number,
                                                        FUN = function(x) { if (is.null(keyword(ff, paste0("$P", x, "V"))[[1]])) {
                                                          NA_real_
                                                        } else {
                                                          as.double(keyword(ff, paste0("$P", x, "V"))[[1]])
                                                        }},
                                                        FUN.VALUE = vector(mode = "double", length = 1))
                      param.table$sample_id <- keyword(ff, "$FIL")[[1]]
                      param.table[!is.na(param.table$par_voltage), c("sample_id", "par_number", "par_name", "par_voltage")]
                    })

  output <- do.call(rbind, c(output, make.row.names = F))

  if (isFALSE(check_identical)) {
    output
  } else {
    voltage_check <- nest(.data = output, data = -par_name)
    voltage_check$identical_voltage <- vapply(X = voltage_check$data,
                                              FUN = function(x) {ifelse(length(unique(x$par_voltage)) == 1L, T, F)},
                                              FUN.VALUE = vector(mode = "logical", length = 1))
    output <- list(voltage_table = output,
                   voltage_check = voltage_check)
  }
}
