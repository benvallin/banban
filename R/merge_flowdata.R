#' Merge flow cytometry data
#'
#' merge_flowdata() merges the flow cytometry data from different FCS files into a single flowFrame.
#'
#' @param filenames a character vector providing the path and name of the FCS files to merge into a flowFrame.
#' @param giv_indiv a logical indicating if flowFrames corresponding to each FCS file should also be provided.
#'
#' @return a flowFrame corresponding to the combined FCS files. If giv_indiv = TRUE, merge_flowdata() returns a list composed of the combined_ff and the indiv_ff.
#'
#' @export
#'
#' @importFrom flowCore read.FCS
#' @importFrom flowCore keyword
#' @importFrom flowCore exprs
#' @importFrom tidyr nest
#'
#'
#' @examples
#'
merge_flowdata <- function(filenames, give_indiv = FALSE) {

  list.ff <- lapply(X = filenames,
                    FUN = function(fcs.file) {read.FCS(filename = fcs.file)})

  names(list.ff) <- vapply(X = list.ff,
                           FUN = function(nm) { keyword(nm, "$FIL")[[1]] },
                           FUN.VALUE = vector(mode = "character", length = 1))

  matrix.values <- do.call(rbind,
                           lapply(X = list.ff, FUN = function(ff) { exprs(ff)}))

  exprs(list.ff[[1]]) <- matrix.values

  combined_ff <- list.ff[[1]]

  if (isFALSE(give_indiv)) {
    combined_ff
  } else {
    list(combined_ff = combined_ff,
         indiv_ff = list.ff)
  }
}
