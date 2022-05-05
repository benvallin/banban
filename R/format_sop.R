#' Format qPCR results from Step One Plus software
#'
#'format_sop() tidies the "Results" sheet provided by StepOne Software v2.3, and makes it ready for processing by ddct.compute().
#'
#' @param file a string specifying the path to the .csv file containing qPCR results.
#'
#' @return A dataframe summarising results in 3 columns: sample.id, tar.nm, ct.
#'
#' @export
#'
#' @examples
#'
format_sop <- function(file) {

  # Read and clean raw data

  data <- read.csv(file = file, skip = 7L)
  data <- data[c(1L:(nrow(data)-5L)), c(2, 3, 10)]
  names(data) <- c("sample.id", "tar.nm", "ct")

  data$ct <- vapply(X = data$ct,
                    FUN = function(x) { ifelse(x == "Undetermined", NA_real_, as.double(x)) },
                    FUN.VALUE = vector(mode = "double", length = 1L))

  return(data)

}
