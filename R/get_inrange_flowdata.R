#' Get in range flow cytometry data
#'
#' get_inrange_flowdata() "tidies" each flowFrame of a flowSet by removing all events for which the recorded value of one or more parameter is out of range.
#'
#' @param fs a flowSet to tidy.
#'
#' @return a tidied flowSet.
#'
#' @export
#'
#' @importFrom flowCore fsApply
#' @importFrom flowCore keyword
#' @importFrom flowCore rectangleGate
#' @importFrom flowCore Subset
#'
#' @examples
#'
get_inrange_flowdata <- function(fs) {

  message("\nLet's go!\n")

  output <- fsApply(x = fs,
                    FUN = function(ff) {

                      message(paste0("Processing file ", keyword(ff, "$FIL")[[1]], "...\n"))

                      param.table <- data.frame(par_number = as.character(seq(as.integer(keyword(ff, c("$PAR"))[[1]]))))
                      param.table$par_name <- vapply(X = param.table$par_number,
                                                     FUN = function(x) { keyword(ff, paste0("$P", x, "N"))[[1]] },
                                                     FUN.VALUE = vector(mode = "character", length = 1))
                      param.table$par_maximum <- vapply(X = param.table$par_number,
                                                        FUN = function(x) { as.double(keyword(ff, paste0("$P", x, "R"))[[1]]) },
                                                        FUN.VALUE = vector(mode = "double", length = 1))

                      count <- 1L
                      while(count < length(param.table[[1]])) {
                        ff <- Subset(x = ff,
                                     subset = !rectangleGate(.gate = setNames(object = list(c(param.table$par_maximum[[count]]-1,
                                                                                              Inf)),
                                                                              nm = param.table$par_name[[count]])))
                        count <- count + 1L
                      }
                      ff
                    },
                    simplify = T)

  message("All done!\n")

  output
}
