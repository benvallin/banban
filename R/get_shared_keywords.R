#' Get shared keyword values in a flowSet/GatingSet
#'
#' get_shared_keywords() summarizes the keyword values shared across elements of a flowSet or GatingSet. If requested, it can return the keyword values which are not shared instead.
#'
#' @param flow.obj a flowSet or GatingSet to analyse.
#' @param invert if TRUE, return the keyword values not shared.
#'
#' @return a data frame summarizing the keyword values shared (or not shared) between the elements of the flowSet/GatingSet.
#'
#' @export
#'
#' @importFrom flowCore fsApply
#' @importFrom flowCore keyword
#' @importFrom flowWorkspace lapply
#'
#' @examples
#'
get_shared_keywords <- function(flow.obj, invert = FALSE) {

  flow.obj.class <- class(flow.obj)

  keyword.names <- names(keyword(flow.obj[[1]]))

  if(flow.obj.class == "GatingSet") {
    keyword.values <- lapply(X = flow.obj, FUN = function(x) { keyword(x) })
  } else {
    keyword.values <- fsApply(x = flow.obj, FUN = function(x) { keyword(x) })
  }

  if(invert == FALSE) {

    out <- lapply(X = keyword.names,
                  FUN = function(x) {
                    unique(lapply(X = keyword.values, FUN = function(y) { y[[x]] }))
                  })
    names(out) <- keyword.names

    keyword.length <- vapply(X = out, FUN = function(x) { length(x) }, FUN.VALUE = vector(mode = "integer", length = 1))

    keyword.length <- ifelse(keyword.length == 1L, TRUE, FALSE)

    out <- data.frame(keyword_name = names(out[keyword.length]),
                      keyword_val = unlist(out[keyword.length]),
                      row.names = NULL)
  } else {

    out <- lapply(X = keyword.names,
                  FUN = function(x) {
                    as.character(lapply(X = keyword.values, FUN = function(y) { y[[x]] }))

                  })

    names(out) <- keyword.names

    keyword.length <- vapply(X = out, FUN = function(x) { length(unique(x)) }, FUN.VALUE = vector(mode = "integer", length = 1))

    keyword.length <- ifelse(keyword.length == 1L, FALSE, TRUE)

    out <- data.frame(keyword_name = names(out[keyword.length]),
                      keyword_val = cbind(out[keyword.length]),
                      row.names = NULL)

  }
}
