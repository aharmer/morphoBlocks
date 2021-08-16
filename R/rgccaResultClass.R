#' S4 rgccaResult class
#'
#' Used as output for multi-block rgcca analysis.
#'
#' @import methods
#' @export rgccaResult_out
#' @exportClass rgccaResult
setOldClass("rgcca")
rgccaResult_out = setClass("rgccaResult", slots = c(result = "rgcca", option = "character", block.list = "list", scores = "list", block.loadings = "list", p = "integer", k = "integer", n = "integer"))
