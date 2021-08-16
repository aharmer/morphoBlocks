#' S4 prcompResult class
#'
#' Used as output for multi-block prcomp analysis.
#'
#' @import methods
#' @export prcompResult_out
#' @exportClass prcompResult
setOldClass("prcomp")
prcompResult_out = setClass("prcompResult", slots = c(result = "prcomp", option = "character", block.list = "list", scores = "list", block.loadings = "list", p = "integer", k = "integer", n = "integer"))
