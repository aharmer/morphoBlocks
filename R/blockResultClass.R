#' S4 blockResult class
#'
#' Used as output for multiblock analysis.
#'
#' @name blockResult-class
#' @rdname blockResult-class
#' @import methods
#' @export
blockResult_out = setClass("blockResult", slots = c(result = class("result"), option = "character", block.list = "list", scores = "list", block.loadings = "list", p = "integer", k = "integer", n = "integer"))
