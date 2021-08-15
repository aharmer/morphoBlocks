#' @import methods
blockResult_out = setClass("blockResult", slots = c(result = class(result), option = "character", block.list = "list", scores = "list", block.loadings = "list", p = "integer", k = "integer", n = "integer"))
