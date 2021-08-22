#' S4 blockList class
#'
#' Used as input for multiblock analysis.
#'
#' @name blockList-class
#' @rdname blockList-class
#' @importFrom methods setClass new
#' @export
blockList_out <- setClass("blockList", slots = c(block.list = "list", p = "integer", k = "integer", n = "integer"))
