#' S4 block class
#'
#' Used as input for multi-block analysis.
#'
#' @name block-class
#' @rdname block-class
#' @import methods
#' @export
block_out = setClass("block", slots = c(gpa.3D = "array", gpa.2D = "matrix", raw ="array",  centroid = "numeric", p = "numeric", k = "numeric", n = "numeric", curves = "array", surfaces = "integer"))
