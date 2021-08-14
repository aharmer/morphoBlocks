#' @title Read pseudolandmarks into a data block
#'
#' @description Reads and organises pseudolandmarks from Generalised Procrustes Surface Analysis (GPSA; Pomidor et al. 2015) into a data block. \code{readGPSA} expects a directory containing a single .dat file containing homologised points for one block of pseudolandmark configurations.
#'
#' @param dirpath the directory path where a file containing homologised pseudolandmark points will be found. The homologised points file will have a .dat extension and must be produced using the Generalised Procrustes Surface Analysis software tool (Pomidor et al. 2015).
#'
#' @details \code{readGPSA} reads a matrix of pseudolandmark configurations from a .dat file into the R environment and organises the configurations into an array. Several objects are calculated for downstream analyses including centroid sizes for each pseudolandmark configuration. Note that centroid sizes are not calculated for the original meshes processed using GPSA, but are instead calculated for the homologised pseudolandmark points from each of those meshes after processing with GPSA.
#' \code{readGPSA} is a wrapper function for the \code{fread} function from the \code{data.table} library (Dowle and Srinivasan 2020), the \code{cSize} function from the \code{Morpho} package (Schlager 2017), and the \code{arrayspecs} function from the \code{geomorph} package (Adams and Otárola-Castillo 2013).
#'
#' @return a 'block' class object, used for downstream analyses. The list contains the elements:
#' \itemize{
#'   \item \code{raw} collation of pseudolandmark configurations organised as a 2D matrix
#'   \item \code{gpa.3D} pseudolandmark configurations organised into a 3D array
#'   \item \code{gpa.2D} duplication of \code{raw}. The duplication occurs so that the object structure is consistent with the structure produced by closely-related functions (e.g. \code{readPts})
#'   \item \code{centroid} centroid sizes of pseudolandmark configurations
#'   \item \code{p} number of pseudolandmarks that each configuration within the data block has
#'   \item \code{k} number of dimensions that each pseudolandmark within the data block has
#'   \item \code{n} number of pseudolandmark configurations included in the data block
#'   }
#'
#' @examples
#' # For this example to work a directory (/...) containing .dat file must first be prepared.
#' dirpath = "/..."
#' block1 = readGPSA(dirpath)
#' block1@p
#' block1@k
#' block1@n
#'
#' @import geomorph
#' @import Morpho
#' @import data.table
#' @export
readGPSA = function(dirpath) {

  require(geomorph, quietly = TRUE, warn.conflicts = FALSE)
  require(Morpho, quietly = TRUE, warn.conflicts = FALSE)
  require(data.table, quietly = TRUE, warn.conflicts = FALSE)

  readpath = list.files(dirpath, pattern = "gpsa_homologized_points.dat", full.names = T)
  if (length(readpath) < 1) {
    stop("No homologized points file found. Remember to specify the path of just the directory containing the file, and not the path of the file itself")
  }
  step1 = fread(readpath)
  gpa.2D = matrix(as.numeric(unlist(step1)), nrow = dim(step1)[1])
  gpa.3D = arrayspecs(gpa.2D, p = dim(step1)[2] / 3, k = 3)

  centroid = c()
  for(i in 1: dim(gpa.2D)[1]) {
    centroid[i] = cSize(gpa.2D[i, ])
  }

  block_out = setClass("block", slots = c(gpa.3D = "array", gpa.2D = "matrix", raw = "array", centroid = "numeric", p = "numeric", k = "numeric", n = "numeric"))
  output = block_out(gpa.3D = gpa.3D, gpa.2D = gpa.2D, raw = gpa.2D, centroid = centroid, p = dim(gpa.3D)[1], k = dim(gpa.3D)[2], n = dim(gpa.3D)[3])

  return(output)

}
