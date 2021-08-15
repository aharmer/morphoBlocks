#' @title Read and combine multiple Landmark Editor files into a single data block
#'
#' @description Reads, collates and transforms landmark configurations from multiple specimens. \code{readPts} expects a directory containing two or more files with .pts extensions (i.e. landmark configurations exported from the Landmark Editor, Wiley et al. 2005).
#'
#' @param dirpath the directory path where two or more landmark configurations with  .pts extensions will be found.
#' @param landmarkRM a vector of landmarks to be excluded from the data block.
#' @param gpa a logical value indicating whether generalised Procrustes analyses should be performed. Default is \code{TRUE}.
#'
#' @details \code{readPts} reads landmark configurations from .pts files into the R environment, organises the configurations into a single array, and performs generalised Procrustes analysis on the array if required. Several objects are calculated for downstream analyses including centroid sizes for each landmark configuration.
#' \code{readPts} is a wrapper function for the \code{read.pts} and \code{cSize} functions from the \code{Morpho} package (Schlager 2017), and the \code{gpagen}, \code{two.d.array} and \code{arrayspecs} functions from the \code{geomorph} package (Adams and Otárola-Castillo 2013).
#' Landmarks are identified by their sequence within the configuration, which is an important consideration when excluding landmarks from the data block. For example, \code{landmarkRM = c(1,3)} would remove the first and third landmarks from all configurations once they were read into the R environment, and thus the data block would not include landmark 1 and landmark 3. The \code{landmarkRM} term might be used for analyses that want to test the sensitivity of dataset covariation on one or more landmarks.
#'
#' @return a 'block' class object, used for downstream analyses. The list contains the elements:
#' \itemize{
#'   \item \code{raw} collation of landmark configurations without Procrustes transformation
#'   \item \code{gpa.3D} landmark configurations after Procrustes transformation organised into a 3D array
#'   \item \code{gpa.2D} landmark configurations after Procrustes transformation organised into a 2D matrix
#'   \item \code{centroid} centroid sizes of landmark configurations without Procrustes transformation
#'   \item \code{p} number of landmarks that each configuration within the data block has
#'   \item \code{k} number of dimensions that each landmark within the data block has
#'   \item \code{n} number of landmark configurations included in the data block
#'   \item \code{curves} matrix for correctly calculating position of curve semilandmarks (see \code{geomorph::gpagen} for more detail)
#'   \item \code{surfaces} vector passed for defining positions of semilandmarks on surfaces (see \code{geomorph::gpagen} for more detail)
#'   }
#'
#' @examples
#' # Example 1
#' # For this example to work a directory (/...) containing .pts files must first be prepared.
#' dirpath = "/..."
#' block1 = readPts(dirpath)
#' block1@p
#' block1@k
#' block1@n
#'
#' # Example 2
#' # Exclude the first and third landmarks from the data block
#' block2 = readPts(dirpath, landmarkRM = c(1,3))
#' block2@p
#' block2@k
#' block2@n
#'
#' @references
#' Adams DC, Otárola-Castillo E. 2013. geomorph: an R package for the collection and analysis of geometric morphometric shape data. Methods in Ecology and Evolution 4:393–399 https://doi.org/10.1111/2041-210X.12035
#' Schlager S. 2017. Morpho and Rvcg–shape analysis in R. In Zheng G, Li S, Székely (eds.) Statistical shape and deformation analysis. Academic Press, London. Pp. 217–256.
#' Wiley DF, Amenta N, Alcantara DA, Ghosh D, Kil YJ, Delson E, Harcourt-Smith W, Rohlf FJ, St. John K, Hamann B. 2005. Evolutionary morphing. Proceedings of the IEEE Visualization 2005 (VIS’05), 431–438.
#'
#' @import geomorph
#' @import Morpho
#' @export
readPts = function (dirpath, landmarkRM = c(), gpa = T) {

  # require(geomorph, quietly = TRUE, warn.conflicts = FALSE)
  # require(Morpho, quietly = TRUE, warn.conflicts = FALSE)

  file.list = list.files(dirpath, full.names = T, pattern = ".pts")
  if (length(file.list) < 1) {
    stop("No .pts files found")
  }

  LM.read = as.matrix(read.pts(file.list[1]))
  LM.rm = c()
  for (i in 1: length(landmarkRM)) {
    LM.x = which(substr(rownames(LM.read), 0, 4) == landmarkRM[i])
    LM.rm = c(LM.rm, LM.x)

    LM.edit = LM.read
    if (length(landmarkRM) > 0) {
      LM.edit = LM.read[-c(LM.rm), ]
    }

    LM.fixed.ID = which(substring(unlist(dimnames(LM.edit)), 0, 1) == "S")
    LM.curve.ID = which(substring(unlist(dimnames(LM.edit)), 0, 1) == "C")
    LM.patch.ID = which(substring(unlist(dimnames(LM.edit)), 0, 1) == "P")

    LM.matrix.all = as.matrix(LM.edit)[1: nrow(LM.edit), ]
    for (i in 2: length(file.list)) {
      LM.read = as.matrix(read.pts(file.list[i]))
      LM.edit = LM.read
      if (length(landmarkRM) > 0) {
        LM.edit = LM.read[-c(LM.rm), ]
      }
      LM.matrix.all = rbind(LM.matrix.all, LM.edit[1: nrow(LM.edit), ])
    }
    LM.all = arrayspecs(unname(LM.matrix.all), nrow(LM.edit), 3)

    LM.curve.count = length(which(substring(rownames(LM.read), 0, 4) == "C000"))
    curves.map = c()
    for (i in ((1: (length(LM.curve.ID) / LM.curve.count)) * LM.curve.count) - (LM.curve.count - 1)) {
      m.curve = define.sliders(c(LM.curve.ID[i]: (LM.curve.ID[i] + (LM.curve.count - 1))))
      curves.map = rbind(curves.map, m.curve)
    }
    LM.fixed.ID = LM.fixed.ID
    LM.curve.ID = LM.curve.ID
    LM.patch.ID = LM.patch.ID

    if(gpa == T){
       proc.all = gpagen(LM.all, curves = curves.map, surfaces = LM.patch.ID, ProcD = F)
       gpa.2D = two.d.array(proc.all$coords)
       gpa.3D = arrayspecs(gpa.2D, p = dim(LM.all)[1], k = 3)
       } else {
       gpa.2D = two.d.array(LM.all)
       gpa.3D = arrayspecs(gpa.2D, p = dim(LM.all)[1], k = 3)
       }
  }
  centroid = c()
  for (i in 1: dim(LM.all)[3]) {
    centroid[i] = cSize(LM.all[, , i])
  }

  # block_out = setClass("block", slots = c(gpa.3D = "array", gpa.2D = "matrix", raw ="array",  centroid = "numeric", p = "numeric", k = "numeric", n = "numeric", curves = "array", surfaces = "integer"))
  output = block_out(gpa.3D = gpa.3D, gpa.2D = gpa.2D, raw = LM.all, centroid = centroid, p = dim(gpa.3D)[1], k = dim(gpa.3D)[2], n = dim(gpa.3D)[3], curves = curves.map, surfaces = LM.patch.ID)

  return (output)

}
