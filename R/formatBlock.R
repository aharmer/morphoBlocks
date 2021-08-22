#' @title Format configurations into a data block
#'
#' @description Organises data representing a block of landmark configurations (or pseudolandmark configurations) into a single data block. \code{formatBlock} is useful for formatting a data block of configurations from sources other than Landmark Editor (Wiley et al. 2005) or Generalised Procrustes Surface Analysis (Pomidor et al. 2016).
#'
#' @param block two-dimensional matrix or three-dimensional array of configurations to be formatted into a data block.
#' @param curves optional matrix passed to \code{gpagen} for correctly calculating position of curve semilandmarks (see \code{geomorph::gpagen} for more detail).
#' @param surfaces optional vector passed to \code{gpagen} for defining positions of semilandmarks on surfaces (see \code{geomorph::gpagen} for more detail).
#' @param cs optional vector of centroid sizes.
#' @param k number of dimensions that each point within the data block has. Required for configurations presented as a two-dimensional matrix. \code{k = 2} by default, but is updated when data are presented as a three-dimensional array.
#' @param gpa a logical value indicating whether generalised Procrustes analyses should be performed. Default is \code{TRUE}.
#'
#' @details \code{formatBlock} expects configurations to be presented as either a two-dimensional matrix or a three-dimensional array. If data are presented as a two-dimensional matrix then the number of \emph{k} dimensions for the landmarks (or pseudolandmarks) must be specified (default \code{k = 2}), each configuration should be a separate row in the matrix, and columns should contain corresponding points and dimensions in each configuration. If data are presented as a three-dimensional array then each of the \emph{n} configurations should be a separate matrix in the array (i.e. \code{array[p, k, n]}), and each matrix should have the same number of rows and columns where rows represent \emph{p} points and columns represent \emph{k} dimensions of each point.
#' By default \code{formatBlock} will perform generalised Procustes analysis (gpa) on the block using \code{gpagen} from the \code{geomorph} package (Adams and Otárola-Castillo 2013). Curve and surface information must be supplied if they are required. Centroid size will be calculated from the original configurations before Procrustes transformation. If gpa is not required then set \code{gpa = FALSE}. Centroid sizes of the original configurations should be supplied using the \code{cs} argument if \code{gpa = FALSE}, or else centroid sizes will be calculated from the (presumably Procrustes-transformed) configurations that are supplied.
#' \code{formatBlock} is a wrapper function for the \code{cSize} function from the \code{Morpho} package (Schlager 2017), and the \code{gpagen}, \code{two.d.array} and \code{arrayspecs} functions from the \code{geomorph} package (Adams and Otárola-Castillo 2013).
#'
#' @return a 'block' object, used for downstream analyses. The list contains the elements:
#' @return \item{raw}{collation of original configuration data}
#' @return \item{gpa.3D}{if \code{gpa = TRUE}, configurations after Procrustes transformation organised into a 3D array. If \code{gpa = FALSE}, original configuration data organised into a 3D array (may duplicate \code{raw})}
#' @return \item{gpa.2D}{if \code{gpa = TRUE}, configurations after Procrustes transformation organised into a 2D matrix. If \code{gpa = FALSE}, original configuration data organised into a 2D matrix (may duplicate \code{raw})}
#' @return \item{centroid}{centroid sizes either calculated from original configurations or supplied using the \code{cs} argument}
#' @return \item{p}{number of points in the configurations of each data block}
#' @return \item{k}{number of dimensions that the points in each configuration has}
#' @return \item{n}{number of configurations included in each data block}
#'
#' @references Adams DC, Otárola-Castillo E. 2013. geomorph: an \R package for the collection and analysis of geometric morphometric shape data. Methods in Ecology and Evolution 4:393–399 https://doi.org/10.1111/2041-210X.12035
#' @references Schlager S. 2017. Morpho and Rvcg–shape analysis in \R. In Zheng G, Li S, Székely (eds.) Statistical shape and deformation analysis. Academic Press, London. Pp. 217–256.
#' @references Pomidor BJ, Makedonska J, Slice DE. 2016. A landmark-free method for three-dimensional shape analysis. PLoS One 11: e0150368 https://doi.org/10.1371/journal.pone.0150368
#' @references Wiley DF, Amenta N, Alcantara DA, Ghosh D, Kil YJ, Delson E, Harcourt-Smith W, Rohlf FJ, St. John K, Hamann B. 2005. Evolutionary morphing. Proceedings of the IEEE Visualization 2005 (VIS’05), 431–438.
#'
#' @examples
#' # Format the head and tail data from the Plethodon dataset in the geomorph library
#' library(geomorph)
#' data(larvalMorph)
#' block1 <- formatBlock(block = larvalMorph$headcoords, curves = larvalMorph$head.sliders)
#' block2 <- formatBlock(block = larvalMorph$tailcoords, curves = larvalMorph$tail.sliders)
#'
#' @importFrom geomorph arrayspecs gpagen two.d.array
#' @importFrom Morpho cSize
#' @export
formatBlock <- function(block, curves = NULL, surfaces = NULL, cs = NULL, k = 2, gpa = TRUE) {
  # require(geomorph, quietly = TRUE, warn.conflicts = FALSE)
  # require(Morpho, quietly = TRUE, warn.conflicts = FALSE)

  if (gpa == TRUE) {
    if (length(dim(block)) == 3) {
      p <- dim(block)[1]
      k <- dim(block)[2]
      n <- dim(block)[3]
      proc.all <- gpagen(block, curves = curves, surfaces = surfaces, ProcD = FALSE, print.progress = FALSE)
      gpa.2D <- two.d.array(proc.all$coords)
      gpa.3D <- arrayspecs(gpa.2D, p, k)

      centroid <- c()
      for (i in 1:n) {
        centroid[i] <- cSize(block[, , i])
      }
    }
    if (length(dim(block)) == 2) {
      p <- dim(block)[2] / k
      k <- k
      n <- dim(block)[1]

      block <- arrayspecs(block, p, k)
      proc.all <- gpagen(block, curves = curves, surfaces = surfaces, ProcD = FALSE, print.progress = FALSE)
      gpa.2D <- two.d.array(proc.all$coords)
      gpa.3D <- arrayspecs(gpa.2D, p, k)

      centroid <- c()
      for (i in 1:n) {
        centroid[i] <- cSize(block[, , i])
      }
    }
  } else {
    if (length(dim(block)) == 3) {
      p <- dim(block)[1]
      k <- dim(block)[2]
      n <- dim(block)[3]
      gpa.3D <- block
      gpa.2D <- two.d.array(block)
      centroid <- cs
    }
    if (length(dim(block)) == 2) {
      p <- dim(block)[2] / k
      k <- k
      n <- dim(block)[1]

      block <- arrayspecs(block, p, k)
      gpa.3D <- block
      gpa.2D <- two.d.array(block)
      centroid <- cs
    }
  }

  output <- block_out(gpa.3D = gpa.3D, gpa.2D = gpa.2D, raw = block, centroid = centroid, p = dim(gpa.3D)[1], k = dim(gpa.3D)[2], n = dim(gpa.3D)[3])

  return(output)
}
