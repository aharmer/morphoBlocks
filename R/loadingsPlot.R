#' @title Plot component loadings from an analysis of multiple data blocks
#'
#' @description Two- or three-dimensional scatterplots showing the loadings from an analysis of data blocks containing shape information
#'
#' @param result result produced by the analyseBlocks function.
#' @param comp the component selected to be shown in the loadings plots. Default is component 1. The selected component must be within the range of components calculated by \code{analyseBlocks}.
#' @param cex.2d value for specifying point size when plotting loadings from analysis of two-dimensional shape configurations. Passed to cex term of \code{plot}. Default is 2.
#' @param cex.3d value for specifying point size when plotting loadings from analysis of three-dimensional shape configurations. Passed to cex term of \code{plot3d}. Default is 5.
#'
#' @details \code{loadingsPlot} helps to visualise the result from the \code{analyseBlocks} function by using the component loadings to colour the mean position of each point in the consensus space (for \code{option = "rcpca"} in \code{analyseBlocks}) or the concatenated superblock (for \code{option = "pca"} in \code{analyseBlocks}). Points are coloured along a gradient from orange to blue. Points that are more orange have relatively large loadings values (i.e. larger within-block covariation), and points that are more blue have relatively small loadings values (i.e. smaller within-block covariation).
#' \code{loadingsPlot} will present a two-dimensional scatterplot if the list of data blocks that were analysed had \code{k = 2} dimensions. Likewise, \code{loadingsPlot} will present a three-dimensional scatterplot if the list of data blocks had \code{k = 3} dimensions.
#' Data from the individual blocks that contributed to either the consensus space or concatenated superblock are shown as separate scatterplots for visual clarity. Separate plots are produced because the points in each block have the same origin after Procrustes transformation and would plot on top of one another if they were presented along the same set of Cartesian axes.
#' \code{loadingsPlot} uses the \code{arrayspecs} and \code{mshape} functions from the \code{geomorph} library (Adams and Otárola-Castillo 2013).
#'
#' @return a two- or three-dimensional scatterplot
#'
#' @examples
#' block1 = dodecBlock()
#' block2 = dodecBlock()
#' blocklist = combineBlocks(blocks = c(block1, block2))
#' result1 = analyseBlocks(blocklist)
#' result2 = analyseBlocks(blocklist, option = "pca", ncomp = 10)
#' loadingsPlot(result1)
#' loadingsPlot(result2, comp = 2)
#'
#' @references
#' \itemize{
#'   \item Adams DC, Otárola-Castillo E. 2013. geomorph: an \R package for the collection and analysis of geometric morphometric shape data. Methods in Ecology and Evolution 4:393–399 https://doi.org/10.1111/2041-210X.12035
#'   }
#'
#' @import geomorph
#' @export
loadingsPlot = function(result, comp = 1, cex.2d = 2, cex.3d = 5) {

  require(geomorph, quietly = TRUE, warn.conflicts = FALSE)

  if(class(result) != "blockResult") {
    stop("Object is not the expected format. Use analyseBlock function to first analyse the data")
  }

  n = result@n[1]
  k = result@k[1]
  block.list = result@block.list
  block.loadings = result@block.loadings
  nvert = c(result@p * k, (sum(result@p * k)))
  J = (length(block.list)) - 1

  if(result@option == "pca") {
    superblock = result@block.list[length(result@block.list)]

    if(k == 3) {
      cl = round(sqrt(J))
      rw = ceiling(J / (round(sqrt(J))))
      mfrow3d(rw, cl)
      pos.st = cumsum(c(1, nvert))
      pos.en = cumsum(nvert)

      for(i in 1: J) {
        A = arrayspecs(result@block.list[[i]], p = result@p[i], k = result@k[1])
        ref = mshape(A)

        block.loadings.bone = matrix(block.loadings[pos.st[i]: pos.en[i], ], ncol = ncol(block.loadings), byrow = F)
        block.loadings.PC = matrix(block.loadings.bone[, comp], ncol = k, byrow = T)
        block.loadings.len = sqrt(((block.loadings.PC[, 1] ^ 2) + (block.loadings.PC[, 3] ^ 2)) + (block.loadings.PC[, 2] ^ 2))
        block.loadings.minmax = (block.loadings.len - min(block.loadings.len)) / (max(block.loadings.len) - min(block.loadings.len))

        plot3d(ref, add = F, col = rgb(block.loadings.minmax, 0.5, 1 - block.loadings.minmax), aspect = "iso", axes = F, xlab = "x", ylab = "y", zlab = "z", size = cex.3d)

        }

    }
    if(k == 2) {
      cl = round(sqrt(J))
      rw = ceiling(J / (round(sqrt(J))))
      layout(matrix(seq(1: (cl * rw)), nrow = rw, ncol = cl))
      pos.st = cumsum(c(1, nvert))
      pos.en = cumsum(nvert)
      for(i in 1: J) {
        A = arrayspecs(result@block.list[[i]], p = result@p[i], k = result@k[1])
        ref = mshape(A)

        block.loadings.bone = matrix(block.loadings[pos.st[i]: pos.en[i], ], ncol = ncol(block.loadings), byrow = F)
        block.loadings.PC = matrix(block.loadings.bone[, comp], ncol = k, byrow = T)
        block.loadings.PC = cbind(block.loadings.PC, rep(0, result@p[i]))
        block.loadings.len = sqrt(((block.loadings.PC[, 1] ^ 2) + (block.loadings.PC[, 3] ^ 2)) + (block.loadings.PC[, 2] ^ 2))
        block.loadings.minmax = (block.loadings.len - min(block.loadings.len)) / (max(block.loadings.len) - min(block.loadings.len))
        plot(x = ref[, 1], y = ref[, 2], col = rgb(block.loadings.minmax, 0.5, 1 - block.loadings.minmax), xlab = "x", ylab = "y", pch = 16, cex = cex.2d) # #
      }

    }
  }

  if(result@option == "rcpca") {
    #block.loadings = result@result$astar

    if(k == 3) {
      cl = round(sqrt(J))
      rw = ceiling(J / (round(sqrt(J))))
      mfrow3d(rw, cl)
      pos.st = cumsum(c(1, nvert))
      pos.en = cumsum(nvert)

      for(i in 1: J) {
        A = arrayspecs(result@block.list[[i]], p = result@p[i], k = result@k[1])
        ref = mshape(A)

        block.loadings.bone = matrix(block.loadings[[J + 1]][pos.st[i]: pos.en[i], ], ncol = ncol(block.loadings[[J + 1]]), byrow = F)
        block.loadings.PC = matrix(block.loadings.bone[, comp], ncol = k, byrow = T)
        block.loadings.len = sqrt(((block.loadings.PC[, 1] ^ 2) + (block.loadings.PC[, 3] ^ 2)) + (block.loadings.PC[, 2] ^ 2))
        block.loadings.minmax = (block.loadings.len - min(block.loadings.len)) / (max(block.loadings.len) - min(block.loadings.len))
        plot3d(ref, add = F, col = rgb(block.loadings.minmax, 0.5, 1 - block.loadings.minmax), aspect = "iso", axes = F, xlab = "x", ylab = "y", zlab = "z", size = cex.3d)

        }

    }
    if(k == 2) {
      cl = round(sqrt(J))
      rw = ceiling(J / (round(sqrt(J))))
      layout(matrix(seq(1: (cl * rw)), nrow = rw, ncol = cl))
      pos.st = cumsum(c(1, nvert))
      pos.en = cumsum(nvert)
      for(i in 1: J) {
        A = arrayspecs(result@block.list[[i]], p = result@p[i], k = result@k[1])
        ref = mshape(A)

        block.loadings.bone = matrix(block.loadings[[J + 1]][pos.st[i]: pos.en[i], ], ncol = ncol(block.loadings[[J + 1]]), byrow = F)
        block.loadings.PC = matrix(block.loadings.bone[, comp], ncol = k, byrow = T)
        block.loadings.PC = cbind(block.loadings.PC, rep(0, result@p[i]))
        block.loadings.len = sqrt(((block.loadings.PC[, 1] ^ 2) + (block.loadings.PC[, 3] ^ 2)) + (block.loadings.PC[, 2] ^ 2))
        block.loadings.minmax = (block.loadings.len - min(block.loadings.len)) / (max(block.loadings.len) - min(block.loadings.len))
        plot(x = ref[, 1], y = ref[, 2], col = rgb(block.loadings.minmax, 0.5, 1 - block.loadings.minmax), xlab = "x", ylab = "y", pch = 16, cex = cex.2d) # # #layout
      }

    }
  }
}
