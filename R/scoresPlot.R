#' @title Plot component scores from an analysis of multiple data blocks
#'
#' @description Bivariate scatterplot (or plots) of the score values from an analysis of data blocks containing shape information
#'
#' @param result result produced by the analyseBlocks function.
#' @param comp the components selected to be shown along the scatter plot axes. Default is component one and component two. The first selected component will be shown along the horizontal axis and the second selected component will be shown along the vertical axis. The selected components must be within the range of components calculated by \code{analyseBlocks}.
#' @param pcol optional colour value (integer, hex code, colour name) or vector of colour values to be applied to the points in the scatterplot. If no value is specified then points will be coloured in a gradient from black to green according to their sequence in the data blocks. If a single integer is supplied (e.g. \code{pcol = 1} or \code{pcol = "red"} or \code{pcol = "#ffffff"}) then all points will have the same colour. If a vector of length \code{n} is supplied (e.g. \code{pcol = 1:10}) then each point will be coloured by a value corresponding to its position in the vector.
#' @param plabels optional user-supplied vector of labels for labelling points in scores plot.
#' @param consensus.only a logical value (default value \code{FALSE}), relevant only if analyses were performed using regularised consensus principal component analysis (i.e. \code{analyseBlocks}, \code{option = 'rcpca'}). If \code{TRUE}, only plot the global scores and the consensus space.
#'
#' @details \code{scoresPlot} helps to visualise the result from the \code{analyseBlocks} function and gives a different result depending on whether 1) \code{option = "rcpca"} or 2) \code{option = "pca"}) was used for \code{analyseBlocks}.
#' 1) If \code{option = "pca"} was used for \code{analyseBlocks} then a scatterplot will be produced for the selected components from each individual block and for the consensus space. Axes will display the average variance explained by the components of the individual blocks and the consensus. Average variance explained by a selected component is presented as a proportion of the total variance. For detail about average variance explained see Tenenhaus and Guillemot (2017) and Tenenhaus et al. (2017).
#' 2) If \code{option = "pca"} was selected for \code{analyseBlocks} then a scatterplot will be produced for the selected components from the analysis of the superblock. Axes will display the variance explained by the components of the superblock. Variance explained by a selected component is presented as a proportion of the total variance explained by all components.
#'
#' @return a two-dimensional scatterplot
#'
#' @examples
#' block1 = dodecBlock()
#' block2 = dodecBlock()
#' blocklist = combineBlocks(blocks = c(block1, block2))
#' result1 = analyseBlocks(blocklist)
#' result2 = analyseBlocks(blocklist, option = "pca", ncomp = 10)
#' scoresPlot(result1)
#' dev.off()
#' scoresPlot(result2, comp = c(2,3), pcol = colorRampPalette(c(rgb(1, 0.7, 0,1), rgb(0, 0, 1, 1)), alpha = TRUE)(result2@n[1]))
#'
#' @references
#' \itemize{
#'   \item Tenenhaus A, Guillemot V. 2017. RGCCA: Regularized and sparse generalized canonical correlation analysis for multiblock data 2.1.2. https://CRAN.R-project.org/package=RGCCA.
#'   \item Tenenhaus M, Tenenhaus A, Groenen PJF. 2017. Regularized generalized canonical correlation analysis: A framework for sequential multiblock component methods. Psychometrika 82: 737-777 https://doi.org/10.1007/s11336-017-9573-x
#'   }
#'
#' @export
scoresPlot = function (result, comp = c(1, 2), pcol = NULL, plabels = NULL, consensus.only = F) {
  if (class(result) != "blockResult") {
    stop("Object is not the expected format. Use analyseBlock function to first analyse the data")
  }

  n = result@n[1]
  scores = result@scores
  compx = comp[1]
  compy = comp[2]

  if (length(pcol) == 0) {
    pcol = rep(0, n)
    pcol = colorRampPalette(c(rgb(0, 1, 0, 1), rgb(0, 0, 0, 1)), alpha = TRUE)(n)
  }
  if (length(pcol) == 1) {
    pcol = rep(pcol, n)
  }
  if (length(pcol) > 1) {
    if (length(pcol) < n) {
      stop("pcol has fewer than expected values")
    }
    pcol = pcol
  }

  if (result@option == "rcpca") {
    if (consensus.only == F) {
      J = length(result @result$Y)
      cl = round(sqrt(J))
      rw = ceiling(J / (round(sqrt(J))))
      layout(matrix(1: (rw * cl), nrow = rw, ncol = cl, byrow = T))

      AVE = result@result$AVE

      for (i in 1: (J - 1)) {
        plot(scores[[i]][, compx], scores[[i]][, compy], xlab = paste("Component ", compx, " (", round(AVE[[1]][[i]][compx], 3), " average variance explained)", sep = ""),
          ylab = paste("Component ", compy, " (", round(AVE[[1]][[i]][compy], 3), "  average variance explained)", sep = ""),
          main = paste("Block", LETTERS[i]), col = "black", pch = 21, bg = pcol, cex = 2)
        if (length(plabels) > 0) {
          text(scores[[i]][, compx], scores[[i]][, compy], labels = plabels, pos = 2)
        }
      }
      plot(scores[[J]][, compx], scores[[J]][, compy], xlab = paste("Global component ", compx, " (", round(AVE[[1]][[i]][compx], 3), " average variance explained)", sep = ""),
        ylab = paste("Global component ", compy, " (", round(AVE[[1]][[J]][compy], 3), "  average variance explained)", sep = ""),
        main = "Consensus", col = "black", pch = 21, bg = pcol, cex = 2)
      if (length(plabels) > 0) {
        text(scores[[J]][, compx], scores[[J]][, compy], labels = plabels, pos = 2)
      }
    }

    if (consensus.only == T) {
      J = length(result@result$Y)
      AVE = result@result$AVE

      plot(scores[[J]][, compx], scores[[J]][, compy], xlab = paste("Global component ", compx, " (", round(AVE[[1]][[i]][compx], 3), " average variance explained)", sep = ""),
        ylab = paste("Global component ", compy, " (", round(AVE[[1]][[J]][compy], 3), "  average variance explained)", sep = ""),
        main = "Consensus", col = "black", pch = 21, bg = pcol, cex = 2)
      if (length(plabels) > 0) {
        text(scores[[J]][, compx], scores[[J]][, compy], labels = plabels, pos = 2)
      }
    }
  }

  if (result@option == "pca") {
    PCx.exp = (round(result@result$sdev ^ 2 / (sum(result@result$sdev ^ 2)), 3))[compx]
    PCy.exp = (round(result@result$sdev ^ 2 / (sum(result@result$sdev ^ 2)), 3))[compy]

    plot(scores[, compx], scores[, compy], ann = FALSE, col = "black", pch = 21, bg = pcol, cex = 2)
    title(main = "Superblock")
    title(xlab = paste("Principal component ", compx, " (", PCx.exp, " variance explained)", sep = ""))
    title(ylab = paste("Principal component ", compy, " (", PCy.exp, " variance explained)", sep = ""))
    if (length(plabels) > 0) {
      text(scores[, compx], scores[, compy], labels = plabels, pos = 2)
    }

  }
}
