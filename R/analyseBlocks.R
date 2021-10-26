#' @title Analyse a list of data blocks containing shape information
#'
#' @description Performs either Regularized Consensus Principal Component Analysis on a list of data blocks containing shape information, or performs principal component analysis on a superblock produced by the column-wise concatenation of the individual data blocks.

#' @param blockList list of 'block' objects produced by the \code{dodecBlock}, \code{formatBlock}, \code{readPts} or \code{readGPSA} functions.
#' @param option either \code{option = "rcpca"} (default) for Regularized Consensus Principal Component Analysis in mode 2 applied to the entire block list, or \code{option = "pca"} for principal component analysis applied to the superblock item in the block list.
#' @param ncomp an integer specifying how many components should be calculated (default is 3)

#' @details \code{analyseBlocks} is applied to an object of class "blockList" produced by the \code{combineBlocks} function and has two options: 1) \code{option = "rcpca"} and 2) \code{option = "pca"}. The \code{option = "rcpca"} will perform Regularized Consensus Principal Component Analysis using the \code{rgcca} function from the \code{RGCCA} package (Tenenhaus and Guillemot 2017), and is the default option for \code{analyseBlocks}. The \code{rgcca} function itself has many options that each perform a different type of analysis. Here the \code{analyseBlocks} function is specifically calling the Regularized Consensus Principal Component Analysis in mode 2 option with scaling applied. For further detail see Tenenhaus and Guillemot (2017) and Tenenhaus et al. (2017). \code{option = "pca"} will perform principal component analysis on the superblock item in the block list using the \code{prcomp} function from \code{base} \R.
#'
#' @return A list object containing output from the Regularized Consensus Principal Component Analysis or principal component analysis. The list contains the elements:
#' @return \item{result}{output from the Regularized Consensus Principal Component Analysis in mode 2 produced by the \code{rgcca} function from the \code{RGCCA} package, or output from principal component analysis produced by the \code{prcomp} function in \code{base} \R.}
#' @return \item{option}{either "rcpca" or "pca".}
#' @return \item{block.list}{a list containing the data blocks and a concatenated superblock. Inherited from the supplied blockList object and retained for downstream analyses.}
#' @return \item{scores}{component score values (for individual blocks and the consensus if \code{option = "rpca"}; for the superblock if \code{option = "pca"}) (see \code{scoresPlot} for more detail).}
#' @return \item{block.loadings}{component loadings (for individual blocks and the consensus if \code{option = "rpca"}; for the superblock if \code{option = "pca"}) (see \code{loadingsPlot} for more detail).}
#' @return \item{p}{number of points in the configurations of each data block. Inherited from the supplied blockList object and retained for downstream analyses.}
#' @return \item{k}{number of dimensions that the points in each configuration has. Inherited from the supplied blockList object and retained for downstream analyses.}
#' @return \item{n}{number of configurations included in each data block. Inherited from the supplied blockList object and retained for downstream analyses.}
#'
#' @references Tenenhaus A, Guillemot V. 2017. RGCCA: Regularized and Sparse Generalized Canonical Correlation Analysis for multiblock data 2.1.2. https://CRAN.R-project.org/package=RGCCA.
#' @references Tenenhaus M, Tenenhaus A, Groenen PJF. 2017. Regularized Generalized Canonical Correlation Analysis: A framework for sequential multiblock component methods. Psychometrika 82: 737-777 https://doi.org/10.1007/s11336-017-9573-x
#'
#' @examples
#' block1 <- dodecBlock()
#' block2 <- dodecBlock()
#' blocklist <- combineBlocks(blocks = c(block1, block2))
#' result1 <- analyseBlocks(blocklist)
#' result2 <- analyseBlocks(blocklist, option = "pca", ncomp = 10)
#'
#' @importFrom RGCCA rgcca
#' @importFrom stats prcomp
#' @export
analyseBlocks <- function(blockList, option = "rcpca", ncomp = 3) {
  if (class(blockList) != "blockList") {
    stop("Object is not the expected format. Use combineBlocks function to format data into a block list")
  }

  p <- blockList@p
  k <- blockList@k
  n <- blockList@n
  blockList <- blockList@block.list

  if (option == "pca") {
    result <- prcomp(blockList[[length(blockList)]])
    scores <- result$x
    block.loadings <- result$rotation

    output <- list(result = result, option = option, block.list = blockList, scores = scores, block.loadings = block.loadings, p = p, k = k, n = n)
  }

  if (option == "rcpca") {
    # require(RGCCA, quietly = TRUE, warn.conflicts = FALSE)
    J <- length(names(blockList)) - 1

    C <- matrix(rep(0, (J + 1) * (J + 1)), ncol = (J + 1))
    C[(J + 1), c(1:J)] <- 1
    C[c(1:J), (J + 1)] <- 1

    tau <- rep(1, length(blockList))

    result <- rgcca(blockList, C = C, tau = tau, ncomp = rep(ncomp, J + 1), scheme = function(x) x^2, scale = TRUE, verbose = FALSE)
    scores <- result$Y
    block.loadings <- result$astar

    output <- list(result = result, option = option, block.list = blockList, scores = scores, block.loadings = block.loadings, p = p, k = k, n = n)
  }

  return(output)
}
