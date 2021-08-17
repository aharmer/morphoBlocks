#' @title Analyse a list of data blocks containing shape information
#'
#' @description Performs either regularised consensus principal component analysis on a list of data blocks containing shape information, or performs principal component analysis on a superblock produced by the column-wise concatenation of the individual data blocks.

#' @param blockList list of ‘block’ objects produced by the \code{dodecBlock}, \code{formatBlock}, \code{readPts} or \code{readGPSA} functions.
#' @param option either \code{option = "rcpca"} (default) for regularised consensus principal component analysis in mode 2 applied to the entire block list, or \code{option = "pca"} for principal component analysis applied to the superblock item in the block list.
#' @param ncomp an integer specifying how many components should be calculated (default is 3)

#' @details \code{analyseBlocks} is applied to an object of class "blockList" produced by the \code{combineBlocks} function and has two options: 1) \code{option = "rcpca"} and 2) \code{option = "pca"}. \code{option = "rcpca"} will perform regularised consensus principal component analysis using the \code{rgcca} function from the \code{RGCCA} package (Tenenhaus and Guillemot 2017), and is the default option for \code{analyseBlocks}. The \code{rgcca} function itself has many options that each perform a different type of analysis. Here the \code{analyseBlocks} function is specifically calling the regularised consensus principal component analysis in mode 2 option with scaling applied. For further detail see Tenenhaus and Guillemot (2017) and Tenenhaus et al. (2017). \code{option = "pca"} will perform principal component analysis on the superblock item in the block list using the \code{prcomp} function from \code{base} \R.
#'
#'
#' @return A blockResult object containing output from the regularised consensus principal component analysis or principal component analysis. The list contains the elements:
#' \itemize{
#'   \item \code{result} output from the regularised consensus principal component analysis in mode 2 produced by the \code{rgcca} function from the \code{RGCCA} package, or output from principal component analysis produced by the \code{prcomp} function in \code{base} \R.
#'   \item \code{option} either "rcpca" or "pca".
#'   \item \code{block.list} a list containing the data blocks and a concatenated superblock. Inherited from the supplied blockList object and retained for downstream analyses.
#'   \item \code{scores} component score values (for individual blocks and the consensus if \code{option = "rpca"}; for the superblock if \code{option = "pca"}) (see \code{scoresPlot} for more detail).
#'   \item \code{block.loadings} component loadings (for individual blocks and the consensus if \code{option = "rpca"}; for the superblock if \code{option = "pca"}) (see \code{loadingsPlot} for more detail).
#'   \item \code{p} number of points in the configurations of each data block. Inherited from the supplied blockList object and retained for downstream analyses.
#'   \item \code{k} number of dimensions that the points in each configuration has. Inherited from the supplied blockList object and retained for downstream analyses.
#'   \item \code{n} number of configurations included in each data block. Inherited from the supplied blockList object and retained for downstream analyses.
#'   }
#'
#' @examples
#' block1 = dodecBlock()
#' block2 = dodecBlock()
#' blocklist = combineBlocks(blocks = c(block1, block2))
#' result1 = analyseBlocks(blocklist)
#' result2 = analyseBlocks(blocklist, option = "pca", ncomp = 10)
#'
#' @references
#' \itemize{
#'   \item Tenenhaus A, Guillemot V. 2017. RGCCA: Regularized and sparse generalized canonical correlation analysis for multiblock data 2.1.2. https://CRAN.R-project.org/package=RGCCA.
#'   \item Tenenhaus M, Tenenhaus A, Groenen PJF. 2017. Regularized generalized canonical correlation analysis: A framework for sequential multiblock component methods. Psychometrika 82: 737-777 https://doi.org/10.1007/s11336-017-9573-x
#'  }
#'
#' @import RGCCA
#' @export
analyseBlocks = function(blockList, option = "rcpca", ncomp = 3) {
  if(class(blockList) != "blockList") {
    stop("Object is not the expected format. Use combineBlocks function to format data into a block list")
  }

  p = blockList@p
  k = blockList@k
  n = blockList@n
  blockList = blockList@block.list

  if(option == "pca") {

    result = prcomp(blockList[[length(blockList)]])
    scores = result$x
    block.loadings = result$rotation
    # setClass("prcomp")

    # blockResult_out = setClass("blockResult", slots = c(result = class(result), option = "character", block.list = "list", scores = "matrix", block.loadings = "matrix", p = "integer", k = "integer", n = "integer"))

    output = list(result = result, option = option, block.list = blockList, scores = scores, block.loadings = block.loadings, p = p, k = k, n = n)
  }

  if(option == "rcpca") {
    require(RGCCA, quietly = TRUE, warn.conflicts = FALSE)
    #Determine the number of blocks that have been included in the block list
    J = length(names(blockList)) - 1

    #Build the design matrix expected by the rgcca() function for RCPCA with m = 2
    C = matrix(rep(0, (J + 1) * (J + 1)), ncol = (J + 1))
    C[(J + 1), c(1: J)] = 1
    C[c(1: J), (J + 1)] = 1

    # Set the shrinkage parameters expected by the rgcca() function for RCPCA with m = 2
    tau = rep(1, length(blockList))

    # Perform RCPCA mode 2
    result = rgcca(blockList, C = C, tau = tau, ncomp = rep(ncomp, J + 1), scheme = function (x) x ^ 2, scale = TRUE, verbose = FALSE)
    scores = result$Y
    block.loadings = result$astar
    # setClass("rgcca")

    # blockResult_out = setClass("blockResult", slots = c(result = class(result), option = "character", block.list = "list", scores = "list", block.loadings = "list", p = "integer", k = "integer", n = "integer"))

    output = list(result = result, option = option, block.list = blockList, scores = scores, block.loadings = block.loadings, p = p, k = k, n = n)
  }

  # output = blockResult_out(result = result, option = option, block.list = blockList, scores = scores, block.loadings = block.loadings, p = p, k = k, n = n)

  return (output)
}
