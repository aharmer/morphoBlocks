#' @title Combine and scale blocks of shape configuration data
#'
#' @description Combines two or more blocks of shape configuration data into a block list. Scales all data blocks using the normalised centroid size method from Profico et al. (2019) (and further described by Collyer et al. 2020). Such scaling occurs by default but is optional and can be prevented. A superblock is produced by column-wise concatenation of the individual blocks.
#'
#' @param blocks vector of 'block' object names produced by the \code{dodecBlock}, \code{formatBlock}, \code{readPts} or \code{readGPSA} functions. e.g. blocks = c(block1, block2, block3)
#' @param cent.scale a logical value indicating if the blocks should be scaled using the normalised centroid size method adapted from \code{geomorph::combine.subsets}. Default value \code{TRUE}.
#'
#' @details The block data is scaled using the normalised centroid size method from Profico et al. (2019) and Collyer et al. (2020) if \code{cent.scale = TRUE}. Dimensions of the blocks are retained for downstream analyses. \code{combineBlocks} adapts the normalised centroid size method from the \code{combine.subsets} function in the \code{geomorph} package (Adams and Ot?rola-Castillo 2013).
#'
#' @return a 'blockList' object, used for downstream analyses, that contains the data blocks and a superblock formed from the column-wise concatenation of all data blocks. The list contains the elements:
#' @return \item{block.list}{a list containing the data blocks and a concatenated superblockblock}
#' @return \item{p}{number of points in the configurations of each data block}
#' @return \item{k}{number of dimensions that the points in each configuration has}
#' @return \item{n}{number of configurations included in each data block}
#'
#' @references Adams DC, Ot?rola-Castillo E. 2013. geomorph: an \R package for the collection and analysis of geometric morphometric shape data. Methods in Ecology and Evolution 4:393-399 https://doi.org/10.1111/2041-210X.12035
#' @references Collyer ML, Davis MA, Adams DC. 2020. Making heads or tails of combined landmark configurations in geometric morphometric data. Evolutionary Biology 47: 193-205 https://doi.org/10.1007/s11692-020-09503-z
#' @references Profico A, Piras P, Buzi C, Del Bove A, Melchionna M, Senczuk G,  Varano V, Veneziano A, Raia P, Manzi G. 2019. Seeing the wood through the trees. Combining shape information from different landmark configurations. Hystrix, the Italian Journal of Mammalogy, 30: 157-165 https://doi.org/10.4404/hystrix-00206-2019
#'
#' @examples
#' block1 <- dodecBlock()
#' block2 <- dodecBlock()
#' combineBlocks(blocks = c(block1, block2))
#'
#' @export
combineBlocks <- function(blocks, cent.scale = TRUE) {
  p <- c()
  k <- c()
  n <- c()
  block.list <- c()

  if (length(blocks) == 0) {
    stop("No files were found")
  }
  if (length(blocks) == 1) {
    stop("Analysis requires two or more more blocks of data")
  }

  for (i in 1:length(blocks)) {
    if (class(blocks[[i]]) != "block") {
      stop("Objects are not the expected format. Use dodecBlock, formatBlock, readPts or readGPSA to format data into a block")
    }
    p[i] <- blocks[[i]]@p
    k[i] <- blocks[[i]]@k
    n[i] <- blocks[[i]]@n
  }
  if (length(unique(n)) > 1) {
    stop("Blocks contain different numbers of samples or need to be transposed")
  }

  if (sum(as.numeric(k < 2)) < 0) {
    stop("Items in list need to be 2- or 3-dimensional arrays")
  }

  if (cent.scale == TRUE) {
    step1 <- c()

    J <- length(blocks)
    for (j in 1:J) {
      step1 <- rbind(step1, blocks[[j]]@centroid / sqrt(blocks[[j]]@p))
    }

    step2 <- matrix(ncol = blocks[[1]]@n, nrow = J)
    for (j in 1:J) {
      for (i in 1:blocks[[1]]@n) {
        step2[j, i] <- sqrt(step1[j, i]^2 / (sum(step1[, i]^2)))
      }
      block.norm <- list()
      for (j in 1:J) {
        step3 <- blocks[[j]]@gpa.2D
        for (i in 1:blocks[[1]]@n) {
          step3[i, ] <- blocks[[j]]@gpa.2D[i, ] * step2[j, i]
        }
        block.norm[[j]] <- step3
      }
    }
    for (i in 1:J) {
      block.list[[i]] <- block.norm[[i]]
    }
    block.list[[J + 1]] <- do.call(cbind, block.list)
    names(block.list) <- c(paste(rep("block_"), LETTERS[1:J], sep = ""), "superblock")
  } else {
    J <- length(blocks)
    for (i in 1:J) {
      block.list[[i]] <- blocks[[i]]@gpa.2D
    }
    block.list[[J + 1]] <- do.call(cbind, block.list)
    names(block.list) <- c(paste(rep("block_"), LETTERS[1:J], sep = ""), "superblock")
  }

  output <- blockList_out(block.list = block.list, p = p, k = k, n = n)

  return(output)
}
