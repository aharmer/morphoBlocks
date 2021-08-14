#' @title Plot a phylomorphospace for a component from an analysis of multiple data blocks
#'
#' @description A plot to visualise the evolutionary relationships between specimens that are represented by values projected into a morphospace (e.g. Sidlauskas, 2008). Here the projected values are from a single component from an analysis of multiple data blocks. Rather than generating a two- or three-dimensional morphospace, this function plots consensus or principal component values against tree distance to show changes in the component value across the evolutionary history of the clade. The component values that are visualised using this function will either be from 1) the consensus space of an analysis performed with regularised consensus principal component analysis, or 2) will be from an analysis of a superblock (i.e. column-wise concatenation of individual data blocks) if principal component analysis was performed instead.
#'
#' @param result result produced by the \code{analyseBlocks} function.
#' @param phy a phylogenetic tree (class “phylo”) with \code{n} tips, where \code{n} is the number of samples in each block. Tip labels must be included and must correspond to the names of the samples in each block.
#' @param n.names a vector of character strings where entries in the vector represent names of the \code{n} samples present in each block. IMPORTANT NOTE: The sequence of character strings in \code{n.names} must match the sequence that samples are named in each block, and the names must exactly match the tip labels in the phylogenetic tree. For example, consider a block containing data from five samples where the first, second and third rows of the block contain data from samples "A", "B" and "C", and the fourth and fifth rows contain data from samples "D" and "E". Consequently, 1) \code{n.names = c("A", "B", "C", "D", "E")}, 2) \code{phy} should only include five tips, and 3) \code{phy$tip.labels} should be "A", "B", "C", "D" and "E" (in any order).
#' @param  comp the component selected to be shown in the phylomorphospace plot. Default is component one. The selected component must be within the range of components calculated by \code{analyseBlocks}.
#' @param pcol optional colour value (integer, hex code, colour name) or vector of colour values to be applied to the points in the phylomorphospace plot. If no value is specified then points will be coloured in a gradient from black to green according to their sequence in the data blocks. If a single integer is supplied (e.g. \code{pcol = 1} or \code{pcol = "red"} or \code{pcol = "#ffffff"} then all points will have the same colour. If a vector of length \code{n} is supplied (e.g. \code{pcol = 1:10}) then each point will be coloured by a value corresponding to its position in the vector.
#' @param xlab character string displayed along the horizontal axis of the phlyomorphospace plot. Default is "Time (Ma)".
#' @param label a logical value indicating if the points in the phylomorphospace should be labelled using \code{n.names}. Default value is \code{off} (i.e. names are not shown).
#' @param fsize size of the \code{n.names} labels if they are included in the plot. Default is 0.8.
#'
#' @details \code{phylomsPlot} helps to visualise the result from the \code{analyseBlocks} function by presenting an evolutionary context for score values from a consensus space (for \code{option = "rcpca"} in \code{analyseBlocks}) or from a concatenated superblock (for \code{option = "pca"} in \code{analyseBlocks}).
#' \code{phylomsPlot} is a wrapper for \code{phylomorphospace} from the \code{phytools} package (Revell, 2012), which is in turn based on the projection of a phylogenetic tree into a morphospace by Sidlauskas (2008). The function \code{phylomsPlot} uses the \code{distRoot} function from the \code{adephylo} package (Jombart et al. 2010). \code{phylomsPlot} does not specifically use functions from the \code{ape} package (Paradis et al. 2004), but this package may be called to read \code{phy} into R.
#'
#' @examples
#' # Simulate a phylogenetic tree with 20 tips
#' library(phytools)
#' phy = pbtree(b = 0.1, n = 20)
#'
#' # Simulate, combine and analyse two data blocks with 20 samples each
#' block1 = dodecBlock(n = 20)
#' block2 = dodecBlock(n = 20)
#' blocklist = combineBlocks(blocks = c(block1, block2))
#' result = analyseBlocks(blocklist)
#'
#' # Simulate a vector of names (here just use the tip labels)
#' n.names = phy$tip.label
#'
#' # Generate phylomorphospace plot
#' phylomsPlot(result, phy, n.names, comp = 1, xlab = "Tip to root distance")
#'
#' @references
#' Jombart T, Balloux F, Dray S. 2010. adephylo: new tools for investigating the phylogenetic signal in biological traits. Bioinformatics 26: 1907–1909 https://doi.org/10.1093/bioinformatics/btq292
#' Paradis E, Claude J, Strimmer K. 2004. APE: Analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289–290 https://doi.org/10.1093/bioinformatics/btg412
#' Revell LJ. 2012. phytools: An R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution 3: 217–223 https://doi.org/10.1111/j.2041-210X.2011.00169.x
#' Sidlauskas B. 2008. Continuous and arrested morphological diversification in sister clades of characiform fishes: a phylomorphospace approach. Evolution 62: 3135–3156 https://doi.org/10.1111/j.1558-5646.2008.00519.x
#'
#' @import phytools
#' @import adephylo
#' @export
phylomsPlot = function(result, phy, n.names, comp = 1, pcol = NULL, xlab = "Time (Ma)", label="off", fsize=0.8) {

  if(class(result) != "blockResult") {
    stop("Object is not the expected format. Use analyseBlock function to first analyse the data")
  }

  require(phytools, quietly = TRUE, warn.conflicts = FALSE)
  require(adephylo, quietly = TRUE, warn.conflicts = FALSE)

  n = result@n[1]
  scores = result@scores

  if(result@option == "rcpca") {
   comp.scores = scores[[length(scores)]][,comp]
   }

  if(result@option == "pca") {
    comp.scores = scores[,comp]
    }

  if(length(pcol) == 0) {
    pcol = rep(0, n)
    pcol = colorRampPalette(c(rgb(0, 1, 0, 1), rgb(0, 0, 0, 1)), alpha = TRUE)(n)
  }
  if(length(pcol) == 1) {
    pcol = rep(pcol, n)
  }
  if(length(pcol) > 1) {
    if(length(pcol) < n) {
      stop("pcol has fewer than expected values")
    }
    pcol = pcol
  }

  dat = matrix(ncol=2,nrow=n)
  rownames(dat) = phy$tip.label
  dat.sort=match(phy$tip.label,n.names)
  dat[,1] = max(as.numeric(distRoot(phy)))-as.numeric(distRoot(phy))
  dat[,2] = comp.scores[dat.sort]

  tip.cols = pcol[dat.sort]
  names(tip.cols) = phy$tip.label
  cols = c(tip.cols,rep("black", phy$Nnode))
  names(cols) = 1:(length(phy$tip) + phy$Nnode)

  phylomorphospace(phy, dat,label = label, fsize = fsize, control=list(col.node=cols), node.size = c(0,1.8), xlab = xlab , ylab = paste("Component",comp))
}
