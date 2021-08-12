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