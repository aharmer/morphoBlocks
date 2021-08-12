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