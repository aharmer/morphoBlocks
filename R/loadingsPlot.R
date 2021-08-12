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