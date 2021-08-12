formatBlock = function (block, curves = NULL, surfaces = NULL, cs = NULL, k = 2, gpa = T) {

  require(geomorph, quietly = TRUE, warn.conflicts = FALSE)
  require(Morpho, quietly = TRUE, warn.conflicts = FALSE)

  if (gpa == T) {

    if (length(dim(block)) == 3) {
      p = dim(block)[1]
      k = dim(block)[2]
      n = dim(block)[3]
      proc.all = gpagen(block, curves = curves, surfaces = surfaces, ProcD = F)
      gpa.2D = two.d.array(proc.all$coords)
      gpa.3D = arrayspecs(gpa.2D, p, k)

      centroid = c()
      for (i in 1: n) {
        centroid[i] = cSize(block[, , i])
      }

    }
    if (length(dim(block)) == 2) {
      p = dim(block)[2] / k
      k = k
      n = dim(block)[1]

      block = arrayspecs(block, p, k)
      proc.all = gpagen(block, curves = curves, surfaces = surfaces, ProcD = F)
      gpa.2D = two.d.array(proc.all$coords)
      gpa.3D = arrayspecs(gpa.2D, p, k)

      centroid = c()
      for (i in 1: n) {
        centroid[i] = cSize(block[, , i])
      }

    }

  } else {

    if (length(dim(block)) == 3) {
      p = dim(block)[1]
      k = dim(block)[2]
      n = dim(block)[3]
      gpa.3D = block
      gpa.2D = two.d.array(block)
      centroid = cs
    }
    if (length(dim(block)) == 2) {
      p = dim(block)[2] / k
      k = k
      n = dim(block)[1]

      block = arrayspecs(block, p, k)
      gpa.3D = block
      gpa.2D = two.d.array(block)
      centroid = cs

    }

  }

  block_out = setClass("block", slots = c(gpa.3D = "array", gpa.2D = "matrix", raw = "array", centroid = "numeric", p = "numeric", k = "numeric", n = "numeric"))
  output = block_out(gpa.3D = gpa.3D, gpa.2D = gpa.2D, raw = block, centroid = centroid, p = dim(gpa.3D)[1], k = dim(gpa.3D)[2], n = dim(gpa.3D)[3])

  return (output)

}