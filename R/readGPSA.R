readGPSA = function(dirpath) {

  require(geomorph, quietly = TRUE, warn.conflicts = FALSE)
  require(Morpho, quietly = TRUE, warn.conflicts = FALSE)
  require(data.table, quietly = TRUE, warn.conflicts = FALSE)

  readpath = list.files(dirpath, pattern = "gpsa_homologized_points.dat", full.names = T)
  if (length(readpath) < 1) {
    stop("No homologized points file found. Remember to specify the path of just the directory containing the file, and not the path of the file itself")
  }
  step1 = fread(readpath)
  gpa.2D = matrix(as.numeric(unlist(step1)), nrow = dim(step1)[1])
  gpa.3D = arrayspecs(gpa.2D, p = dim(step1)[2] / 3, k = 3)

  centroid = c()
  for(i in 1: dim(gpa.2D)[1]) {
    centroid[i] = cSize(gpa.2D[i, ])
  }

  block_out = setClass("block", slots = c(gpa.3D = "array", gpa.2D = "matrix", raw = "array", centroid = "numeric", p = "numeric", k = "numeric", n = "numeric"))
  output = block_out(gpa.3D = gpa.3D, gpa.2D = gpa.2D, raw = gpa.2D, centroid = centroid, p = dim(gpa.3D)[1], k = dim(gpa.3D)[2], n = dim(gpa.3D)[3])

  return(output)

}