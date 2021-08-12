readPts = function (dirpath, landmarkRM = c(), gpa = T) {

  require(geomorph, quietly = TRUE, warn.conflicts = FALSE)
  require(Morpho, quietly = TRUE, warn.conflicts = FALSE)

  file.list = list.files(dirpath, full.names = T, pattern = ".pts")
  if (length(file.list) < 1) {
    stop("No .pts files found")
  }
  
  LM.read = as.matrix(read.pts(file.list[1]))
  LM.rm = c()
  for (i in 1: length(landmarkRM)) {
    LM.x = which(substr(rownames(LM.read), 0, 4) == landmarkRM[i])
    LM.rm = c(LM.rm, LM.x)

    LM.edit = LM.read
    if (length(landmarkRM) > 0) {
      LM.edit = LM.read[-c(LM.rm), ]
    }

    LM.fixed.ID = which(substring(unlist(dimnames(LM.edit)), 0, 1) == "S")
    LM.curve.ID = which(substring(unlist(dimnames(LM.edit)), 0, 1) == "C")
    LM.patch.ID = which(substring(unlist(dimnames(LM.edit)), 0, 1) == "P")

    LM.matrix.all = as.matrix(LM.edit)[1: nrow(LM.edit), ]
    for (i in 2: length(file.list)) {
      LM.read = as.matrix(read.pts(file.list[i]))
      LM.edit = LM.read
      if (length(landmarkRM) > 0) {
        LM.edit = LM.read[-c(LM.rm), ]
      }
      LM.matrix.all = rbind(LM.matrix.all, LM.edit[1: nrow(LM.edit), ])
    }
    LM.all = arrayspecs(unname(LM.matrix.all), nrow(LM.edit), 3)

    LM.curve.count = length(which(substring(rownames(LM.read), 0, 4) == "C000"))
    curves.map = c()
    for (i in ((1: (length(LM.curve.ID) / LM.curve.count)) * LM.curve.count) - (LM.curve.count - 1)) {
      m.curve = define.sliders(c(LM.curve.ID[i]: (LM.curve.ID[i] + (LM.curve.count - 1))))
      curves.map = rbind(curves.map, m.curve)
    }
    LM.fixed.ID = LM.fixed.ID
    LM.curve.ID = LM.curve.ID
    LM.patch.ID = LM.patch.ID

    if(gpa == T){
       proc.all = gpagen(LM.all, curves = curves.map, surfaces = LM.patch.ID, ProcD = F)
       gpa.2D = two.d.array(proc.all$coords)
       gpa.3D = arrayspecs(gpa.2D, p = dim(LM.all)[1], k = 3)
       } else {
       gpa.2D = two.d.array(LM.all)
       gpa.3D = arrayspecs(gpa.2D, p = dim(LM.all)[1], k = 3)
       }
  }
  centroid = c()
  for (i in 1: dim(LM.all)[3]) {
    centroid[i] = cSize(LM.all[, , i])
  }
  
  block_out = setClass("block", slots = c(gpa.3D = "array", gpa.2D = "matrix", raw ="array",  centroid = "numeric", p = "numeric", k = "numeric", n = "numeric", curves = "array", surfaces = "integer"))
  output = block_out(gpa.3D = gpa.3D, gpa.2D = gpa.2D, raw = LM.all, centroid = centroid, p = dim(gpa.3D)[1], k = dim(gpa.3D)[2], n = dim(gpa.3D)[3], curves = curves.map, surfaces = LM.patch.ID)

  return (output)

}