reMesh = function(dirpath, scale = 1) {

  #require(Morpho, quietly = TRUE, warn.conflicts = FALSE)
  #require(rgl, quietly = TRUE, warn.conflicts = FALSE)
  require(Rvcg, quietly = TRUE, warn.conflicts = FALSE)

  files = list.files(dirpath, full.names = TRUE, pattern = "\\.obj|ply|stl$")
  if (length(files) >= 1) {
    mesh.list = lapply(files, vcgImport)
  } else {
    stop("No mesh files detected")
  }

  vert = c()
  for(i in 1: length(mesh.list)) {
    vert[i] = dim(mesh.list[[i]] $vb)[2]
  }
  vertscale = min(vert) / vert

  dec.list = list()
  for(i in 1: length(mesh.list)) {
    dec.list[[i]] = vcgQEdecim(mesh.list[[i]], percent = vertscale[i], topo = FALSE)
  }

  dir.create(paste(dirpath, "/reMesh", sep = ""), showWarnings = F)

  rem = lapply(dec.list, function(mesh, scale) {
    edge.length = as.numeric(vcgMeshres(mesh)[1])
    #new.mesh = vcgUniformRemesh(mesh, voxelSize = edge.length * (1 / sqrt(scale)), multiSample = TRUE)
    new.mesh = vcgUniformRemesh(mesh, voxelSize = edge.length * scale, multiSample = TRUE)
    smooth.mesh = vcgSmooth(new.mesh, type = "surfPreserveLaplace", iteration = 1, delta = 0.3)
  }, scale)

  plynames = gsub(".obj|.stl", ".ply", list.files(dirpath, full.names = F, pattern = "\\.obj|ply|stl$"))
  for (i in 1: length(rem)) {
    vcgPlyWrite(rem[[i]], paste(dirpath, "/reMesh/", plynames[i], sep = ""))
  }
}