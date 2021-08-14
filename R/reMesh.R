#' @title Produce new versions of 3D meshes with similar numbers of evenly distributed vertices
#'
#' @description Generates new versions of 3D meshes (.obj, .ply or .stl) in a selected directory. The vertex count in the new meshes is reduced to become similar to the smallest vertex count amongst the selected meshes. Vertices in the new meshes are uniformly distributed. Meshes are automatically exported in polygon file format (.ply). \code{reMesh} is useful for studies that want to apply whole-mesh analyses (e.g. Generalized Procrustes Surface Analysis from Pomidor et al. 2016) to a set of 3D meshes that have very large vertex counts, or have large disparity in vertex counts between meshes.
#'
#' @param dirpath the directory path where two or 3D meshes with .obj, .ply, or .stl extensions will be found.
#' @param scale a proportion for specifying how the vertex counts in the new versions of the 3D meshes should be scaled relative to the smallest vertex count amongst the 3D meshes. For \code{scale = 1} (default), all new versions of meshes will have vertex counts that are similar to the vertex count of the smallest mesh.
#'
#' @details \code{reMesh} can be used for 3D mesh preprocessing before computationally demanding analyses.
#' \code{reMesh} is a wrapper for \code{vcgImport}, \code{vcgMeshres}, \code{vcgPlyWrite} and \code{vcgQEdecim} from the \code{Rvcg} package (Schlager, 2017).
#'
#' @references
#' Pomidor BJ, Makedonska J, Slice DE. 2016. A landmark-free method for three-dimensional shape analysis. PLoS One 11: e0150368 https://doi.org/10.1371/journal.pone.0150368
#' Schlager S. 2017. Morpho and Rvcg–shape analysis in R. In Zheng G, Li S, Székely (eds.) Statistical shape and deformation analysis. Academic Press, London. Pp. 217–256.
#' @import Rvcg
#' @export
reMesh = function(dirpath, scale = 1) {

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
