#' @title Produce new decimated versions of a set of 3D meshes that all have the same number of evenly distributed vertices
#'
#' @description Generates new versions of 3D meshes (.obj, .ply or .stl) in a selected directory. The vertex count in the new meshes is reduced to the smallest vertex count amongst the selected meshes. Users may scale the vertex count down further. Vertices in the new meshes are approximately uniformly distributed. Meshes are automatically exported in polygon file format (.ply). \code{decMesh} is useful for studies that want to reduce the vertex count of 3D meshes before applying whole-mesh analyses (e.g. Generalized Procrustes Surface Analysis from Pomidor et al. 2016).
#' @param dirpath	the directory path where 3D meshes with .obj, .ply, or .stl extensions will be found
#' @param scale	a proportion for specifying how the vertex counts in the new versions of the 3D meshes should be scaled relative to the smallest vertex count amongst the 3D meshes. For scale = 1 (default), all new versions of meshes will have vertex counts that are similar to the vertex count of the smallest mesh. \code{decMesh} expects a scale value to be between 0 and 1
#' @param vsize	a constant used to start the remesh process. Default is 0.25. Should only need to be adjusted during function troubleshooting
#'
#' @details \code{decMesh} can be used for 3D mesh preprocessing before computationally demanding analyses. \code{decMesh} is a wrapper for \code{pcAlign} from the \code{Morpho} package (Schlager, 2017) and \code{vcgImport}, \code{vcgPlyWrite}, \code{vcgQEdecim} and \code{vcgUniformRemesh} from the \code{Rvcg} package (Schlager, 2017). Users may consider Instant Meshes as an alternative to \code{decMesh} (Wenzel et al. 2015).
#'
#' @note Function involves several iterative steps and can therefore be slow.
#'
#' @examples
#' # Any directory path where 3D meshes with .obj, .ply, or .stl extensions will be found
#' path = "/."
#'
#' # Example 1
#' # Vertex count in new meshes becomes the same as the smallest vertex count amongst meshes
#' decMesh(path)
#'
#' # Example 2
#' # Vertex count in new meshes becomes 10 percent of smallest vertex count amongst meshes
#' decMesh(path, scale = 0.1)
#'
#' @references
#' \itemize{
#'   \item Pomidor BJ, Makedonska J, Slice DE. 2016. A landmark-free method for three-dimensional shape analysis. PLoS One 11: e0150368 https://doi.org/10.1371/journal.pone.0150368
#'   \item Schlager S. 2017. Morpho and Rvcg-shape analysis in \R. In Zheng G, Li S, Sz?kely (eds.) Statistical shape and deformation analysis. Academic Press, London. Pp. 217-256.
#'   \item Wenzel J, Tarini M, Panozzo D, Sorkine-Hormung O. 2015. Instant field-alighed meshes. ACM Transactions on Graphics (Proceedings of SIGGRAPH Asia 2015) 34: 189-181.
#'   }
#'
#' @import Morpho
#' @import Rvcg
#' @export
decMesh = function(dirpath, scale = 1, vsize = 0.25) {

  require(Rvcg, quietly = TRUE, warn.conflicts = FALSE)
  require(Morpho, quietly = TRUE, warn.conflicts = FALSE)

  if(scale > 1){
    stop("Scale should be between 0 and 1")
  }

  files = list.files(dirpath, full.names = TRUE, pattern = "\\.obj|ply|stl$")
  if(length(files) >= 1) {
    mesh.list = lapply(files, vcgImport)
  } else {
    stop("No mesh files detected")
  }

  vert = c()
  for (i in 1: length(mesh.list)) {
    vert[i] = dim(mesh.list[[i]] $vb)[2]
  }
  vertscale = min(vert) / vert

  ref = mesh.list[[which(vertscale == max(vertscale))]]

  dec = c()
  dec.list = list()
  rem.list = list()
  for (j in 1: length(mesh.list)) {

    iter = 3
    align = mesh.list[[j]]

    for (i in 1:iter) {
      align = pcAlign(align, ref, optim = TRUE, subsample = 1000, iterations = 10)
      zscale = as.numeric(na.omit((colMeans(meshcube(ref)[c(3, 4, 7, 8), ]) - colMeans(meshcube(ref)[c(1, 2, 5, 6), ])) / (colMeans(meshcube(align)[c(3, 4, 7, 8), ]) - colMeans(meshcube(align)[c(1, 2, 5, 6), ]))))
      align = scalemesh(align, zscale, center = "bbox")
    }

    vsize = vsize
    sc = 1.1
    while (sc > 1) {
      rem = vcgUniformRemesh(align, voxelSize = vsize, multiSample = TRUE)
      rem.list[[j]] = rem
      sc = min(vert) / dim(rem$vb)[2]
      vsize = vsize - vsize * 0.25
    }

    dec = vcgQEdecim(rem, percent = sc * scale, topo = TRUE, quality = TRUE)
    dec.list[[j]] = dec

  }

  dir.create(paste(dirpath, "/decMesh", sep = ""), showWarnings = FALSE)

  plynames = gsub(".obj|.stl", ".ply", list.files(dirpath, full.names = F, pattern = "\\.obj|ply|stl$"))
  for(i in 1: length(dec.list)) {
    vcgPlyWrite(dec.list[[i]], paste(dirpath, "/decMesh/", plynames[i], sep = ""))
  }

}
