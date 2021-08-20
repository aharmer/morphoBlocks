#' @title Simulate a data block of landmark configurations
#'
#' @description Generates a block of \emph{n} dodecahedra where the vertices of each dodecahedron represent a landmark configuration with \emph{p = 20} points in \emph{k = 3} dimensions. The data block is transformed using generalised Procrustes analysis.
#'
#' @param n	the number of configurations to generate for the data block. Default is 10. Must be 2 or greater.
#' @param dist_x,dist_y,dist_z optional integer or vector to specify a distance along the stated axis (i.e. \code{dist_x} for the x axis) for translating the configurations. If no value is specified then configurations will not be translated. If a single integer is supplied (e.g. \code{dist_x = 1}) then all configurations will be translated along an axis by the same amount. If a vector of length \emph{n} is supplied (e.g. \code{dist_x = 1:10}) then each configuration will be translated by a value corresponding to its position in the vector.
#' @param theta_x,theta_y,theta_z optional integer or vector to specify an angle for rotating the configuration about its origin in the stated axis (i.e. \code{theta_x} for the x axis). If no value is specified then configurations will not be rotated. If a single integer is supplied (e.g. \code{theta_x = 30}) then all configurations will be rotated by the same amount. If a vector of length \emph{n} is supplied (e.g. \code{theta_x = seq(from = 0, to = 210, by = 30)}) then each configuration will be rotated by a value corresponding to its position in the vector.
#' @param size optional integer or vector for scaling the configurations. If no value is specified then configurations will not be scaled relative to one another. If a single integer is supplied (e.g. \code{size = 1}) then all configurations will be scaled by the same amount. If a vector of length \emph{n} is supplied (e.g. \code{size = 1:10}) then each configuration will be scaled by a value corresponding to its position in the vector.
#' @param noise	optional integer for introducing noise to the configurations. If no value is specified then \code{runif(1, min = 0, max = (1 / 100))} is added to each of the \code{x}, \code{y} and \code{z} positions of each point in each configuration. This small amount of noise is added so that no configuration has identical coordinates. The amount of noise that is introduced to the data block can be increased by specifying a noise value that is greater than 1.
#' @param vertex_shift optional integer or vector to specify a distance along the z axis for translating the first vertex in each configuration. This argument can be used to introduce variation into the data block. If no value is specified then vertex 1 in each configuration will not be translated. If a single integer is supplied (e.g. \code{vertex_shift = 1}) then vertex 1 in all configurations will be translated along the z axis by the same amount. If a vector of length \emph{n} is supplied (e.g. \code{vertex_shift = 1:10}) then vertex 1 in each configuration will be translated by a value corresponding to its position in the vector. Things get really fun when you use both \code{vertex_shift} and \code{theta_z} together (see example 2 below).
#' @param plot a logical value indicating whether the configurations should be plotted. Default is FALSE.
#'
#' @details \code{dodecBlock} builds a data block of configurations that have been transformed using generalised Procrustes analysis. The resulting data block is a three dimensional array with \emph{n = 2} or greater configurations, with each configuration represented by \emph{p = 20} points that have \emph{k = 3} dimensions. Several objects are calculated for downstream analyses including centroid sizes for each configuration.
#' \code{formatBlock} uses the \code{cSize} function from the \code{Morpho} package (Schlager 2017), and the \code{gpagen}, \code{two.d.array} and \code{arrayspecs} functions from the \code{geomorph} package (Adams and Ot?rola-Castillo 2013).
#'
#' @return a 'blockList' object of \emph{n} dodecahedra, used for downstream analyses. The list contains the elements:
#' \itemize{
#'   \item \code{raw} configurations without Procrustes transformation
#'   \item \code{gpa.3D} configurations after Procrustes transformation organised into a 3D array
#'   \item \code{gpa.2D} configurations after Procrustes transformation organised into a 2D matrix
#'   \item \code{centroid} centroid sizes of the configurations
#'   \item \code{p} number of points in the configurations of each data block
#'   \item \code{k} number of dimensions that the points in each configuration has
#'   \item \code{n} number of configurations included in each data block
#'   }
#'
#' @examples
#' # Example 1: Generate a block with ten configurations
#' block = dodecBlock()
#' # Example 2: Generate a block with five configurations. Rotate the configurations and translate the first vertex in each configuration. Plot the result.
#' block2 = dodecBlock(n = 5, theta_z = seq(from = 0, to = 270, by = 60), vertex_shift = 0:4, plot = TRUE)
#'
#' @references
#' \itemize{
#'   \item Adams DC, Ot?rola-Castillo E. 2013. geomorph: an \R package for the collection and analysis of geometric morphometric shape data. Methods in Ecology and Evolution 4:393-399 https://doi.org/10.1111/2041-210X.12035
#'   \item Schlager S. 2017. Morpho and Rvcg-shape analysis in \R. In Zheng G, Li S, Sz?kely (eds.) Statistical shape and deformation analysis. Academic Press, London. Pp. 217-256.
#'   }
#'
#' @import Morpho
#' @import Rvcg
#' @export
dodecBlock = function(n = 10, dist_x = NULL, dist_y = NULL, dist_z = NULL, theta_x = NULL, theta_y = NULL, theta_z = NULL, size = NULL, noise = NULL, vertex_shift = NULL, plot = FALSE) {

  require(geomorph, quietly = TRUE, warn.conflicts = FALSE)
  require(Morpho, quietly = TRUE, warn.conflicts = FALSE)

  if(n < 2) {
  stop("n must be 2 or greater")
  }


  #Build the configuration base model
  p = 20
  k = 3
  dodec = matrix(nrow = p, ncol = k)
  colnames(dodec) = c("x", "y", "z")
  phi = (1 + sqrt(5)) / 2
  dodec[, 1] = c(phi, -phi, -phi, phi, 1 / phi, 1 / phi, -1 / phi, -1 / phi, 0, 0, 0, 0, 1, 1, -1, -1, -1, 1, 1, -1)
  dodec[, 2] = c(0, 0, 0, 0, phi, -phi, -phi, phi, 1 / phi, 1 / phi, -1 / phi, -1 / phi, 1, -1, -1, 1, 1, 1, -1, -1)
  dodec[, 3] = c(1 / phi, 1 / phi, -1 / phi, -1 / phi, 0, 0, 0, 0, phi, -phi, -phi, phi, 1, 1, 1, 1, -1, -1, -1, -1)

  dodec.block = array(dodec, dim = c(p, k, n))

  # Offset configurations relative to one another
  if(length(dist_x) == 0) {
    dist_x = rep(0, n)
  }
  if(length(dist_x) == 1) {
    dist_x = rep(dist_x, n)
  }
  if(length(dist_x) > 1) {
    if(length(dist_x) < n) {
      stop("dist_x has fewer than expected values")
    }
    dist_x = dist_x
  }
  for(i in 1: n) {
    dodec.block[, 1, i] = dodec.block[, 1, i] + dist_x[i]
  }

  if(length(dist_y) == 0) {
    dist_y = rep(0, n)
  }
  if(length(dist_y) == 1) {
    dist_y = rep(dist_y, n)
  }
  if(length(dist_y) > 1) {
    if(length(dist_y) < n) {
      stop("dist_y has fewer than expected values")
    }
    dist_y = dist_y
  }
  for(i in 1: n) {
    dodec.block[, 1, i] = dodec.block[, 1, i] + dist_y[i]
  }

  if(length(dist_z) == 0) {
    dist_z = rep(0, n)
  }
  if(length(dist_z) == 1) {
    dist_z = rep(dist_z, n)
  }
  if(length(dist_z) > 1) {
    if(length(dist_z) < n) {
      stop("dist_z has fewer than expected values")
    }
    dist_z = dist_z
  }
  for(i in 1: n) {
    dodec.block[, 1, i] = dodec.block[, 1, i] + dist_z[i]
  }

  # Scale the configurations
  if(length(size) == 0) {
    size = rep(1, n)
  }
  if(length(size) == 1) {
    size = rep(size, n)
  }
  if(length(size) > 1) {
    if(length(size) < n) {
      stop("size has fewer than expected values")
    }
    size = size
  }
  if(length(which(size == 0)) > 0) {
    stop("size is a multiplier applied to all coordinates in the configuration and cannot be 0")
  }

  for(i in 1: n) {
    dodec.block[, , i] = dodec.block[, , i] * size[i]
  }

  # Establish the level of noise added to each vertex in the configuration
  if(length(noise) == 0) {
    noise = rep(1, n)
  }
  if(length(noise) == 1) {
    noise = rep(noise, n)
  }
  if(length(noise) > 1) {
    if(length(noise) < n) {
      stop("noise has fewer than expected values")
    }
    noise = noise
  }

  for(i in 1: n) {
    dodec.block[, , i] = dodec.block[, , i] + matrix(runif((p * k), min = 0, max = (noise / 100)), ncol = k)
  }

  #Shift the z value of vertex 1
  if(length(vertex_shift) == 0) {
    vertex_shift = rep(1, n)
  }
  if(length(vertex_shift) == 1) {
    vertex_shift = rep(vertex_shift, n)
  }
  if(length(vertex_shift) > 1) {
    if(length(vertex_shift) < n) {
      stop("vertex_shift has fewer than expected values")
    }
    vertex_shift = vertex_shift
  }

  for(i in 1: n) {
    dodec.block[1, 3, i] = dodec.block[1, 3, i] + vertex_shift[i]
  }

  # Clockwise rotate configurations about the centre
  if(length(theta_x) == 0) {
    theta_x = rep(0, n)
  }
  if(length(theta_x) == 1) {
    theta_x = rep(theta_x, n)
  }
  if(length(theta_x) > 1) {
    if(length(theta_x) < n) {
      stop("theta_x has fewer than expected values")
    }
    theta_x = theta_x
  }
  if(length(theta_y) == 0) {
    theta_y = rep(0, n)
  }
  if(length(theta_y) == 1) {
    theta_y = rep(theta_y, n)
  }
  if(length(theta_y) > 1) {
    if(length(theta_y) < n) {
      stop("theta_y has fewer than expected values")
    }
    theta_y = theta_y
  }

  if(length(theta_z) == 0) {
    theta_z = rep(0, n)
  }
  if(length(theta_z) == 1) {
    theta_z = rep(theta_z, n)
  }
  if(length(theta_z) > 1) {
    if(length(theta_z) < n) {
      stop("theta_z has fewer than expected values")
    }
    theta_z = theta_z
  }
  r.step1 = dodec.block
  r.step2 = r.step1

  for(i in 1: n) {
    r.step2[, 2, i] = mean(r.step1[, 2, i]) + (cos(theta_x[i] * pi / 180) * (r.step1[, 2, i] - mean(r.step1[, 2, i])) + sin(theta_x[i] * pi / 180) * (r.step1[, 3, i] - mean(r.step1[, 3, i])))
    r.step2[, 3, i] = mean(r.step1[, 3, i]) + (-sin(theta_x[i] * pi / 180) * (r.step1[, 2, i] - mean(r.step1[, 2, i])) + cos(theta_x[i] * pi / 180) * (r.step1[, 3, i] - mean(r.step1[, 3, i])))
  }

  r.step3 = r.step2
  for(i in 1: n) {
    r.step3[, 1, i] = mean(r.step2[, 1, i]) + (cos(theta_y[i] * pi / 180) * (r.step2[, 1, i] - mean(r.step2[, 1, i])) + sin(theta_y[i] * pi / 180) * (r.step2[, 3, i] - mean(r.step2[, 3, i])))
    r.step3[, 3, i] = mean(r.step2[, 3, i]) + (-sin(theta_y[i] * pi / 180) * (r.step2[, 1, i] - mean(r.step2[, 1, i])) + cos(theta_y[i] * pi / 180) * (r.step2[, 3, i] - mean(r.step2[, 3, i])))
  }

  r.step4 = r.step3
  for(i in 1: n) {
    r.step4[, 1, i] = mean(r.step3[, 1, i]) + (cos(theta_z[i] * pi / 180) * (r.step3[, 1, i] - mean(r.step3[, 1, i])) + sin(theta_z[i] * pi / 180) * (r.step3[, 2, i] - mean(r.step3[, 2, i])))
    r.step4[, 2, i] = mean(r.step3[, 2, i]) + (-sin(theta_z[i] * pi / 180) * (r.step3[, 1, i] - mean(r.step3[, 1, i])) + cos(theta_z[i] * pi / 180) * (r.step3[, 2, i] - mean(r.step3[, 2, i])))
  }

  # Perform generalised Procrustes analysis
  proc.all = gpagen(r.step4, ProcD = FALSE, print.progress = FALSE)
  gpa.2D = two.d.array(proc.all$coords)
  gpa.3D = arrayspecs(gpa.2D, p, k)

  centroid = c()
  for(i in 1: n) {
    centroid[i] = cSize(r.step4[, , i])
  }

  if(plot == T) {
    plot3d(gpa.3D[, , 1], size = 8, col = "red")
    for (i in 2: n) {
      points3d(gpa.3D[, , i], size = 8, col = i, add = T)
    }
  }

  # block_out = setClass("block", slots = c(gpa.3D = "array", gpa.2D = "matrix", raw = "array", centroid = "numeric", p = "numeric", k = "numeric", n = "numeric"))
  output = block_out(gpa.3D = gpa.3D, gpa.2D = gpa.2D, raw = r.step4, centroid = centroid, p = dim(gpa.3D)[1], k = dim(gpa.3D)[2], n = dim(gpa.3D)[3])

  return(output)

}
