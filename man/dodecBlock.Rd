% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dodecBlock.R
\name{dodecBlock}
\alias{dodecBlock}
\title{Simulate a data block of landmark configurations}
\usage{
dodecBlock(
  n = 10,
  dist_x = NULL,
  dist_y = NULL,
  dist_z = NULL,
  theta_x = NULL,
  theta_y = NULL,
  theta_z = NULL,
  size = NULL,
  noise = NULL,
  vertex_shift = NULL,
  plot = FALSE
)
}
\arguments{
\item{n}{the number of configurations to generate for the data block. Default is 10. Must be 2 or greater.}

\item{dist_x, dist_y, dist_z}{optional integer or vector to specify a distance along the stated axis (i.e. \code{dist_x} for the x axis) for translating the configurations. If no value is specified then configurations will not be translated. If a single integer is supplied (e.g. \code{dist_x = 1}) then all configurations will be translated along an axis by the same amount. If a vector of length \emph{n} is supplied (e.g. \code{dist_x = 1:10}) then each configuration will be translated by a value corresponding to its position in the vector.}

\item{theta_x, theta_y, theta_z}{optional integer or vector to specify an angle for rotating the configuration about its origin in the stated axis (i.e. \code{theta_x} for the x axis). If no value is specified then configurations will not be rotated. If a single integer is supplied (e.g. \code{theta_x = 30}) then all configurations will be rotated by the same amount. If a vector of length \emph{n} is supplied (e.g. \code{theta_x = seq(from = 0, to = 210, by = 30)}) then each configuration will be rotated by a value corresponding to its position in the vector.}

\item{size}{optional integer or vector for scaling the configurations. If no value is specified then configurations will not be scaled relative to one another. If a single integer is supplied (e.g. \code{size = 1}) then all configurations will be scaled by the same amount. If a vector of length \emph{n} is supplied (e.g. \code{size = 1:10}) then each configuration will be scaled by a value corresponding to its position in the vector.}

\item{noise}{optional integer for introducing noise to the configurations. If no value is specified then \code{runif(1, min = 0, max = (1 / 100))} is added to each of the \code{x}, \code{y} and \code{z} positions of each point in each configuration. This small amount of noise is added so that no configuration has identical coordinates. The amount of noise that is introduced to the data block can be increased by specifying a noise value that is greater than 1.}

\item{vertex_shift}{optional integer or vector to specify a distance along the z axis for translating the first vertex in each configuration. This argument can be used to introduce variation into the data block. If no value is specified then vertex 1 in each configuration will not be translated. If a single integer is supplied (e.g. \code{vertex_shift = 1}) then vertex 1 in all configurations will be translated along the z axis by the same amount. If a vector of length \emph{n} is supplied (e.g. \code{vertex_shift = 1:10}) then vertex 1 in each configuration will be translated by a value corresponding to its position in the vector. Things get really fun when you use both \code{vertex_shift} and \code{theta_z} together (see example 2 below).}

\item{plot}{a logical value indicating whether the configurations should be plotted. Default is FALSE.}
}
\value{
a 'blockList' object of \emph{n} dodecahedra, used for downstream analyses. The list contains the elements:

\item{raw}{configurations without Procrustes transformation}

\item{gpa.3D}{configurations after Procrustes transformation organised into a 3D array}

\item{gpa.2D}{configurations after Procrustes transformation organised into a 2D matrix}

\item{centroid}{centroid sizes of the configurations}

\item{p}{number of points in the configurations of each data block}

\item{k}{number of dimensions that the points in each configuration has}

\item{n}{number of configurations included in each data block}
}
\description{
Generates a block of \emph{n} dodecahedra where the vertices of each dodecahedron represent a landmark configuration with \emph{p = 20} points in \emph{k = 3} dimensions. The data block is transformed using generalized Procrustes analysis.
}
\details{
\code{dodecBlock} builds a data block of configurations that have been transformed using generalized Procrustes analysis. The resulting data block is a three dimensional array with \emph{n = 2} or greater configurations, with each configuration represented by \emph{p = 20} points that have \emph{k = 3} dimensions. Several objects are calculated for downstream analyses including centroid sizes for each configuration.
\code{formatBlock} uses the \code{cSize} function from the \code{Morpho} package (Schlager 2017), and the \code{gpagen}, \code{two.d.array} and \code{arrayspecs} functions from the \code{geomorph} package (Adams and Ot?rola-Castillo 2013).
}
\examples{
# Example 1: Generate a block with ten configurations
block <- dodecBlock()

# Example 2: Generate a block with five configurations.
# Rotate the configurations and translate the first vertex in each configuration.
# Plot the result.
block2 <- dodecBlock(n = 5, theta_z = seq(from = 0, to = 270, by = 60), vertex_shift = 0:4, plot = TRUE)

}
\references{
Adams DC, Ot?rola-Castillo E. 2013. geomorph: an \R package for the collection and analysis of geometric morphometric shape data. Methods in Ecology and Evolution 4:393-399 https://doi.org/10.1111/2041-210X.12035

Schlager S. 2017. Morpho and Rvcg-shape analysis in \R. In Zheng G, Li S, Sz?kely (eds.) Statistical shape and deformation analysis. Academic Press, London. Pp. 217-256.
}
