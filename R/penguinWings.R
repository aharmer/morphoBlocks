#' @title Penguin partial wing skeleton landmarks dataset
#'
#' @description The multi-part objects in this example are partial wing skeletons comprised of humerus, radius and ulna. These partial wing skeletons are from 15 extant species of penguin and five fossil species of penguin, which together constitute a dataset of 60 wing bones. Shape data from the multi-part objects are provided as landmark configurations. Three sets of landmarks were produced for each of the digital replicas so that the replicates could be averaged to mitigate effects of placement error. The nine sets of landmark configurations (i.e. three sets from each of the three bones) were separately read into \code{R} using \code{readPts} with \code{gpa = FALSE} (i.e. generalised Procrustes transformation not performed). Metadata for each bone are available on github (https://github.com/aharmer/morphoBlocks).
#'
#' @docType data
#' @keywords datasets
#' @name penguinWings
#' @usage data(penguinWings)
#' @format A set of 'block' class objects.
"hum1"
#' @rdname penguinWings
"hum2"
#' @rdname penguinWings
"hum3"
#' @rdname penguinWings
"rad1"
#' @rdname penguinWings
"rad2"
#' @rdname penguinWings
"rad3"
#' @rdname penguinWings
"uln1"
#' @rdname penguinWings
"uln2"
#' @rdname penguinWings
"uln3"
