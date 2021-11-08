#' umap
#'
#' @export
umap <- NULL
np <- NULL

.onLoad <- function(libname, pkgname) {
	umap <<- reticulate::import("umap", delay_load = list(
		priority = 5
#		environment = "r411"
	))
	np <<- reticulate::import("numpy", delay_load = list(
		priority = 5
#		environment = "r411"
	))
}
