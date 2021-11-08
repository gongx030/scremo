#' Combine multiple SingleCellExperiment objects into a SingleCellExperimentList object
#'
#' @exportMethod coerce
#'
setAs('list', "SingleCellExperimentList", function(from) {
	out <- as(from, 'SimpleList')
	new('SingleCellExperimentList', out)
})
