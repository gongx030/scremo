setClassUnion('listOrNULL', members = c('list', 'NULL'))
setClassUnion('matrixOrNULL', members = c('matrix', 'NULL'))
setOldClass('kerastools.model.RModel')
setOldClass('tf_dataset')
setOldClass('tensorflow.tensor')
setOldClass('tensorflow.python.framework.sparse_tensor.SparseTensor')
setOldClass('tensorflow.python.training.tracking.data_structures.ListWrapper')

#' Model
#'
setClass('Model', slot = c(model = 'kerastools.model.RModel'))

#' SingleCellExperimentList
#' 
#' @export
#'
setClass(
	'SingleCellExperimentList',
	contains = 'SimpleList'
)

