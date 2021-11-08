#' prepare_data
#'
#' The generic function of prepare_data
#'
#' @param model a model object
#' @param x a data object 
#' @param ... Other arguments
#'
setGeneric('prepare_data', function(model, x, ...) standardGeneric('prepare_data'))

#' fit
#'
#' The generic function of fit
#'
#' @param model a model object
#' @param x a data object 
#' @param ... Other arguments
#'
setGeneric('fit', function(model, x, ...) standardGeneric('fit'))

#' predict
#'
#' The generic function of predict
#'
#' @param model a model object
#' @param x a data object that
#' @param ... Other arguments
#'
setGeneric('predict', function(model, x, ...) standardGeneric('predict'))
