#' RecurrentMoEModel
#'
setClass('RecurrentMoEModel', contains = 'Model')


#' RecurrentMoEModel
#' 
#' A batch VAE model for scRNA-seq
#'
#' @export
#'
RecurrentMoEModel <- function(
	latent_dim = 10L,
	n_batches = 1000L,
	modality = NULL,
	n_genes = NULL,
	hidden1 = 128L,
	hidden2 = 256L,
	normalize = FALSE,
	rate = 0.1,
	iter = 2L,
	mask = FALSE,
	batch_modality = FALSE,
	n_cell_types = 1L,
	name = NULL
){

	keras_model_custom(name = name, function(self){

		self$latent_dim <- latent_dim
		self$n_batches <- n_batches
		self$modality <- modality
		self$n_genes <- n_genes
		self$n_modality <- length(modality)
		self$iter <- iter
		self$mask <- mask
		self$batch_modality <- batch_modality
		self$n_cell_types <- n_cell_types

		self$dense_input <- lapply(1:length(modality), function(i){
			tf$keras$layers$Dense(hidden2)
		})

		self$encoder <- tf$keras$Sequential(list(
			tf$keras$layers$BatchNormalization(),
			tf$keras$layers$Dropout(rate),
			tf$keras$layers$Activation('relu'),
			tf$keras$layers$Dense(hidden1),
			tf$keras$layers$BatchNormalization(),
			tf$keras$layers$Dropout(rate),
			tf$keras$layers$Activation('relu')
		))

		self$representation <- tf$keras$layers$Dense(latent_dim)

		if (batch_modality){
			self$embedding_1 <- tf$keras$layers$Embedding(n_batches, self$latent_dim)
			self$embedding_2 <- tf$keras$layers$Embedding(n_batches, self$latent_dim)
		}else{
			self$embedding <- tf$keras$layers$Embedding(n_batches, self$latent_dim)
		}

		self$mod1 <- tf$keras$layers$Embedding(2L, self$latent_dim)
		self$mod2 <- tf$keras$layers$Embedding(2L, self$latent_dim)

		self$decoder <- tf$keras$Sequential(list(
			tf$keras$layers$Dense(hidden1),
			tf$keras$layers$BatchNormalization(),
			tf$keras$layers$Dropout(rate),
			tf$keras$layers$Activation('relu'),
			tf$keras$layers$Dense(hidden2),
			tf$keras$layers$BatchNormalization(),
			tf$keras$layers$Dropout(rate),
			tf$keras$layers$Activation('relu')
		))

		self$dense_output <- lapply(1:length(modality), function(i){
			tf$keras$layers$Dense(n_genes[[i]], activation = 'relu')
		})	

		if (self$n_cell_types > 1L){
			self$cell_type_class <- tf$keras$layers$Dense(n_cell_types)
		}

		function(x, ..., training = TRUE){

			if (self$batch_modality){
				b1 <- x$batch %>% self$embedding_1()
				b2 <- x$batch %>% self$embedding_2()
			}else{
				b1 <- x$batch %>% self$embedding()
				b2 <- x$batch %>% self$embedding()
			}

			m1 <- tf$zeros_like(x$batch) %>% self$mod1()
			m2 <- tf$zeros_like(x$batch) %>% self$mod2()

			x1 <- x$x1
			x2 <- x$x2

			if (self$mask && training){
				w <- tf$random$uniform(shape(x1$shape[[1]]), maxval = self$n_modality, dtype = tf$int64) %>%
					tf$one_hot(self$n_modality)
				maskable <- x$w %>% tf$reduce_sum(axis = 1L, keepdims = TRUE)
				maskable <- tf$cast(maskable > 1L, tf$float32)	# maskable samples  (e.g. with data from two mods)
				w <- x$w * (1 - maskable) + w * maskable
			}else{
				w <- x$w
			}

			for (i in 1:self$iter){

				z1 <- x1 %>%
					self$dense_input[[0L]]() %>%
					self$encoder() %>%
					self$representation()

				z2 <- x2 %>%
					self$dense_input[[1L]]() %>%
					self$encoder() %>%
					self$representation()

				if (normalize){
					z1 <- z1 %>% tf$math$l2_normalize(-1L)
					z2 <- z2 %>% tf$math$l2_normalize(-1L)
				}

				# use the mod with data
				z <- (z1 * w[, 1L, drop = FALSE] + z2 * w[, 2L, drop = FALSE]) / tf$reduce_sum(w, 1L, keepdims = TRUE)

				x1_pred <- (z + b1 + m1) %>%
					self$decoder() %>%
					self$dense_output[[0L]]()

				x2_pred <- (z + b2 + m2) %>%
					self$decoder() %>%
					self$dense_output[[1L]]()

				if (i < self$iter){
					x1 <- x1 * w[, 1L, drop = FALSE] + x1_pred * (1 - w[, 1L, drop = FALSE])
					x2 <- x2 * w[, 2L, drop = FALSE] + x2_pred * (1 - w[, 2L, drop = FALSE])
				}
			}

			if (self$n_cell_types > 1L){
				cell_type <- z %>% self$cell_type_class()
			}else{
				cell_type <- tf$zeros(x$batch$shape[[1L]])
			}

			list(
				z = z,
				z1 = z1,
				z2 = z2,
				x1 = x1_pred,
				x2 = x2_pred,
				cell_type = cell_type
			)
		}
	})
}

#' prepare_data
#'
#' Prepare dataset for training a model
#' 
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'prepare_data',
	signature(
		model = 'RecurrentMoEModel',
		x = 'SingleCellExperimentList'
	),
	function(
		model,
		x,
		batch_field = c('batch', 'batch'),
		assay_field = c('X', 'X'),
		weight = NULL,
		cell_type = NULL,
		...
	){

		stopifnot(length(x) == model@model$n_modality)
		stopifnot(length(batch_field) == length(x))

		for (i in 1:length(x)){
			stopifnot(!is.null(assays(x[[i]])[[assay_field[i]]]))
			stopifnot(!is.null(colData(x[[i]])[[batch_field[i]]]))
		}

		stopifnot(length(unique(sapply(x, ncol))) == 1)

		if (!is.null(cell_type)){
			for (i in 1:length(x)){
				stopifnot(!is.null(colData(x[[i]])[[cell_type]]))
			}
			for (i in 2:length(x) ){
				stopifnot(all(colData(x[[i - 1]])[[cell_type]] == colData(x[[i]])[[cell_type]]))
			}
		}

		if (is.null(weight)){
			weight <- matrix(TRUE, nrow = ncol(x[[1]]), ncol = length(x))
		}

		d <- lapply(1:length(x), function(i){
			assays(x[[i]])[[assay_field[[i]]]] %>%
				t() %>%
				as.matrix() %>% 
				tf$cast(tf$float32) 
		})
		names(d) <- sprintf('x%d', 1:length(x))

		batch <- colData(x[[1]])[[batch_field[1]]] %>% as.numeric()
		batch[is.na(batch)] <- 0L
		batch <- batch + 1L
		batch <- batch %>% tf$cast(tf$int64)

		w <- lapply(1:length(x), function(i) weight[, i] %>% as.numeric() %>% tf$cast(tf$float32) %>% tf$expand_dims(1L)) %>%
			tf$concat(axis = 1L)

		if (!is.null(cell_type)){
			ct <- colData(x[[1]])[[cell_type]] %>% as.numeric()  %>% tf$cast(tf$int64)
			ct <- ct - 1L
		}else{
			ct <- rep(0, ncol(x[[1]])) %>% tf$cast(tf$int64)
		}

		c(d, w = w, batch = batch, cell_type = ct)

	}
)

#' fit
#'
#' Fit a VaeModel
#'
#' @param model a VaeModel object, initialized by `new('VaeModel', model = VaeModel(...))`
#' @param x a tf_dataset object
#' @param batch_size Batch size (default: 256L)
#' @param epochs Number of training epochs (default: 500L)
#' @param learning_rate Learning rate (default: 1e-4)
#' @param beta a weight for KL loss (default: 5e-5)
#' @param compile Whether or not compile the tensorflow model (default: TRUE)
#' @export
#'
#' @return a VaeModel
#'
setMethod(
	'fit',
	signature(
		model = 'RecurrentMoEModel',
		x = 'tf_dataset'
	),
	function(
		model,
		x,
		batch_size = 128L,
		epochs = 100L,
		learning_rate = 1e-3,
		match_loss_type = 'mse',
		compile = TRUE,
		training = TRUE
	){

		x <- x %>%
			dataset_batch(batch_size)

		optimizer <- tf$keras$optimizers$Adam(learning_rate)
		mse_loss <- tf$keras$losses$MeanSquaredError(reduction = 'none')
		cce_loss <- tf$keras$losses$SparseCategoricalCrossentropy(from_logits = TRUE)

		train_step <- function(batch){
			n <- batch$x1$shape[[1]]
			with(tf$GradientTape(persistent = TRUE) %as% tape, {

				res <- model@model(batch, training = training)

				loss_hsic_1 <- (hsic_loss(res$z1, batch$batch %>% tf$one_hot(model@model$n_batches)))
				loss_hsic_2 <- (hsic_loss(res$z2, batch$batch %>% tf$one_hot(model@model$n_batches)))

				loss_1 <- (tf$reduce_sum(batch$w[, 1L] * mse_loss(batch$x1, res$x1)) / (1 + tf$reduce_sum(batch$w[, 1L])))
				loss_2 <- (tf$reduce_sum(batch$w[, 2L] * mse_loss(batch$x2, res$x2)) / (1 + tf$reduce_sum(batch$w[, 2L])))

				if (match_loss_type == 'contrastive'){
					temperature <- 1
					labels <- tf$one_hot(tf$range(n), n * 2L)
					logits_aa <- tf$matmul(res$z1, res$z1, transpose_b = TRUE) / temperature
					logits_bb <- tf$matmul(res$z2, res$z2, transpose_b = TRUE) / temperature
					logits_ab <- tf$matmul(res$z1, res$z2, transpose_b = TRUE) / temperature
					logits_ba <- tf$matmul(res$z2, res$z1, transpose_b = TRUE) / temperature
					loss_a <- tf$nn$softmax_cross_entropy_with_logits(labels, tf$concat(list(logits_ab, logits_aa), 1L))
					loss_b <- tf$nn$softmax_cross_entropy_with_logits(labels, tf$concat(list(logits_ba, logits_bb), 1L))
					loss_match <- tf$reduce_mean(loss_a + loss_b) / n
				}else if (match_loss_type == 'mse'){
					loss_match <- (mse_loss(res$z1, res$z2) %>% tf$reduce_mean())
				}else if (match_loss_type == 'task2'){
					s <- tf$matmul(res$z1, res$z2, transpose_b = TRUE)
					labels <- tf$one_hot(tf$range(n), n)
					loss_match <- tf$nn$softmax_cross_entropy_with_logits(labels, s) %>% tf$reduce_mean()
					loss_match 
				}

				if (model@model$n_cell_types > 1){
					loss_class <- cce_loss(batch$cell_type, res$cell_type)
				}else{
					loss_class <- 0
				}

				loss <- loss_1 + loss_2 + loss_hsic_1 + loss_hsic_2 + loss_match + loss_class
			})

			gradients <- tape$gradient(loss, model@model$trainable_variables)
			list(gradients, model@model$trainable_variables) %>%
				purrr::transpose() %>%
				optimizer$apply_gradients()
			list(
				loss_1 = loss_1,
				loss_2 = loss_2,
				loss_hsic_1 = loss_hsic_1,
				loss_hsic_2 = loss_hsic_2,
				loss_match = loss_match,
				loss_class = loss_class,
				loss = loss
			)
		}

		if (compile){
			train_step <- tf_function(train_step) # convert to graph mode
		}

		for (epoch in seq_len(epochs)){
			loss <- NULL
			iter <- x %>%
				dataset_shuffle(1000L) %>%
				make_iterator_one_shot()
			res <- until_out_of_range({
				batch <- iterator_get_next(iter)
				res <- train_step(batch)
				loss <- rbind(loss, sapply(res, as.numeric))
			})

			loss <- colMeans(loss)
			sprintf('epoch=%6.d/%6.d | %s', epoch, epochs, paste(sapply(1:length(loss), function(i) sprintf('%s=%13.7f', names(loss)[i], loss[i])), collapse = ' | ')) %>%
				message()
		}
		model
	}
)


#' predict 
#'
#' Predict counts
#'
#' @param model a trained model
#' @param x a SingleCellExperiment object
#' @param batch_size Batch size (default: 256L)
#'
#' @export
#' @author Wuming Gong (gongx030@umn.edu)
#'
setMethod(
	'predict',
	signature(
		model = 'RecurrentMoEModel',
		x = 'SingleCellExperimentList'
	),
	function(
		model,
		x,
		batch_size = 256L, 
		assay_field = c('X', 'X'),
		output = FALSE,
		...
	){

		d <- model %>%
			prepare_data(x, assay_field = assay_field, ...) %>%
			tensor_slices_dataset() %>%
			dataset_batch(batch_size)

		iter <- d %>% make_iterator_one_shot()

		latent_combined <- NULL
		latent_1 <- NULL
		latent_2 <- NULL

		if (output){
			x1 <- NULL
			x2 <- NULL
		}

		res <- until_out_of_range({
			batch <- iterator_get_next(iter)
			res <- model@model(batch, training = FALSE)
			latent_combined <- c(latent_combined, res$z)
			latent_1 <- c(latent_1, res$z1)
			latent_2 <- c(latent_2, res$z2)
			if (output){
				x1 <- c(x1, res$x1)
				x2 <- c(x2, res$x2)
			}
		})

		latent_combined <- latent_combined %>% tf$concat(axis = 0L)
		latent_1 <- latent_1 %>% tf$concat(axis = 0L)
		latent_2 <- latent_2 %>% tf$concat(axis = 0L)
		reducedDim(x[[1]], 'latent') <- as.matrix(latent_1)
		reducedDim(x[[2]], 'latent') <- as.matrix(latent_2)
		reducedDim(x[[1]], 'latent_combined') <- as.matrix(latent_combined)
		reducedDim(x[[2]], 'latent_combined') <- as.matrix(latent_combined)

		if (output){
			x1 <- x1 %>% tf$concat(axis = 0L)
			x2 <- x2 %>% tf$concat(axis = 0L)
			assay_field <- sprintf('%s_predicted', assay_field)
			assays(x[[1]], withDimnames = FALSE)[[assay_field[1]]] <- as.matrix(x1) %>% as('dgCMatrix') %>% t()
			assays(x[[2]], withDimnames = FALSE)[[assay_field[2]]] <- as.matrix(x2) %>% as('dgCMatrix') %>% t()
		}

		x
	}
)	
