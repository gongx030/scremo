% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recurrent_moe.r
\name{RecurrentMoEModel}
\alias{RecurrentMoEModel}
\title{RecurrentMoEModel}
\usage{
RecurrentMoEModel(
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
)
}
\description{
A batch VAE model for scRNA-seq
}
