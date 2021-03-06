% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recurrent_moe.r
\name{predict,RecurrentMoEModel,SingleCellExperimentList-method}
\alias{predict,RecurrentMoEModel,SingleCellExperimentList-method}
\title{predict}
\usage{
\S4method{predict}{RecurrentMoEModel,SingleCellExperimentList}(
  model,
  x,
  batch_size = 256L,
  assay_field = c("X", "X"),
  output = FALSE,
  ...
)
}
\arguments{
\item{model}{a trained model}

\item{x}{a SingleCellExperiment object}

\item{batch_size}{Batch size (default: 256L)}
}
\description{
Predict counts
}
\author{
Wuming Gong (gongx030@umn.edu)
}
