#' Hilbert-Schmidt Independence Criteria (HSIC)
#'
#' Compute HSIC loss beteween two matrices
#'
#' @param x input data
#' @param y input data
#' @return HSIC loss
#'
#' @references Tsai, Y.-H. H., Bai, S., Morency, L.-P. & Salakhutdinov, R.
#' A Note on Connecting Barlow Twins with Negative-Sample-Free Contrastive Learning. Arxiv (2021).
#'
#' @references  Weinberger, E. & Lee, S.-I. Transferable representations of single-cell transcriptomic data. Biorxiv 2021.04.13.439707 (2021) doi:10.1101/2021.04.13.439707.
#'
hsic_loss <- function(x, y){

  n <- x$shape[[1]]
  h <- tf$linalg$diag(tf$ones(n)) - tf$ones(shape(n, n)) / n  # the centering matrix

  kx <- tf$matmul(x, x, transpose_b = TRUE)
  kx <- tf$matmul(kx, h)

  ky <- tf$matmul(y, y, transpose_b = TRUE)
  ky <- tf$matmul(ky, h)

  loss <- tf$matmul(kx, ky) %>% tf$linalg$trace()
  loss <- loss / n^2
  loss
}
