#' Given a dataframe, calculate the standard deviation of values in
#' a len x len submatrix, sampled randomly with replacement from the larger
#' dataframe n_iter times.
bootstrap_sd_from_random_submatrices <- function(df,len,n_iter){
  boot.means <- rowMeans(matrix(sample(df,len*len*n_iter,replace=T),n_iter,len*len))
  return(sd(boot.means))
}
