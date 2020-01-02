#' Given a square dataframe, calculate the expected variance of values in a
#' square submatrix <= the size of the input dataframe, and return these values
#' in a vector in order of increasing submatrix size. These values can be used
#' to determine the signficance of correlation for a given square submatrix
#' relative to the baseline expected correlation, calculated from the total
#' dataset.
bootstrap_baselines <- function(df,n_iter){
  # df is a square input dataframe containing gene-gene correlation coefficients
  n_points <- 50
  stdevs <- rep(0,n_points)
  for(i in 1:n_points){
    stdevs[i] <- bootstrap_sd_from_random_submatrices(df,i,n_iter)
  }
  s <- stdevs

  fit <- nls(y~a/(x),data=data.frame(x=1:length(s),y=s[1:length(s)]), start=list(a=1))
  final.stdevs <- predict(fit,list(x=1:dim(df)[1]))
  return(final.stdevs)
}
