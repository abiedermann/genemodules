#' Converts Pearson correlation coefficients into Z-scores. Sets diagonal
#' entries to zero to avoid Inf values. Note: This function should only be
#' applied to square correlation matrices.
apply_fischer_transform <- function(this.matrix){
  # transforming from correlation matrix to z-score matrix with Fischer transform
  for(i in 1:length(this.matrix[,1])){
      if(i==1){
          this.matrix[i,i] <- mean(this.matrix[i+1,i],
                                   this.matrix[i,i+1])
          next
      }
      if(i==length(this.matrix[,1])){
          this.matrix[i,i] <- mean(this.matrix[i-1,i],
                                   this.matrix[i,i-1])
          next
      }
      this.matrix[i,i] <- mean(this.matrix[i-1,i],
                               this.matrix[i+1,i],
                               this.matrix[i,i-1],
                               this.matrix[i,i+1])
  }
  this.matrix <- atanh(this.matrix)
  return(this.matrix)
}
