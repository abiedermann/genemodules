#' Assign cluster numbers
get_module_interaction_df <- function(cluster.nums,this.cor){
  # cluster.nums is the assign cluster numbers output.
  # Note: All clusters that are not deemed significant are grouped into cluster 0 together.
  x <- max(as.numeric(names(lu))) # length of interaction matrix
  imat <- matrix(0,x+1,x+1) # empty x by x matrix (adjusting for 0 index)

  for(i in 0:x){
      for(j in 0:x){
          n_idx <- which(rownames(this.cor) %in% 
                         cluster.nums[which(names(cluster.nums)==i)])
          m_idx <- which(rownames(this.cor) %in% 
                         cluster.nums[which(names(cluster.nums)==j)])
          imat[i+1,j+1] = mean(this.cor[n_idx,m_idx])
      }
  }
  imat <- data.frame(imat)
  rownames(imat) <- 0:x
  colnames(imat) <- 0:x
  return(imat)
}



