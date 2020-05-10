#' Calculate significance of subtree
get_subtree_zscore <- function(this.node,df,this.sd=bootstrap_baselines(df,200)){
  # this.node is the root of the subtree of interest
  # df is the gene-gene correlation matrix
  this.mean <- mean(df)
  # get leaves in subnode
  leaves <- get_leaves(this.node)
  sign <- mean(df[leaves,leaves])-this.mean
  sign <- 1#sign/abs(sign)
  return(sign*(mean(df[leaves,leaves])-this.mean)^2/this.sd[length(leaves)])
}
