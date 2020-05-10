#' Assign cluster numbers
assign_cluster_numbers <- function(this.gm,this.cor,alpha=0.05){
  # this.gm is the gene module output (a list containing lists of genesets, with the z-score used to name each list)
  # alpha is the level of significance (by default 0.05, corresponding to a z-score cutoff of 1.645)
  # Note: All clusters that are not deemed significant are grouped into cluster 0 together.
  # Output: a vector of gene names, named according to cluster
  z.threshold <- qnorm(1-alpha)

  this.genes <- unlist(this.gm)
  #names(this.genes) <- rep(0,length(this.genes))

  cluster.nums <- c()
  for(i in 1:length(this.gm[names(this.gm)>=z.threshold])){
    if(length(this.gm[[i]])==1){
      cluster.nums <- c(cluster.nums,0)
      next
    }
    cluster.nums <- c(cluster.nums,rep(i,length(this.gm[[i]])))
  }

  # Note: this assumes that this.gm is ordered according to decreasing z-score genesets
  cluster.nums <- c(cluster.nums,rep(0,length(this.genes)-length(cluster.nums)))
  names(this.genes) <- cluster.nums

  # Calculating correlation of each cluster num to cluster 1 (the most
  # significant cluster)
  initial.nums <- unique(names(this.genes))
  initial.nums <- initial.nums[which(!(initial.nums %in% c(0)))]
  final.nums <- rep(0,length(initial.nums))
  for(i in 1:length(initial.nums)){
      final.nums[i] <- mean(this.cor[this.genes[names(this.genes)==1],this.genes[names(this.genes)==initial.nums[i]]])
  }
  
  # Renaming clusters
  this.names <- names(this.genes)
  new.names <- this.names
  for(i in 1:length(initial.nums)){
      new.names[this.names==initial.nums[order(final.nums,decreasing=T)[i]]] = i
  }
  names(this.genes) <- new.names

  this.genes <- this.genes[order(names(this.genes))]

  return(this.genes)
}
