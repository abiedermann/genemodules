#' Gene module discovery algorithm (fast implementation)
get_gene_modules <- function(this.matrix,this.tree=as.phylo(hclust(parDist(this.matrix)))){
  library(data.tree)
  library(ape)
# this.matrix is a gene correlation matrix
  
  # Note: the fast implementation avoids recalculating all significance scores on each iterations
  # and only updates parents of the trimmed node
  
  # Accounting for the fact that data.trees seems to replace _ with spaces
  mod.list <- list()
  mod.zscores <- c()
  
  this.tree <- as.Node(this.tree)
  this.sd <- bootstrap_baselines(this.matrix,200)

  test.out <- this.tree$Get(function(x) get_subtree_zscore(x,this.matrix,this.sd),filterFun=function(x) isNotLeaf(x))
  
  while(TRUE){
    # Find node to delete
    node.name <- FindNode(this.tree,names(which(test.out==max(test.out))))
    
    # Isolate module genes, add their names and z-scores to lists
    this.mod <- get_leaves(node.name)
    
    mod.list <- append(mod.list,list(this.mod))
    mod.zscores <- c(mod.zscores,test.out[node.name$name])
    
    # get parent node of the node we are about to prune
    this.parent <- node.name$parent
    
    n.del <- Prune(this.tree,pruneFun=function(x) x$name!=node.name$name)
    
    # break the loop when the tree has no more leaves
    if((this.tree$leafCount==0) | (n.del==0)){
      break
    }
    
    # drop leaves which are no longer in test.out
    test.out <- test.out[names(this.tree$Get("name",filterFun=function(x) !isLeaf(x)))]
    #test.out <- test.out[get_leaves(this.tree)]
    
    # Finally, updating the root node
    while(!isRoot(this.parent)){
      # update subtree score
      test.out[this.parent$name] <- get_subtree_zscore(this.parent,this.matrix,this.sd)
      this.parent <- this.parent$parent
    }
    test.out[this.parent$name] <- get_subtree_zscore(this.parent,this.matrix,this.sd)
    
    # print remaining number of leaves
    print(this.tree$leafCount)
  }
  
  names(mod.list) <- as.numeric(mod.zscores)
  
  return(mod.list)
}
