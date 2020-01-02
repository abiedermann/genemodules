library(data.tree)
#' Get leaves in the subtree of a given node.
get_leaves <- function(this.node){
  return(gsub(" ","_",as.character(this.node$Get("name",filterFun=function(x) isLeaf(x)))))
}
