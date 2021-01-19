#' Runs standard gene module analysis workflow
#'
#' Given a correlation matrix (or any square matrix), this function returns
#' cluster numbering suitable for plotting with ComplexHeatmap, as well as
#' interaction scores, clustered.
run_gene_module_analysis <- function(this.cor,alpha=0.05,cluster_abs=F,apply_fischer=F,
                                     fontsize=12,plot_non_sig=T,
                                     this.tree=as.phylo(hclust(parDist(this.cor)))){
    library(genemodules)
    library(loveseq)
    this.gm <- get_gene_modules(this.cor,this.tree)
    clust.nums <- assign_cluster_numbers(this.gm,this.cor,alpha=alpha)
    plot_gene_modules(clust.nums,this.cor,fontsize=fontsize,plot_non_sig=plot_non_sig)
    imat <- get_module_interaction_df(clust.nums,this.cor)
    abs.imat <- cluster_columns_and_rows(abs(imat))
    imat <- imat[rownames(abs.imat),colnames(abs.imat)]

    return(list(clust.nums,imat))
}


