#' Plots gene modules
plot_gene_modules <- function(gene.clusters,gene.cors,fontsize=12,plot_non_sig=T){
    library(ComplexHeatmap)
    library(circlize)
    n=0
    if(!plot_non_sig){
	    n=1
    }
    col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    idx <- gene.clusters[names(gene.clusters)>=n]
    names(idx) <- as.factor(as.numeric(names(idx)))
    print(Heatmap(gene.cors[idx,idx],cluster_rows=F,cluster_columns=F,
        column_split = as.factor(as.numeric(names(idx))), 
        row_split=as.factor(as.numeric(names(idx))),
	row_names_gp = gpar(fontsize = fontsize),column_names_gp = gpar(fontsize = fontsize)))

}
