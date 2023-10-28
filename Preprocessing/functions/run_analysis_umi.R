run_analysis_umi <- function(sce , known_labels , number_of_genes) {
  sce$known_labels <- known_labels
  stats <- perCellQCMetrics(sce, subsets=list(Mito=which(rownames(sce) == "MT- ")))
  high.sum <- isOutlier(stats$sum, nmads = 3, type="lower")
  high.total <- isOutlier(stats$total, nmads = 3, type="lower")
  high.mito <- isOutlier(stats$subsets_Mito_percent , nmads = 3 , type = "higher")
  sce <- sce[,!(high.sum | high.total | high.mito)]
  clusters <- quickCluster(sce)
  sce <- computeSumFactors(sce, cluster=clusters)
  sce <- logNormCounts(sce)
  dec <- modelGeneVar(sce)
  top <- getTopHVGs(dec , n = number_of_genes)
  sce <- denoisePCA(sce, subset.row=top, technical=dec)
  number_of_pcs <- ncol(reducedDim(sce))
  denoised_pca <- getDenoisedPCs(logcounts(sce) , subset.row = top , technical = dec)
  denoised_pca_perc <- denoised_pca[["percent.var"]]
  explained_var_perc <- sum(denoised_pca_perc[1:number_of_pcs])
  pc.dist <- pc.dist <- dist(reducedDim(sce , "PCA"))
  algorithms <- list(igraph::cluster_walktrap , igraph::cluster_leiden , igraph::cluster_louvain)
  algorihm_names <- c("Walktrap","Leiden","Louvain")
  clusters <- matrix(nrow = 50,ncol = length(algorithms))
  ARI <- matrix(nrow = 50,ncol = length(algorithms))
  for (j in 1:length(algorithms)) {
    for (i in 1:50) {
      clust <- algorithms[[j]](buildSNNGraph(sce, k=i, use.dimred = 'PCA'))$membership
      clusters[i,j] <- max(clust)
      ARI[i,j] <- ARI(clust , sce$known_labels)
    }
  }
  clusters <- as.data.frame(clusters)
  colnames(clusters) <- paste("Clusters",algorihm_names,sep="_")
  ARI <- as.data.frame(ARI)
  colnames(ARI) <- paste("ARI",algorihm_names,sep="_")
  all_results <- cbind(clusters,ARI,rep(number_of_pcs , 50),rep(explained_var_perc,50))
  return(all_results)
}
