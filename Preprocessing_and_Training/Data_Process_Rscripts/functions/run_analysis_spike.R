run_analysis_spike <- function(sce , known_labels , number_of_genes) {
  sce$known_labels <- known_labels
  is.spike <- grepl("^ERCC", rownames(sce))
  sce <- splitAltExps(sce, ifelse(is.spike, "ERCC", "gene"))
  stats <- perCellQCMetrics(sce)
  high.sum <- isOutlier(stats$sum, nmads = 3, type="lower")
  high.total <- isOutlier(stats$total, nmads = 3, type="lower")
  spike.drop <- isOutlier(stats$altexps_ERCC_percent, nmads=3, type="higher")
  sce <- sce[,!(high.sum | high.total | spike.drop)]
  num.cells <- nexprs(sce, byrow=TRUE)
  to.keep <- num.cells > 0
  sce <- sce[to.keep,]
  clusters <- quickCluster(sce)
  sce <- computeSumFactors(sce, cluster=clusters)
  sce <- logNormCounts(sce)
  dec <- modelGeneVarWithSpikes(sce , "ERCC")
  top <- getTopHVGs(dec ,  n = number_of_genes)
  sce <- denoisePCA(sce, subset.row=top, technical=dec)
  number_of_pcs <- ncol(reducedDim(sce))
  pc.dist <- pc.dist <- dist(reducedDim(sce , "PCA"))
  algorithms <- list(igraph::cluster_walktrap , igraph::cluster_leiden , igraph::cluster_louvain)
  algorihm_names <- c("Walktrap","Leiden","Louvain")
  clusters <- matrix(nrow = 50,ncol = length(algorithms))
  AMI <- matrix(nrow = 50,ncol = length(algorithms))
  ARI <- matrix(nrow = 50,ncol = length(algorithms))
  for (j in 1:length(algorithms)) {
    for (i in 1:50) {
      clust <- algorithms[[j]](buildSNNGraph(sce, k=i, use.dimred = 'PCA'))$membership
      clusters[i,j] <- max(clust)
      AMI[i,j] <- AMI(clust , sce$known_labels)
      ARI[i,j] <- ARI(clust , sce$known_labels)
    }
  }
  clusters <- as.data.frame(clusters)
  colnames(clusters) <- paste("Clusters",algorihm_names,sep="_")
  AMI <- as.data.frame(AMI)
  colnames(AMI) <- paste("AMI",algorihm_names,sep="_")
  ARI <- as.data.frame(ARI)
  colnames(ARI) <- paste("ARI",algorihm_names,sep="_")
  all_results <- cbind(clusters,AMI,ARI,rep(number_of_pcs , 50))
  return(all_results)
}
