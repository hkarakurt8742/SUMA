library(SingleCellExperiment)
library(randomForest)
library(Seurat)
library(scran)
library(ConsensusClusterPlus)
library(clValid)
library(igraph)
library(easystats)

count_matrix <- read.table("C:/Users/pc/Documents/R/studio/shiny/sc_clust/counts_matrix_droplet.csv" , header = T , row.names = 1 , sep = ",")
sce <- SingleCellExperiment(assays = list(counts=count_matrix))
mito_gene <- "mt- "
number_of_genes <- 2000
rf_model <- readRDS("rf_model2.RDS")

stats <- perCellQCMetrics(sce, subsets=list(Mito=which(rownames(sce) == mito_gene)))
high.sum <- isOutlier(stats$sum, nmads = 3, type="both")
high.total <- isOutlier(stats$total, nmads = 3, type="both")
high.mito <- isOutlier(stats$subsets_Mito_percent , nmads = 3 , type = "higher")
discard <- high.sum | high.total | high.mito
sce$discard <- discard
sce <- sce[,!(high.sum | high.total | high.mito)]
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clusters)
sce <- logNormCounts(sce)
dec <- modelGeneVar(sce)
top <- getTopHVGs(dec , n = number_of_genes)
sce <- denoisePCA(sce, subset.row=top, technical=dec)
number_of_pcs <- ncol(reducedDim(sce))
### CLUSTERING VIA MODEL ###
denoised_pca <- getDenoisedPCs(logcounts(sce) , subset.row = top , technical = dec)
denoised_pca_perc <- denoised_pca[["percent.var"]]
explained_var_perc <- sum(denoised_pca_perc[1:number_of_pcs])
algorithms <- c(rep(1,50),rep(2,50),rep(3,50))
model_data <- as.data.frame(cbind(rep(ncol(sce), 150),
                                  rep(1:50,3),
                                  rep(number_of_pcs,150),
                                  rep(number_of_genes,150),
                                  algorithms,
                                  rep(explained_var_perc,150),
                                  rep(1,150)))
colnames(model_data) <- c("Number_of_Cells" , "Number_of_Neighbours" , 
                          "Number_of_PCs" , "Number_of_HVGs" , "Algorithm" , 
                          "var_perc" , "Exp_Type")
model_data$Number_of_Cells <- as.integer(model_data$Number_of_Cells)
model_data$Number_of_Neighbours <- as.integer(model_data$Number_of_Neighbours)
model_data$Number_of_PCs <- as.integer(model_data$Number_of_PCs)
model_data$Number_of_HVGs <- as.integer(model_data$Number_of_HVGs)
model_data$Algorithm <- as.factor(model_data$Algorithm)
model_data$Exp_Type <- as.integer(model_data$Exp_Type)
rf_pred <- predict(rf_model , newdata = model_data)
selected_k <- as.numeric(names(sort(rf_pred , decreasing = TRUE)))[1]
if (model_data[as.numeric(names(sort(rf_pred , decreasing = TRUE)))[1],]$Algorithm == 1) {
  model_clusters <- cluster_leiden(buildSNNGraph(sce , k = selected_k))$membership
}
if (model_data[as.numeric(names(sort(rf_pred , decreasing = TRUE)))[1],]$Algorithm == 2) {
  model_clusters <- cluster_louvain(buildSNNGraph(sce , k = selected_k))$membership
}
if (model_data[as.numeric(names(sort(rf_pred , decreasing = TRUE)))[1],]$Algorithm == 3) {
  model_clusters <- cluster_walktrap(buildSNNGraph(sce , k = selected_k))$membership
}

model_clusters <- as.factor(model_clusters)

### CLUSTERING VIA CONSENSUSCLUSTERING ###
pca_res <- as.matrix(reducedDim(sce , "PCA"))
cc_res <- ConsensusClusterPlus(as.matrix(counts(sce)) ,  maxK = 25)
Kvec = 2:25
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
for(i in Kvec){
  M = cc_res[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}#end for i
# The optimal K
optK = Kvec[which.min(PAC)]
consensus_clusters <- as.factor(cc_res[[optK]]$consensusClass)

### Clustering Based on Dunn Index (Louvain) ###
dunn_indexes_louv <- as.numeric()

for (i in 2:50) {
  dunn_indexes_louv[i] <- dunn(clusters = cluster_louvain(buildSNNGraph(sce , k = i))$membership
                          , Data = pca_res)
}

dunn_indexes_louv[is.na(dunn_indexes_louv)] <- 0
names(dunn_indexes_louv) <- c(1:50)
dunn_selected_k <- as.numeric(names(sort(dunn_indexes_louv , decreasing = TRUE))[1])
dunn_clusters_louv <- cluster_louvain(buildSNNGraph(sce , k = dunn_selected_k))$membership
dunn_clusters_louv <- as.factor(dunn_clusters_louv)


### Clustering Based on Dunn Index (Louvain) ###
dunn_indexes_leid <- as.numeric()

for (i in 2:50) {
  dunn_indexes_leid[i] <- dunn(clusters = cluster_leiden(buildSNNGraph(sce , k = i))$membership
                               , Data = pca_res)
}

dunn_indexes_leid[is.na(dunn_indexes_leid)] <- 0
names(dunn_indexes_leid) <- c(1:50)
dunn_selected_k <- as.numeric(names(sort(dunn_indexes_leid , decreasing = TRUE))[1])
dunn_clusters_leid <- cluster_leiden(buildSNNGraph(sce , k = dunn_selected_k))$membership
dunn_clusters_leid <- as.factor(dunn_clusters_leid)

######################### SEURAT CLUSTERING ###############################

sce.seurat <- as.Seurat(sce , counts = "counts" , data = "logcounts")
sce.seurat <- FindVariableFeatures(sce.seurat, 
                                   selection.method = "vst", nfeatures = number_of_genes)
sce.seurat <- ScaleData(sce.seurat)
sce.seurat <- RunPCA(sce.seurat , npcs = number_of_pcs)
sce.seurat <- RunUMAP(sce.seurat , dims = 1:number_of_pcs)
sce.seurat <- FindNeighbors(sce.seurat, dims = 1:number_of_pcs)
sce.seurat <- FindClusters(sce.seurat)

##################### MERGING CLUSTERS ####################

seurat_clusters <- sce.seurat@active.ident
sce.seurat@meta.data$consensus_clusters <- consensus_clusters
sce.seurat@meta.data$dunn_clusters_louv <- dunn_clusters_louv
sce.seurat@meta.data$dunn_clusters_leid <- dunn_clusters_leid
sce.seurat@meta.data$model_clusters <- model_clusters


##################################################

UMAPPlot(sce.seurat , group.by = "consensus_clusters")
UMAPPlot(sce.seurat , group.by = "dunn_clusters_louv")
UMAPPlot(sce.seurat , group.by = "dunn_clusters_leid")
UMAPPlot(sce.seurat , group.by = "model_clusters")
UMAPPlot(sce.seurat)

cluster_results <- list(seurat_clusters , consensus_clusters , 
                        dunn_clusters_louv , dunn_clusters_leid , model_clusters)

meta_clusters <- cluster_meta(cluster_results)

heatmap(meta_clusters)

sce.seurat@meta.data$model_clusters <- meta_clusters_2

my_group <- as.numeric(as.factor(substr(rownames(meta_clusters), 1 , 1)))
colSide <- brewer.pal(9, "Set1")[my_group]
colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
heatmap(meta_clusters, scale="column", 
        col= cm.colors(256) , main="Meta-Clustering Heatmap",RowSideColors=colSide)



