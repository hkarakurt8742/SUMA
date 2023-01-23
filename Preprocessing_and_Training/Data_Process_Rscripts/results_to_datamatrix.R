# THIS SCRIPT MANIPULATES THE RESULTS
# RESULTS MUST BE MANIPULATED ONE BY ONE (EACH DATASET GROUP) AND SAVED


load("zhengmix_results.RData")
load("droplet_results.RData")
load("spike_results.RData")


sc_clust_data_manipulator <- function(algorithm_name , number_of_variants , datasets , results) {
  cells_in_datasets <- c()
  for (i in 1:length(datasets)) {
    cells_in_datasets[i] <- ncol(datasets[[i]])
  }
  cells_in_datasets <- rep(cells_in_datasets , each = 50)
  number_of_neighbours <- c(1:50)
  number_of_neighbours <- rep(number_of_neighbours ,length(datasets))
  number_of_pcs <- c(results[[1]][["rep(number_of_pcs, 50)"]],
                     results[[2]][["rep(number_of_pcs, 50)"]])
  for (i in 3:length(datasets)) {
    number_of_pcs <- c(number_of_pcs , results[[i]][["rep(number_of_pcs, 50)"]])
  }
  number_of_variant_genes <- rep(number_of_variants , (length(datasets) * 50))
  aris <- c(get(paste("results_spike_",number_of_variants,sep=""))[[1]][[paste("ARI",algorithm_name,sep = "_")]],
            get(paste("results_spike_",number_of_variants,sep=""))[[2]][[paste("ARI",algorithm_name,sep = "_")]])
  
  for (i in 3:length(datasets)) {
    aris <- c(aris , get(paste("results_spike_",number_of_variants,sep=""))[[i]][[paste("ARI",algorithm_name,sep = "_")]])
  }
  
  cluster_results <- as.data.frame(cbind(number_of_neighbours, number_of_pcs, cells_in_datasets,
                                         number_of_variant_genes , aris , rep(algorithm_name , (length(datasets) * 50))))
  cluster_results
}

algorithms <- c("Walktrap" , "Louvain" , "Leiden")
number_of_variants <- c(500 , 1000 , 1500 , 2000)


for (i in 1:3) {
  for (j in 1:4) {
    assign(paste(algorithms[i],number_of_variants[j],sep="_") , 
           sc_clust_data_manipulator(algorithms[i] , number_of_variants[j] , datasets_spike , results_spike_500))
  }
}

all_results <- rbind(Walktrap_500 , Walktrap_1000 , Walktrap_1500 , Walktrap_2000,
                     Louvain_500 , Louvain_1000 , Louvain_1500 , Louvain_2000,
                     Leiden_500 , Leiden_1000, Leiden_1500 , Leiden_2000)


all_results$number_of_neighbours <- as.numeric(all_results$number_of_neighbours)
all_results$number_of_pcs <- as.numeric(all_results$number_of_pcs)
all_results$cells_in_datasets <- as.numeric(all_results$cells_in_datasets)
all_results$number_of_variant_genes <- as.numeric(all_results$number_of_variant_genes)
all_results$aris <- round(as.numeric(all_results$aris) , digits = 5)

colnames(all_results) <- c("Number_of_Neighbours","Number_of_PCs",
                           "Number_of_Cells","Number_of_HVGs",
                           "ARI","Algorithm")

