result_extractions_w <- function(results_data) {
  number_of_pcs_used <- results_data$`rep(number_of_pcs, 50)`[1]
  highest_ami <- max(results_data$AMI_Walktrap)
  highest_silhouette <- max(results_data$Silhouette_Walktrap)
  best_k_with_highest_ami <- which.max(results_data$AMI_Walktrap)
  best_k_with_highest_silhouette <- which.max(results_data$Silhouette_Walktrap)
  number_of_cluster_with_highest_ami <- results_data$Clusters_Walktrap[best_k_with_highest_ami]
  number_of_cluster_with_highest_silhouette <- results_data$Clusters_Walktrap[best_k_with_highest_silhouette]
  comb_results <- (cbind(number_of_pcs_used , highest_ami , highest_silhouette , best_k_with_highest_ami , 
                         best_k_with_highest_silhouette , number_of_cluster_with_highest_ami,
                         number_of_cluster_with_highest_silhouette))
  print(as.numeric(comb_results))
}


result_extractions_louv <- function(results_data) {
  number_of_pcs_used <- results_data$`rep(number_of_pcs, 50)`[1]
  highest_ami <- max(results_data$AMI_Louvain)
  highest_silhouette <- max(results_data$Silhouette_Louvain)
  best_k_with_highest_ami <- which.max(results_data$AMI_Louvain)
  best_k_with_highest_silhouette <- which.max(results_data$Silhouette_Louvain)
  number_of_cluster_with_highest_ami <- results_data$Clusters_Louvain[best_k_with_highest_ami]
  number_of_cluster_with_highest_silhouette <- results_data$Clusters_Louvain[best_k_with_highest_silhouette]
  comb_results <- (cbind(number_of_pcs_used , highest_ami ,highest_silhouette, best_k_with_highest_ami , 
                         best_k_with_highest_silhouette , number_of_cluster_with_highest_ami,
                         number_of_cluster_with_highest_silhouette))
  print(as.numeric(comb_results))
}


result_extractions_leid <- function(results_data) {
  number_of_pcs_used <- results_data$`rep(number_of_pcs, 50)`[1]
  highest_ami <- max(results_data$AMI_Leiden)
  highest_silhouette <- max(results_data$Silhouette_Leiden)
  best_k_with_highest_ami <- which.max(results_data$AMI_Leiden)
  best_k_with_highest_silhouette <- which.max(results_data$Silhouette_Leiden)
  number_of_cluster_with_highest_ami <- results_data$Clusters_Leiden[best_k_with_highest_ami]
  number_of_cluster_with_highest_silhouette <- results_data$Clusters_Leiden[best_k_with_highest_silhouette]
  comb_results <- (cbind(number_of_pcs_used , highest_ami ,highest_silhouette, best_k_with_highest_ami , 
                         best_k_with_highest_silhouette , number_of_cluster_with_highest_ami,
                         number_of_cluster_with_highest_silhouette))
  print(as.numeric(comb_results))
}

result_extractions_info <- function(results_data) {
  number_of_pcs_used <- results_data$`rep(number_of_pcs, 50)`[1]
  highest_ami <- max(results_data$AMI_Infomap)
  highest_silhouette <- max(results_data$Silhouette_Infomap)
  best_k_with_highest_ami <- which.max(results_data$AMI_Infomap)
  best_k_with_highest_silhouette <- which.max(results_data$Silhouette_Infomap)
  number_of_cluster_with_highest_ami <- results_data$Clusters_Infomap[best_k_with_highest_ami]
  number_of_cluster_with_highest_silhouette <- results_data$Clusters_Infomap[best_k_with_highest_silhouette]
  comb_results <- (cbind(number_of_pcs_used , highest_ami ,highest_silhouette, best_k_with_highest_ami , 
                         best_k_with_highest_silhouette , number_of_cluster_with_highest_ami,
                         number_of_cluster_with_highest_silhouette))
  print(as.numeric(comb_results))
}