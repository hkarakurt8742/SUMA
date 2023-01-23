### THIS SCRIPT USES THE FUNCTIONS (FROM FUNCTIONS FOLDERS) TO DO CLUSTERING
### ON EACH DATASET WITH DIFFERENT HVGs, NEIGHBOUR AND CLUSTERING ALGORITHM

##REQUIRED PACKAGES
library(SingleCellExperiment)
library(scran)
library(scater)
library(aricode)
library(clValid)

##FUNCTION IMPORT
source("functions/run_analysis_umi.R")
source("functions/run_analysis_spike.R")

########################## TABULA MURIS SPIKE #####################


load("spikedata.RData")

#############################################################################

for (i in 1:length(datasets_spike)) {
  assign(paste(data_names_spike[i] , "_500_results" , sep = ""),
         run_analysis_spike(datasets_spike[[i]] , labels_spike[[i]] , 500))
}

results_spike_500 <- list()

for (i in 1:length(ls(pattern = "_500_results"))){
  results_spike_500[[i]] <- get(ls(pattern = "_500_results")[i])
}

rm(list = ls(pattern = "_500_results"))

##########################################################################

for (i in 1:length(datasets_spike)) {
  assign(paste(data_names_spike[i] , "_1000_results" , sep = ""),
         run_analysis_spike(datasets_spike[[i]] , labels_spike[[i]] , 1000))
}

results_spike_1000 <- list()

for (i in 1:length(ls(pattern = "_1000_results"))){
  results_spike_1000[[i]] <- get(ls(pattern = "_1000_results")[i])
}

rm(list = ls(pattern = "_1000_results"))

############################################################################

for (i in 1:length(datasets_spike)) {
  assign(paste(data_names_spike[i] , "_1500_results" , sep = ""),
         run_analysis_spike(datasets_spike[[i]] , labels_spike[[i]] , 1500))
}

results_spike_1500 <- list()

for (i in 1:length(ls(pattern = "_1500_results"))){
  results_spike_1500[[i]] <- get(ls(pattern = "_1500_results")[i])
}

rm(list = ls(pattern = "_1500_results"))


############################################################################

for (i in 1:length(datasets_spike)) {
  assign(paste(data_names_spike[i] , "_2000_results" , sep = ""),
         run_analysis_spike(datasets_spike[[i]] , labels_spike[[i]] , 2000))
}

results_spike_2000 <- list()

for (i in 1:length(ls(pattern = "_2000_results"))){
  results_spike_2000[[i]] <- get(ls(pattern = "_2000_results")[i])
}

rm(list = ls(pattern = "_2000_results"))

save(results_spike_500 , results_spike_1000 , 
     results_spike_1500 , results_spike_2000 , 
     data_names_spike,datasets_spike , file = "spike_results.RData")

rm(list = ls())

#############################################################################

########################## TABULA MURIS DROPLET #####################


load("dropletdata.RData")

for (i in 1:length(datasets_droplet)) {
  assign(paste(data_names_droplet[i] , "_500_results" , sep = ""),
         run_analysis_umi(datasets_droplet[[i]] , labels_droplet[[i]] , 500))
}

results_droplet_500 <- list()

for (i in 1:length(ls(pattern = "_500_results"))){
  results_droplet_500[[i]] <- get(ls(pattern = "_500_results")[i])
}

rm(list = ls(pattern = "_500_results"))

##########################################################################

for (i in 1:length(datasets_droplet)) {
  assign(paste(data_names_droplet[i] , "_1000_results" , sep = ""),
         run_analysis_umi(datasets_droplet[[i]] , labels_droplet[[i]] , 1000))
}

results_droplet_1000 <- list()

for (i in 1:length(ls(pattern = "_1000_results"))){
  results_droplet_1000[[i]] <- get(ls(pattern = "_1000_results")[i])
}

rm(list = ls(pattern = "_1000_results"))

############################################################################

for (i in 1:length(datasets_droplet)) {
  assign(paste(data_names_droplet[i] , "_1500_results" , sep = ""),
         run_analysis_umi(datasets_droplet[[i]] , labels_droplet[[i]] , 1500))
}

results_droplet_1500 <- list()

for (i in 1:length(ls(pattern = "_1500_results"))){
  results_droplet_1500[[i]] <- get(ls(pattern = "_1500_results")[i])
}

rm(list = ls(pattern = "_1500_results"))


##################################################################################

for (i in 1:length(datasets_droplet)) {
  assign(paste(data_names_droplet[i] , "_2000_results" , sep = ""),
         run_analysis_umi(datasets_droplet[[i]] , labels_droplet[[i]] , 2000))
}

results_droplet_2000 <- list()

for (i in 1:length(ls(pattern = "_2000_results"))){
  results_droplet_2000[[i]] <- get(ls(pattern = "_2000_results")[i])
}

rm(list = ls(pattern = "_2000_results"))


save(results_droplet_500 , results_droplet_1000 , 
     results_droplet_1500 , results_droplet_2000 , 
     data_names_droplet,datasets_droplet , file = "droplet_results.RData")

rm(list = ls())
####################################################################################

########################## ZHENGMIX DATASETS #####################


load("zhengmixdata.RData")

for (i in 1:length(datasets_zhengmix)) {
  assign(paste(data_names_zheng[i] , "_500_results" , sep = ""),
         run_analysis_umi(datasets_zhengmix[[i]] , labels_zhengmix[[i]] , 500))
}

results_zhengmix_500 <- list()

for (i in 1:length(ls(pattern = "_500_results"))){
  results_zhengmix_500[[i]] <- get(ls(pattern = "_500_results")[i])
}

rm(list = ls(pattern = "_500_results"))

##########################################################################

for (i in 1:length(datasets_zhengmix)) {
  assign(paste(data_names_zheng[i] , "_1000_results" , sep = ""),
         run_analysis_umi(datasets_zhengmix[[i]] , labels_zhengmix[[i]] , 1000))
}

results_zhengmix_1000 <- list()

for (i in 1:length(ls(pattern = "_1000_results"))){
  results_zhengmix_1000[[i]] <- get(ls(pattern = "_1000_results")[i])
}

rm(list = ls(pattern = "_1000_results"))

##########################################################################

for (i in 1:length(datasets_zhengmix)) {
  assign(paste(data_names_zheng[i] , "_1500_results" , sep = ""),
         run_analysis_umi(datasets_zhengmix[[i]] , labels_zhengmix[[i]] , 1500))
}

results_zhengmix_1500 <- list()

for (i in 1:length(ls(pattern = "_1500_results"))){
  results_zhengmix_1500[[i]] <- get(ls(pattern = "_1500_results")[i])
}

rm(list = ls(pattern = "_1500_results"))

##########################################################################

for (i in 1:length(datasets_zhengmix)) {
  assign(paste(data_names_zheng[i] , "_2000_results" , sep = ""),
         run_analysis_umi(datasets_zhengmix[[i]] , labels_zhengmix[[i]] , 2000))
}

results_zhengmix_2000 <- list()

for (i in 1:length(ls(pattern = "_2000_results"))){
  results_zhengmix_2000[[i]] <- get(ls(pattern = "_2000_results")[i])
}

rm(list = ls(pattern = "_2000_results"))


save(results_zhengmix_500 , results_zhengmix_1000 , 
     results_zhengmix_1500 , results_zhengmix_2000 , 
     data_names_zheng , datasets_zhengmix , file = "zhengmix_results.RData")
