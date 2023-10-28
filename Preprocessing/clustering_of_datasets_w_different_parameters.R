### THIS SCRIPT USES THE FUNCTIONS (FROM FUNCTIONS FOLDERS) TO DO CLUSTERING
### ON EACH DATASET WITH DIFFERENT HVGs, NEIGHBOUR AND CLUSTERING ALGORITHM

##REQUIRED PACKAGES
library(SingleCellExperiment)
library(scran)
library(scater)
library(aricode)
library(clValid)
library(foreach)
library(doParallel)

##FUNCTION IMPORT
source("functions/run_analysis_umi.R")
source("functions/run_analysis_spike.R")

########################## TABULA MURIS SPIKE #####################


load("spikedata.RData")
packages <- c('SingleCellExperiment' , 'scran' , 'scater' ,
              'aricode' , 'clValid')

#############################################################################

cl <- makeCluster(6)
registerDoParallel(cl)

results_spike_500 <- list()


results_spike_500 <- foreach(i = 1:length(datasets_spike) , .packages=packages) %dopar% {
  run_analysis_spike(datasets_spike[[i]] , labels_spike[[i]] , 500)
}

stopCluster(cl)


gc()

##########################################################################

cl <- makeCluster(6)
registerDoParallel(cl)


results_spike_1000 <- list()


results_spike_1000 <- foreach(i = 1:length(datasets_spike) , .packages=packages) %dopar% {
  run_analysis_spike(datasets_spike[[i]] , labels_spike[[i]] , 1000)
}

stopCluster(cl)


gc()

############################################################################

cl <- makeCluster(8)
registerDoParallel(cl)


results_spike_1500 <- list()


results_spike_1500 <- foreach(i = 1:length(datasets_spike) , .packages=packages) %dopar% {
  run_analysis_spike(datasets_spike[[i]] , labels_spike[[i]] , 1500)
}

stopCluster(cl)


gc()

############################################################################

cl <- makeCluster(8)
registerDoParallel(cl)


results_spike_2000 <- list()


results_spike_2000 <- foreach(i = 1:length(datasets_spike) , .packages=packages) %dopar% {
  run_analysis_spike(datasets_spike[[i]] , labels_spike[[i]] , 2000)
}

stopCluster(cl)


gc()

##############################################################################

cl <- makeCluster(8)
registerDoParallel(cl)


results_spike_2500 <- list()


results_spike_2500 <- foreach(i = 1:length(datasets_spike) , .packages=packages) %dopar% {
  run_analysis_spike(datasets_spike[[i]] , labels_spike[[i]] , 2500)
}

stopCluster(cl)


gc()

#############################################################################

cl <- makeCluster(8)
registerDoParallel(cl)


results_spike_3000 <- list()


results_spike_3000 <- foreach(i = 1:length(datasets_spike) , .packages=packages) %dopar% {
  run_analysis_spike(datasets_spike[[i]] , labels_spike[[i]] , 3000)
}

stopCluster(cl)


gc()


save(results_spike_500 , results_spike_1000 , 
     results_spike_1500 , results_spike_2000 , 
     results_spike_2500 , results_spike_3000 ,
     data_names_spike , datasets_spike , file = "spike_results.RData")



#############################################################################

########################## TABULA MURIS DROPLET #####################


load("dropletdata.RData")
packages <- c('SingleCellExperiment' , 'scran' , 'scater' ,
              'aricode' , 'clValid')

cl <- makeCluster(8)
registerDoParallel(cl)


results_spike_500 <- list()


results_spike_500 <- foreach(i = 1:length(datasets_droplet) , .packages=packages) %dopar% {
  run_analysis_spike(datasets_droplet[[i]] , labels_droplet[[i]] , 500)
}

stopCluster(cl)


gc()

##########################################################################

cl <- makeCluster(8)
registerDoParallel(cl)


results_spike_1000 <- list()


results_spike_1000 <- foreach(i = 1:length(datasets_droplet) , .packages=packages) %dopar% {
  run_analysis_spike(datasets_droplet[[i]] , labels_droplet[[i]] , 1000)
}

stopCluster(cl)


gc()

############################################################################

cl <- makeCluster(8)
registerDoParallel(cl)


results_spike_1500 <- list()


results_spike_1500 <- foreach(i = 1:length(datasets_droplet) , .packages=packages) %dopar% {
  run_analysis_spike(datasets_droplet[[i]] , labels_droplet[[i]] , 1500)
}

stopCluster(cl)


gc()

############################################################################

cl <- makeCluster(8)
registerDoParallel(cl)


results_spike_2000 <- list()


results_spike_2000 <- foreach(i = 1:length(datasets_droplet) , .packages=packages) %dopar% {
  run_analysis_spike(datasets_droplet[[i]] , labels_droplet[[i]] , 2000)
}

stopCluster(cl)


gc()

##############################################################################

cl <- makeCluster(8)
registerDoParallel(cl)


results_spike_2500 <- list()


results_spike_2500 <- foreach(i = 1:length(datasets_droplet) , .packages=packages) %dopar% {
  run_analysis_spike(datasets_droplet[[i]] , labels_droplet[[i]] , 2500)
}

stopCluster(cl)


gc()

#############################################################################

cl <- makeCluster(8)
registerDoParallel(cl)


results_spike_3000 <- list()


results_spike_3000 <- foreach(i = 1:length(datasets_droplet) , .packages=packages) %dopar% {
  run_analysis_spike(datasets_droplet[[i]] , labels_droplet[[i]] , 3000)
}

stopCluster(cl)


gc()


save(results_droplet_500 , results_droplet_1000 , 
     results_droplet_1500 , results_droplet_2000 , 
     results_droplet_2500 , results_droplet_3000 ,
     data_names_droplet , datasets_droplet , file = "droplet_results.RData")

####################################################################################

########################## ZHENGMIX DATASETS #####################


load("zhengmixdata.RData")

cl <- makeCluster(8)
registerDoParallel(cl)


results_spike_500 <- list()


results_spike_500 <- foreach(i = 1:length(datasets_zhengmix) , .packages=packages) %dopar% {
  run_analysis_spike(datasets_zhengmix[[i]] , labels_zhengmix[[i]] , 500)
}

stopCluster(cl)


gc()

##########################################################################

cl <- makeCluster(8)
registerDoParallel(cl)


results_spike_1000 <- list()


results_spike_1000 <- foreach(i = 1:length(datasets_zhengmix) , .packages=packages) %dopar% {
  run_analysis_spike(datasets_zhengmix[[i]] , labels_zhengmix[[i]] , 1000)
}

stopCluster(cl)


gc()

############################################################################

cl <- makeCluster(8)
registerDoParallel(cl)


results_spike_1500 <- list()


results_spike_1500 <- foreach(i = 1:length(datasets_zhengmix) , .packages=packages) %dopar% {
  run_analysis_spike(datasets_zhengmix[[i]] , labels_zhengmix[[i]] , 1500)
}

stopCluster(cl)


gc()

############################################################################

cl <- makeCluster(8)
registerDoParallel(cl)


results_spike_2000 <- list()


results_spike_2000 <- foreach(i = 1:length(datasets_zhengmix) , .packages=packages) %dopar% {
  run_analysis_spike(datasets_zhengmix[[i]] , labels_zhengmix[[i]] , 2000)
}

stopCluster(cl)


gc()

##############################################################################

cl <- makeCluster(8)
registerDoParallel(cl)


results_spike_2500 <- list()


results_spike_2500 <- foreach(i = 1:length(datasets_zhengmix) , .packages=packages) %dopar% {
  run_analysis_spike(datasets_zhengmix[[i]] , labels_zhengmix[[i]] , 2500)
}

stopCluster(cl)


gc()

#############################################################################

cl <- makeCluster(8)
registerDoParallel(cl)


results_spike_3000 <- list()


results_spike_3000 <- foreach(i = 1:length(datasets_zhengmix) , .packages=packages) %dopar% {
  run_analysis_spike(datasets_zhengmix[[i]] , labels_zhengmix[[i]] , 3000)
}

stopCluster(cl)


gc()

save(results_zhengmix_500 , results_zhengmix_1000 , 
     results_zhengmix_1500 , results_zhengmix_2000 , 
     results_zhengmix_2500 , results_zhengmix_3000 ,
     data_names_zhengmix , datasets_zhengmix , file = "zhengmix_results.RData")
