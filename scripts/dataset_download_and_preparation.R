### THIS SCRIPT DOWNLOADS THE REQUIRED DATASETS FOR THE ANALYSIS
### EACH DATASET IS SAVED AS SEPERATE RDATA OBJECTS

## REQUIRED PACKAGES

library(TabulaMurisData)
library(DuoClustering2018)
library(ExperimentHub)
library(SingleCellExperiment)
library(scran)
library(scater)
library(aricode)
library(clValid)


################### ZHENGMIX DATA PREPREPARE ################## 


zhengmix4eq <- sce_full_Zhengmix4eq(metadata = FALSE)
zheng4mixeq.sce.labels <- zhengmix4eq@colData@listData[["phenoid"]]
zhengmix4eq_counts <- zhengmix4eq@assays@data@listData[["counts"]]
rownames(zhengmix4eq_counts) <- zhengmix4eq@rowRanges@elementMetadata@listData[["symbol"]]
zhengmix4eq_sce <- SingleCellExperiment(assays = list(counts = zhengmix4eq_counts))
rm(zhengmix4eq_counts , zhengmix4eq)

zhengmix4uneq <- sce_full_Zhengmix4uneq(metadata = FALSE)
zheng4mixuneq.sce.labels <- zhengmix4uneq@colData@listData[["phenoid"]]
zhengmix4uneq_counts <- zhengmix4uneq@assays@data@listData[["counts"]]
rownames(zhengmix4uneq_counts) <- zhengmix4uneq@rowRanges@elementMetadata@listData[["symbol"]]
zhengmix4uneq_sce <- SingleCellExperiment(assays = list(counts = zhengmix4uneq_counts))
rm(zhengmix4uneq_counts , zhengmix4uneq)

zhengmix8eq <- sce_full_Zhengmix8eq(metadata = FALSE)
zheng8mixeq.sce.labels <- zhengmix8eq@colData@listData[["phenoid"]]
zhengmix8eq_counts <- zhengmix8eq@assays@data@listData[["counts"]]
rownames(zhengmix8eq_counts) <- zhengmix8eq@rowRanges@elementMetadata@listData[["symbol"]]
zhengmix8eq_sce <- SingleCellExperiment(assays = list(counts = zhengmix8eq_counts))
rm(zhengmix8eq_counts , zhengmix8eq)

################### ZHENGMIX DATASETS AND LABELS ###############

labels_zhengmix <- list()

for (i in 1:length(ls(pattern = ".sce.labels"))){
  labels_zhengmix[[i]] <- get(ls(pattern = ".sce.labels")[i])
}

datasets_zhengmix <- list()

for (i in 1:length(ls(pattern = "_sce"))){
  datasets_zhengmix[[i]] <- get(ls(pattern = "_sce")[i])
}

data_names_zhengmix <- ls(pattern = "_sce")

rm(list = ls(pattern = ".sce.labels"))
rm(list = ls(pattern = "_sce") , i)

save(datasets_zhengmix , labels_zhengmix , data_names_zhengmix , file = "zhengmixdata.RData")

rm(list = ls())

#################################################################

################### TABULA MURIS DROPLET ################## 

eh <- ExperimentHub()
tabula_drop<- eh[["EH1617"]]

tabula_drop<- tabula_drop[,!(is.na(tabula_drop$cell_ontology_class))]
num.cells <- nexprs(tabula_drop, byrow=TRUE)
to.keep <- num.cells > 0
tabula_drop<- tabula_drop[to.keep,]

rm(num.cells , to.keep)

tissues <- unique(tabula_drop$tissue)

for (i in 1:length(tissues)) {
  assign(paste(tissues[i],"_droplet",sep = "") ,
         SingleCellExperiment(assays = 
                                list(counts = counts(tabula_drop
                                                     [,tabula_drop$tissue == tissues[i]]))))
  
  assign(paste(tissues[i],"_labels",sep = "") ,
         tabula_drop[,tabula_drop$tissue == tissues[i]]@colData@listData[["cell_ontology_class"]])
}


###################3 TABULA MURIS DROPLET DATASETS AND LABELS ###############

labels <- list()

for (i in 1:length(ls(pattern = "_labels"))){
  labels[[i]] <- get(ls(pattern = "_labels")[i])
}

datasets <- list()

for (i in 1:length(ls(pattern = "_droplet"))){
  datasets[[i]] <- get(ls(pattern = "_droplet")[i])
}

data_names <- ls(pattern = "_droplet")
rm(list = ls(pattern = "_labels"))
rm(list = ls(pattern = "_droplet"))
rm(eh , i , tissues)

datasets_droplet <- datasets
labels_droplet <- labels
data_names_droplet <- data_names

rm(datasets , data_names , labels , tabula_drop)

merged_data_1_drop <- cbind(counts(datasets_droplet[[3]]) , counts(datasets_droplet[[7]]))
merged_data_1_drop_labels <- c(labels_droplet[[3]] , labels_droplet[[7]])

merged_data_2_drop <- cbind(counts(datasets_droplet[[4]]) , counts(datasets_droplet[[6]]))
merged_data_2_drop_labels <- c(labels_droplet[[4]] , labels_droplet[[6]])

merged_data_3_drop <- cbind(counts(datasets_droplet[[9]]) , counts(datasets_droplet[[11]]))
merged_data_3_drop_labels <- c(labels_droplet[[9]] , labels_droplet[[11]])

merged_data_4_drop <- cbind(counts(datasets_droplet[[7]]) , counts(datasets_droplet[[9]]))
merged_data_4_drop_labels <- c(labels_droplet[[7]] , labels_droplet[[9]])

merged_data_1_drop <- SingleCellExperiment(assays = list(counts = merged_data_1_drop))
merged_data_2_drop <- SingleCellExperiment(assays = list(counts = merged_data_2_drop))
merged_data_3_drop <- SingleCellExperiment(assays = list(counts = merged_data_3_drop))
merged_data_4_drop <- SingleCellExperiment(assays = list(counts = merged_data_4_drop))
datasets_droplet[[13]] <- merged_data_1_drop
datasets_droplet[[14]] <- merged_data_2_drop
datasets_droplet[[15]] <- merged_data_3_drop
datasets_droplet[[16]] <- merged_data_4_drop
labels_droplet[[13]] <- merged_data_1_drop_labels
labels_droplet[[14]] <- merged_data_2_drop_labels
labels_droplet[[15]] <- merged_data_3_drop_labels
labels_droplet[[16]] <- merged_data_4_drop_labels
data_names_droplet <- c(data_names_droplet , "merged1" , "merged2" , "merged3" , "merged4")



save(datasets_droplet , labels_droplet , data_names_droplet , file = "dropletdata.RData")

rm(list = ls())

###############################################################################3

################### TABULA MURIS SPIKE ################## 

eh <- ExperimentHub()
tabula_smart <- eh[["EH1618"]]

tabula_smart <- tabula_smart[,!(is.na(tabula_smart$cell_ontology_class))]
num.cells <- nexprs(tabula_smart, byrow=TRUE)
to.keep <- num.cells > 0
tabula_smart <- tabula_smart[to.keep,]

tissues_smart <- unique(tabula_smart$tissue)

for (i in 1:length(tissues_smart)) {
  assign(paste(tissues_smart[i],"_spike",sep = "") ,
         SingleCellExperiment(assays = 
                                list(counts = counts(tabula_smart
                                                     [,tabula_smart$tissue == tissues_smart[i]]))))
  
  assign(paste(tissues_smart[i],".spike.labels",sep = "") ,
         tabula_smart[,tabula_smart$tissue == tissues_smart[i]]@colData@listData[["cell_ontology_class"]])
}



labels <- list()

for (i in 1:length(ls(pattern = ".spike.labels"))){
  labels[[i]] <- get(ls(pattern = ".spike.labels")[i])
}

datasets <- list()

for (i in 1:length(ls(pattern = "_spike"))){
  datasets[[i]] <- get(ls(pattern = "_spike")[i])
}

####################### TABULA MURIS SPIKE DATASET AND LABELS ###################

data_names <- ls(pattern = "_spike")
rm(list = ls(pattern = ".spike.labels"))
rm(list = ls(pattern = "_spike"))
rm(eh , i , tissues_smart ,num.cells , to.keep)

datasets_spike <- datasets
labels_spike <- labels
data_names_spike <- data_names

rm(datasets , data_names , labels , tabula_smart)

merged_data_1_spike <- cbind(counts(datasets_spike[[3]]) , counts(datasets_spike[[7]]))
merged_data_1_spike_labels <- c(labels_spike[[3]] , labels_spike[[7]])
merged_data_2_spike <- cbind(counts(datasets_spike[[4]]) , counts(datasets_spike[[6]]))
merged_data_2_spike_labels <- c(labels_spike[[4]] , labels_spike[[6]])
merged_data_3_spike <- cbind(counts(datasets_spike[[9]]) , counts(datasets_spike[[11]]))
merged_data_3_spike_labels <- c(labels_spike[[9]] , labels_spike[[11]])
merged_data_4_spike <- cbind(counts(datasets_spike[[7]]) , counts(datasets_spike[[9]]))
merged_data_4_spike_labels <- c(labels_spike[[7]] , labels_spike[[9]])
merged_data_1_spike <- SingleCellExperiment(assays = list(counts = merged_data_1_spike))
merged_data_2_spike <- SingleCellExperiment(assays = list(counts = merged_data_2_spike))
merged_data_3_spike <- SingleCellExperiment(assays = list(counts = merged_data_3_spike))
merged_data_4_spike <- SingleCellExperiment(assays = list(counts = merged_data_4_spike))
datasets_spike[[19]] <- merged_data_1_spike
datasets_spike[[20]] <- merged_data_2_spike
datasets_spike[[21]] <- merged_data_3_spike
datasets_spike[[22]] <- merged_data_4_spike
labels_spike[[19]] <- merged_data_1_spike_labels
labels_spike[[20]] <- merged_data_2_spike_labels
labels_spike[[21]] <- merged_data_3_spike_labels
labels_spike[[22]] <- merged_data_4_spike_labels
data_names_spike <- c(data_names_spike , "merged1" , "merged2" , "merged3" , "merged4")

save(datasets_spike , labels_spike , data_names_spike , file = "spikedata.RData")

rm(list = ls())


