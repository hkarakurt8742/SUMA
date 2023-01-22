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

###################3 ZHENGMIX DATASETS AND LABELS ###############

labels_zhengmix <- list()

for (i in 1:length(ls(pattern = ".sce.labels"))){
  labels_zhengmix[[i]] <- get(ls(pattern = ".sce.labels")[i])
}

datasets_zhengmix <- list()

for (i in 1:length(ls(pattern = "_sce"))){
  datasets_zhengmix[[i]] <- get(ls(pattern = "_sce")[i])
}

data_names_zheng <- ls(pattern = "_sce")

rm(list = ls(pattern = ".sce.labels"))
rm(list = ls(pattern = "_sce") , i)

save(datasets_zhengmix , labels_zhengmix , data_names_zheng , file = "zhengmixdata.RData")

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

save(datasets_droplet , labels_droplet , data_names_droplet , file = "dropletdata.RData")

rm(list = ls())

###################################################################


###################### SINCELL DATA PREPROCESS ########################


load("~/Documents/R/sc_clustering/analyses/sincell_with_class.RData")

sce.10x.qc.labels.sincell <- colData(sce_sc_10x_qc)["cell_line"][,1]
sce_10x_qc_counts <- sce_sc_10x_qc@assays[["data"]]@listData[["counts"]]
sce_10x_qc <- SingleCellExperiment(assays = list(counts = sce_10x_qc_counts))
rm(sce_10x_qc_counts , sce_sc_10x_qc)

sce.celseq.qc.labels.sincell <- colData(sce_sc_CELseq2_qc)["cell_line"][,1]
sce_celseq_qc_counts <- sce_sc_CELseq2_qc@assays[["data"]]@listData[["counts"]]
sce_celseq_qc <- SingleCellExperiment(assays = list(counts = sce_celseq_qc_counts))
rm(sce_celseq_qc_counts , sce_sc_CELseq2_qc)

sce.dropseq.qc.labels.sincell <- colData(sce_sc_Dropseq_qc)["cell_line"][,1]
sce_dropseq_qc_counts <- sce_sc_Dropseq_qc@assays[["data"]]@listData[["counts"]]
sce_dropseq_qc <- SingleCellExperiment(assays = list(counts = sce_dropseq_qc_counts))
rm(sce_dropseq_qc_counts , sce_sc_Dropseq_qc)


################### SINCELL DATASETS AND LABELS ###############

labels <- list()

for (i in 1:length(ls(pattern = ".labels.sincell"))){
  labels[[i]] <- get(ls(pattern = ".labels.sincell")[i])
}

datasets <- list()

for (i in 1:length(ls(pattern = "sce_"))){
  datasets[[i]] <- get(ls(pattern = "sce_")[i])
}

data_names <- ls(pattern = "sce_")
rm(list = ls(pattern = ".labels.sincell"))
rm(list = ls(pattern = "sce_"))
rm(i)

datasets_sincell <- datasets
labels_sincell <- labels
data_names_sincell <- data_names

rm(datasets , data_names , labels)

save(datasets_sincell , labels_sincell , data_names_sincell , file = "sincelldata.RData")

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

rm(datasets , data_names , labels , tabula_drop)

save(datasets_spike , labels_spike , data_names_spike , file = "spikedata.RData")

rm(list = ls())


