load("~/R/only_10x/results_2_10x.RData") # Results Files

walktrap_max_500 <- rbind(results_500[[1]][which.max(results_500[[1]]$AMI_Walktrap),],
                          results_500[[2]][which.max(results_500[[2]]$AMI_Walktrap),])

for (i in 3:15) {
  walktrap_max_500 <- rbind(walktrap_max_500,results_500[[i]][which.max(results_500[[i]]$AMI_Walktrap),])
}

walktrap_max_500$Number_of_Neighbours <- c()

for (i in 1:15) {
  walktrap_max_500$Number_of_Neighbours[i] <- which.max(results_500[[i]]$AMI_Walktrap)
}


##############################
  
walktrap_max_1000 <- rbind(results_1000[[1]][which.max(results_1000[[1]]$AMI_Walktrap),],
                          results_1000[[2]][which.max(results_1000[[2]]$AMI_Walktrap),])

for (i in 3:15) {
  walktrap_max_1000 <- rbind(walktrap_max_1000,results_1000[[i]][which.max(results_1000[[i]]$AMI_Walktrap),])
}

walktrap_max_1000$Number_of_Neighbours <- c()

for (i in 1:15) {
  walktrap_max_1000$Number_of_Neighbours[i] <- which.max(results_1000[[i]]$AMI_Walktrap)
}

###############################

walktrap_max_2000 <- rbind(results_2000[[1]][which.max(results_2000[[1]]$AMI_Walktrap),],
                           results_2000[[2]][which.max(results_2000[[2]]$AMI_Walktrap),])

for (i in 3:15) {
  walktrap_max_2000 <- rbind(walktrap_max_2000,results_2000[[i]][which.max(results_2000[[i]]$AMI_Walktrap),])
}


walktrap_max_2000$Number_of_Neighbours <- c()

for (i in 1:15) {
  walktrap_max_2000$Number_of_Neighbours[i] <- which.max(results_2000[[i]]$AMI_Walktrap)
}

###############################

Louvain_max_500 <- rbind(results_500[[1]][which.max(results_500[[1]]$AMI_Louvain),],
                         results_500[[2]][which.max(results_500[[2]]$AMI_Louvain),])

for (i in 3:15) {
  Louvain_max_500 <- rbind(Louvain_max_500,results_500[[i]][which.max(results_500[[i]]$AMI_Louvain),])
}

Louvain_max_500$Number_of_Neighbours <- c()

for (i in 1:15) {
  Louvain_max_500$Number_of_Neighbours[i] <- which.max(results_500[[i]]$AMI_Louvain)
}

#####################################

Louvain_max_1000 <- rbind(results_1000[[1]][which.max(results_1000[[1]]$AMI_Louvain),],
                          results_1000[[2]][which.max(results_1000[[2]]$AMI_Louvain),])


for (i in 3:15) {
  Louvain_max_1000 <- rbind(Louvain_max_1000,results_1000[[i]][which.max(results_1000[[i]]$AMI_Louvain),])
}

Louvain_max_1000$Number_of_Neighbours <- c()

for (i in 1:15) {
  Louvain_max_1000$Number_of_Neighbours[i] <- which.max(results_1000[[i]]$AMI_Louvain)
}

#########################################

Louvain_max_2000 <- rbind(results_2000[[1]][which.max(results_2000[[1]]$AMI_Louvain),],
                          results_2000[[2]][which.max(results_2000[[2]]$AMI_Louvain),])

for (i in 3:15) {
  Louvain_max_2000 <- rbind(Louvain_max_2000,results_2000[[i]][which.max(results_2000[[i]]$AMI_Louvain),])
}

Louvain_max_2000$Number_of_Neighbours <- c()

for (i in 1:15) {
  Louvain_max_2000$Number_of_Neighbours[i] <- which.max(results_2000[[i]]$AMI_Louvain)
}

##########################################

Leiden_max_500 <- rbind(results_500[[1]][which.max(results_500[[1]]$AMI_Leiden),],
                        results_500[[2]][which.max(results_500[[2]]$AMI_Leiden),])

for (i in 3:15) {
  Leiden_max_500 <- rbind(Leiden_max_500,results_500[[i]][which.max(results_500[[i]]$AMI_Leiden),])
}

Leiden_max_500$Number_of_Neighbours <- c()

for (i in 1:15) {
  Leiden_max_500$Number_of_Neighbours[i] <- which.max(results_500[[i]]$AMI_Leiden)
}

########################################

Leiden_max_1000 <- rbind(results_1000[[1]][which.max(results_1000[[1]]$AMI_Leiden),],
                         results_1000[[2]][which.max(results_1000[[2]]$AMI_Leiden),])

for (i in 3:15) {
  Leiden_max_1000 <- rbind(Leiden_max_1000,results_1000[[i]][which.max(results_1000[[i]]$AMI_Leiden),])
}

Leiden_max_1000$Number_of_Neighbours <- c()

for (i in 1:15) {
  Leiden_max_1000$Number_of_Neighbours[i] <- which.max(results_1000[[i]]$AMI_Leiden)
}

##########################################

Leiden_max_2000 <- rbind(results_2000[[1]][which.max(results_2000[[1]]$AMI_Leiden),],
                         results_2000[[2]][which.max(results_2000[[2]]$AMI_Leiden),])

for (i in 3:15) {
  Leiden_max_2000 <- rbind(Leiden_max_2000,results_2000[[i]][which.max(results_2000[[i]]$AMI_Leiden),])
}

Leiden_max_2000$Number_of_Neighbours <- c()

for (i in 1:15) {
  Leiden_max_2000$Number_of_Neighbours[i] <- which.max(results_2000[[i]]$AMI_Leiden)
}

##############################################

cells_in_datasets <- c()
for (i in 1:15) {
  cells_in_datasets[i] <- ncol(datasets[[i]])
}


max_values_walktrap <- as.data.frame(rbind(cbind(walktrap_max_500$Number_of_Neighbours,walktrap_max_500$`rep(number_of_pcs, 50)`,
                             cells_in_datasets,rep(500,15),walktrap_max_500$AMI_Walktrap),
                             cbind(walktrap_max_1000$Number_of_Neighbours,walktrap_max_1000$`rep(number_of_pcs, 50)`,
                             cells_in_datasets,rep(1000,15),walktrap_max_1000$AMI_Walktrap),
                             cbind(walktrap_max_2000$Number_of_Neighbours,walktrap_max_2000$`rep(number_of_pcs, 50)`,
                             cells_in_datasets,rep(2000,15),walktrap_max_2000$AMI_Walktrap)))

max_values_walktrap$Algorithm <- rep("Walktrap" , 45)


##########################################################

max_values_louvain <- as.data.frame(rbind(cbind(Louvain_max_500$Number_of_Neighbours,Louvain_max_500$`rep(number_of_pcs, 50)`,
                                                 cells_in_datasets,rep(500,15),Louvain_max_500$AMI_Louvain),
                                           cbind(Louvain_max_1000$Number_of_Neighbours,Louvain_max_1000$`rep(number_of_pcs, 50)`,
                                                 cells_in_datasets,rep(1000,15),Louvain_max_1000$AMI_Louvain),
                                           cbind(Louvain_max_2000$Number_of_Neighbours,Louvain_max_2000$`rep(number_of_pcs, 50)`,
                                                 cells_in_datasets,rep(2000,15),Louvain_max_2000$AMI_Louvain)))

max_values_louvain$Algorithm <- rep("Louvain" , 45)


#######################################################

max_values_leiden <- as.data.frame(rbind(cbind(Leiden_max_500$Number_of_Neighbours,Leiden_max_500$`rep(number_of_pcs, 50)`,
                                                cells_in_datasets,rep(500,15),Leiden_max_500$AMI_Leiden),
                                          cbind(Leiden_max_500$Number_of_Neighbours,Leiden_max_1000$`rep(number_of_pcs, 50)`,
                                                cells_in_datasets,rep(1000,15),Leiden_max_1000$AMI_Leiden),
                                          cbind(Leiden_max_500$Number_of_Neighbours,Leiden_max_2000$`rep(number_of_pcs, 50)`,
                                                cells_in_datasets,rep(2000,15),Leiden_max_2000$AMI_Leiden)))

max_values_leiden$Algorithm <- rep("Leiden" , 45)


#############################################


colnames <-  c("Number_of_Neighbours" , "Number_of_PCs" , "Number_of_Cells" ,
               "Number_of_HVGs" , "AMI" , "Algorithm")

colnames(max_values_walktrap) <- colnames
colnames(max_values_louvain) <- colnames
colnames(max_values_leiden) <- colnames

max_values <- rbind(max_values_walktrap , max_values_louvain , max_values_leiden)

library(caret)

max_values$Number_of_PCs <- as.numeric(max_values$Number_of_PCs)
max_values$Number_of_Cells <- as.numeric(max_values$Number_of_Cells)
max_values$Number_of_HVGs <- as.numeric(max_values$Number_of_HVGs)
max_values$AMI <- as.numeric(max_values$AMI)

dummy <- dummyVars(" ~ .", data=max_values)
oh_results <- data.frame(predict(dummy, newdata=max_values))

write.table(oh_results , file = "oh_encoded_results.csv" , sep = "," , quote = F , 
            row.names = F , col.names = T)

oh_results_max_values <- oh_results

save(oh_results_max_values , file = "oh_encoded_max_results.RData")
