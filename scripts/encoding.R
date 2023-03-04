load("zhengmix_final_results.RData")
load("droplet_final_results.RData")
load("spike_final_results.RData")

all_results_droplet$Exp_Type <- rep("Droplet" , 14400)
all_results_spike$Exp_Type <- rep("Spike" , 19800)
all_results_zhengmix$Exp_Type <- rep("Droplet" , 2700)


all_results <- rbind(all_results_droplet , all_results_spike , all_results_zhengmix)
all_results$Algorithm <- as.factor(all_results$Algorithm)
all_results$Exp_Type <- as.factor(all_results$Exp_Type)

library(caret)

dummy <- dummyVars(" ~ .", data=all_results)
all_results <- data.frame(predict(dummy, newdata=all_results))

all_results$ARI <- round(all_results$ARI , digits = 3)

all_results$acceptable_ari <- all_results$ARI >= 0.65
all_results$ARI <- NULL
all_results$Exp_Type.Spike <- NULL
all_results$acceptable_ari <- as.integer(as.logical(all_results$acceptable_ari))

library(corrplot)
correlations <- cor(all_results[,1:8])
corrplot(correlations, method="circle")

write.table(all_results , file = "oh_encoded_classification_data_65.txt" , sep = "\t" , col.names = T,
            row.names = F , quote = F)

