load("final_results_all.RData")

all_results_droplet$Exp_Type <- rep("Droplet" , 7200)
all_results_spike$Exp_Type <- rep("Spike" , 10800)
all_results_zhengmix$Exp_Type <- rep("Droplet" , 1800)

all_results <- rbind(all_results_droplet , all_results_spike , all_results_zhengmix)

all_results$acceptable_ari <- all_results$ARI >= 0.65
all_results$ARI <- NULL
all_results$acceptable_ari <- as.numeric(all_results$acceptable_ari)

library(caret)

dummy <- dummyVars(" ~ .", data=all_results)
encoded_results <- data.frame(predict(dummy, newdata=all_results))
encoded_results$Exp_TypeSpike <- NULL

write.table(encoded_results , file = "encoded_dataset.txt" , quote = F , sep = "\t",
            col.names = T , row.names = F)
