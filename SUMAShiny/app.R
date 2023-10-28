library(shiny)
library(DT)
library(SingleCellExperiment)
library(randomForest)
library(Seurat)
library(scran)
library(ConsensusClusterPlus)
library(clValid)
library(igraph)
library(easystats)
library(gridExtra)
library(shinybusy)
library(gridExtra)
library(scater)
library(ggplot2)
library(RColorBrewer)
options(shiny.maxRequestSize=100*1024^2)
# Load other necessary libraries if required.

load("suma_r.RData")

# UI Section
ui <- fluidPage(
  tags$style(type="text/css", "
             body {
               overflow-y: scroll;
             }
             "),
  titlePanel("SUMAShiny: A Multi-Algorithm scRNA-Seq Clustering App"),
  # Somewhere in UI
  add_busy_gif(src = "https://media2.giphy.com/media/7aqTC87afdWGXF3o8T/giphy.gif?cid=ecf05e475vgr4jui6uha1neqjhfbpsqwsa6fwklob3kzcns1&ep=v1_gifs_related&rid=giphy.gif&ct=s",
               height = 150 , width = 150 , position = "full-page"),
  sidebarLayout(
    sidebarPanel(
      fileInput('file1', 'Choose CSV File', accept=c('.csv')),
      radioButtons("method", label = "Method:", choices = c("UMI", "Spike"), selected = "UMI"),
      textInput("mito_gene", "Mitochondrial Gene Symbol", value = "^mt-"),
      textInput("num_genes", "Number of Genes", value = "1000"),
      actionButton("start_button", "Start Process"),
      downloadButton("downloadCSV", "Download Results as CSV"),
      actionButton("exampleBtn", "Run Example"),
      tags$br(),  # Adds a line break
      "You can just upload your count matrix (Columns are cells, rows are features) as
      CSV file and click 'Start Process' button. Based on your data size it may take 
      up to 30 minutes.
      This tool may require lots of memory (> 8GB) based on your data. It may crash
      in ShinyApps server due to this reason.
      If you want to see how the tool works and the outputs, you can just click
      'Run Example' button."
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data Quality",
                 h3("Data Summary"),
                 verbatimTextOutput("summary"),  # For displaying text
                 h3("Data Statistics"),
                 DT::DTOutput("stats_table"),
                 plotOutput("data_quality_plots" , height="600px")),        # For displaying plot
        
        tabPanel("Random Forest Based Clustering",
                 verbatimTextOutput("modelselect"),
                 plotOutput("model_plot" , height="600px"),
                 plotOutput("cell_freq_model", height="600px" )),
        
        tabPanel("Clustering Results", 
                 selectInput("clustering",
                             "Select Clustering Method", choices = c("Seurat Based", "Dunn Index Based Louvain",
                                                                     "Dunn Index Based Leiden",
                                                                     "Consensus Clustering (Auto)")),
                 plotOutput("cluster_plot" , height="600px"),
                 plotOutput("cell_freq_clust", height="600px" ))
      )
    )
  )
)

# Server Section
server <- function(input, output) {
  observeEvent(input$start_button, {
    inFile <- input$file1
    
    if (!is.null(inFile)) {
      # Step 1
      data_matrix <- as.matrix(read.csv(inFile$datapath, row.names=1))
      
      # Step 2
      sce <- SingleCellExperiment(assays = list(counts=data_matrix))
      output$summary <- renderPrint({
        print(sce)
      })
      number_of_genes <- as.integer(input$num_genes)
      mito_gene <- input$mito_gene
      
      if (input$method == "UMI") {
        # UMI method steps...
        stats <- perCellQCMetrics(sce, subsets=list(Mito=which(rownames(sce) == mito_gene)))
        stats_table <- cbind(summary(stats$sum) , summary(stats$detected) , 
                             summary(stats$subsets_Mito_percent))
        colnames(stats_table) <- c("Sum" , "Detected" , "Mitochondrial Gene Percentage")
        
        output$stats_table <- DT::renderDT({
          print(stats_table)
        })
        
        high.sum <- isOutlier(stats$sum, nmads = 3, type="both")
        high.total <- isOutlier(stats$total, nmads = 3, type="both")
        high.mito <- isOutlier(stats$subsets_Mito_percent , nmads = 3 , type = "both")
        discard <- high.sum | high.total | high.mito
        sce$discard <- discard
        is.mito <- grepl(mito_gene , rownames(sce))
        sce <- addPerCellQCMetrics(sce, subsets=list(Mito=is.mito))
        
        output$data_quality_plots <- renderPlot({
          grid.arrange(
            plotColData(sce,  y="sum", colour_by="discard") + 
              scale_y_log10() + ggtitle("Total count"),
            plotColData(sce,  y="detected", colour_by="discard") + 
              scale_y_log10() + ggtitle("Detected Features"),
            plotColData(sce,  y="subsets_Mito_percent", 
                        colour_by="discard") +
              ggtitle("Mito percent"),
            ncol=1
          )
        })
        
        sce <- sce[,!(high.sum | high.total | high.mito)]
        is.mito <- grepl(mito_gene , rownames(sce))
        clusters <- quickCluster(sce)
        sce <- computeSumFactors(sce, cluster=clusters)
        sce <- logNormCounts(sce)
        dec <- modelGeneVar(sce)
        top <- getTopHVGs(dec , n = number_of_genes)
        sce <- denoisePCA(sce, subset.row=top, technical=dec)
        number_of_pcs <- ncol(reducedDim(sce))
      } else {
        # Not UMI method steps...
        is.spike <- grepl("^ERCC", rownames(sce))
        sce <- splitAltExps(sce, ifelse(is.spike, "ERCC", "gene"))
        stats <- perCellQCMetrics(sce)
        stats_table <- cbind(summary(stats$sum) , summary(stats$detected) , 
                             summary(stats$altexps_ERCC_percent))
        colnames(stats_table) <- c("Sum" , "Detected" , "ERCC Expression Percentage")
        
        output$stats_table <- DT::renderDT({
          print(stats_table)
        })
        high.sum <- isOutlier(stats$sum, nmads = 3, type="both")
        high.total <- isOutlier(stats$total, nmads = 3, type="both")
        spike.drop <- isOutlier(stats$altexps_ERCC_percent, nmads=3, type="both")
        discard <- high.sum | high.total | spike.drop
        sce <- addPerCellQCMetrics(sce, subsets=list(Spike=is.spike))
        output$data_quality_plots <- renderPlot({
          grid.arrange(
            plotColData(sce,  y="sum", colour_by="discard") + 
              scale_y_log10() + ggtitle("Total count"),
            plotColData(sce,  y="detected", colour_by="discard") + 
              scale_y_log10() + ggtitle("Detected Features"),
            plotColData(sce,  y="subsets_Spike_percent", 
                        colour_by="discard") +
              ggtitle("ERCC Spike Percent"),
            ncol=1
          )
        })
        
        sce$discard <- discard
        sce <- sce[,!(high.sum | high.total | spike.drop)]
        clusters <- quickCluster(sce)
        sce <- computeSumFactors(sce, cluster=clusters)
        sce <- logNormCounts(sce)
        dec <- modelGeneVarWithSpikes(sce , "ERCC")
        top <- getTopHVGs(dec ,  n = number_of_genes)
        sce <- denoisePCA(sce, subset.row=top, technical=dec)
        number_of_pcs <- ncol(reducedDim(sce))
      }
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
      model_data$Algorithm <- as.integer(model_data$Algorithm)
      model_data$Exp_Type <- as.integer(model_data$Exp_Type)
      model_data$var_perc <- as.numeric(model_data$var_perc)
      rf_pred <- rf_pred <- predict(suma_r , newdata = model_data , type = "prob")[,2]
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
      model_selected_algorithm <- c()
      if (model_data[as.numeric(names(sort(rf_pred , decreasing = TRUE)))[1],]$Algorithm == 1) {
        model_selected_algorithm <- "Leiden"
      }
      if (model_data[as.numeric(names(sort(rf_pred , decreasing = TRUE)))[1],]$Algorithm == 2) {
        model_selected_algorithm <- "Louvain"
      }
      if (model_data[as.numeric(names(sort(rf_pred , decreasing = TRUE)))[1],]$Algorithm == 3) {
        model_selected_algorithm <- "Walktrap"
      }
      
      model_res_text <- paste0("RandomForest Model Selected " , model_selected_algorithm, " as Clustering Algorithm and ",
                               as.character(selected_k), " as Number of Neighbours as the best parameters.")
      
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
      
      ### Clustering Based on Dunn Index (Leiden) ###
      
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
      sce.seurat@meta.data$model_clusters <- model_clusters
      sce.seurat@meta.data$dunn_clusters_leid <- dunn_clusters_leid
      
      
      output$model_plot <- renderPlot({
        UMAPPlot(sce.seurat , group.by = "model_clusters") + ggtitle("Clustering Based On Random Forest Model Selection")
      })
      
      output$cell_freq_model <- renderPlot({
        ggplot(as.data.frame(table(model_clusters)), aes(x=model_clusters, y=Freq)) +
          geom_bar(stat="identity", fill="steelblue", color="black") + 
          theme_minimal() + 
          labs(title="Number of Cells in Each Cluster (Model Clusters)", x="Cluster", y="Count") + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_text(aes(label=Freq), vjust=-0.5)
      })
      
      output$modelselect <- renderText({
        print(model_res_text)
      })
      
      output$cluster_plot <- renderPlot({
          if (input$clustering == "Dunn Index Based Louvain") {
          UMAPPlot(sce.seurat , group.by = "dunn_clusters_louv") + ggtitle("Dunn Index Based Louvain Clustering")
        } else if (input$clustering == "Dunn Index Based Leiden") {
          UMAPPlot(sce.seurat , group.by = "dunn_clusters_leid") + ggtitle("Dunn Index Based Leiden Clustering")
        } else if (input$clustering == "Seurat Based") {
          UMAPPlot(sce.seurat) + ggtitle("Seurat Clustering (Default Parameters)")
        } else if (input$clustering == "Consensus Clustering (Auto)") {
          UMAPPlot(sce.seurat , group.by = "consensus_clusters") + ggtitle("Consensus Clustering (Automatic K selection)")
        }
      })
      
      output$cell_freq_clust <- renderPlot({
        if (input$clustering == "Dunn Index Based Louvain") {
          ggplot(as.data.frame(table(dunn_clusters_louv)), aes(x=dunn_clusters_louv, y=Freq)) +
            geom_bar(stat="identity", fill="steelblue", color="black") + 
            theme_minimal() + 
            labs(title="Number of Cells in Each Cluster (Dunn Index Based Louvain Clustering)", x="Cluster", y="Count") + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_text(aes(label=Freq), vjust=-0.5)
        } else if (input$clustering == "Dunn Index Based Leiden") {
          ggplot(as.data.frame(table(dunn_clusters_leid)), aes(x=dunn_clusters_leid, y=Freq)) +
            geom_bar(stat="identity", fill="steelblue", color="black") + 
            theme_minimal() + 
            labs(title="Number of Cells in Each Cluster (Dunn Index Based Leiden Clustering)", x="Cluster", y="Count") + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_text(aes(label=Freq), vjust=-0.5)
        } else if (input$clustering == "Seurat Based") {
          ggplot(as.data.frame(table(seurat_clusters)), aes(x=seurat_clusters, y=Freq)) +
            geom_bar(stat="identity", fill="steelblue", color="black") + 
            theme_minimal() + 
            labs(title="Number of Cells in Each Cluster (Seurat Clustering)", x="Cluster", y="Count") + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_text(aes(label=Freq), vjust=-0.5)
        } else if (input$clustering == "Consensus Clustering (Auto)") {
          ggplot(as.data.frame(table(consensus_clusters)), aes(x=consensus_clusters, y=Freq)) +
            geom_bar(stat="identity", fill="steelblue", color="black") + 
            theme_minimal() + 
            labs(title="Number of Cells in Each Cluster (Consensus Clustering)", x="Cluster", y="Count") + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_text(aes(label=Freq), vjust=-0.5)
        }
      })
      
      results_matrix <- cbind(colnames(sce.seurat) , as.numeric(model_clusters) , as.numeric(dunn_clusters_louv) , 
                              as.numeric(dunn_clusters_leid) , as.numeric(sce.seurat@active.ident) , as.numeric(consensus_clusters))
      results_matrix <- as.data.frame(results_matrix)
      colnames(results_matrix) <- c("Model Clusters" , "Louvain Clusters" , "Leiden Clusters" ,
                                    "Seurat Clusters" , "Consensus Clusters")
      
      output$downloadCSV <- downloadHandler(
        filename = function() {
          paste("data-", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          write.table(results_matrix , file , row.names = FALSE , sep = "," , col.names = T , quote = F)
        }
      )
    }
  })
  
  observeEvent(input$exampleBtn, {
    # Read the RDS file
    example_data <- readRDS("example.RDS")
    sce <- example_data$sce
    sce.seurat <- example_data$sce.seurat
    stats_table <- example_data$stats_table
    model_clusters <- example_data$model_clusters
    consensus_clusters <- example_data$consensus_clusters
    dunn_clusters_leid <- example_data$dunn_clusters_leid
    dunn_clusters_louv <- example_data$dunn_clusters_louv
    model_selected_algorithm <- example_data$model_selected_algorithm
    selected_k <- example_data$selected_k
    seurat_clusters <- sce.seurat@active.ident
    rm(example_data)
    
    model_res_text <- paste0("RandomForest Model Selected " , model_selected_algorithm, " as Clustering Algorithm and ",
                             as.character(selected_k), " as Number of Neighbours as the best parameters.")
    
    output$stats_table <- DT::renderDT({
      print(stats_table)
    })
    
    output$data_quality_plots <- renderPlot({
      grid.arrange(
        plotColData(sce,  y="sum", colour_by="discard") + 
          scale_y_log10() + ggtitle("Total count"),
        plotColData(sce,  y="detected", colour_by="discard") + 
          scale_y_log10() + ggtitle("Detected Features"),
        plotColData(sce,  y="subsets_Mito_percent", 
                    colour_by="discard") +
          ggtitle("Mito percent"),
        ncol=1
      )
    })
    
    output$model_plot <- renderPlot({
      UMAPPlot(sce.seurat , group.by = "model_clusters") + ggtitle("Clustering Based On Random Forest Model Selection")
    })
    
    output$cell_freq_model <- renderPlot({
      ggplot(as.data.frame(table(model_clusters)), aes(x=model_clusters, y=Freq)) +
        geom_bar(stat="identity", fill="steelblue", color="black") + 
        theme_minimal() + 
        labs(title="Number of Cells in Each Cluster (Model Clusters)", x="Cluster", y="Count") + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_text(aes(label=Freq), vjust=-0.5)
    })
    
    output$modelselect <- renderText({
      print(model_res_text)
    })
    
    output$cluster_plot <- renderPlot({
      if (input$clustering == "Dunn Index Based Louvain") {
        UMAPPlot(sce.seurat , group.by = "dunn_clusters_louv") + ggtitle("Dunn Index Based Louvain Clustering")
      } else if (input$clustering == "Dunn Index Based Leiden") {
        UMAPPlot(sce.seurat , group.by = "dunn_clusters_leid") + ggtitle("Dunn Index Based Leiden Clustering")
      } else if (input$clustering == "Seurat Based") {
        UMAPPlot(sce.seurat) + ggtitle("Seurat Clustering (Default Parameters)")
      } else if (input$clustering == "Consensus Clustering (Auto)") {
        UMAPPlot(sce.seurat , group.by = "consensus_clusters") + ggtitle("Consensus Clustering (Automatic K selection)")
      }
    })
    
    output$cell_freq_clust <- renderPlot({
      if (input$clustering == "Dunn Index Based Louvain") {
        ggplot(as.data.frame(table(dunn_clusters_louv)), aes(x=dunn_clusters_louv, y=Freq)) +
          geom_bar(stat="identity", fill="steelblue", color="black") + 
          theme_minimal() + 
          labs(title="Number of Cells in Each Cluster (Dunn Index Based Louvain Clustering)", x="Cluster", y="Count") + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_text(aes(label=Freq), vjust=-0.5)
      } else if (input$clustering == "Dunn Index Based Leiden") {
        ggplot(as.data.frame(table(dunn_clusters_leid)), aes(x=dunn_clusters_leid, y=Freq)) +
          geom_bar(stat="identity", fill="steelblue", color="black") + 
          theme_minimal() + 
          labs(title="Number of Cells in Each Cluster (Dunn Index Based Leiden Clustering)", x="Cluster", y="Count") + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_text(aes(label=Freq), vjust=-0.5)
      } else if (input$clustering == "Seurat Based") {
        ggplot(as.data.frame(table(seurat_clusters)), aes(x=seurat_clusters, y=Freq)) +
          geom_bar(stat="identity", fill="steelblue", color="black") + 
          theme_minimal() + 
          labs(title="Number of Cells in Each Cluster (Seurat Clustering)", x="Cluster", y="Count") + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_text(aes(label=Freq), vjust=-0.5)
      } else if (input$clustering == "Consensus Clustering (Auto)") {
        ggplot(as.data.frame(table(consensus_clusters)), aes(x=consensus_clusters, y=Freq)) +
          geom_bar(stat="identity", fill="steelblue", color="black") + 
          theme_minimal() + 
          labs(title="Number of Cells in Each Cluster (Consensus Clustering)", x="Cluster", y="Count") + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_text(aes(label=Freq), vjust=-0.5)
      }
    })
    
    results_matrix <- cbind(colnames(sce.seurat) , as.numeric(model_clusters) , as.numeric(dunn_clusters_louv) , 
                            as.numeric(dunn_clusters_leid) , as.numeric(sce.seurat@active.ident) , as.numeric(consensus_clusters))
    results_matrix <- as.data.frame(results_matrix)
    colnames(results_matrix) <- c("Model Clusters" , "Louvain Clusters" , "Leiden Clusters" ,
                                  "Seurat Clusters" , "Consensus Clusters")
    
    output$downloadCSV <- downloadHandler(
      filename = function() {
        paste("example_results-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.table(results_matrix , file , row.names = FALSE , sep = "," , col.names = T , quote = F)
      }
    )
    
  })
  
}

# Run the Shiny App with the above server function and the provided UI
shinyApp(ui = ui, server = server)
