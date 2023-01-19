TestAllProcesses <- function(object, assay = DefaultAssay(object = object), UMIs = "nCount_RNA",
                             TFs = "All_TFs", stats = t.test){
  TestProcesses_1 <- function(object, cluster){
    if(median(object@meta.data[UMIs][,1])>median(cluster@meta.data[UMIs][,1]) &
       median(object@meta.data[TFs][,1])>median(cluster@meta.data[TFs][,1])){
      sprintf("cluster %s is probably cell processes", i)
    }
    else{
      sprintf("0")}}
  TestProcesses_Outcome <- c()
  a <- c()
  b <- c()
  message("Initiating...")
  for(i in levels(object)){
    level_PC <- which(levels(object) == i)
    cluster <- subset(object, idents = i)
    p = stats(processes_cluster@meta.data[UMIs][,1],object@meta.data[UMIs][,1])$p.value
    a <- append(a, p)
    q = stats(processes_cluster@meta.data[TFs][,1],object@meta.data[TFs][,1])$p.value
    b <- append(b, q)}
  a <- p.adjust(a, method = "fdr")
  b <- p.adjust(b, method = "fdr")
  for(i in levels(object)){
    level_PC <- which(levels(object) == i)
    message(sprintf("Testing Cluster %s", i))
    if(a[level_PC] < 0.005 & b[level_PC] < 0.005){
      cluster <- subset(object, idents = i)
      if(print(TestProcesses_1(object, cluster)) != 0){
        TestProcesses_Outcome <- rbind(TestProcesses_Outcome, i)
        colnames(TestProcesses_Outcome) <- c("Cluster")
        .GlobalEnv$TestProcesses_Outcome <- TestProcesses_Outcome
        #write.table(TestProcesses_Outcome, file = "", 
        #           sep = "\t", row.names = F)
      }
    }
  }
}




TestAllProcesses_Correlation <- function(object, assay = DefaultAssay(object = object), UMIs = "nCount_RNA",
                                         TFs = "All_TFs", stats = t.test, correlation = 0.9, 
                                         correlation_method = NULL, features = VariableFeatures(object)){
  avg  <- AverageExpression(object, assays = assay, features = features)
  cor_matrix <- cor(as.data.frame(avg[1]), method = correlation_method)
  TestProcesses_1 <- function(object, processes = NULL, cell_body = NULL){
    if(median(processes_cluster@meta.data[UMIs][,1])<median(cell_body_cluster@meta.data[UMIs][,1]) &
       median(processes_cluster@meta.data[TFs][,1])<median(cell_body_cluster@meta.data[TFs][,1])){
      sprintf("cluster %s is probably cell processes of cluster %s", processes, cell_body)
    }
    else{
      sprintf("0")}}
  TestProcesses_Outcome <- c()
  for(i in levels(object)){
    level_PC <- which(levels(object) == i)
    message(sprintf("Testing Cluster %s", i))
    processes_cluster <- subset(object, idents = i)
    a <- c()
    b <- c()
    for(j in levels(object)){
      cell_body_cluster <- subset(object, idents = j)
      p = stats(processes_cluster@meta.data[UMIs][,1],cell_body_cluster@meta.data[UMIs][,1])$p.value
      a <- append(a, p)
      q = stats(processes_cluster@meta.data[TFs][,1],cell_body_cluster@meta.data[TFs][,1])$p.value
      b <- append(b, q)}
    a <- p.adjust(a, method = "fdr")
    b <- p.adjust(b, method = "fdr")
    for(j in levels(object)){
      level_CB <- which(levels(object) == j)
      if(a[level_CB] < 0.005 & b[level_CB] < 0.005){
        if(cor_matrix[level_PC,level_CB] >= correlation & cor_matrix[level_PC,level_CB] < 1){
          cell_body_cluster = subset(object, idents = j)
          if(print(TestProcesses_1(object, i, j)) != 0){
            TestProcesses_Outcome <- rbind(TestProcesses_Outcome,c(i,j, cor_matrix[level_PC,level_CB]))
            colnames(TestProcesses_Outcome) <- c("Process", "Cell_Body", "Correlation")
            .GlobalEnv$TestProcesses_Outcome <- TestProcesses_Outcome
            #write.table(TestProcesses_Outcome, file = "", 
            #           sep = "\t", row.names = F)
          }
        }
      }
    }
  }
}


