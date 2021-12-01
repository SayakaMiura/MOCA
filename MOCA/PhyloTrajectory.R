PhyloTrajectory <- function(expression_matrix,ancestry_annotation,titles){
  gene_scale <- c(200, 300, 400, 500, 600, 700, 800, 900, 1000) #Can edit to test on different numbers of genes
  dominant_cutoff <- 0.8 #The percent of cells belonging to a single expression state or genetic ancestry for it to be considered dominant.
  colnames(ancestry_annotation)<- c('Cell', 'GeneticAncestry')
  cluster_matrix <- data.frame(matrix(nrow = nrow(ancestry_annotation), ncol = length(gene_scale)+1))
  cluster_matrix[,1] <- as.numeric(ancestry_annotation[,2])
  genes <- data.frame(row.names(expression_matrix))
  row.names(genes) <- genes[,1] 
  colnames(genes)<-  "gene_short_name"
  library('monocle')
  geneticData <- new("AnnotatedDataFrame", data = ancestry_annotation)
  genes <- new("AnnotatedDataFrame", data = genes)
  
  #Create a CellDataSet for Monocle Trajectory Analysis
  Input_Matrix <- newCellDataSet(as.matrix(expression_matrix),
                                 phenoData = geneticData,
                                 featureData = genes,
                                 expressionFamily=uninormal())
  Input_Matrix <- detectGenes(Input_Matrix, min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(Input_Matrix),
                                      num_cells_expressed >= 10))
  fData(Input_Matrix)$use_for_ordering <-
    fData(Input_Matrix)$num_cells_expressed > 0.05 * ncol(Input_Matrix)
  Input_Matrix_expressed_genes <-  row.names(subset(fData(Input_Matrix),
                                                    num_cells_expressed >= 10))
  clustering_DEG_genes <-
    differentialGeneTest(Input_Matrix[Input_Matrix_expressed_genes,],
                         fullModelFormulaStr = '~GeneticAncestry',
                         cores = 1)
  ancestries <- as.character(sort(unique(pData(Input_Matrix)$GeneticAncestry)))
  final_results <- data.frame(matrix(nrow = length(gene_scale), ncol = length(ancestries)))
  final_results[is.na(final_results)] <- 0
  gg1 <- data.frame(matrix(nrow = 2, ncol =  length(gene_scale)))
  for(a in 1:length(gene_scale)){
    Input_Matrix_ordering_genes <-
      row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:gene_scale[a]]
    Input_Matrix_Filtered <- Input_Matrix[Input_Matrix_ordering_genes,]
    Input_Matrix_Filtered <- reduceDimension(Input_Matrix_Filtered,
                                             max_components = 2,
                                             norm_method = 'none',
                                             num_dim = 3,
                                             reduction_method = 'tSNE',
                                             verbose = T, 
                                             perplexity = 2)
    Input_Matrix_Filtered <- clusterCells(Input_Matrix_Filtered, verbose = F, num_clusters = 4)
    fileName <- paste0(gene_scale[a], "Clusters.png")
    png(filename = fileName)
    print(plot_cell_clusters(Input_Matrix_Filtered))
    dev.off()
    Input_Matrix_Filtered <-
      setOrderingFilter(Input_Matrix_Filtered,
                        ordering_genes = Input_Matrix_ordering_genes)
    Input_Matrix_Filtered <-
      reduceDimension(Input_Matrix_Filtered, method = 'DDRTree', norm_method = 'none')
    
    Input_Matrix_Filtered <-
      orderCells(Input_Matrix_Filtered)
    
    #Populates a table containing pseudotime values, genetic ancestry classification, and expression state annotation.
    
    results <- data.frame(pData(Input_Matrix_Filtered))
    ancestries <- as.character(sort(unique(results$GeneticAncestry)))
    states <- as.character(sort(unique(results$State)))
    Ancestry_Colors <- sample(colors(distinct = TRUE), length(ancestries))
    State_Colors <- sample(colors(distinct = TRUE), length(states))

    
    fileName <- paste0(gene_scale[a], "GeneticAncestryTrajectory.png")
    png(filename = fileName)
    print(plot_cell_trajectory(Input_Matrix_Filtered, color_by = "GeneticAncestry", cell_size = 2) + 
            labs(title  = paste(titles, gene_scale[a], sep = " ")) +
            scale_color_manual(values = c(Ancestry_Colors))+ theme(text = element_text(size = 15)))
    dev.off()
    fileName <- paste0(gene_scale[a], "PseudotimeTrajectory.png")
    png(filename = fileName)
    print(plot_cell_trajectory(Input_Matrix_Filtered, color_by = "Pseudotime", cell_size = 2)+ theme(text = element_text(size = 15))+
            labs(title  = paste(titles, gene_scale[a], sep = " ")))
    dev.off()
    fileName <- paste0(gene_scale[a], "StateTrajectory.png")
    png(filename = fileName)
    print(plot_cell_trajectory(Input_Matrix_Filtered, color_by = "State", cell_size = 2) + 
            labs(title  = paste(titles, gene_scale[a], sep = " ")) +
            scale_color_manual(values = c(State_Colors))+ theme(text = element_text(size = 15)))
    dev.off()
    plot_cell_trajectory(Input_Matrix_Filtered, color_by = "State", cell_size = 2) + 
      scale_color_manual(values = c(State_Colors))+ theme(text = element_text(size = 15))
    
    
    bargraph <- data.frame(matrix(nrow = length(ancestries)*length(states), ncol = 3))
    colnames(bargraph)<-c("GeneticAncestry", "NumberCells", "ExpressionState")
    z=1
    for(x in 1:length(states)){
      for(y in 1:length(ancestries)){
        bargraph[z,1] <- ancestries[y]
        bargraph[z,2] <-length(which(results$GeneticAncestry == ancestries[y] & results$State == states[x]))
        bargraph[z,3] <- states[x]
        z=z+1
      }
      y = 1 
    }
    fileName <- paste0(gene_scale[a], "GeneticAncestryHeterogeneityGeneScale.png")
    png(filename = fileName)
    print(ggplot(data = bargraph, aes(x = GeneticAncestry, y = NumberCells, fill = ExpressionState)) + 
            geom_col() + 
            labs(title = paste("Genetic Ancestry Heterogeneity Gene Scale:", gene_scale[a], sep = " "), subtitle = titles) + 
            scale_fill_manual(values = c(State_Colors))+ theme(text = element_text(size = 15)))
    dev.off()
    fileName <- paste0(gene_scale[a], "ExpressionStateHeterogeneityGeneScale.png")
    png(filename = fileName)
    print(ggplot(data = bargraph, aes(x = ExpressionState, y = NumberCells, fill = GeneticAncestry)) + 
            geom_col() + 
            labs(title = paste("Expression State Heterogeneity Gene Scale", gene_scale[a], sep = " "), subtitle = titles) +
            scale_fill_manual(values = c(Ancestry_Colors))+ theme(text = element_text(size = 15)))
    dev.off()
    fileName <- paste0(gene_scale[a], "GeneticAncestryCompositionGeneScale.png")
    png(filename = fileName)
    print(ggplot(data = bargraph, aes(x = GeneticAncestry, y = NumberCells, fill = ExpressionState)) + 
            geom_bar(position="fill", stat="identity") + 
            labs(title = paste("Genetic Ancestry Composition Gene Scale", gene_scale[a] , sep = " "), subtitle = titles) + 
            scale_y_continuous(labels = scales::percent_format()) + 
            scale_fill_manual(values = c(State_Colors)) + 
            ylab("Percent of Cells")+ theme(text = element_text(size = 15)))
    dev.off()
    fileName <- paste0(gene_scale[a], "ExpressionStateCompGeneScale.png")
    png(filename = fileName)
    print(ggplot(data = bargraph, aes(x = ExpressionState, y = NumberCells, fill = GeneticAncestry)) + 
            geom_bar(position="fill", stat="identity") + 
            labs(title = paste("Expression State Composition Gene Scale:", gene_scale[a], sep = " "), subtitle = titles) + 
            scale_y_continuous(labels = scales::percent_format()) +
            scale_fill_manual(values = c(Ancestry_Colors)) +
            ylab("Percent of Cells")+ theme(text = element_text(size = 15)))
    dev.off()
    boxplot2 <-data.frame(as.factor(results$GeneticAncestry), as.numeric(results$Pseudotime))
    colnames(boxplot2)<-c("GeneticAncestry", "Pseudotime")
    fileName <- paste0(gene_scale[a], "ExpressionStateVariance.png")
    png(filename = fileName)
    print(ggplot(boxplot2, aes(x = GeneticAncestry, y = Pseudotime, fill = GeneticAncestry)) + 
            geom_boxplot() + 
            labs(title = paste("Expression Variance Within Genetic Ancestries Gene Scale:", gene_scale[a], sep =" "), subtitle = titles) +
            scale_fill_manual(values = c(Ancestry_Colors)))
    dev.off()
    cluster_matrix[,a+1] <- as.numeric(results$Cluster)
    colnames(cluster_matrix)[a+1] <- paste(gene_scale[a])
    concordance <- data.frame(matrix(ncol = length(ancestries)+1, nrow = length(states)+2))
    concordance[is.na(concordance) == TRUE] <- 0
    colnames(concordance) <- c(ancestries, "Major Genetic Group")
    rownames(concordance) <- c(states, "Major Expression State", "Sequence Type Expression Standard Deviation")
    for(x in 1:length(states)){
      for(y in 1:length(ancestries)){
        concordance[x,y] <-length(which(results$GeneticAncestry == ancestries[y] & results$State == states[x]))
      }
    }
    for (e in 1:length(states)) {
      total <-sum(concordance[e,1:length(ancestries)])
      greatest <- max(concordance[e,1:length(ancestries)])
      group <- which(concordance[e,]==greatest)
      percent <- greatest / total
      concordance[e, ncol(concordance)] <- percent
    }
    for (c in 1:length(ancestries)) {
      total <-sum(concordance[1:length(states),c])
      greatest <- max(concordance[1:length(states),c])
      state <- which(concordance[,c]==greatest)
      percent <- greatest / total
      concordance[(nrow(concordance)-1), c] <- percent
    }
    for (b in 1:length(ancestries)) {
      concordance[nrow(concordance),b]<-sd(boxplot2$Pseudotime[which(boxplot2$GeneticAncestry == ancestries[b])])
    }
    cth <- length(which(concordance[1:length(states),ncol(concordance)] >= dominant_cutoff))
    gg1[2,a] <- cth
    gg1[1,a] <- length(states) - cth
    row.names(gg1) <- c("States Without", "States With")
    paths <- paste0(gene_scale[a], "concordance_stats.txt")
    colnames(final_results) <- ancestries
    row.names(final_results) <- gene_scale
    cc1 <- which(concordance[,ncol(concordance)] >= dominant_cutoff)
    if(length(cc1) > 0){
      for(h in 1:length(cc1)){
        dom_group <- ancestries[which(concordance[cc1[h],] == max(concordance[cc1[h],]))]
        final_results[a,which(colnames(final_results) == dom_group)] <- final_results[a,which(colnames(final_results) == dom_group)] +1
      }
    }
    colnames(concordance) <- c(ancestries, "Major Genetic Group")
    for(p in 1:length(ancestries) - 1){
      colnames(concordance)[p] <- paste0("Ancestry ", ancestries[p])
    }
    rownames(concordance) <- c(states, "Major Expression State", "Sequence Type Expression Standard Deviation")
    for(q in 1:length(states)){
      row.names(concordance)[q] <- paste0("State ", states[q])
    }
    write.table(concordance, paths, row.names = TRUE, col.names = TRUE)
  }
  fileName <- paste0(gene_scale[a], "ProportionStates.png")
  png(filename = fileName)
  
  barplot(as.matrix(gg1),  main = "States With and Without Dominant Ancestry", xlab = "Gene Scale", ylab = "Number of States", names.arg = gene_scale, legend = row.names(gg1), col = c("orange", "purple"))
  dev.off()
  fileName <- paste0(gene_scale[a], "DominantStates.png")
  png(filename = fileName)
  barplot(as.matrix(final_results), beside = TRUE, xlab = "Genetic Ancestries", ylab = "Number of Dominant States", legend.text = row.names(final_results), col=colors()[c(23,89,12,99,110,1,13,76,81)])
  dev.off()
  row.names(cluster_matrix) <- row.names(ancestry_annotation)
  colnames(cluster_matrix)[1] <- "AncestryAnnotation"
  cluster_matrix <- arrange(cluster_matrix, cluster_matrix$AncestryAnnotation)
  annotations <- cluster_matrix[,1]
  cluster_matrix <- cluster_matrix[,-1]
  colside <- Ancestry_Colors[as.character(annotations)]
  fileName <- paste0("ExpressionClustering_Annotation.png")
  png(filename = fileName)
  heatmap(as.matrix(cluster_matrix),Colv = NA, Rowv = NA, scale = "none", col = terrain.colors(max(cluster_matrix)), xlab = "Gene Scales", RowSideColors =as.character(annotations), main = "Expression Clustering Compared to Ancestry Annotation")
  dev.off()
  
  #Computes the overall statistics across gene scales
  num_new <- length(gene_scale)
  newTable <- paste0(gene_scale[1], "concordance_stats.txt")
  lengths <- data.frame(matrix(nrow = 3, ncol = num_new))
  tables <- data.frame(matrix(nrow = num_new, ncol = 1))
  tables[1,1] <- newTable
  newTable <- read.table(newTable)
  x2 <- newTable[,ncol(newTable)]
  x2 <- x2[-which(x2 ==0)]
  x3 <- newTable[nrow(newTable),]
  x3 <- x3[-which(x3 == 0)]
  x4 <- newTable[nrow(newTable) - 1,]
  x4 <- x4[-which(x4 == 0)]
  lengths[1,1] <- length(x2)
  lengths[2,1] <- length(x3)
  lengths[3,1] <- length(x4)
  for(a in 2:num_new){
    stats1 <-  paste0(gene_scale[a], "concordance_stats.txt")
    newTable <- read.table(stats1)
    tables[a,1] <- stats1
    x2 <- newTable[,ncol(newTable)]
    x2 <- x2[-which(x2 ==0)]
    x3 <- newTable[nrow(newTable),]
    x3 <- x3[-which(x3 == 0)]
    x4 <- newTable[nrow(newTable) - 1,]
    x4 <- x4[-which(x4 == 0)]
    lengths[1,a] <- length(x2)
    lengths[2,a] <- length(x3)
    lengths[3,a] <- length(x4)
  }
  num_states <- max(lengths[1,])
  time_index <- max(lengths[2,])
  num_groups <- max(lengths[3,])
  most <- which(lengths[1,] == num_states)
  state_indices <- read.table(paste0(tables[most[1],1]))
  states_indices <- data.frame(matrix(nrow = num_states, ncol = num_new +1))
  states_indices[,most[1]] <- state_indices[1:num_states,ncol(state_indices)]
  w <- which(1:ncol(lengths) != most[1])
  while(length(w) != 0){
    state_indices2 <- read.table(paste0(tables[w[1],1]))
    state_vector <- state_indices2[1:(nrow(state_indices2)-2),ncol(state_indices2)]
    state_vector <- c(state_vector, rep(0, num_states-length(state_vector)))
    states_indices[,w[1]] <- state_vector
    w <- w[-1]
  }
  most2 <- which(lengths[2,] == time_index)
  time_indices <- read.table(paste0(tables[most2[1],1]))
  times_indices <- data.frame(matrix(nrow = time_index, ncol = num_new +1))
  t_vector <- time_indices[nrow(time_indices),]
  t_vector <- t_vector[which(t_vector != 0)]
  t_vector <- t(t_vector)
  times_indices[1:length(t_vector),most2[1]] <- t_vector
  w <- which(1:ncol(lengths) != most2[1])
  while(length(w) != 0){
    time_indices2 <- read.table(paste0(tables[w[1],1]))
    t_vector <- time_indices2[nrow(time_indices2),]
    t_vector <- t_vector[which(t_vector != 0)]
    t_vector <- t(t_vector)
    t_vector <- c(t_vector, rep(0,(time_index -length(t_vector))))
    times_indices[,w[1]] <- t_vector
    w <- w[-1]
  }
  most3 <- which(lengths[3,] == num_groups)
  type_indices <- read.table(paste0(tables[most3[1],1]))
  types_indices <- data.frame(matrix(nrow = num_groups, ncol = num_new +1))
  type_vector <- type_indices[nrow(type_indices)-1,1:(ncol(type_indices)-1)]
  type_vector <- t(type_vector)
  types_indices[,most3[1]] <- type_vector
  w <- which(1:ncol(lengths) != most3[1])
  while(length(w) != 0){
    type_indices2 <- read.table(paste0(tables[w[1],1]))
    type_vector <- type_indices2[nrow(type_indices2)-1,1:(ncol(type_indices2)-1)]
    type_vector <- t(type_vector)
    t_vector <- c(type_vector, rep(0,num_groups -length(type_vector)))
    types_indices[,w[1]] <- type_vector
    w <- w[-1]
  }
  types_indices <- data.frame(t(types_indices))
  colnames(types_indices) <- colnames(type_indices)[1:ncol(types_indices)]
  states_indices <- data.frame(t(states_indices))
  colnames(states_indices) <- row.names(state_indices)[1:ncol(states_indices)]
  times_indices <- data.frame(t(times_indices))
  colnames(times_indices) <- colnames(time_indices)[1:ncol(times_indices)]
  
  #Generates random colors.
  colors2 <- sample(colors(distinct = TRUE), num_groups)
  
  ends <- str_locate_all(tables[,1], "concordance_stats")
  starts <- str_locate_all(tables[,1], "/")
  legend_titles <- data.frame(matrix(ncol = 1, nrow = nrow(tables)))
  for(v in 1:nrow(tables)){
    ends <- min(str_locate_all(tables[v,1], "concordance_stats")[[1]])-1
    starts <- max(str_locate_all(tables[v,1], "/")[[1]])+1
    legend_titles[v,1] <- substring(paste0(tables[v,1]),starts,ends)
  }
  fileName <- paste0("DominantAncestry1.png")
  png(filename = fileName)
  barplot(as.matrix(types_indices), beside = TRUE, main = "Genetic Ancestry Concordance Index", col = colors2, xlab = "Genetic Ancestries", ylab = "Proportion of Dominant Expression State", legend.text = legend_titles[,1])
  dev.off()
  fileName <- paste0("DominantState1.png")
  png(filename = fileName)
  barplot(as.matrix(states_indices), beside = TRUE, main = "Expression State Concordance Index", col = Ancestry_Colors, xlab = "Expression States", ylab = "Proportion of Dominant Genetic Ancestry", legend.text = legend_titles[,1])
  dev.off()
  fileName <- paste0("Pseudotime1.png")
  png(filename = fileName)
  barplot(as.matrix(times_indices), beside = TRUE, main = "Pseudotemporal Variance Index", col = colors2, xlab = "Genetic Ancestries", ylab = "Standard Deviation of Pseudotime", legend.text = legend_titles[,1])
  dev.off()
  
  #Writes the overall concordance index to a file.
  write.table(final_results,"Res_Table.txt")

  #Performs trajectory analysis using genes most differential in expression evolution.
  genes_of_interest <- data.frame(matrix(nrow =100, ncol = 3))
  colnames(genes_of_interest)<- c("Genetic", "Transcriptomic", "Shared")
  genes_of_interest$Genetic <- Input_Matrix_ordering_genes[1:100]
  
  Input_transcriptomic <- newCellDataSet(as.matrix(expression_matrix),
                                         phenoData = geneticData,
                                         featureData = genes,
                                         expressionFamily=uninormal())
  Input_transcriptomic <- detectGenes(Input_transcriptomic, min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(Input_transcriptomic),
                                      num_cells_expressed >= 10))
  fData(Input_transcriptomic)$use_for_ordering <-
    fData(Input_transcriptomic)$num_cells_expressed > 0.05 * ncol(Input_transcriptomic)
  Input_expressed_genes <-  row.names(subset(fData(Input_transcriptomic),
                                             num_cells_expressed >= 10))
  Input_ordering_genes <-
    Input_expressed_genes
  Input_transcriptomic <- reduceDimension(Input_transcriptomic, max_components = 2, num_dim = 6,
                                          reduction_method = 'tSNE', verbose = T, norm_method = 'none', perplexity = 2)
  Input_transcriptomic <-
    reduceDimension(Input_transcriptomic, method = 'DDRTree', norm_method = 'none')
  Input_transcriptomic <-
    orderCells(Input_transcriptomic)
  diff_test_res <- differentialGeneTest(Input_transcriptomic,
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)")
  genes_of_interest$Transcriptomic <-
    row.names(diff_test_res)[order(diff_test_res$qval)][1:100]
  
  
  
  #Identify and record genes which follow both genetic and transcriptomic evolution
  shared_genes <- intersect(genes_of_interest$Genetic, genes_of_interest$Transcriptomic)
  shared_genes <- c(shared_genes, rep(NA, 100 - length(shared_genes)))
  genes_of_interest$Shared<- shared_genes
  
  
  
  #Write out one file containing genes of interest.
  write.table(genes_of_interest, "genes_of_interest.txt")
  
  #Returns the overall concordance index table.
  return(final_results)
  
}
