#Annotates genetic ancestries from a newick tree file. 
#As input it requires a tree in R, and an integer value representing the number of ancestries to identify.

UnbalancedAnnotation <- function(tree, num_groups){
  tree<-makeNodeLabel(tree, method = "number", prefix = "")
  tree2<-makeNodeLabel(tree, method = "number", prefix = "")
  edges <- data.frame(tree$edge)
  edges<- edges[which(tree$edge.length == 0),]
  annotations <- data.frame(matrix(nrow = length(tree2$tip.label), ncol = 2))
  tot_cells <- length(tree$tip.label)
  cells_per_clade <- tot_cells / num_groups #Computes target number of cells for each ancestry
  finalClades <- data.frame(matrix(nrow = num_groups, ncol = 1))
  
  #Selects clades containing amount of cells closest to the target number.
  for(x in 1:num_groups){
    clade <- data.frame(matrix(nrow = tree$Nnode, ncol = 2))
    if(tree$Nnode < 2){
      break
    }
    for(i in 2:tree$Nnode){
      clade[i,2] <- tree$node.label[i]
      clade[i,1] <- length(extract.clade(tree, paste(clade[i,2]))$tip.label)
      clade[i,1] <- abs(cells_per_clade - clade[i,1])
    }
    clade <- na.omit(clade)
    clade_select <- clade[which(clade[,1] == min(clade[,1])),2]
    if(length(clade_select) > 1){
      sizes <-data.frame(matrix(nrow = length(clade_select), ncol = 1))
      for(a in 1:length(clade_select)){
        sizes[a,1]<- sum(extract.clade(tree, paste(clade_select[a]))$edge.length)
      }
      clade_select <- clade_select[which(sizes[,1] == min(sizes[,1]))]
    }
    old_cells <- extract.clade(tree, paste(clade_select))$tip.label
    clade_1 <- as.integer(clade_select)
    for(e in 1:clade_1){
      test <- as.integer(clade_1)-e
      if(test == 0){
        clade_select <- clade_1
        print("The clades have no branch length seperating them")
        break
      }
      new_sum <- sum(extract.clade(tree2, paste(test))$edge.length)
      new_cells <- extract.clade(tree2, paste(test))$tip.label
      difference <-length(setdiff(old_cells, new_cells))
      if(new_sum == 0 & difference == 0){
        clade_select <- as.character(test)
      }
      else{
        break
      }
    }
    tree <- drop.tip(tree, extract.clade(tree, paste(clade_select))$tip.label)
    finalClades[x,1] <- clade_select
    
  }
  
  finalClades <- na.omit(finalClades)
  annotations <- data.frame(matrix(nrow = length(tree2$tip.label), ncol = 2))
  colors2 <- data.frame(matrix(nrow = length(tree2$tip.label), ncol = 1))
  annotations[,1]<-tree2$tip.label
  for(k in 1:nrow(finalClades)){
    if(k == 1){
      annotations[which(annotations[,1] %in% extract.clade(tree2, paste(finalClades[k,1]))$tip.label),2] <- "1" 
    }else{
      annotations[which(annotations[,1] %in% extract.clade(tree2, paste(finalClades[k,1]))$tip.label ==TRUE & is.na(annotations[,2]) == TRUE),2] <- paste(k)
    }
  }
  annotations[which(annotations[,1] == "Normal"),2] <- "Normal"
  annotations[which(is.na(annotations[,2]) == TRUE),2] <- "None"
  row.names(annotations) <- annotations[,1]
  
  #Generates random colors.
  colors3 <- sample(colors(distinct = TRUE), length(unique(annotations[,2])))
  
  for(p in 1:nrow(finalClades)){
    if(p == 1){
      colors2[which(annotations[,2] == paste(p)),1] <- paste0(colors3[p]) 
    }else{
      colors2[which(annotations[,2] == paste(p) & is.na(colors2[,1]) == TRUE),1] <- paste0(colors3[p])
    }
  }
  
  colors2[which(annotations[,2] == "None"),1] <- paste0(colors3[p+1])
  colors2[which(annotations[,1] == "Normal"),1] <- paste0(colors3[length(colors3)])
  
  #Plots the tree with the cells identified by their ancestry.
  print(plot.phylo(tree2, tip.color = colors2[,1]))
  legend("topright",sort(unique(annotations[,2])),fill=colors3)
  manual_edit <- as.integer(readline(prompt <- "Enter 1 to manually edit Ancestry IDs, otherwise Enter 0: "))
  
  #Allows users to edit the name of the ancestries.
  if(manual_edit == 1){
    annotations <- userDefineClades(annotations)
    for(p in 1:length(unique(annotations[,2]))){
      colors2[which(annotations[,2] == sort(unique(annotations[,2]))[p]),1] <- colors3[p]
    }
    print(plot.phylo(tree2, tip.color = colors2[,1]))
    legend("topright",sort(unique(annotations[,2])),fill=colors3)
  }
  
  #Produces file which can be used for MOCA
  colnames(annotations) <-c("Cell", "GeneticAncesty")
  row.names(annotations) <- annotations[,1]
  write.csv(annotations, 'UnbalancedAnnotation.txt', sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
  return(annotations)
}
