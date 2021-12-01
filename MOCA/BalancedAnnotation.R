#Function which identifies and annotates genetic ancestries from a newick tree.
#Requires a tree in R, a decimal value representing the minimum percent of the total cells to be included in the first two ancestries, a decimal representing the percent of the total cells to be included in each of first two ancestries, and an integer for the number of ancestries to identify.

BalancedAnnotation <- function(tree, min_percent_included, min_percent_included_group, num_groups){
  id_clades5 <- function(p){
    add_clades <- data.frame(matrix(nrow = length(tree$tip.label), ncol = 1))
    old_clades <- data.frame(matrix(nrow = length(tree$tip.label), ncol = 1))
    add_cells <- data.frame(matrix(nrow = length(tree$tip.label), ncol = 1))
    o1 <- edges[which(edges[,1] == p),2]
    o2 <- edges[which(edges[,2] == p),1]
    o3 <- c(o1,o2)
    if(length(o3) == 0){
      return(add_cells)
    }else{
      for(h in 1:length(o3)){
        h1 <- as.integer(o3[h])
        if(h1 <= length(tree$tip.label)){
          add_cells[min(which(is.na(add_cells[,1])==TRUE)),1] <- tree$tip.label[o3[h]]
        }else{
          add_clades[min(which(is.na(add_clades[,1])==TRUE)),1] <- o3[h]
        }
      }
      while(length(which(is.na(add_clades)==FALSE)) > 0){
        z <- add_clades[1,1]
        o1 <- edges[which(edges[,1] == z),2]
        o2 <- edges[which(edges[,2] == z),1]
        o3 <- c(o1,o2)
        old_clades[min(which(is.na(old_clades[,1])==TRUE)),1] <- z
        add_clades<-data.frame(add_clades[-1,])
        for(b in 1:length(o3)){
          if(o3[b] <= length(tree$tip.label)){
            add_cells[min(which(is.na(add_cells[,1])==TRUE)),1] <- tree$tip.label[o3[b]]
          }else{
            if(o3[b] %in% old_clades[,1] | o3[b] %in% add_clades[,1]){
              next
            }else{
              add_clades[min(which(is.na(add_clades[,1])==TRUE)),1] <- o3[b]
            }
          }
        }
      }
      add_cells <- na.omit(add_cells)
      return(add_cells)
    }
  }
  
  clade_combos <- data.frame(matrix(nrow = length(tree$tip.label)*length(tree$tip.label), ncol = 2))

  tree<-makeNodeLabel(tree, method = "number", prefix = "")

  z = 0

  for (clade1 in 1:tree$Nnode) {
    tree3 <- extract.clade(tree,paste(clade1))
    for (clade2 in 1:tree$Nnode) {
      tree4 <- extract.clade(tree,paste(clade2))
      combined_cells <- c(tree3$tip.label, tree4$tip.label)
      if (length(intersect(tree3$tip.label, tree4$tip.label)) == 0 & length(combined_cells) >= min_percent_included*length(tree$tip.label) & length(tree3$tip.label) >= min_percent_included_group*length(tree$tip.label) & length(tree4$tip.label) >= min_percent_included_group*length(tree$tip.label)) {
        z = z+1
        clade_combos[z,1]<- clade1
        clade_combos[z,2]<- clade2
      }
    }
  }
  clade_combos<-na.omit(clade_combos)
  if(nrow(clade_combos) == 0){
    exit <- function() { invokeRestart("abort") }
    print('No combos found, please adjust parameters or consider using UnbalancedAnnotation')
    exit()
    }
  clade_combos <- clade_combos[1:(nrow(clade_combos)/2),]
  branch_size <- data.frame(matrix(nrow = nrow(clade_combos), ncol = 1))
  
  #Determines which of the qualified combinations contains two groups with smallest internal branch length.
  for(a in 1:nrow(clade_combos)){
    sum1 <- sum(extract.clade(tree,paste(clade_combos[a,1]))$edge.length)
    sum2 <- sum(extract.clade(tree,paste(clade_combos[a,2]))$edge.length)
    branch_size[a,1] <- sum1 + sum2
  }
  clade_combos2<-clade_combos[which(branch_size[,1] == min(branch_size[,1])),]
  finalClades<-data.frame(matrix(nrow = num_groups, ncol = 1))
  finalClades[1,1] <- clade_combos2[1,1]
  finalClades[2,1] <- clade_combos2[1,2]
  edges <- data.frame(tree$edge)
  branches <- data.frame(tree$edge.length)
  edges <- edges[which(branches[,1] == 0),]
  v1 <-length(tree$tip.label) + as.integer(finalClades[1,1])
  v2 <-length(tree$tip.label) + as.integer(finalClades[1,1])
  c1 <- id_clades5(v1)
  c1 <- na.omit(c1)
  c2 <- id_clades5(v2)
  c2 <- na.omit(c2)
  cells1_zero <- c1[which(c1[,1] %in% extract.clade(tree, as.character(finalClades[1,1]))$tip.label == FALSE ),1]
  cells2_zero <- c2[which(c2[,1] %in% extract.clade(tree, as.character(finalClades[2,1]))$tip.label == FALSE ),1]
  cells <- c(extract.clade(tree, as.character(finalClades[1,1]))$tip.label, extract.clade(tree, as.character(finalClades[2,1]))$tip.label, cells1_zero, cells2_zero)
  
  #Selects largest remaining clades which are completely unique
  if(num_groups > 2){
    f <- 3
    nodeLabels <- tree$node.label
    nodeLabels <- nodeLabels[-which(nodeLabels == as.character(clade_combos2[1,1]))]
    nodeLabels <- nodeLabels[-which(nodeLabels == as.character(clade_combos2[1,2]))]
    while(f <= num_groups){
      newClades <- data.frame(matrix(nrow = length(nodeLabels), ncol = 1))
      sizes <- data.frame(matrix(nrow = length(nodeLabels), ncol =2))
      sizes[,1] <- nodeLabels
      for (s in 1:length(nodeLabels)) {
        tree5 <- extract.clade(tree, nodeLabels[s])
        if(length(intersect(cells, tree5$tip.label)) == 0){
          newClades[s,1]<- as.character(nodeLabels[s])
          sizes[which(nodeLabels[s] == as.character(sizes[,1])),2] <- length(tree5$tip.label)
        }
      }
      newClades <- na.omit(newClades)
      sizes <- na.omit(sizes)
      #If user specifies more groups than are possible the script automatically does the greatest number possible.
      if(nrow(sizes) == 0){
        print("Number of groups input is not possible. Annotated phylogeny will include greatest number possible")
        finalClades <- na.omit(finalClades)
        break
      }
      newClade <- sizes[which(as.integer(sizes[,2]) == max(as.integer(sizes[,2]))),1]
      
      #If two or more clades are the same size, script selects which has the smallest internal branch lengths.
      if(length(newClade) > 1){
        cladeLength <- data.frame(matrix(nrow = length(newClade), ncol = 1))
        for(h in 1:length(newClade)){
          tree6 <- extract.clade(tree, as.character(newClade[h]))
          cladeLength[h,1] <- sum(tree6$edge.length)
        }
        newOne <- newClade[which(cladeLength[,1] == min(cladeLength[,1]))]
        if(length(newOne) > 0){
          finalClades[f,1] <- newOne[1]
          nodeLabels <- nodeLabels[-which(nodeLabels == newOne[1])]
          cells <- c(cells, extract.clade(tree, newOne[1])$tip.label)
          rm(newOne)
        }else{
          finalClades[f,1] <- newOne
          nodeLabels <- nodeLabels[-which(nodeLabels == newOne)]
          cells <- c(cells, extract.clade(tree, newOne)$tip.label)
          rm(newOne)
        }
      } else {
        finalClades[f,1] <- newClade
        nodeLabels <- nodeLabels[-which(nodeLabels == newClade)]
        cells <- c(cells, extract.clade(tree, newClade)$tip.label)
      }
      f = f + 1
    }
  }
  annotations <- data.frame(matrix(nrow = length(tree$tip.label), ncol = 2))
  colors2 <- data.frame(matrix(nrow = length(tree$tip.label), ncol = 1))
  annotations[,1]<-tree$tip.label
  for(k in 1:nrow(finalClades)){
    annotations[which(annotations[,1] %in% extract.clade(tree,paste(finalClades[k,1]))$tip.label ==TRUE ),2]<-paste(k)
    v1 <-length(tree$tip.label) + as.integer(finalClades[k,1])
    c1 <- id_clades5(v1)
    c1 <- na.omit(c1)
    cells1_zero <- c1[which(c1[,1] %in% extract.clade(tree, as.character(finalClades[1,1]))$tip.label == FALSE ),1]
    annotations[which(annotations[,1] %in% cells1_zero ==TRUE),2] <- paste(k)
  }
  annotations[which(annotations[,1] == "Normal"),2] <- "Normal"
  annotations[which(is.na(annotations[,2]) == TRUE),2] <- "None"
  row.names(annotations) <- annotations[,1]

  #Colors are randomly generated based on the number of clusters. This allows the script to always run without error, however user should specify their color choice for better visuals.
  colors3 <- sample(colors(distinct = FALSE), length(unique(annotations[,2])))

  for(p in 1:nrow(finalClades)){
    colors2[which(annotations[,2] == paste(p)),1] <- colors3[p]
  }
  if ("None" %in% annotations[,2]){
    colors2[which((annotations[,2]) == "None"),1] <- colors3[p+1]
  }
  colors2[which(annotations[,2] == "Normal"),1] <- colors3[length(colors3)]
  print(plot.phylo(tree, tip.color = colors2[,1]))
  legend("topright",sort(unique(annotations[,2])),fill=colors3)
  manual_edit <- as.integer(readline(prompt <- "Enter 1 to manually edit Ancestry IDs, otherwise Enter 0: "))
  if(manual_edit == 1){

    annotations1 <- userDefineClades(annotations)
    for(p in 1:length(unique(annotations1[,2]))){
      colors2[which(annotations1[,2] == sort(unique(annotations1[,2]))[p]),1] <- colors3[p]
    }
    print(plot.phylo(tree, tip.color = colors2[,1]))
    legend("topright",sort(unique(annotations1[,2])),fill=colors3)
  }

  #Produces Annotation File

  if(manual_edit == 1){
    colnames(annotations1) <-c("Cell", "GeneticAncesty")
    write.csv(annotations1, 'BalancedAnnotations.txt', sep = " ", dec = ".",
              row.names = TRUE, col.names = TRUE)
    return(annotations1)
  }else{
    colnames(annotations) <-c("Cell", "GeneticAncesty")
    write.csv(annotations, 'BalancedAnnotations.txt', sep = " ", dec = ".",
              row.names = TRUE, col.names = TRUE)
    return(annotations)
  }
  
}
