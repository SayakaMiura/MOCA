#Compares the annotations produced by different methods, trees, cutoff values, etc.
#Requires as input a string containing the path to the reference annotations. A data frame containing the reference annotations, and an integer containing the number of alternate analyses to run.

AncestryComparison <- function(reference_tree ,reference_annotations,alternates){
  x <- length(alternates)
  num_alternates <- length(alternates)
  file_type <-substring(reference_tree,nchar(reference_tree) -2)
  removes1 <- which(reference_annotations[,ncol(reference_annotations)] == "None")
  if(file_type == "nwk"){
    ref_tree <-phytools::read.newick(reference_tree)
    removes1 <- which(reference_annotations[,ncol(reference_annotations)] == "None")
    tips1 <- which(ref_tree$tip.label %in% row.names(reference_annotations)[removes1])
    ref_tree <- drop.tip(ref_tree, tips1) #Removes tips from the tree which have no ancestry designation.
  }
  if(length(removes1) > 0){
    reference_annotations <- reference_annotations[-removes1,]
  }

  heat <- data.frame(matrix(nrow = nrow(reference_annotations), ncol = num_alternates + 1))
  row.names(heat) <- row.names(reference_annotations)
  heat[,1] <- reference_annotations[,2]
  for (c in 1:num_alternates) {
    annotations2 <- alternates[c]
    anno_type2 <-substring(annotations2,nchar(annotations2) -2)
    if(anno_type2 == "nwk"){
      moca3 <-phytools::read.newick(annotations2)
      if(exists("ref_tree") == FALSE){
        ref_tree <- moca3
        tips1 <- which(ref_tree$tip.label %in% row.names(reference_annotations) == FALSE)
        ref_tree <- drop.tip(ref_tree, tips1) #Removes tips from the tree which have no ancestry designation.
        reference_annotations <- reference_annotations[which(row.names(reference_annotations) %in% ref_tree$tip.label),]
        }
      
      #If user inputs a path to a tree file, the function calls MOCA's ancestry annotation tools.
      balanced_tree_res <- TreeBalance(moca3)
      balanced_tree_request <- as.integer(readline(prompt<- paste0("Suggested is ", balanced_tree_res, ". Enter 1 for balanced annotation or 2 for unbalanced annotation. ")))
      if(balanced_tree_request == 1){
        num_clades <- as.integer(readline(prompt<- paste0("Enter how many ancestries you would like to identify: ")))
        min_percent_included <- as.numeric(readline(prompt<- paste0("Enter A Decimal Value for Percent of Total Cells to Be in First Two Genetic Ancestries :  ")))
        min_percent_included_group <- as.numeric(readline(prompt<- paste0("Enter A Decimal Value for Percent of Total Cells Required in Each of First Two Genetic Ancestries :  ")))
        annotate_file2 <- BalancedAnnotation(moca3, min_percent_included, min_percent_included_group, num_clades)
      }else{
        num_clades <- as.integer(readline(prompt<- paste0("Enter how many ancestries you would like to identify: ")))
        annotate_file2  <- UnbalancedAnnotation(moca3, num_clades)
      }
      #Removes cells from alternate annotations which are absent in the reference annotation.
      annotate_file2 <- annotate_file2[which(annotate_file2[,1] %in% reference_annotations[,1] == TRUE), 1:2]
      redo <- match(annotate_file2[,1], row.names(heat))
      for(g in 1:length(redo)){
        heat[redo[g],c+1] <- annotate_file2[g,2]
      }
      #Cells which are in the reference, but not in the alternate are labeled as 'removed'.
      heat[which(is.na(heat[,c+1]) == TRUE),c+1] <- "Removed"
      begin <- data.frame(str_locate_all(annotations2,"/"))
      if(nrow(begin) == 0){
        begin[nrow(begin)+1, ncol(begin)] <- 0
      }
      anno_name <- substring(annotations2,begin[nrow(begin), ncol(begin)] + 1,nchar(annotations2)-4)
      colnames(heat)[c+1] <- paste0(anno_name, c)
    }else{
      anno_table <- read.table(paste0(annotations2))
      anno_table <- anno_table[which(anno_table[,1] %in% reference_annotations[,1] == TRUE), 1:2]
      redo <- match(anno_table[,1], row.names(heat))
      anno_table2 <- data.frame(matrix(nrow = nrow(heat), ncol = 1))
      for(g in 1:nrow(heat)){
        if(g > length(redo)){
          break
        }
        anno_table2[redo[g],1] <- anno_table[g,2]
      }
      heat[,c+1]<- anno_table2[,1]
      heat[which(is.na(heat[,c+1]) == TRUE),c+1] <- "Removed"
      begin <- data.frame(str_locate_all(annotations2,"/"))
      if(nrow(begin) == 0){
        begin[nrow(begin)+1, ncol(begin)] <- 0
      }
      anno_name <- substring(annotations2,begin[nrow(begin), ncol(begin)] + 1,nchar(annotations2)-4)
      colnames(heat)[c+1] <- paste0(anno_name)
    }
  }

  groups2 <- as.vector(unique(heat[,1]))
  for(y in 2:ncol(heat)){
    groups <- as.vector(unique(heat[,y]))
    groups2 <- unique(union(groups2,groups))
  }
  colors3 <- sample(colors(distinct = TRUE), length(groups2)+1)
  color_code <- data.frame(matrix(ncol = 2, nrow = length(colors3)))
  color_code[,1] <- colors3
  color_code[,2] <- c(groups2, "Unknown")
  
  colors2 <- data.frame(matrix(nrow = nrow(reference_annotations), ncol = 1))
  for(a in 1:nrow(colors2)){
    group <- heat[a,1]
    rowCode <- which(color_code[,2] == group)
    colors2[a,1]<- color_code[rowCode,1]
  }
  color_code <- color_code[order(color_code[,2]),]
  for(d in 1:ncol(heat)){
    if(d ==1){
      colnames(heat)[1] <- "Reference"
    }else{
      colnames(heat)[d] <- paste0("Annotation ", d)
    }
  }
  heatNum <- data.frame(matrix(nrow = nrow(heat), ncol = ncol(heat)))
  for(d in 1:ncol(heat)){
    if(d ==1){
      colnames(heatNum)[1] <- "Reference"
    }else{
      colnames(heatNum)[d] <- paste0("Annotation ", d)
    }
    for(l in 1:nrow(color_code)){
      heatNum[which(heat[,d] == color_code[l,2]),d] <- l
    }
  }
  row.names(heatNum) <- row.names(heat)
  
  if(exists("ref_tree")){
    px <- ggtree(ref_tree) + geom_tippoint(color = colors2[,1], size = 2)+scale_color_manual(values=c(sort(colors3)))
    print(gheatmap(px, heat, font.size = 2)+scale_fill_manual(values=c(color_code[,1]), name="GeneticAncestry"))
  }else{
    heatmap(as.matrix(heatNum), Rowv = NA, Colv = NA)+scale_fill_manual(values=c(colors3), name="GeneticAncestry")
  }
  #If users enter in more than one alternate annotation, the function computes an aggregate annotation.
  if(num_alternates > 1){
    heat2 <- data.frame(matrix(nrow = nrow(heat), ncol = ncol(heat)))
    for(d in 1:ncol(heat)){
      for(e in 1:nrow(color_code)){
        heat2[which(heat[,d] == color_code[e,2]),d] <- as.integer(e)
      }
    }
    agg_tree <- data.frame(matrix(nrow = nrow(heat2), ncol = 1))
    for(b in 1:nrow(heat2)){
      combined <- heat2[b,]
      combined <- as.integer(as.vector(combined[1,]))
      tab1 <- data.frame(table(combined))
      consensusAncestry <- tab1[which(max(tab1$Freq) == tab1$Freq),1]
      if(length(consensusAncestry) > 1){
        agg_tree[b,1] <- as.integer(which(color_code[,2] == "Unknown"));
      }else{
        consensusAncestry <- droplevels(consensusAncestry)
        agg_tree[b,1] <- as.integer(levels(consensusAncestry))
      }
    }
    
    heat3 <- heat2
    heat3[,ncol(heat3)+1] <- agg_tree[,1]
    colnames(heat3)[ncol(heat3)] <- "Consensus"
    heat4 <- data.frame(matrix(nrow = nrow(heat3), ncol = ncol(heat3)))
    row.names(heat4) <- row.names(heat3)
    
    for(d in 1:ncol(heat4)){
      if(d ==1){
        colnames(heat3)[1] <- "Reference"
      }else if(d < ncol(heat3)){
        colnames(heat3)[d] <- paste0("Annotation ", d)
      }else{
        colnames(heat3)[d] <- "Consensus"
      }
      for(l in 1:nrow(color_code)){
        heat4[which(heat3[,d] == l),d] <- color_code[l,2]
      }
      heat4[,d] <- as.factor(heat4[,d])
    }
    
    #If the annotation is ambiguous, the consensus annotation for the cell is labeled 'Unknown'.

    row.names(heat3)<- row.names(heat)
    row.names(heat4)<- row.names(heat)
    colnames(heat4) <- colnames(heat3)
    if(exists("ref_tree")){
      px <- ggtree(ref_tree) + geom_tippoint(color = colors2[,1], size = 2)+scale_color_manual(values=c(sort(colors3)))
      print(gheatmap(px, heat4, font.size = 2)+scale_fill_manual(values=c(color_code[,1]), name="GeneticAncestry"))
    }else{
      heatmap(as.matrix(heat3), Rowv = NA, Colv = NA)+scale_fill_manual(values=c(colors3), name="GeneticAncestry")
    }
    return(heat4)
  }else{
    return(heat)
  }
}
