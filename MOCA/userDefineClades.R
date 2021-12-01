#As input takes the annotation files in the form of a dataframe. Allows users to edit the names of ancestries and returns an updated annotation dataframe.

userDefineClades <- function(anno_table){
  revised_table <- data.frame(matrix(nrow = nrow(anno_table), ncol = ncol(anno_table)))
  revised_table[,1]<- anno_table[,1]
  row.names(revised_table) <- row.names(anno_table)
  for(f in 1:length(unique(anno_table[,2]))){
    new_name <- readline(prompt <- paste("Enter new name or enter 0 to make no change to ancestry ",  unique(anno_table[,2])[f], " :",   sep = " "))
    if(new_name == 0){
      revised_table[which(anno_table[,2] == unique(anno_table[,2])[f]),2] <- paste0(unique(anno_table[,2])[f])
    }else{
      revised_table[which(anno_table[,2] == unique(anno_table[,2])[f]),2] <- new_name
    }
  }
  return(revised_table)
}
