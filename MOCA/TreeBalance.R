#Requires as input a phylogenetic tree. Determines a measure for the tree's balance and returns a classification for the tree.

TreeBalance <- function(x){
  phylogeny <- multi2di(x)
  phylogeny <- as.treeshape(phylogeny, "pda")
  phylogeny_results <- maxlik.betasplit(phylogeny, up = 1, remove.outgroup = TRUE, confidence.interval = "none", conf.level = 0.95, size.bootstrap = 1000)
  balance_value <- phylogeny_results$max_lik
  balance_value <- round(balance_value, digits = 1)
  if(balance_value <= -1.5){
    return("Unbalanced")
  }else{
    return("Balanced")
  }
}
