# Intra-Condition Scatter Plot Function
intra_condition_input = function(sample_table,cts) {
  
  # Establish List of Conditions
  intra_conditions = levels(as.factor(sample_table[,1]))
  
  # Intra-Condition Combinations Preparation
  intra_combinations = NULL
  for (i in intra_conditions) {
    intra_combinations = cbind(intra_combinations,
                               combn(rownames(sample_table)[which(sample_table[,1] == i)],2))
  }
  
  # Gridded Plots Generation
  par(mfrow = c(length(intra_conditions),ceiling(length(intra_combinations)/(2*length(intra_conditions)))),
      mar = c(3,1,1,1),mgp = c(2,0.5,0),pty = "s")
  for (i in 1:length(intra_combinations[1,])){
    plot(log2(replace(cts,cts < 1,1))[,c(intra_combinations[1,i],intra_combinations[2,i])],
         col="lightskyblue3",pch=15,cex=0.2)
  }
}