# Tiny RNA Matrix Preparation
tiny_rna_matrix_prep = function(tiny_rna_matrix) {
  
  # Separate Common Names
  tiny_rna_matrix[,1] = tiny_rna_matrix[,2]
  for (i in 1:length(tiny_rna_matrix[,1])) {
    if (!is.na(strsplit(tiny_rna_matrix[i,1],", ")[[1]][2])) {
      tiny_rna_matrix[i,2] = strsplit(tiny_rna_matrix[i,1],", ")[[1]][2]
    } else {
      tiny_rna_matrix[i,2] = strsplit(tiny_rna_matrix[i,1],", ")[[1]][1]
    }
    tiny_rna_matrix[i,1] = strsplit(tiny_rna_matrix[i,1],", ")[[1]][1]
  }
  
  # Reconcile Repeated Values
  tiny_rna_matrix[,1] = repeated_values_solver(tiny_rna_matrix[,1])
  tiny_rna_matrix[,2] = repeated_values_solver(tiny_rna_matrix[,2])
  
  # Convert First Column to Row Names
  rownames(tiny_rna_matrix) = tiny_rna_matrix[,1]
  tiny_rna_matrix = tiny_rna_matrix[,-1]
  
  # Round Numeric Values
  tiny_rna_matrix[,-c(1,2)] = round(tiny_rna_matrix[,-c(1,2)])
  
  # Convert Numeric Values to Integers
  for (i in 3:length(tiny_rna_matrix[1,])) {
    tiny_rna_matrix[,i] = as.integer(tiny_rna_matrix[,i])
  }
  
  # Object Return
  return(tiny_rna_matrix)
  
}