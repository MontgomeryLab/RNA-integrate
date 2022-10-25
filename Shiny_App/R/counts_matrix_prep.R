# Counts Matrix Preparation
counts_matrix_prep = function(counts_matrix) {
  
  # Reconcile Repeated Values
  counts_matrix[,1] = repeated_values_solver(counts_matrix[,1])
  
  # Convert First Column to Row Names
  rownames(counts_matrix) = counts_matrix[,1]
  counts_matrix = counts_matrix[,-1]
  
  # Round Numeric Values
  counts_matrix = round(counts_matrix)
  
  # Convert Numeric Values to Integers
  for (i in 1:length(counts_matrix[1,])) {
    counts_matrix[,i] = as.integer(counts_matrix[,i])
  }
  
  # Object Return
  return(counts_matrix)
  
}