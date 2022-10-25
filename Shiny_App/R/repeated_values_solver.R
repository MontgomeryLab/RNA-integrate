# Repeated Values Solver for Gene Table and Count Matrices/tiny RNA Matrices
repeated_values_solver = function(ids) {
  
  # Check for Repeated Values, Then Add Subscripts
  if (length(ids) == length(unique(ids))) {
    return(ids)
    
  } else {
    
    # Produce List of Subscripts
    ids_reps = rep(NA,length(ids))
    for (i in 1:length(ids)) {
      if (length(which(ids[i] == ids[1:i])) > 1) {
        ids_reps[i] = paste0("_",as.character(length(which(ids[i] == ids[1:i]))))
      } else {
        ids_reps[i] = ""
      }
    }
    
    # Add Subscripts to Repeated Values
    for (i in 1:length(ids)) {
      ids[i] = paste0(ids[i],ids_reps[i])
    }
    
    # Object Return
    return(ids)
    
  }
}