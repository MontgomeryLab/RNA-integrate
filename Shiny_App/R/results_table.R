# Results Table Function
results_table = function(sample_table,cts,common_names,class_names,res,group1,group2) {
  
  # Results Object Ordering by P-Value
  res_ordered = res[order(res$pvalue),]
  res_ordered = data.frame(res_ordered)
  
  # Isolating Fold Change and Adjusted P-Value Columns
  res_ordered = res_ordered[,c(2,6)]
  
  # Add Counts Columns for the Given Contrast
  res_ordered = cbind(cts[rownames(res_ordered),
                          rownames(sample_table)[which(sample_table == as.character(group1))]],
                      cts[rownames(res_ordered),
                          rownames(sample_table)[which(sample_table == as.character(group2))]],
                      res_ordered)
  
  # Fold Change Column Transformation
  res_ordered$log2FoldChange[which(res_ordered$log2FoldChange > 0)] = 
    2^res_ordered$log2FoldChange[which(res_ordered$log2FoldChange > 0)]
  res_ordered$log2FoldChange[which(res_ordered$log2FoldChange < 0)] = 
    -1/(2^res_ordered$log2FoldChange[which(res_ordered$log2FoldChange < 0)])
  
  # Fold Change Column Renaming
  colnames(res_ordered)[length(colnames(res_ordered))-1] = "Fold_Change"
  
  # Values Rounding to Three Decimal Places
  res_ordered = round(res_ordered,digits = 3)
  
  # Conditional Addition of Class Names Column
  if (!is.null(class_names)) {
    res_ordered = cbind(class_names[rownames(res_ordered)],res_ordered)
    colnames(res_ordered)[1] = "Class"
  }
  
  # Conditional Addition of Common Names Column
  if (!is.null(common_names)) {
    res_ordered = cbind(common_names[rownames(res_ordered)],res_ordered)
    colnames(res_ordered)[1] = "Common_Name"
  }
  
  # Object Return
  return(res_ordered)
}