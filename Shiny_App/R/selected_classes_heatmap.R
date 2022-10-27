# Selected Classes Heatmap Function
selected_classes_heatmap = function(cts_heat,cts_names,class_names,heatmap_selected_classes) {
  
  heatmap_full = NULL
  
  # Iteration for Each Selected Class
  for (i in 1:length(heatmap_selected_classes)) {
    
    # Isolation of Class-Based Matrix
    cts_heat_class = cts_heat[intersect(names(class_names[which(class_names == heatmap_selected_classes[i])]),rownames(cts_heat)),]
    
    if (length(cts_heat_class) > 0) {
      cts_heat_class = matrix(cts_heat_class,nrow = length(cts_heat_class)/length(cts_heat[1,]),byrow = FALSE,
                              dimnames = list(rownames(cts_heat_class),colnames(cts_heat_class)))
      
      # Common Names Conversion
      rownames(cts_heat_class) = cts_names[rownames(cts_heat_class)]
      
      # Clustering and Re-Ordering for Class-Based Matrices of >1 Row
      if (length(cts_heat_class) > length(cts_heat[1,])) {
        d = dist(cts_heat_class,method = "euclidean")
        h = hclust(d,method = "complete")
        cts_heat_class = cts_heat_class[h$order,]
      }
      
      # Heatmaply Call
      hm_class = heatmaply(
        cts_heat_class,
        colors = c("#020532","#03095B","#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F","#000000"),
        limits = c(floor(min(cts_heat)),ceiling(max(cts_heat))),
        showticklabels = c(TRUE,FALSE),
        Colv = FALSE,Rowv = FALSE,
        ylab = heatmap_selected_classes[i],
        main = "Selected Classes Heatmap",
        margins = c(50,50,50,50)
      )
      
      # Addition of Heatmaply Object to Subplots List
      heatmap_full = append(heatmap_full,list(hm_class))
    }
  }
  
  # Subplot Compilation
  return(subplot(heatmap_full,nrows = length(heatmap_selected_classes),shareX = TRUE,titleY = TRUE,margin = 0.005))
}