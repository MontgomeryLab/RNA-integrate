# All Classes Heatmap Function
all_classes_heatmap = function(cts_heat,cts_names,class_names) {
  
  # Class-separated Matrix and Class Text Matrix Setup
  cts_heat_full = NULL
  cts_heat_text = NULL
  
  # Iteration for Each Class
  for (i in 1:length(levels(as.factor(class_names)))) {
    
    # Isolation of Class-Based Matrix
    cts_heat_class = cts_heat[intersect(names(class_names[which(class_names == levels(as.factor(class_names))[i])]),rownames(cts_heat)),]
    
    # Check for Null Matrices Before Integration Into Heatmap
    if (length(cts_heat_class) > 0) {
      
      # Produce Named Class-Based Matrix
      cts_heat_class = matrix(
        cts_heat_class,
        nrow = length(cts_heat_class)/length(cts_heat[1,]),
        byrow = FALSE,
        dimnames = list(rownames(cts_heat_class),colnames(cts_heat_class))
      )
      
      # Addition of Class Text Entries to Overall Text Matrix
      cts_heat_text = rbind(
        cts_heat_text,
        matrix(
          rep(paste0("Class: ",class_names[rownames(cts_heat_class)]),dim(cts_heat_class)[2]),
          nrow = dim(cts_heat_class)[1],
          ncol = dim(cts_heat_class)[2]
        )
      )
      
      # Common Names Conversion
      rownames(cts_heat_class) = cts_names[rownames(cts_heat_class)]
      
      # Clustering and Re-Ordering for Class-Based Matrices of >1 Row
      if (dim(cts_heat_class)[1] > 1) {
        d = dist(cts_heat_class,method = "euclidean")
        h = hclust(d,method = "complete")
        cts_heat_class = cts_heat_class[h$order,]
      }
      
      # Addition of Class-Based Matrix to Overall Matrix
      cts_heat_full = rbind(cts_heat_full,cts_heat_class)
    }
  }
  
  # Heatmaply Call
  return(
    heatmaply(
      cts_heat_full,
      colors = c("#020532","#03095B","#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F","#000000"),
      limits = c(floor(min(cts_heat)),ceiling(max(cts_heat))),
      showticklabels = c(TRUE,FALSE),
      Colv = FALSE,Rowv = FALSE,
      main = "All Classes Heatmap",
      margins = c(50,50,50,50),
      custom_hovertext = cts_heat_text
    )
  )
}