# Complete Heatmap Function
complete_heatmap = function(cts_heat,cts_names) {
  
  # Common Names Conversion
  rownames(cts_heat) = cts_names[rownames(cts_heat)]
  
  # Clustering and Re-Ordering
  d = dist(cts_heat,method = "euclidean")
  h = hclust(d,method = "complete")
  cts_heat = cts_heat[h$order,]
  
  # Heatmaply Call
  return(
    heatmaply(
      cts_heat,
      colors = c("#020532","#03095B","#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F","#000000"),
      limits = c(floor(min(cts_heat)),ceiling(max(cts_heat))),
      showticklabels = c(TRUE,FALSE),
      Colv = FALSE,Rowv = FALSE,
      main = "Complete Heatmap"
    )
  )
}