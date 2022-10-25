# PCA Plot Function
pca_input = function(dds) {
  
  # Create PCA Plot
  pca = plotPCA(vst(dds,blind = FALSE),intgroup="condition")
  
  # Add Plot Title
  pca = ggplot_add(ggtitle("PCA Plot"),pca)
  
  # Revert to Classic Theme
  pca = ggplot_add(theme_classic() + theme(plot.title = element_text(hjust = 0.5)),pca)
  
  # Object Return
  return(pca)
}