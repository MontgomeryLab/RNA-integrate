# Mean Reads Scatter Plot Function
mean_reads_input = function(
    p_value_threshold,
    fold_change_threshold,
    lower_transparency,
    upper_transparency,
    cts_avg,
    cts_names,
    class_names,
    customize_by_significance,
    point_color,
    point_size,
    res,
    group1,
    group2) {
  
  # Significant Genes Classification
  res = res[rownames(cts_avg),]
  res_sig = as.data.frame(subset(res,padj < p_value_threshold & abs(log2FoldChange) > log2(fold_change_threshold)))
  
  # GGplot2 Aesthetic Variables
  x = cts_avg[,group1]
  y = cts_avg[,group2]
  n = cts_names[rownames(cts_avg)]
  
  # Point Significance Variable
  point_sig = rownames(cts_avg) %in% rownames(res_sig)
  
  # Conditional Determination of Size and Color Factor Variables
  
  # Point Colors Assignment
  if (is.null(point_color)) {
    
    # Default Point Colors
    color_class = paste0(c("p < 0.05","p > 0.05")[2-point_sig])
    color_levels = c("#009ACD","#D3D3D3")
    
  } else {
    
    if (customize_by_significance == TRUE) {
      
      # Significance-Customized Point Colors
      color_class = paste0(c("p < 0.05","p > 0.05")[2-point_sig])
      color_class[color_class == "p < 0.05"] = class_names[names(n)[color_class == "p < 0.05"]]
      color_levels = c("p > 0.05" = "#D3D3D3",point_color)
      
    } else {
      
      # Fully Customized Point Colors
      color_class = class_names[names(n)]
      color_levels = point_color
      
    }
    
  }
  
  # Point Sizes Assignment
  if (is.null(point_size)) {
    
    # Default Point Sizes
    size_class = paste0(c("p < 0.05","p > 0.05")[2-point_sig])
    size_levels = c(0.7,0.5)
    
  } else {
    
    # Customized Point Sizes
    size_class = class_names[names(n)]
    size_levels = point_size
    
  }
  
  # Formatting Gene Name Text Variable
  for (i in 1:length(n)) {
    n[i] = paste0("Gene: ",n[i],"<br>","Class: ",class_names[names(n)[i]])
  }
  
  # Final Data Frame Assembly for GGplot
  points_data = data.frame(x,y,color_class,size_class,point_sig,row.names = n)
  
  # Major and Minor Tick Mark Locations
  tick_sep = seq(-4,20,0.8)
  tick_sep = 2^tick_sep
  tick_sep[1:6] = seq(tick_sep[1],tick_sep[6],length.out = 6)
  tick_sep[6:11] = seq(tick_sep[6],tick_sep[11],length.out = 6)
  tick_sep[11:16] = seq(tick_sep[11],tick_sep[16],length.out = 6)
  tick_sep[16:21] = seq(tick_sep[16],tick_sep[21],length.out = 6)
  tick_sep[21:26] = seq(tick_sep[21],tick_sep[26],length.out = 6)
  tick_sep[26:31] = seq(tick_sep[26],tick_sep[31],length.out = 6)
  
  # Tick Mark Labels
  tick_labs = rep("",31)
  tick_labs[seq(1,31,5)] = prettyNum(log2(tick_sep[seq(1,31,5)]),digits=2,format="f")
  
  # GGplott Object Generation
  mean_reads = ggplot(points_data,aes(x,y)) + 
    
    # Point Plotting, Sizing, Coloring, and Text Addition
    geom_point(
      aes(
        color = factor(color_class),
        size = factor(size_class),
        alpha = factor(point_sig),
        shape = factor(20),
        text = paste0(group1,": ",round(x,digits = 3),"<br>",group2,": ",round(y,digits = 3),"<br>",n)
      )
    ) + 
    
    # Color, Size, and Transparency Levels
    scale_color_manual(values = color_levels,name = "Feature Class") + 
    scale_size_manual(values = size_levels,guide = "none") +
    scale_alpha_manual(values = c(lower_transparency,upper_transparency),name = "Significance",labels = c("p > 0.05","p < 0.05")) +
    scale_shape(guide = "none") + 
    
    # Guide Lines
    geom_abline(intercept = 0,slope = 1,color = "grey60") + 
    geom_abline(intercept = 1,slope = 1,color = "grey60") +
    geom_abline(intercept = -1,slope = 1,color = "grey60") + 
    geom_hline(yintercept = log2(10),color = "grey60",linetype = "dashed") + 
    geom_vline(xintercept = log2(10),color = "grey60",linetype = "dashed") + 
    
    # Tick Marks
    scale_y_continuous(
      breaks = log2(tick_sep),
      labels = tick_labs,
      limits = c(-4,max(cts_avg[,c(group1,group2)]))
    ) + 
    
    scale_x_continuous(
      breaks = log2(tick_sep),
      labels = tick_labs,
      limits = c(-4,max(cts_avg[,c(group1,group2)]))
    ) + 
    
    # Title
    ggtitle(paste0(group2,' vs ',group1)) +
    
    # Axis Labeling
    xlab(paste0('Log2 Average reads in ',group1)) +
    ylab(paste0('Log2 Average reads in ',group2)) +
    
    # Theme Specification
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(legend.position = "right")
  
  # Object Return
  return(mean_reads)
}