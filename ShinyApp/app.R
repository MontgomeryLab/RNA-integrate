# DESeq2 Graphics App - Montgomery Lab at Colorado State University
# Dr. Taiowa Montgomery and Spencer Kuhn

## Environment Setup

# Package Installation (if Necessary)
if (require(shiny) == FALSE) {
  install.packages("shiny")
}
if (require(shinycssloaders) == FALSE) {
  install.packages("shinycssloaders")
}
if (require(DT) == FALSE) {
  install.packages("DT")
}
if (require(heatmaply) == FALSE) {
  install.packages("heatmaply")
}
if (require(ggplot2) == FALSE) {
  install.packages("ggplot2")
}
if (require(BiocManager) == FALSE) {
  install.packages("BiocManager")
}
if (require(tximport) == FALSE) {
  BiocManager::install("tximport",force = TRUE,update = FALSE)
}
if (require(DESeq2) == FALSE) {
  BiocManager::install("DESeq2",force = TRUE,update = FALSE)
}

# Load Installed Packages
library(shiny)
library(shinycssloaders)
library(DESeq2)
library(tximport)
library(DT)
library(heatmaply)
library(ggplot2)

# Data-Processing Functions

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

# Table/Plot-Generating Functions

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

# Mean Reads Scatter Plot Function
mean_reads_input = function(p_value_threshold,fold_change_threshold,
                            lower_transparency,upper_transparency,
                            cts_avg,cts_names,class_names,
                            customize_by_significance,
                            point_color,point_size,
                            res,group1,group2) {
  
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
    geom_point(aes(color = factor(color_class),size = factor(size_class),alpha = factor(point_sig),shape = factor(20),
                   text = paste0(group1,": ",round(x,digits = 3),"<br>",
                                 group2,": ",round(y,digits = 3),"<br>",n))) + 
    
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
    scale_y_continuous(breaks = log2(tick_sep),labels = tick_labs,
                       limits = c(-4,max(cts_avg[,c(group1,group2)]))) + 
    scale_x_continuous(breaks = log2(tick_sep),labels = tick_labs,
                       limits = c(-4,max(cts_avg[,c(group1,group2)]))) + 
    
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

# MA Plot Function
ma_input = function(res,group1,group2) {
  return(plotMA(res,ylim=c(-4,4),main = paste0(group2,' vs ',group1)))
}

# Complete Heatmap Function
complete_heatmap = function(cts_heat,cts_names) {
  
  # Common Names Conversion
  rownames(cts_heat) = cts_names[rownames(cts_heat)]
  
  # Clustering and Re-Ordering
  d = dist(cts_heat,method = "euclidean")
  h = hclust(d,method = "complete")
  cts_heat = cts_heat[h$order,]
  
  # Heatmaply Call
  return(heatmaply(cts_heat,
                   colors = c("#020532","#03095B","#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0",
                              "#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F","#000000"),
                   limits = c(floor(min(cts_heat)),ceiling(max(cts_heat))),
                   showticklabels = c(TRUE,FALSE),Colv = FALSE,Rowv = FALSE,
                   main = "Complete Heatmap"))
}

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
      cts_heat_class = matrix(cts_heat_class,nrow = length(cts_heat_class)/length(cts_heat[1,]),byrow = FALSE,
                              dimnames = list(rownames(cts_heat_class),colnames(cts_heat_class)))
      
      # Addition of Class Text Entries to Overall Text Matrix
      cts_heat_text = rbind(cts_heat_text,
                            matrix(rep(paste0("Class: ",class_names[rownames(cts_heat_class)]),dim(cts_heat_class)[2]),
                                   nrow = dim(cts_heat_class)[1],ncol = dim(cts_heat_class)[2]))
      
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
  return(heatmaply(cts_heat_full,
                   colors = c("#020532","#03095B","#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0",
                              "#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F","#000000"),
                   limits = c(floor(min(cts_heat)),ceiling(max(cts_heat))),
                   showticklabels = c(TRUE,FALSE),Colv = FALSE,Rowv = FALSE,
                   main = "All Classes Heatmap",margins = c(50,50,50,50),
                   custom_hovertext = cts_heat_text))
}

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
      hm_class = heatmaply(cts_heat_class,
                           colors = c("#020532","#03095B","#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0",
                                      "#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F","#000000"),
                           limits = c(floor(min(cts_heat)),ceiling(max(cts_heat))),
                           showticklabels = c(TRUE,FALSE),Colv = FALSE,Rowv = FALSE,
                           ylab = heatmap_selected_classes[i],main = "Selected Classes Heatmap",
                           margins = c(50,50,50,50))
      
      # Addition of Heatmaply Object to Subplots List
      heatmap_full = append(heatmap_full,list(hm_class))
    }
  }
  
  # Subplot Compilation
  return(subplot(heatmap_full,nrows = length(heatmap_selected_classes),shareX = TRUE,titleY = TRUE,margin = 0.005))
}

# Define UI for application that draws a histogram
ui = fluidPage(
  
  # Title
  titlePanel("DESeq2 Graphics App, Montgomery Lab at Colorado State University"),
  
  # Body
  mainPanel(
    tabsetPanel(type = "tabs",
                
                # Parameters Tab
                tabPanel("Parameters",
                         
                         # Header
                         h2("Parameters"),
                         br(),
                         
                         # Introduction
                         HTML("<p>Welcome to the Montgomery DESeq2 Graphics and Analysis App, developed by the Montgomery Lab at Colorado State University. 
                              This app is designed to convert tabulated counts files from mRNA or small RNA sequencing Differential Gene Expression (DGE) experiments into user-friendly results tables and plots. 
                              The app does not assist with upstream processes of the RNA-sequencing pipeline, such as filtering or aligment. 
                              Compatible tabulation software includes Salmon, Kallisto, RSEM, and HTSeq. Furthermore, standard count matrices, 
                              (such as the kind produced by featureCounts), as well as tiny RNA matrices (developed by the Montgomery Lab) are acceptable input types.</p>"),
                         br(),
                         
                         # Experiment ID
                         h2("Step One: Experiment ID"),
                         
                         HTML("<p>Enter a brief identifier (no spaces or special characters) for this experimental run (ex. <strong>elegans1</strong>).</p> 
                              <p>The experiment ID will precede file names of saved tables and plots.</p>"),
                         
                         textInput("experiment_id",NULL),
                         
                         # Counts Files Parameters
                         h2("Step Two: Counts Files Parameters"),
                         
                         HTML("<p>First, select the tabulation software used to produce counts files. Options include <strong>RSEM, Salmon, Kallisto, and HTSeq</strong>.
                              for FeatureCounts matrices and general counts matrices for which samples are separated by column, select <strong>Counts Matrix</strong>.</p>
                              <p>Next, if <strong>Counts Matrix</strong> or <strong>tiny RNA</strong> is selected, choose the appropriate matrix file from the working directory.</p>"),
                         br(),
                         
                         fluidRow(
                           column(6,
                                  
                                  # Tabulation Software
                                  selectInput("software_method", "Tabulation Software",
                                              c("RSEM" = "rsem","Salmon" = "salmon","Kallisto" = "kallisto","HTSeq" = "htseq",
                                                "Counts Matrix" = "counts_matrix","tiny RNA" = "tiny_rna"),
                                              selected = "rsem",multiple = FALSE)
                                  
                           ),
                           
                           column(6,
                                  
                                  # Conditional Counts Matrix Selection from Working Directory
                                  conditionalPanel(condition = "input.software_method == 'counts_matrix' | input.software_method == 'tiny_rna'",
                                                   selectInput("selected_counts_matrix","Select counts matrix for this experiment",
                                                               c("select csv file" = "",dir()),selected = character(0),multiple = FALSE)),
                                  
                           )
                         ),
                         
                         # Metadata and Gene Table Files
                         h2("Step Three: Metadata and Gene Table"),
                         
                         HTML("<p>First, select the metadata file corresponding to the experimental run, or choose to generate a new metadata file within the app.
                              To generate a new metadata file, select <strong>generate new</strong> as the metadata method and then follow the prompts that appear.
                              Otherwise, choose <strong>working directory</strong> and choose the appropriate csv file.</p>
                              <p>Next, choose what type of gene table, if any, should be uploaded for the purposes of customizing tables and plots.
                              an uploaded gene table csv can include gene IDs and their <strong>Common Names</strong> (public names), <strong>Class Names</strong> (for small RNA or other classifications), or <strong>Both</strong>.
                              To proceed withouth using a gene table, select <strong>No Table</strong>.</p>"),
                         br(),
                         
                         fluidRow(
                           column(6,
                                  
                                  # Header Text
                                  h2("Metadata Parameters"),
                                  
                                  # Metadata Method
                                  selectInput("metadata_method", "Use existing csv or generate a new csv",
                                              c("select method" = "","working directory" = "working_directory","generate new" = "generate"),selected = character(0)),
                                  
                                  # Conditional Csv Selection from Working Directory
                                  conditionalPanel(condition = "input.metadata_method == 'working_directory'",
                                                   selectInput("selected_metadata","Select csv for this experiment",
                                                               c("select csv file" = "",dir()),selected = character(0),multiple = FALSE)),
                                  
                                  # Condition Number Selection
                                  conditionalPanel(condition = "input.metadata_method == 'generate'",
                                                   radioButtons("condition_number","Select number of experimental conditions",2:10)),
                                  
                                  # Condition Parameter Entries UI
                                  uiOutput("condition_ui")
                                  
                           ),
                           
                           column(6,
                                  
                                  # Header Text
                                  h2("Gene Table Parameters"),
                                  
                                  # Gene Table Method
                                  selectInput("gene_table_method","Choose the type of gene table to import",
                                              c("select method" = "","Full Table" = "full_table","Common Names Only" = "common_names_only",
                                                "Gene Class Only" = "gene_class_only","No Table" = "no_table"),
                                              selected = character(0),multiple = FALSE),
                                  
                                  # Conditional Gene Table Selection from Working Directory
                                  conditionalPanel(condition = "input.gene_table_method != 'no_table' & input.gene_table_method != ''",
                                                   selectInput("selected_gene_table","Select gene table csv for this Experiment",
                                                               c("select csv file" = "",dir()),selected = character(0),multiple = FALSE)),
                                  
                           )
                         ),
                         
                         # Plot Selection
                         h2("Step Four: Plot Selection"),
                         
                         HTML("<p>Choose which plots to render and save for this experiment.</p>"),
                         
                         checkboxInput("results_tables_render","Results Tables",value = TRUE),
                         checkboxInput("pca_render","PCA Plot",value = TRUE),
                         checkboxInput("intra_condition_render","Intra-Condition Scatter Plot",value = TRUE),
                         checkboxInput("mean_reads_render","Mean Reads Scatter Plots",value = TRUE),
                         checkboxInput("ma_render","MA Plots",value = TRUE),
                         checkboxInput("heatmap_render","Heatmap",value = TRUE),
                         
                         # Plot Customization
                         h2("Step Five: Plot Customization"),
                         
                         HTML("<p>Follow the prompts to customize mean reads scatter plots and/or the heatmap. 
                              Without importing a gene table, only the p-value, fold change, and transparency thresholds of the mean reads scatter plots will be customizable.
                              With a gene table imported, both the mean reads scatter plots and the heatmap can be customized by the classes included in the gene table.</p> 
                              <p>A <strong>class parameters</strong> csv can be imported from the working directory or generated within the app to color or size points by class.
                              Furthermore, genes with insignificant p-values can be colored grey regardless of their classification color.</p>
                              <p>Three heatmap types are available if a gene table is imported. The <strong>complete</strong> heatmap clusters all genes in the experiment regardless of class.
                              The <strong>All Classes</strong> heatmap clusters genes within each class and compiles the separate clusters into a single plot.
                              The <strong>Selected Classes</strong> heatmap isolates user-chosen classes into a series of individually-clustered subplots.</p>"),
                         br(),
                         
                         fluidRow(
                           column(7,
                                  
                                  # Mean Reads Scatter Plot Parameters
                                  conditionalPanel(condition = "input.mean_reads_render",
                                                   
                                                   # Header Text
                                                   h2("Mean Reads Scatter Plot Parameters"),
                                                   
                                                   # P-Value Threshold
                                                   sliderInput("p_value_threshold","select significant gene p-value threshold",
                                                               value = 0.05,min = 0.01,max = 0.1,step = 0.01),
                                                   
                                                   # Fold Change Threshold
                                                   sliderInput("fold_change_threshold","select significant gene fold change threshold",
                                                               value = 1.3,min = 1.1,max = 2.0,step = 0.05),
                                                   
                                                   # Lower Transparency Threshold
                                                   sliderInput("lower_transparency","transparency level for insignificant genes",
                                                               value = 0.2,min = 0.1,max = 1.0,step = 0.05),
                                                   
                                                   # Fold Change Threshold
                                                   sliderInput("upper_transparency","transparency level for significant genes",
                                                               value = 0.9,min = 0.1,max = 1.0,step = 0.05)
                                                   
                                                   ),
                                  
                                  # Conditional Class Parameters Method Selection
                                  conditionalPanel(condition = "input.selected_gene_table != '' & input.mean_reads_render",
                                                   checkboxInput("customize_by_class","customize plots by class?",value = FALSE),
                                                   checkboxInput("customize_by_significance","color insignificant genes grey?",value = FALSE)),
                                  
                                  # Conditional Class Parameters Generation Option
                                  conditionalPanel(condition = "input.customize_by_class",
                                                   checkboxInput("generate_class_parameters","generate new class parameters csv?",value = FALSE)),       
                           
                                  # Conditional Csv Selection from Working Directory
                                  conditionalPanel(condition = "input.customize_by_class & !input.generate_class_parameters",
                                                   selectInput("select_class_parameters","Select class parameters csv or leave blank for default csv",
                                                               c("select csv file" = "",dir()),selected = character(0),multiple = FALSE)),
                                  
                                  # Conditional Plot Parameters Class Selection
                                  uiOutput("class_parameters_ui"),
                                  
                                  # Condition Parameter Entries UI
                                  uiOutput("class_ui")
                                  
                           ),
                           
                           column(5,
                                  
                                  # Heatmap Parameters
                                  conditionalPanel(condition = "input.selected_gene_table != '' & input.heatmap_render",
                                                   
                                                   # Header Text
                                                   h2("Heatmap Parameters"),
                                                   
                                                   # Heatmap Type Selection
                                                   selectInput("heatmap_type", "choose heatmap type",
                                                               c("complete" = "complete","class-separated: all classes" = "all_classes",
                                                                 "class-separated: selected classes" = "selected_classes"),selected = "complete")),
                                  
                                  # Conditional Heatmap Class Selection
                                  uiOutput("heatmap_classes_ui"),
                                  
                           )
                         ),
                         
                         # Run Button
                         h2("Step Six: Run Analysis"),
                         
                         HTML("<p>Running a full analysis may take several minutes depending on the size of the data and the level of customization.</p>
                              <p>A <strong>progress bar</strong> will appear in the lower right corner once plot saving has started.</p>"),
                         br(),
                         
                         radioButtons("run","",c("Start Analysis","Set Parameters/Pause Analysis"),selected = "Set Parameters/Pause Analysis",inline = TRUE),
                         br(),
                         fluidRow(column(3,withSpinner(textOutput("all_save"),type = 7,size = 0.5,proxy.height = "25px"))),
                         br()
                         ),
                
                tabPanel("Results Tables",
                         
                         h2("Results Tables"),
                         HTML("Results tables display the gene-wise counts from each replicate in a given contrast, followed by the fold change value of each gene and the associated p-value (adjusted) of the negative binomial hypothesis test conducted by DESeq2. 
                              Lower p-values indicate a lower probability of the null hypothesis that counts between the two conditions are derived from the same distributional parameters. 
                              Furthermore, any common gene names and/or classes from an uploaded gene table will be included in columns beside the Gene ID column. 
                              In addition to csv outputs, these tables are rendered within the app as html widgets with <strong>sorting, searching, and page size customization</strong> features."),
                         br(),
                         uiOutput("results_tables_ui"),
                         br()
                         
                         ),
                
                tabPanel("PCA Plot",
                         
                         h2("PCA Plot"),
                         HTML("<p>The PCA plot displays loadings of the first two principal components for each sample/biological replicate in the experimental design. 
                              Colored by experimental condition, the points of the PCA plot provide a visualization of clustering amongst the samples, both within conditions and across conditions.</p>"),
                         br(),
                         uiOutput("pca_ui"),
                         br()
                         
                         ),
                
                tabPanel("Intra-Condition Scatter Plot",
                         
                         h2("Intra-Condition Scatter Plot"),
                         HTML("<p>The Intra-Condition Scatter Plots display log2 counts between pairs of biological replicates within each condition of the experimental design. 
                              Any counts in the data below 1 are replaced with a value of 1 to simplify the log2 transformation. All replicate pairs in each condition are displayed.</p>"),
                         br(),
                         uiOutput("intra_condition_ui"),
                         br()
                         
                         ),
                
                tabPanel("Mean Reads Scatter Plots",
                         
                         h2("Mean Reads Scatter Plots"),
                         HTML("<p>The Mean Reads Scatter Plots display average log2 counts across biological replicates of experimental condition pairs for a provided contrast. 
                              Average counts with a value of 0 are assigned the value -4 following the log2 transformation. Guidelines are added to assist with visualization. 
                              Statistically significant genes are less transparent than insignificant genes, and in default plots, 
                              they are colored <strong>blue</strong> and sized with a value of <strong>0.5</strong> (insignificant genes are sized with a value of <strong>0.3</strong>). 
                              Axes are scaled to more closely resemble the log2 transformation.Unlabeled tick marks, therefore, do not always represent whole number intervals between labeled tick marks. 
                              In addition to pdf outputs, Mean Reads Scatter Plots can be saved as html widgets with <strong>draggable zooming, panning, and hover text</strong> features using the download button. 
                              To reset the scatter plot axes within the widget, the <strong>home</strong> button in the top right corner of the widget can be pressed.</p>
                              <p>The mean reads scatter plots can be customized by changing the <strong>p-value</strong> and <strong>fold change</strong> thresholds for distinguishing statistically significant genes. 
                              Furthermore, <strong>upper and lower transparency</strong> thresholds can be set for distinguishing statistical significance. 
                              The lower threshold corresponds to insignificant genes. The points of the plot can also be <strong>colored and sized</strong> according to their gene classes as specified by the <strong>gene table</strong>. 
                              If no gene table has been selected/uploaded, only the p-value, fold change, and transparency thresholds will be customizable. 
                              If a gene table has been selected/uploaded, a class parameters csv can be selected from the working directory or built within the app. If used, only points/genes corresponding to classes in the table will be plotted. 
                              If no class parameters csv is selected, every class will be automatically colored from a list of 15 colors, and every point will be sized with a cex value of 0.3. 
                              Insignificant genes will be colored grey if the corresponding checkbox option is selected.</p>"),
                         br(),
                         uiOutput("mean_reads_ui"),
                         br()
                         
                ),
                
                tabPanel("MA Plots",
                         
                         h2("MA Plots"),
                         HTML("<p>The Standard MA Plots are built from the DESeq2 package and display the fold change of a gene over its mean counts value (normalized) for a provided contrast between experimental conditions. 
                              Genes with a statistically significant p-value (p < 0.05) are colored in <strong>blue</strong>.</p>"),
                         br(),
                         uiOutput("ma_ui"),
                         br()
                         
                ),
                
                tabPanel("Heatmap",
                         
                         h2("Heatmap"),
                         HTML("<p>The <strong>Complete</strong> heatmap shows log2 counts across all samples for any genes above a certain mean count threshold of three (meaning an average of eight counts across all samples). 
                              Using the package <strong>heatmaply</strong>, an interactive html widget is rendered with <strong>draggable zooming, panning, and hover text</strong> features. 
                              To reset the heatmap axes, the home button in the top right corner of the widget can be pressed. Darker blue cells indicate lower log2 counts values, while darker red cells indicate higher log2 counts values. 
                              Furthermore, rows (genes) are clustered using the complete hierarchical clustering method in R (hclust) with eudclidean distances.</p>
                              <p>The <strong>Class-Separated: All Classes</strong> heatmap relies on the gene table to map all genes with an associated class. 
                              Such genes are grouped alphabetically by class and clustered hierarchically within their respective group. Hover text indicates the class assigned to each gene.</p>
                              <p>The <strong>Class-Separated: Selected Classes</strong> heatmap produces subplots for each class listed in the associated yaml parameter entry. 
                              Heatmaps for each listed class are scaled according to the limits of the complete heatmap, although cell sizes are scaled according to the number of genes associated with each class.</p>"),
                         br(),
                         uiOutput("heatmap_ui"),
                         br()
                         
                )
                
                )
  )
)

# Define server logic required to draw a histogram
server = function(input,output) {
  
  # Create Output Directory and Save Parameters YAML
  output_directory = observeEvent(input$run == "Start Analysis",{
    
    # Create Output Directory
    dir.create(paste0(getwd(),"/",format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/"))
    
    # Compile List of Parameters
    
    # Conditional Metadata Entry
    if (input$metadata_method == "generate") {
      metadata_param = paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_metadata_generated.csv")
    } else {
      metadata_param = input$selected_metadata
    }
    
    # Conditional Class Parameters Entry
    if (input$generate_class_parameters) {
      class_parameters_param = paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_class_parameters_generated.csv")
    } else {
      class_parameters_param = input$select_class_parameters
    }
    
    # Parameter List
    params = list("experiment_id" = as.character(input$experiment_id),
                  "software_method" = input$software_method,
                  "counts_matrix" = input$selected_counts_matrix,
                  "metadata" = metadata_param,
                  "gene_table_method" = input$gene_table_method,
                  "gene_table" = input$selected_gene_table,
                  "generate_results_tables" = input$results_tables_render,
                  "generate_pca" = input$pca_render,
                  "generate_intra_condition" = input$intra_condition_render,
                  "generate_mean_reads" = input$mean_reads_render,
                  "save_mean_reads_interactive" = FALSE,
                  "p_value_threshold" = input$p_value_threshold,
                  "fold_change_threshold" = input$fold_change_threshold,
                  "lower_transparency" = input$lower_transparency,
                  "upper_transparency" = input$upper_transparency,
                  "customize_by_class" = input$customize_by_class,
                  "customize_by_significance" = input$customize_by_significance,
                  "class_parameters" = class_parameters_param,
                  "generate_ma" = input$ma_render,
                  "generate_heatmap" = input$heatmap_render,
                  "heatmap_type" = input$heatmap_type,
                  "heatmap_selected_classes" = paste0(input$heatmap_selected_classes,collapse = ","))
    
    # Save Parameters YAML
    yaml::write_yaml(params,paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
                                   format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_params.yml"))
  },ignoreInit = TRUE)
  
  # Experimental Setup and Metadata
  
  # UI for Specifying the Parameters of Each Condition
  output$condition_ui = renderUI({
    
    # Looping Condition Parameters UI
    lapply(1:input$condition_number,function(i) {
      fluidRow(
        
        # Condition Name
        conditionalPanel(condition = "input.metadata_method == 'generate'",
                         textInput(paste0("condition_",i),h3(paste0("Enter name of condition ",i)))),
        
        # Control Group Specification
        conditionalPanel(condition = "input.metadata_method == 'generate'",
                         radioButtons(paste0("control_",i),"Is this condition a control group?",
                                      choiceNames = c("yes","no"),choiceValues = c(TRUE,FALSE),selected = FALSE)),
        
        # Working Directory File Selections
        conditionalPanel(condition = "input.metadata_method == 'generate'",
                         selectInput(paste0("files_",i),"Select files that belong to this condition",
                                     c("select files" = "",dir()),selected = character(0),multiple = TRUE))
      )
    })
  })
  
  # Generated Metadata Function from User Inputs
  metadata_csv_generator = eventReactive(input$run == "Start Analysis",{
    
    # Column Setup
    files_list = NULL
    replicates_list = NULL
    conditions_list = NULL
    control_list = NULL
    
    if (input$metadata_method == "generate") {
      for (i in 1:input$condition_number) {
        
        # Files Column
        files_list = c(files_list,eval(parse(text = paste0("input$files_",i))))
        
        # Replicates Column
        replicates_list = c(replicates_list,
                            paste0(rep(eval(parse(text = paste0("input$condition_",i))),
                                       length(eval(parse(text = paste0("input$files_",i))))),
                                   "_",1:length(eval(parse(text = paste0("input$files_",i))))))
        
        # Conditions Column
        conditions_list = c(conditions_list,
                            rep(eval(parse(text = paste0("input$condition_",i))),
                                length(eval(parse(text = paste0("input$files_",i))))))
        
        # Controls Column
        control_list = c(control_list,
                         rep(eval(parse(text = paste0("input$control_",i))),
                             length(eval(parse(text = paste0("input$files_",i))))))
      }
      
      # Object Return
      return(cbind(files_list,replicates_list,conditions_list,control_list))
    }
  })
  
  # Metadata Saving in Working Directory (if applicable)
  observeEvent(input$run == "Start Analysis",{
    if (input$metadata_method == "generate") {
      write.csv(metadata_csv_generator(),
                paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
                       format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_metadata_generated.csv"),
                row.names = FALSE,quote = TRUE)
    }
  })
  
  # Metadata Reading Function
  deseq2_metadata = eventReactive(input$run == "Start Analysis",{
    
    # Working Directory File Inputs, Working Directory Metadata Inputs
    if (input$metadata_method == "working_directory") {
      metadata = read.csv(input$selected_metadata)
      
      # File Upload Inputs or Generated Metadata Inputs
    } else if (input$metadata_method == "generate") {
      metadata = metadata_csv_generator()
    }
    
    # Object Return
    return(metadata)
  })
  
  # Gene Table Import
  gene_table = eventReactive(input$selected_gene_table != "",{
    
    # Conditional Gene Table Import
    if (input$selected_gene_table == "") {
      gene_table = NULL
    } else {
      gene_table = read.csv(input$selected_gene_table)
    }
    
    # Gene Table Values Preparation
    if (is.null(gene_table)) {
      
      # No Gene Table Option
      common_names = NULL
      class_names = NULL
      
    } else {
      if (input$gene_table_method == "full_table") {
        
        # Full Gene Table Option
        common_names = repeated_values_solver(gene_table[,2])
        names(common_names) = repeated_values_solver(gene_table[,1])
        class_names = gene_table[,3]
        names(class_names) = repeated_values_solver(gene_table[,1])
        
      } else if (input$gene_table_method == "common_names_only") {
        
        # Common Names Only Option
        common_names = repeated_values_solver(gene_table[,2])
        names(common_names) = repeated_values_solver(gene_table[,1])
        class_names = NULL
        
      } else if (input$gene_table_method == "gene_class_only") {
        
        # Gene Class Only Option
        common_names = NULL
        class_names = gene_table[,2]
        names(class_names) = repeated_values_solver(gene_table[,1])
      }
    }
    
    # Object Return
    return(list(common_names,class_names))
  })
  
  output$class_parameters_ui = renderUI({
    conditionalPanel(condition = "input.generate_class_parameters & input.selected_gene_table != ''",
                     checkboxGroupInput("mean_reads_selected_classes","select classes to customize",
                                        levels(as.factor(gene_table()[[2]]))))
  })
  
  # UI for Specifying the Mean Reads Scatter Plot Parameters
  output$class_ui = renderUI({
    
    if (!is.null(input$mean_reads_selected_classes[1])) {
      # Looping Condition Parameters UI
      lapply(1:length(input$mean_reads_selected_classes),function(i) {
        fluidRow(
          
          # Color Specification
          conditionalPanel(condition = "input.generate_class_parameters",
                           selectInput(paste0("color_",i),paste0("Color of genes/points for ",input$mean_reads_selected_classes[i]),
                                       c("select a color" = "","light blue" = "#A6CEE3","blue" = "#1F78B4","teal" = "#1B9E77",
                                         "light green" = "#B2DF8A","green" = "#33A02C","lime green" = "#66A61E",
                                         "light red" = "#FB9A99","red" = "#E31A1C","bright pink" = "#E7298A",
                                         "light orange" = "#FDBF6F","orange" = "#FF7F00","burnt orange" = "#D95F02",
                                         "light purple" = "#CAB2D6","purple" = "#6A3D9A","plum" = "#7570B3",
                                         "gold" = "#E6AB02","bronze" = "#A6761D","dark grey" = "#525252","black" = "#000000"),
                                       selected = character(0))),
          
          # Size Specification
          conditionalPanel(condition = "input.generate_class_parameters",
                           sliderInput(paste0("size_",i),paste0("Size of genes/points for class ",input$mean_reads_selected_classes[i]),
                                       value = 0.5,min = 0.1,max = 0.9,step = 0.05))
        )
      })
    }
  })
  
  # Generated Plot Parameters CSV from User Inputs
  class_parameters_generator = eventReactive(input$run == "Start Analysis",{
    
    # Column Setup
    color_list = NULL
    size_list = NULL
    
    if (input$generate_class_parameters & !is.null(input$mean_reads_selected_classes[1])) {
      for (i in 1:length(input$mean_reads_selected_classes)) {
        
        # Color Column
        color_list = c(color_list,eval(parse(text = paste0("input$color_",i))))
                                
        
        # Size Column
        size_list = c(size_list,eval(parse(text = paste0("input$size_",i))))
      }
      
      # Object Return
      return(data.frame("point_class" = input$mean_reads_selected_classes,
                        "point_colors" = color_list,
                        "point_sizes" = size_list))
    }
  })
  
  # Class Parameters Saving in Working Directory (if applicable)
  observeEvent(input$run == "Start Analysis",{
    if (input$generate_class_parameters) {
      write.csv(class_parameters_generator(),
                paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
                       format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_class_parameters_generated.csv"),
                row.names = FALSE,quote = TRUE)
    }
  })
  
  # Plot Parameters Reading Function
  class_parameters = eventReactive(input$run == "Start Analysis",{
    
    # Conditional Import of Working Directory or Generated Plot Parameters
    if (input$generate_class_parameters) {
      
      # Generated Class Parameters Option
      class_parameters = class_parameters_generator()
      
    } else if (input$select_class_parameters != "") {
      
      # Imported Class Parameters Option
      class_parameters = read.csv(input$select_class_parameters)
      
    } else if (input$select_class_parameters == "" & input$customize_by_class) {
      
      # Default Class Parameters Option
      
      # Default Colors List
      class_colors = c("#1F78B4","#1B9E77","#33A02C","#014421","#E6AB02","#FF7F00","#FF9E84","#8B0000",
                       "#E31A1C","#E7298A","#6A3D9A","#907567","#A6761D","#525252","#000000")
      
      # Default Table Generation
      class_parameters = data.frame(point_class = levels(as.factor(gene_table()[[2]])),
                                    point_colors = rep(class_colors,ceiling(length(levels(as.factor(gene_table()[[2]])))/length(class_colors)))[1:length(levels(as.factor(gene_table()[[2]])))],
                                    point_sizes = rep(0.3,length(levels(as.factor(gene_table()[[2]])))))
      
      # Default Table Saving
      write.csv(class_parameters,paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
                                        format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_class_parameters_generated.csv"),row.names = FALSE)
      
    } else {
      
      # Null Class Parameters
      class_parameters = NULL
      
    }
    
    # Plot Parameter Values Preparation
    if (is.null(class_parameters)) {
      
      # No Plot Parameters Option
      point_color = NULL
      point_size = NULL
      
    } else {
      
      # Full Plot Parameters Option
      point_color = class_parameters[,2]
      names(point_color) = class_parameters[,1]
      point_size = class_parameters[,3]
      names(point_size) = class_parameters[,1]
        
    }
    
    # Object Return
    return(list(point_color,point_size))
  })
  
  # Heatmap Class Selection UI
  output$heatmap_classes_ui = renderUI({
    conditionalPanel(condition = "input.heatmap_type == 'selected_classes'",
                     selectInput("heatmap_selected_classes","select classes to isolate",
                                 c("select classes" = "",levels(as.factor(gene_table()[[2]]))),
                                 selected = character(0),multiple = TRUE))
  })
  
  # Sample Table Preparation
  deseq2_sample_table = eventReactive(input$run == "Start Analysis",{
    
    # Data Import and Table Arrangement
    sample_table = data.frame(condition = as.factor(deseq2_metadata()[,3]))
    rownames(sample_table) = deseq2_metadata()[,2]
    
    # Object Return
    return(sample_table)
  })
  
  # DESeq2 Initialization
  deseq2 = eventReactive(input$run == "Start Analysis",{
    
    # Metadata Import
    metadata = deseq2_metadata()
    
    # Sample Table Import
    sample_table = deseq2_sample_table()
    
    # Conditional DESeq Data Import Calls
    if (input$software_method == "htseq") {
      
      # HTSeq Option
      metadata_htseq = data.frame(sampleName = metadata[,2],fileName = metadata[,1],condition = factor(metadata[,3]))
      dds = DESeqDataSetFromHTSeqCount(metadata_htseq,getwd(),~condition)
      
    } else if (input$software_method == "counts_matrix") {
      
      # Counts Matrix/FeatureCounts Option
      counts_matrix = counts_matrix_prep(read.csv(input$selected_counts_matrix))
      colnames(counts_matrix) = rownames(sample_table)
      dds = DESeqDataSetFromMatrix(counts_matrix,sample_table,~condition)
      
    } else if (input$software_method == "tiny_rna") {
      
      # tiny RNA Option
      tiny_rna_matrix = tiny_rna_matrix_prep(read.csv(input$selected_counts_matrix))[,-c(1,2)]
      colnames(tiny_rna_matrix) = rownames(sample_table)
      dds = DESeqDataSetFromMatrix(tiny_rna_matrix,sample_table,~condition)
      
    } else {
      
      # RSEM, Salmon, and Kallisto Options
      txi = tximport(metadata[,1],type = input$software_method,txIn = FALSE,txOut = FALSE)
      txi$length[txi$length == 0] = 1
      dds = DESeqDataSetFromTximport(txi,sample_table,~condition)
      
    }
    
    # DESeq2 Analysis
    dds = DESeq(dds)
    
    # Object Return
    return(dds)
  })
  
  # Contrast List Generation
  deseq2_contrasts = eventReactive(input$run == "Start Analysis",{
    
    # Metadata Import
    metadata = deseq2_metadata()
    
    # Control Condition Status Preparation
    conditions_set = as.factor(metadata[,3])
    control_status = as.logical(metadata[,4])
    names(control_status) = conditions_set
    control_status = control_status[unique(names(control_status))]
    control_group = names(which(control_status == TRUE))
    
    # Combination Generation
    if (all(!control_status) | sum(control_status > 1)) {
      combinations = combn(levels(conditions_set),2)
    } else {
      combinations = data.frame(rbind(rep(control_group,length(levels(conditions_set))-1),
                                      levels(conditions_set)[-which(levels(conditions_set) == control_group)]))
    }
    
    # Column Naming
    for (i in 1:length(combinations[1,])) {
      colnames(combinations)[i] = paste0(combinations[2,i]," vs ",combinations[1,i])
    }
    
    # Object Return
    return(combinations)
  })
  
  # Results List Generation
  deseq2_results = eventReactive(input$run == "Start Analysis",{
    
    # Data Import
    dds = deseq2()
    combinations = deseq2_contrasts()
    
    # List Generation
    res = list()
    for (i in 1:length(combinations[1,])) {
      res = append(res,list(results(dds,contrast=c("condition",combinations[2,i],combinations[1,i]))))
    }
    
    # Object Return
    return(res)
  })
  
  deseq2_cts_names = eventReactive(input$run == "Start Analysis",{
    
    cts = counts(deseq2(),normalized = TRUE)
    common_names = gene_table()[[1]]
    
    # Replace Gene IDs with Common Names (if Applicable)
    cts_names = rownames(cts)
    if (!is.null(common_names)) {
      for (i in 1:length(cts_names)) {
        if (cts_names[i] %in% names(common_names)) {
          cts_names[i] = common_names[cts_names[i]]
        }
      }
    }
    
    # Reconcile Repeated Values
    cts_names = repeated_values_solver(cts_names)
    
    # Assign Gene IDs to Common Names
    names(cts_names) = rownames(cts)
    
    # Object Return
    return(cts_names)
  })

  # Results Tables
  
  # List Generation
  results_tables = eventReactive(input$run == "Start Analysis",{
    
    # Data Import
    sample_table = deseq2_sample_table()
    cts = counts(deseq2(),normalized = TRUE)
    common_names = gene_table()[[1]]
    class_names = gene_table()[[2]]
    res = deseq2_results()
    combinations = deseq2_contrasts()
    
    # Tables List Setup
    tables = list()
    
    # Conditional Rendering with Progress Bar
    if (input$results_tables_render) {
      for (i in 1:length(combinations[1,])) {
        tables = append(tables,list(results_table(sample_table,cts,common_names,class_names,
                                                  res[[i]],combinations[1,i],combinations[2,i])))
      }
    }
    
    # Object Return
    return(tables)
  })
  
  # Results Tables UI
  output$results_tables_ui = renderUI({
    req(input$run == "Start Analysis")
    if (input$results_tables_render) {
      fluidRow(selectInput("results_tables_contrast","Select contrast",colnames(deseq2_contrasts()),selected = 1),dataTableOutput("results_tables_render"))
    } else {
      HTML("To view results tables, check the corresponding box under <strong>Plot Selection</strong> and run the analysis again.")
    }
  })
  
  # Results Tables Rendering by Contrast
  output$results_tables_render = renderDataTable({
    req(input$run == "Start Analysis")
    return(datatable(results_tables()[[which(colnames(deseq2_contrasts()) == input$results_tables_contrast)]]))
  })
  
  # PCA Plotting
  
  # PCA Plot UI
  output$pca_ui = renderUI({
    req(input$run == "Start Analysis")
    if (input$pca_render) {
      withSpinner(plotOutput("pca_plot_render",height = "600px",width = "800px"),type = 1)
    } else {
      HTML("To view the pca plot, check the corresponding box under <strong>Plot Selection</strong> and run the analysis again.")
    }
  })
  
  # PCA Plot Rendering
  output$pca_plot_render = renderPlot({
    req(input$run == "Start Analysis")
    return(pca_input(deseq2()))
  })
  
  # Intra-Condition Scatter Plotting
  
  # Intra-Condition UI
  output$intra_condition_ui = renderUI({
    req(input$run == "Start Analysis")
    if (input$intra_condition_render) {
      withSpinner(plotOutput("intra_condition_plot_render",height = "800px",width = "800px"),type = 1)
    } else {
      HTML("To view the intra-condition scatter plot, check the corresponding box under <strong>Plot Selection</strong> and run the analysis again.")
    }
  })
  
  # Intra-Condition Rendering
  output$intra_condition_plot_render = renderPlot({
    req(input$run == "Start Analysis")
    return(intra_condition_input(deseq2_sample_table(),counts(deseq2(),normalized = TRUE)))
  })
  
  # Mean Reads Scatter Plotting
  
  # Average Counts Table Assembly
  deseq2_cts_avg = eventReactive(input$run == "Start Analysis",{
    
    # Data Import
    metadata = deseq2_metadata()
    cts = counts(deseq2(),normalized = TRUE)
    
    # Object Setup
    cts_avg = NULL
    for (i in levels(as.factor(metadata[,3]))) {
      cts_avg = cbind(cts_avg,log2(rowMeans(cts[,metadata[,2][metadata[,3] == i]])))
    }
    
    # Replace Infinite Values With -4
    cts_avg = replace(cts_avg,cts_avg=='-Inf',-4)
    
    # Assign Appropriate Column Names
    colnames(cts_avg) = levels(as.factor(metadata[,3]))
    
    # Object Return
    return(cts_avg)
  })
  
  # List Generation
  mean_reads_plots = eventReactive(input$run == "Start Analysis",{
    
    # Data Import
    p_value_threshold = input$p_value_threshold
    fold_change_threshold = input$fold_change_threshold
    lower_transparency = input$lower_transparency
    upper_transparency = input$upper_transparency
    cts_avg = deseq2_cts_avg()
    cts_names = deseq2_cts_names()
    class_names = gene_table()[[2]]
    customize_by_significance = input$customize_by_significance
    point_color = class_parameters()[[1]]
    point_size = class_parameters()[[2]]
    res = deseq2_results()
    combinations = deseq2_contrasts()
    
    # Isolate Specific Classes
    if (!is.null(point_color)) {
      cts_avg = cts_avg[intersect(names(class_names[class_names %in% names(point_color)]),names(cts_names)),]
    }
    
    # Plots List Setup
    plots = list()
    
    # Conditional Rendering with Progress Bar
    if (input$mean_reads_render) {
      for (i in 1:length(combinations[1,])) {
        plots = append(plots,list(mean_reads_input(p_value_threshold,fold_change_threshold,
                                                   lower_transparency,upper_transparency,
                                                   cts_avg,cts_names,class_names,
                                                   customize_by_significance,
                                                   point_color,point_size,
                                                   res[[i]],combinations[1,i],combinations[2,i])))
      }
    }
    
    # Object Return
    return(plots)
  })
  
  # Mean Reads Scatter Plots UI
  output$mean_reads_ui = renderUI({
    req(input$run == "Start Analysis")
    if (input$mean_reads_render) {
      fluidRow(selectInput("mean_reads_contrast","Select contrast",colnames(deseq2_contrasts()),selected = 1),
               withSpinner(plotOutput("mean_reads_render",height = "650px",width = "800px"),type = 1),
               downloadButton("download_mean_reads","Download Interactive Version"))
    } else {
      HTML("To view mean reads scatter plots, check the corresponding box under <strong>Plot Selection</strong> and run the analysis again.")
    }
  })
  
  # Plot Rendering by Contrast
  output$mean_reads_render = renderPlot({
    req(input$run == "Start Analysis")
    return(mean_reads_plots()[[which(colnames(deseq2_contrasts()) == input$mean_reads_contrast)]])
  })
  
  # Interactive Plot Download
  output$download_mean_reads = downloadHandler(
    filename = paste0(format(Sys.Date(),format = "%d_%m_%y_"),"Mean_Reads.html"),
    content = function(file) {
      htmlwidgets::saveWidget(ggplotly(mean_reads_plots()[[which(colnames(deseq2_contrasts()) == input$mean_reads_contrast)]],
                                       tooltip = "text"),file = file)
    }
  )
  
  # MA Plotting
  
  # Mean Reads Scatter Plots UI
  output$ma_ui = renderUI({
    req(input$run == "Start Analysis")
    if (input$ma_render) {
      fluidRow(selectInput("ma_contrast","Select contrast",colnames(deseq2_contrasts()),selected = 1),
               withSpinner(plotOutput("ma_render",height = "650px",width = "800px"),type = 1))
    } else {
      HTML("To view ma plots, check the corresponding box under <strong>Plot Selection</strong> and run the analysis again.")
    }
  })
  
  # Plot Rendering by Contrast
  output$ma_render = renderPlot({
    req(input$run == "Start Analysis")
    
    # Data Import
    combinations = deseq2_contrasts()
    res = deseq2_results()
    
    # Plotting Call
    return(ma_input(res[[which(colnames(combinations) == input$ma_contrast)]],
                    combinations[1,which(colnames(combinations) == input$ma_contrast)],
                    combinations[2,which(colnames(combinations) == input$ma_contrast)]))
  })
  
  # Heatmap Plotting
  
  # Heatmap UI
  output$heatmap_ui = renderUI({
    req(input$run == "Start Analysis")
    if (input$heatmap_render) {
      withSpinner(plotlyOutput("heatmap_plot",height = "800px",width = "1000px"),type = 1)
    } else {
      HTML("To view the heatmap, check the corresponding box under <strong>Plot Selection</strong> and run the analysis again.")
    }
  })
  
  # Heatmap Matrix Generation
  deseq2_cts_heat = eventReactive(input$run == "Start Analysis",{
    return(log2(counts(deseq2(),normalized = TRUE))[which(rowMeans(log2(counts(deseq2(),normalized = TRUE))) >= 3),])
  })
  
  # Heatmap Rendering
  output$heatmap_plot = renderPlotly({
    req(input$run == "Start Analysis")
    
    # Data Import
    cts_heat = deseq2_cts_heat()
    cts_names = deseq2_cts_names()
    class_names = gene_table()[[2]]
    
    # Conditional Plotting Calls
    if (input$heatmap_type == "complete") {
      complete_heatmap(cts_heat,cts_names)
    } else if (input$heatmap_type == "all_classes") {
      all_classes_heatmap(cts_heat,cts_names,class_names)
    } else if (input$heatmap_type == "selected_classes") {
      selected_classes_heatmap(cts_heat,cts_names,class_names,input$heatmap_selected_classes)
    }
  })
  
  # Saving Selected Tables and Plots
  output$all_save = renderText({
    req(input$run == "Start Analysis")
    
    # Data Import
    sample_table = deseq2_sample_table()
    cts = counts(deseq2(),normalized = TRUE)
    cts_names = deseq2_cts_names()
    common_names = gene_table()[[1]]
    class_names = gene_table()[[2]]
    res = deseq2_results()
    combinations = deseq2_contrasts()
    
    # Table and Plot Saving
    withProgress(message = "Saving Plots",{
      if (input$results_tables_render) {
        for (i in 1:length(combinations[1,])){
          incProgress(0.15/length(combinations[1,]),message = "saving results tables")
          write.csv(results_tables()[[i]],
                    paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
                           format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_",
                           combinations[2,i],"_vs_",combinations[1,i],"_Results.csv"))
        }
      } else {
        incProgress(0.15,message = "results tables not selected")
      }
      
      if (input$pca_render) {
        incProgress(0.15,message = "saving pca plot")
        pdf(paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
                   format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_PCA_Plot.pdf"),
            title = paste0(input$experiment_id," PCA Plot"),
            height = 6.5,width = 8)
        print(pca_input(deseq2()))
        dev.off()
      } else {
        incProgress(0.15,message = "pca plot not selected")
      }
      
      if (input$intra_condition_render) {
        incProgress(0.15,message = "saving intra-condition scatter plot")
        pdf(paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
                   format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_Intra_Condition.pdf"),
            title = paste0(input$experiment_id," Intra Condition Scatter Plot"),
            height = 6.5,width = 8)
        print(intra_condition_input(sample_table,cts))
        dev.off()
      } else {
        incProgress(0.15,message = "intra-condition scatter plots not selected")
      }
      
      if (input$mean_reads_render) {
        for (i in 1:length(combinations[1,])){
          incProgress(0.15/length(combinations[1,]),message = "saving mean reads scatter plots")
          pdf(paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
                     format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_",
                     combinations[2,i],"_vs_",combinations[1,i],"_Mean_Reads.pdf"),
              title = paste0(input$experiment_id," ",combinations[2,i]," vs ",combinations[1,i]," Mean Reads Scatter Plot"),
              height = 6.5,width = 8)
          print(mean_reads_plots()[[i]])
          dev.off()
        }
      } else {
        incProgress(0.15,message = "mean reads scatter plots not selected")
      }
      
      if (input$ma_render) {
        for (i in 1:length(combinations[1,])){
          incProgress(0.15/length(combinations[1,]),message = "saving ma plots")
          pdf(paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
                     format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_",
                     combinations[2,i],"_vs_",combinations[1,i],"_MA.pdf"),
              title = paste0(input$experiment_id," ",combinations[2,i]," vs ",combinations[1,i]," MA Plot"),
              height = 6.5,width = 8)
          print(ma_input(res[[i]],combinations[1,i],combinations[2,i]))
          dev.off()
        }
      } else {
        incProgress(0.15,message = "ma plots not selected")
      }
      
      if (input$heatmap_render) {
        incProgress(0.15,message = "saving heatmap")
        if (input$heatmap_type == "complete") {
          
          # Complete Heatmap Plot Saving
          htmlwidgets::saveWidget(complete_heatmap(deseq2_cts_heat(),cts_names),
                                  file = paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
                                                format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_Heatmap.html"))
            
        } else if (input$heatmap_type == "all_classes") {
            
          # All Classes Heatmap Plot Saving
          htmlwidgets::saveWidget(all_classes_heatmap(deseq2_cts_heat(),cts_names,class_names),
                                  file = paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
                                                format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_Heatmap.html"))
          
        } else if (input$heatmap_type == "selected_classes") {
            
          # Selected Classes Heatmap Plot Saving
          htmlwidgets::saveWidget(selected_classes_heatmap(deseq2_cts_heat(),cts_names,class_names,input$heatmap_selected_classes),
                                  file = paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
                                                format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_Heatmap.html"))
          
        }
      } else {
        incProgress(0.15,message = "heatmap not selected")
      }
      
      # Delete Accessory Files
      unlink(paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
                    format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_Heatmap_files/"),recursive = TRUE)
    })
    
    return("All Selected Tables and Plots Saved!")
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
