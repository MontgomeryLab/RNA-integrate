# Montgomery R Shiny App for mRNA and small RNA DGE Experiments

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Configuration](#configuration)
  - [User Inputs](#user-inputs)
  - [Metadata CSV](#metadata-csv)
  - [Gene Table](#gene-table)
  - [Class Parameters CSV](#class-parameters-csv)
- [Output](#output)
  - [Results Tables](#results-tables)
  - [PCA Plot](#pca-plot)
  - [Intra-Condition Scatter Plot](#intra-condition-scatter-plot)
  - [Mean Reads Scatter Plots](#mean-reads-scatter-plots)
  - [MA Plots](#ma-plots)
  - [Heatmap](#heatmap)
- [Authors](#authors)


## Introduction

Welcome to the R Shiny DESeq2 graphics app by the Montgomery Lab at Colorado State University. This pipeline uses tabulated counts files from mRNA or small RNA sequencing experiments to produce differential gene expression (DGE) analysis results and visualizations. In addition to the interactive app interface, the pipeline saves csv, pdf, or html files to an output folder within the local directory off the app script as well as a YAML describing the user input parameters used for that experimental run.

The app can be launched within R studio using the following command:

`shiny::runApp("Montgomery_DESeq2_App.R")`

Or in the command line as such:

`R -e 'shiny::runApp("Montgomery_DESeq2_App.R")'`

## Installation

This app has been developed and tested with R version <strong>4.1.2</strong>. Several packages are required in order to launch and run the app. These packages will be installed if necessary within the markdown code, but note that some packages may be difficult to install if version inconsistencies and dependency issues arise. For troubleshooting, we suggest trying to install any troublesome packages individually within R studio. We recommend using BioConductor for installing <strong>tximport and DESeq2</strong>. To install a package individually, use the following command in R:

`install.packages("Package Name")`

Or use the command line as such:

`R -e 'install.packages("Package Name")'`

The following packages are used within the pipeline:

- shiny: used for producing interactive and reactive app interface
  - Winston Chang, Joe Cheng, JJ Allaire, Carson Sievert, Barret Schloerke, Yihui Xie, Jeff Allen, Jonathan McPherson, Alan Dipert and Barbara Borges (2021). shiny: Web Application Framework for R. R package version 1.7.1. https://CRAN.R-project.org/package=shiny
- shinycssloaders: used to display loading icons while tables and plots are produced
  - Andras Sali and Dean Attali (2020). shinycssloaders: Add Loading Animations to a 'shiny' Output While It's Recalculating. R package version 1.0.0. https://CRAN.R-project.org/package=shinycssloaders
- BiocManager: BioConductor package used to install tximport and DESeq2
  - Martin Morgan (2021). BiocManager: Access the Bioconductor Project Package Repository. R package version 1.30.16. https://CRAN.R-project.org/package=BiocManager
- tximport: used to format tabulated counts files into usable forms for DESeq2
  - Soneson C, Love MI, Robinson MD (2015). “Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences.” F1000Research, 4. doi: 10.12688/f1000research.7563.1.
- DESeq2: negative binomial statistical testing and data analysis for DGE experiments
  - Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)
- DT: produces HTML widgets for results tables
  - Yihui Xie, Joe Cheng and Xianying Tan (2021). DT: A Wrapper of the JavaScript Library 'DataTables'. R package version 0.20. https://CRAN.R-project.org/package=DT
- ggplot2: dynamic graphics package for PCA Plots and Mean Reads Scatter Plots
  - Hadley Wickham (2016). ggplot2: Elegant Graphics for Data Analysis. R package version 3.3.5. https://ggplot2.tidyverse.org
- heatmaply: heatmap package based on plotly, a package for producing interactive graphics
  - Tal Galili, Alan O'Callaghan, Jonathan Sidi and Carson Sievert (2017). heatmaply: an R package for creating interactive cluster heatmaps for online publishing. R package version 1.3.0. https://doi.org/10.1093/bioinformatics/btx657

## Configuration

### User Inputs

A <strong>Parameters</strong> tab, shown upon launching the app, walks users through the input parameters necessary for an experimental run. The following list describes the options for each input:

<strong>Initial Setup</strong>
- Experiment ID: a short string with no spaces or special characters that distinguishes a particular experimental run. The experiment id will be prepended in the names of the output directory and saved files. 
- Tabulation Software: the tabulation software used to produce counts files or a counts matrix. Options include <strong>Counts Matrix</strong> (FeatureCounts), <strong>tiny RNA</strong> (developed by the Montgomery Lab - https://github.com/MontgomeryLab/tinyRNA), <strong>RSEM</strong>, <strong>HTSeq</strong>, <strong>Salmon</strong>, and <strong>Kallisto</strong>. 
  - Counts Matrix: selection from the working directory of the counts matrix file if the <strong>Tabulation Software</strong> selection is <strong>Counts Matrix</strong> or <strong>tiny RNA</strong>.

<strong>Metadata</strong>
- Metadata Method: describes whether to use an existing metadata file in the working directory or generate a new metadata file within the app. See [Metadata CSV](#metadata-csv).
  - If <strong>Working Directory</strong> is chosen as the <strong>Metadata Method</strong>:
    - Selected Metadata: selection from the working directory of the metadata csv.
  - If <strong>Generate New</strong> is chosen as the <strong>Metadata Method</strong>:
    - Condition Number: radio button selection of the number of experimental conditions if <strong>Generate new</strong> is chosen as the <strong>Metadata Method</strong>.
    - Condition Name, Control Status, and Condition Files: prompts users to enter information about each experimental condition if <strong>Generate New</strong> is chosen as the <strong>Metadata Method</strong>.

<strong>Gene Table</strong>
- Gene Table Method: describes the type of gene table, if any, that will be imported for the experiment. options include <strong>Full Table</strong>, <strong>Common Names Only</strong>, <strong>Gene Class Only</strong>, and <strong>No Table</strong>. see [Gene Table](#gene-table).
  - Gene Table: selection from the working directory the gene table, unless <strong>No Table</strong> is chosen as the <strong>Gene Table Method</strong>.

<strong>Plot Parameters</strong>
- Plot Selection: series of checkbox selection describing whether or not to render and save each type of plot, including the [Results Tables](#results-tables), [pca plot](#pca-plot), [intra-condition scatter plot](#intra-condition-scatter-plot), [mean reads scatter plots](#mean-reads-scatter-plots), [ma plots](#ma-plots), and [heatmap](#heatmap)

<strong>Mean Reads Scatter Plots</strong>
- P-Value Threshold: numeric slider input between 0.01 and 0.1 for classifying significant genes by p-value.
- Fold Change Threshold: numeric slider input between 1.1 and 2.0 for classifying significant genes by fold change.
- Lower Transparency: numeric slider input between 0.1 and 1.0 describing the transparency level (alpha level) for insignificant genes. Values closer to 0 produce more transparent points. 
- Upper Transparency: numeric slider input between 0.1 and 1.0 describing the transparency level (alpha level) for significant genes. Values closer to 0 produce more transparent points. 
- Customize By Class: checkbox input describing whether mean reads scatter plots should be customized by their classes, as listed in the gene table. Note that a gene table must be imported in order to customize by class. see [Class Parameters CSV](#class-parameters-csv).
  - If <strong>Customize By Class</strong> is checked:
    - Customize By Significance: checkbox input describing whether points representing insignificant genes should be colored grey. Otherwise, all points will be colored, and their p-value significance will be distinguished only by transparency.
    - Generate Class Parameters: checkbox input describing whether a class parameters csv should be generated within the app.
      - If <strong>Generate Class Parameters</strong> is checked:
        - Class Selection: selection of classes to plot genes for if the <strong>Generate Class Parameters</strong> box is checked. Class options are derived from the factor levels of the class column in the gene table.
        - Point Color and Point Size: prompts users to enter information about each selected class if the <strong>Generate Class Parameters</strong> box is checked.
      - If <strong>Generate Class Parameters</strong> is not checked
        - Select Class Parameters: selection from working directory of the class parameters csv. If left blank, a default class parameters csv will be generated including all classes listed in the gene table.

<strong>Heatmap</strong>
- Heatmap Type: describes what type of heatmap is to be rendered and saved. Options include <strong>Complete</strong>, <strong>All Classes</strong>, and <strong>Selected Slasses</strong>. See [Heatmap](#heatmap). If no gene table is uploaded, a complete heatmap will be rendered, and the input widget will not be shown. 
  - Selected Classes: multiple selction list of classes for which heatmaps will be generated if the heatmap type is <strong>Selected Classes</strong>. Class options are derived from the factor levels of the class column in the gene table.

### Metadata CSV

The metadata csv summarizes the experimental design and provides information to the DESeq2 analysis. It can be imported from the working directory or generated within the app. The first column of the metadata csv should include counts file paths corresponding to each sample. If the <strong>Software Method</strong> being used is <strong>Counts Matrix</strong> or <strong>tiny RNA</strong> then the first column can consist of either arbitrary/empty strings or name of the counts matrix repeated on each row. The second column of the metadata csv describes specific replicates for each sample in the experiment (ex. "WT_1","prg-1_2", etc.). The third column describes the condition associated with each of the samples/replicates in the previous column (ex. simply "WT" or "prg-1"). The fourth column of the csv consists of logical values (TRUE/FALSE) describing whether each sample belongs to a control group/condition. An example of the metadata csv can be viewed [here](metadata.csv), and a preview is shown below.

| files            | replicates | condition | control_condition |
|------------------|------------|-----------|-------------------|
| counts1.results  | N2_1       | N2        | TRUE              |
| counts2.results  | prg-1      | prg-1     | FALSE             |
| ...              | ...        | ...       | ...               |

### Gene Table

The gene table is a key component of customizing the pipeline's data analysis and visualization, although it is optional for individual runs of the markdown pipeline. To disable the gene table, select <strong>no_table</strong> as the gene table method. Choosing <strong>common_names_only</strong> will eliminate the class-based customization, and choosing <strong>class_names_only</strong> will eliminate common name substitution. Choosing <strong>full_table</strong> as the gene table method will maintain both class-based customization and common name substitution. The first column of the gene table should consist of any number of gene IDs identical to corresponding IDs from the tabulated counts files or matrix. The second column of the gene table should consist of any common names that should replace a corresponding gene ID in results tables, mean reads scatter plots, and heatmaps. To avoid replacing the gene ID for a particular gene, the gene ID should be repeated in the second column. The third column of the gene table should consist of the class to which each gene ID belongs. Class names are used to characterize and color visualized genes for the mean reads scatter plots and heatmaps (if <strong>selected_classes</strong> is chosen as the heatmap type). A fourth column of the gene table may be included that lists target gene ids, separated by commas, but this fourth column us only usable in the [integrative pipeline](../Integrative_Pipeline/). They are also included in the results tables. A template gene table can be viewed [here](template_gene_table.csv), and a preview is shown below.

| Gene_ID      | Public_Name   | Class         | Targets        |
|--------------|---------------|---------------|----------------|
| MIMAT0000002 | lin-4-5p      | miRNA         | WBGene00003003 |
| ...          | ...           | ...           | ...            |

### Class Parameters CSV

The class parameters csv is used to customize [Mean Reads Scatter Plots](#mean-reads-scatter-plots). A gene table must be used in the experiment in order for the class parameters csv to be usable. If used, only points/genes corresponding to classes in the table will be plotted. The first column of the csv should consist of any number of class names identical to those found in the gene table. Points (genes) corresponding to these classes can be colored and re-sized in the mean reads scatter plots. The second column of the table should consist of hex color values corresponding to each class (ex. "#D95F02"). The third column of the table should consist of numeric values between 0.1 and 1.0 that describe the cex point sizes for genes of the corresponding class (default size is 0.3). The default class parameters table will color all classes in the gene table according to a standard list of 15 colors, and all classes will be given a size value 0.3. An example plot parameters csv can be seen [here](class_parameters.csv), and a preview is shown below. 

| point_class   | point_colors | point_sizes  |
|---------------|--------------|--------------|
| ALG           | #D95F02      | 0.7          |
| ...           | ...          | ...          |

## Output

### Results Tables

Results tables display the gene-wise counts from each replicate in a given contrast, followed by the fold change value of each gene and the associated p-value (adjusted) of the negative binomial hypothesis test conducted by DESeq2. Lower p-values indicate a lower probability of the null hypothesis that counts between the two conditions are derived from the same distributional parameters. Furthermore, any common gene names and/or classes from an uploaded gene table will be included in columns beside the Gene ID column. In addition to csv outputs, these tables are rendered within the app as html widgets with <strong>sorting, searching, and page size customization</strong> features.

<img src=../Demo/mRNA_Demo_Results/mrna_demo_piwi_vs_wt_Results.png width="800" height="500">

### PCA Plot

The PCA plot displays loadings of the first two principal components for each sample/biological replicate in the experimental design. Colored by experimental condition, the points of the PCA plot provide a visualization of clustering amongst the samples, both within conditions and across conditions.

<img src=../Demo/mRNA_Demo_Results/mrna_demo_PCA_Plot.jpg width="600" height="500">

### Intra-Condition Scatter Plot

The Intra-Condition Scatter Plots display log2 counts between pairs of biological replicates within each condition of the experimental design. Any counts in the data below 1 are replaced with a value of 1 to simplify the log2 transformation. All replicate pairs in each condition are displayed.

<img src=../Demo/mRNA_Demo_Results/mrna_demo_Intra_Condition.jpg width="700" height="600">

### Mean Reads Scatter Plots

The Mean Reads Scatter Plots display average log2 counts across biological replicates of experimental condition pairs for a provided contrast. Average counts with a value of 0 are assigned the value -4 following the log2 transformation. Guidelines are added to assist with visualization. Statistically significant genes are less transparent than insignificant genes, and in default plots, they are colored <strong>blue</strong> and sized with a value of <strong>0.5</strong> (insignificant genes are colored grey and sized with a value of <strong>0.3</strong>). Axes are scaled to more closely resemble the log2 transformation. Unlabeled tick marks, therefore, do not always represent whole number intervals between labeled tick marks. In addition to pdf outputs, Mean Reads Scatter Plots can be downloaded as html widgets with <strong>draggable zooming, panning, and hover text</strong> features. To reset the scatter plot axes within the widget, the <strong>home</strong> button in the top right corner of the widget can be pressed.

The mean reads scatter plots can be customized by changing the <strong>p-value</strong> and <strong>fold change</strong> thresholds for distinguishing statistically significant genes. Furthermore, <strong>upper and lower transparency</strong> thresholds can be set for distinguishing statistical significance. The lower threshold corresponds to insignificant genes. The points of the plot can also be <strong>colored and sized</strong> according to their gene classes as specified by the <strong>gene table</strong>. If no gene table has been selected/uploaded, only the p-value, fold change, and transparency thresholds will be customizable. If a gene table has been selected/uploaded, a class parameters csv can be selected from the working directory or produced within the app. If used, only points/genes corresponding to classes in the table will be plotted. If no class parameters csv is selected, every class will be automatically colored from a default list of 15 colors, and every point will be sized with a cex value of 0.3. Insignificant genes will be colored grey if the <strong>Customize By Significance</strong> box is checked.

<img src=../Demo/mRNA_Demo_Results/mrna_demo_piwi_vs_wt_Mean_Reads.jpg width="700" height="600">

### MA Plots

The Standard MA Plots are built from the DESeq2 package and display the fold change of a gene over its mean counts value (normalized) for a provided contrast between experimental conditions. Genes with a statistically significant p-value (p < 0.05) are colored in <strong>blue</strong>.

<img src=../Demo/mRNA_Demo_Results/mrna_demo_piwi_vs_wt_MA_Plot.jpg width="600" height="600">

### Heatmap

The <strong>Complete</strong> heatmap shows log2 counts across all samples for any genes above a certain mean count threshold of three (meaning an average of eight counts across all samples). Using the package <strong>heatmaply</strong>, an interactive html widget is rendered with <strong>draggable zooming, panning, and hover text</strong> features. To reset the heatmap axes, the home button in the top right corner of the widget can be pressed. Darker blue cells indicate lower log2 counts values, while darker red cells indicate higher log2 counts values. Furthermore, rows (genes) are clustered using the complete hierarchical clustering method in R (hclust) with eudclidean distances.

<img src=../Demo/mRNA_Demo_Results/mrna_demo_complete_heatmap.png width="800" height="500">

The <strong>All Classes</strong> heatmap relies on the gene table to map all genes with an associated class. Such genes are grouped alphabetically by class and clustered hierarchically within their respective group. Hover text indicates the class assigned to each gene.

<img src=../Demo/mRNA_Demo_Results/mrna_demo_all_classes_heatmap.png width="800" height="500">

The <strong>Selected Classes</strong> heatmap produces subplots for each class listed in the <strong>Heatmap Selected Classes</strong> input. Heatmaps for each listed class are scaled according to the limits of the complete heatmap, although cell sizes are scaled according to the number of genes associated with each class. 

<img src=../Demo/mRNA_Demo_Results/mrna_demo_selected_classes_heatmap.png width="800" height="500">

## Authors

* **Dr. Taiowa Montgomery** - 05/2021-present - Colorado State University - [taimontgomery](https://github.com/taimontgomery)
* **Spencer Kuhn** - 06/2021-present - Colorado State University - [smcguirekuhn](https://github.com/smcguirekuhn)
