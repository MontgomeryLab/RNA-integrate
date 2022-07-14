## RNA-Integrate: Integrative Pipeline for Parallel mRNA and small RNA Data

:warning: **Under Development** :warning:

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Configuration Files](#configuration-files)
  - [Integrative Parameters YAML](#integrative-parameters-yaml)
  - [Metadata CSV](#metadata-csv)
  - [Gene Table](#gene-table)
  - [Class Parameters CSV](#class-parameters-csv)
- [Output](#output)
  - [Integrative Results Tables](#integrative-results-tables)
  - [Galaxy Plots](#galaxy-plots)
  - [Slope Plots](#slope-plots)
- [Authors](#authors)


## Introduction

Welcome to the Montgomery DESeq2 Pipeline for small RNA and mRNA integrative experiments in R Markdown. This pipeline allows for the analysis and visualization of data describing the interactions between small RNA and mRNA components of the same experiment. First, separate runs of an R Markdown pipeline for individual experiments are executed, once for the small RNA data and once for the mRNA data. Then, information found in an imported gene table describing targets for each small RNA or mRNA gene is used to produce [integrative results tables](#integrative-results-tables), [galaxy scatter plots](#galaxy-plots), and [slope plots](#slope-plots). To see more information about the individual experiment pipeline, click [here](../Markdown_Pipeline/). For each individual experiment run, as well as for the integrative data, knitted R Markdown reports will be produced in addition to plots, tables, and a copy of the yaml parameters file, all saved in separate directories within the markdown pipeline's local directory. 

The pipeline can be executed within R studio using the following command:

`rmarkdown::render("RNA_Integrate.Rmd",params = yaml::read_yaml("int_params.yml"))`

Or in the command line as such:

`R -e 'rmarkdown::render("RNA_Integrate.Rmd",params = yaml::read_yaml("int_params.yml"))'`

Execution of the pipeline can take up to fifteen minutes depending on the complexity of the experiment and/or the level of customization.

## Installation

This pipeline has been developed and tested with R version <strong>4.1.2</strong>. Several packages are required in order to execute the pipeline. These packages will be installed if necessary within the markdown code, but note that some packages may be difficult to install if version inconsistencies and dependency issues arise. For troubleshooting, we suggest trying to install any troublesome packages individually within R studio. We recommend using BioConductor for installing packages like <strong>tximport and DESeq2</strong>. To install a package individually, use the following command in R:

`install.packages("Package Name")`

Or use the command line as such:

`R -e 'install.packages("Package Name")'`

The following packages are used within the pipeline:

- knitr: used for knitting markdown files
  - Xie Y (2021). knitr: A General-Purpose Package for Dynamic Report Generation in R. R package version 1.37, https://yihui.org/knitr/.
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

Furthermore, RColorBrewer was used to inspire baseline color palettes for the heatmap and the galaxy scatter plots. RColorBrewer is cited below:
- Erich Neuwirth (2014). RColorBrewer: ColorBrewer Palettes. R package version 1.1-2. https://CRAN.R-project.org/package=RColorBrewer

## Configuration Files

### Integrative Parameters YAML

The integrative parameters YAML file contains the details of the experimental configuration, including a number of plot customization options, for the individual mRNA and small RNA runs as well as the integrative run. A template int_params.yml file can be seen [here](template_int_params.yml). The following parameters are read from the YAML file during the pipeline's execution:

#### Setup and Metadata Parameters

- mrna_experiment_id/small_rna_experiment_id/integrative_experiment_id: short strings with no spaces or special characters that distinguish particular experimental runs for the mRNA, small RNA, and integrative data. The experiment id will be prepended to the names of the output directory and saved files. 
- mrna_software_method/small_rna_software_method: the tabulation software used to produce counts files or a counts matrix for the mRNA and small RNA data. Options include "counts_matrix" (FeatureCounts), "tiny_rna" (developed by the Montgomery Lab - https://github.com/MontgomeryLab/tinyRNA), "rsem", "htseq", "salmon", or "kallisto". 
- mrna_counts_matrix/small_rna_counts_matrix: a string describing the file path for the counts matrix or the tiny RNA matrix, if applicable, for the mRNA and/or small RNA data. Leave as an empty string if count matrices were not used. 
- mrna_metadata/small_rna_metadata: a string describing the file path for the [Metadata CSV](#metadata-csv) pertaining to the mRNA or small RNA data.
- mrna_gene_table_method/small_rna_gene_table_method: options include "full_table", "common_names_only", "gene_class_only", and "no_table". see [Gene Table](#gene-table).
- gene_table: a string describing the file path for the gene table. This gene table should include information about all three data sets. 

#### Results Tables Parameters

- generate_results_tables: TRUE or FALSE logical value describing whether Results Tables should be rendered and saved for both mRNA and small RNA data.

#### PCA Plot Parameters

- generate_pca: TRUE or FALSE logical value describing whether the PCA Plot should be rendered and saved for both mRNA and small RNA data.

#### Intra-Condition Scatter Plot Parameters

- generate_intra_condition: TRUE or FALSE logical value describing whether the Intra-Condition Scatter Plot should be rendered and saved for both mRNA and small RNA data.

#### Mean Reads Scatter Plot Parameters

- generate_mean_reads: TRUE or FALSE logical value describing whether the Mean Reads Scatter Plots should be rendered and saved for both mRNA and small RNA data.
- save_mean_reads_interactive: TRUE or FALSE logical value describing whether or not to save interactive html widget versions of mean reads scatter plots
- p_value_threshold: numeric value threshold for classifying significant genes by p-value (between 0.01 and 0.1 for best results).
- fold_change_threshold: numeric value threshold for classifying significant genes by fold change (between 1.1 and 2.0 for best results).
- lower_transparency: numeric value between 0.1 and 1.0 describing the transparency level (alpha level) for insignificant genes. Values closer to 0 produce more transparent points. 
- upper_transparency: numeric value between 0.1 and 1.0 describing the transparency level (alpha level) for significant genes. Values closer to 0 produce more transparent points. 
- customize_by_class: TRUE or FALSE logical value describing whether mean reads scatter plots should be customized by their classes, as listed in the gene table. Note that a gene table must be imported in order to customize by class. see [Class Parameters CSV](#class-parameters-csv).
- customize_by_significance: TRUE or FALSE logical value describing whether points representing insignificant genes should be colored grey. Otherwise, all points will be colored, and their p-value significance will be distinguished only by transparency.
- mrna_class_parameters/small_rna_class_parameters: a string describing the file path for the plot parameters csv, if applicable, for the mRNA and small RNA data. If left as an empty string, a default class parameters csv will be generated including all classes listed in the gene table. 

#### MA Plot Parameters

- generate_ma: TRUE or FALSE logical value describing whether the MA Plots should be rendered and saved for both mRNA and small RNA data.

#### Heatmap Parameters

- generate_heatmap: TRUE or FALSE logical value describing whether the Heatmap should be rendered and saved for both mRNA and small RNA data.
- heatmap_type: options include "complete", "all_classes", and "selected_classes".
- heatmap_selected_classes: comma-separated list of classes for which heatmaps will be generated if the heatmap type is "selected_classes". Classes should be identical to those found in the gene table.

#### Integrative Parameters

- cross_comparisons: Comma-separated List of comparisons for which to produce interactive plots. Comparisons should written in the form 'treatment_vs_control' to compare groups 'treatment' and 'control'. Names of treatment and control groups should be identical to those found in the small RNA and mRNA metadata files. Leave as empty string to produce all comparisons present in both small RNA and mRNA epxerimental runs.


- generate_integrative_results_tables: TRUE or FALSE logical value describing whether or not to generate and save integrative results tables.


- generate_cosmic_plots: TRUE or FALSE logical value describing whether or not the galaxy scatter plots should be rendered and saved.


- generate_slope_plots: TRUE or FALSE logical value describing whether or not to generate and save slope plots.


- slope_plot_classes: comma-separated list of classes for which slope plot lines/combinations will be generated. Classes should be identical to those found in the gene table. Slope Plot classes refer to the classes associated with small RNA genes.

### Metadata CSV

The metadata csv summarizes the experimental design and provides information to the DESeq2 analysis for individual mRNA and small RNA experimental runs. The first column of the metadata csv should include counts file names corresponding to each sample. If the software method being used is <strong>Counts Matrix</strong> or <strong>tiny RNA</strong> then the first column can consist of either arbitrary/empty strings or name of the counts matrix repeated on each row. The second column of the metadata csv describes specific replicates for each sample in the experiment (ex. WT_1,prg-1_2, etc.). The third column describes the condition associated with each of the samples/replicates in the previous column (ex. simply "WT" or "prg-1"). The fourth column of the csv consists of logical values (TRUE/FALSE) describing whether each sample belongs to a control group/condition. A template metadata csv can be viewed [here](template_metadata.csv), and a preview is shown below.

| files            | replicates | condition | control_condition |
|------------------|------------|-----------|-------------------|
| counts1.results  | N2_1       | N2        | TRUE              |
| counts2.results  | prg-1      | prg-1     | FALSE             |
| ...              | ...        | ...       | ...               |

### Gene Table

The gene table is a key component of customizing the pipeline's data analysis and visualization. While the gene table is optional for the mRNA and small RNA experimental runs, it is required for the integrative analysis. To disable the gene table for the mRNA and/or small RNA experimental runs, select <strong>no_table</strong> as the gene table method for the corresponding run type in the [integrative parameters yaml](#integrative-parameters-yaml). Choosing <strong>common_names_only</strong> will eliminate the class-based customization in an individual mRNA or small RNA run, and choosing <strong>class_names_only</strong> will eliminate common name substitution in an individual mRNA or small RNA run. The first column of the gene table should consist of any number of gene IDs identical to corresponding IDs from the tabulated counts files or matrix. The second column of the gene table should consist of any common names that should replace a corresponding gene ID in results tables, mean reads scatter plots, and heatmaps. To avoid replacing the gene ID for a particular gene, the gene ID should be repeated in the second column. The third column of the gene table should consist of the class to which each gene ID belongs. Class names are used to characterize and color visualized genes for the mean reads scatter plots, heatmaps (if <strong>selected_classes</strong> is chosen as the heatmap type), and slope plots. They are also included in both individual and integrative results tables. The fourth column of the gene table should consist of targets (ideally other genes in the table but potentially other genes present in the expeerimental data) to be paired with corresponding gene IDs. Multiple targets for a certain gene ID can be listed and separated by commas. A template gene table can be viewed [here](template_gene_table.csv), and a preview is shown below.

| Gene_ID   | Common_Name | Gene_Class  |
|-----------|-------------|-------------|
| Y53H1C.1  | aat-9       | CSR         |
| ...       | ...         | ...         |

### Class Parameters CSV

The class parameters csv is used to customize Mean Reads Scatter Plots for individual mRNA and small RNA experimental runs. If used, only points/genes corresponding to classes in the table will be plotted. The first column of the csv should consist of any number of class names identical to those found in the gene table. Points (genes) corresponding to these classes can be colored and re-sized in the mean reads scatter plots. The second column of the table should consist of hex color values corresponding to each class (ex. "#D95F02"). The third column of the table should consist of numeric values between 0.1 and 1.0 that describe the cex point sizes for genes of the corresponding class (default size is 0.3). The default class parameters table will color all classes in the gene table according to a standard list of 15 colors, and all classes will be given a size value 0.3. A template plot parameters csv can be seen [here](template_class_parameters.csv), and a preview is shown below. 

| point_class   | point_colors | point_sizes  |
|---------------|--------------|--------------|
| ALG           | #D95F02      | 0.7          |
| ...           | ...          | ...          |

## Output

For the mRNA and small RNA experimental runs, output will consist of results tables, pca plots, intra-condition scatter plots, mean reads scatter plots, and heatmaps, unless otherwise specified. Information about these outputs can be seen [here](../Markdown_Pipeline/)

### Integrative Results Tables

Integrative Results Tables include the class, base mean, fold change, and negative binomial test p-value for every small RNA and mRNA target pairing for which valid data is available (some low count genes have indeterminable fold changes and p-values. These genes have therefore been excluded). Classes are derived from the imported gene table. Base mean, fold change, and p-value entries are derived from the DESeq2 analysis. 

<img src=../Demo/Integrative_Demo_Results/integrative_demo_integrative_data.png width="800" height="500">

### Cosmic Plots

Cosmic Plots utilize ggplot2 techniques for representing five-dimensions of data within a two-dimensional scatter plot. The log2 fold change of the small RNA gene in each small RNA and mRNA pairing is represented on the x-axis, while the log2 fold change of the mRNA target of the corresponding small RNA is represented on the y-axis. Furthermore, pairings with an insignificant small RNA p-value (p > 0.05) are plotted in light grey, and pairings with a significant mRNA p-value (p < 0.05) are given a dark grey border around the point. Points are sized according to their mRNA log2 mean, and they are colored (if the small RNA p-value is significant) on a modified 'blues' color scale from the <strong>RColorBrewer</strong> package according to their small RNA mean (not log2 transformed for better texturing), with darker blues and black representing higher means. Thus, a comprehensive view of every small RNA and mRNA pairing across the two individual experiments can be visualized amongst the plot's four quadrants. Using the plotly package, these scatter plots offer zooming, panning, significance group isolation, and hover text features. Plots for each experimental contrast (provided the contrast is present in both experiments and specified in the associated yaml parameter) are saved as html widgets, which can be printed within a browser window to produce high-quality pdf images.

<img src=../Demo/Integrative_Demo_Results/integrative_demo_cosmic_plot.jpg width="800" height="500">

### Slope Plots

Slope Plots provide another way of visualizing fold change relationships between small RNA and mRNA pairings, namely by tracing a straight line between the fold change level of the small RNA gene on the left side of the plot and the fold change level of the corresponding mRNA target on the right side of the plot. Slope plots rely on ggplot2 plotting techniques. Lines are colored by their class as specified in the gene table (colors are built from a high-contrast, accessible palette developed by the Montgomery Lab), and it is recommended that classes representing smaller numbers of genes be singled out for plotting using the associated yaml parameter. This plot may be most useful for analyzing micro RNA (miRNA) and mRNA target relationships from one experimental condition to the other. Furthermore, line opacity for a certain class is inversely proportional to the number of genes in that class according to the function `opacity = (1/(1 + n/100))`, where n is the number of genes for a class. Thus, classes with an abundance of genes will appear more transparent. Plots for each experimental contrast (provided the contrast is present in both experiments and specified in the associated yaml parameter) are saved as pdf files. 

<img src=../Demo/Integrative_Demo_Results/integrative_demo_slope_plot.jpg width="800" height="500">

## Authors

* **Dr. Taiowa Montgomery** - 05/2021-present - Colorado State University - [taimontgomery](https://github.com/taimontgomery)
* **Spencer Kuhn** - 06/2021-present - Colorado State University - [smcguirekuhn](https://github.com/smcguirekuhn)
