# Montgomery Lab DESeq2 Pipeline

Welcome to the Montgomery DESeq2 Pipeline in R Markdown and R Shiny, developed by the Montgomery Lab at Colorado State University. Our pipelines use tabulated counts files from a variety of software types to produce data analysis and visualization results for differential gene expression experiments. Our R Shiny app provides a user-friendly and interactive interface for first-time users, while our markdown pipeline provides rapid, malleable analysis for more experienced users. The Integrative Markdown Pipeline not only executes parallel runs of mRNA and small RNA data analysis for the same experimental design, but also produces integrated tables and plots for the analysis of small RNA and mRNA target combinations. 

The quick start guides below illustrate our basic execution process for each version of the pipeline. 

## Shiny App Pipeline

1. Tabulate Counts Files from RSEM, Salmon, Kallisto, HTSeq, or FeatureCounts
2. (optional) Prepare [Gene Table CSV](Shiny_App/gene_table.csv) for common names conversion and plot customization
3. Launch [App](Shiny_App/Montgomery_DESeq2_App.R) from RStudio or Command Line
4. Generate [Metadata CSV](Shiny_App/metadata.csv) within app or import an existing metadata file to specify experimental design
5. (optional) Import Gene Table CSV within the app
6. (optional) Generate [Class Parameters CSV](Shiny_App/class_parameters.csv) within the app or import an existing csv to customize mean reads scatter plots by gene class
7. Run Analysis (5-10 minute runtime) to produce, save, and interact with a variety of tables and plots

See the [Shiny App Directory](Shiny_App/) for more details.

## Markdown Pipeline
1. Tabulate Counts Files from RSEM, Salmon, Kallisto, HTSeq, or FeatureCounts
2. Prepare [Metadata CSV](Markdown_Pipeline/metadata.csv) to specify experimental design
3. (optional) Prepare [Gene Table CSV](Markdown_Pipeline/gene_table.csv) for common names conversion and plot customization
4. (optional) Prepare [Class Parameters CSV](Markdown_Pipeline/class_parameters.csv) to customize mean reads scatter plots by gene class
5. Prepare [Parameters YAML](Markdown_Pipeline/params.yml) to specify table and plot options
6. Execute [Markdown Script](Markdown_Pipeline/Montgomery_DESeq2_Pipeline.Rmd) from RStudio or Command Line (3-5 minute runtime) to produce and save a variety of tables and plots, including interactive html files

See the [Markdown Pipeline Directory](Markdown_Pipeline/) for more details.

## Integrative Pipeline
1. Tabulate mRNA Counts Files from RSEM, Salmon, Kallisto, HTSeq, or FeatureCounts
2. Tabulate small RNA Counts Files from RSEM, Salmon, Kallisto, HTSeq, or FeatureCounts
3. Prepare [Gene Table CSV](Integrative_Pipeline/gene_table.csv) for common names conversion and plot customization
4. Prepare [Metadata CSV](Integrative_Pipeline/metadata.csv) files for mRNA and small RNA experiments to specify experimental design
5. (optional) Prepare [Class Parameters CSV](Integrative_Pipeline/class_parameters.csv) for mRNA and small RNA experiments to customize mean reads scatter plots by gene class
6. Prepare Integrative [Parameters YAML](Integrative_Pipeline/int_params.yml) to specify table and plot options
7. Execute [Markdown Script](Integrative_Pipeline/Montgomery_Integrative_DESeq2_Pipeline.Rmd) from RStudio or Command Line (10-15 minute runtime) to produce and save a variety of tables and plots, including interactive html files

See the [Integrative Pipeline Directory](Integrative_Pipeline/) for more details.

## Authors

* **Dr. Taiowa Montgomery** - 05/2021-present - Colorado State University - [taimontgomery](https://github.com/taimontgomery)
* **Spencer Kuhn** - 06/2021-present - Colorado State University - [smcguirekuhn](https://github.com/smcguirekuhn)
