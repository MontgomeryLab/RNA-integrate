:construction: Under Construction! :construction:

# RNA-integrate (RNA-i): integrative analysis of small RNA and mRNA high-throughput sequencing data

RNA-i performs differential small RNA and mRNA expression analysis using DESeq2 and integrates the results to identify small RNA-mRNA target relationships. The application can be run as an R script or as an R Shiny app. The Shiny app provides a user-friendly and interactive interface for first-time users, while the markdown provides a rapid and highly configurable workflow for more experienced users. Counts tables and a file with information on small RNA-mRNA target pairs are input and integrate differential gene expresion tables and plots are output. 

The workflow can be run individually for mRNA or small RNA data analysis or in the integrative mode. The quick start guides below illustrate our basic execution process for each version of the workflow. 

## Shiny App Pipeline

1. Tabulate Counts Files from RSEM, Salmon, Kallisto, HTSeq, or FeatureCounts
2. (optional) Prepare [Gene Table CSV](Shiny_App/template_gene_table.csv) for common names conversion and plot customization
3. Launch [App](Shiny_App/Montgomery_DESeq2_App.R) from RStudio or Command Line
4. Generate [Metadata CSV](Shiny_App/template_metadata.csv) within app or import an existing metadata file to specify experimental design
5. (optional) Import Gene Table CSV within the app
6. (optional) Generate [Class Parameters CSV](Shiny_App/template_class_parameters.csv) within the app or import an existing csv to customize mean reads scatter plots by gene class
7. Run Analysis (5-10 minute runtime) to produce, save, and interact with a variety of tables and plots

See the [Shiny App Directory](Shiny_App/) for more details.

## Markdown Pipeline
1. Tabulate Counts Files from RSEM, Salmon, Kallisto, HTSeq, or FeatureCounts
2. Prepare [Metadata CSV](Markdown_Pipeline/template_metadata.csv) to specify experimental design
3. (optional) Prepare [Gene Table CSV](Markdown_Pipeline/template_gene_table.csv) for common names conversion and plot customization
4. (optional) Prepare [Class Parameters CSV](Markdown_Pipeline/template_class_parameters.csv) to customize mean reads scatter plots by gene class
5. Prepare [Parameters YAML](Markdown_Pipeline/template_params.yml) to specify table and plot options
6. Execute [Markdown Script](Markdown_Pipeline/Montgomery_DESeq2_Pipeline.Rmd) from RStudio or Command Line (3-5 minute runtime) to produce and save a variety of tables and plots, including interactive html files

See the [Markdown Pipeline Directory](Markdown_Pipeline/) for more details.

## Integrative Pipeline
1. Tabulate mRNA Counts Files from RSEM, Salmon, Kallisto, HTSeq, or FeatureCounts
2. Tabulate small RNA Counts Files from RSEM, Salmon, Kallisto, HTSeq, or FeatureCounts
3. Prepare [Gene Table CSV](Integrative_Pipeline/template_gene_table.csv) for common names conversion and plot customization
4. Prepare [Metadata CSV](Integrative_Pipeline/template_metadata.csv) files for mRNA and small RNA experiments to specify experimental design
5. (optional) Prepare [Class Parameters CSV](Integrative_Pipeline/template_class_parameters.csv) for mRNA and small RNA experiments to customize mean reads scatter plots by gene class
6. Prepare Integrative [Parameters YAML](Integrative_Pipeline/template_int_params.yml) to specify table and plot options
7. Execute [Markdown Script](Integrative_Pipeline/RNA_Integrate.Rmd) from RStudio or Command Line (10-15 minute runtime) to produce and save a variety of tables and plots, including interactive html files

See the [Integrative Pipeline Directory](Integrative_Pipeline/) for more details.

## Authors

* **Taiowa Montgomery** - 05/2021-present - Colorado State University - [taimontgomery](https://github.com/taimontgomery)
* **Spencer Kuhn** - 06/2021-present - Colorado State University - [smcguirekuhn](https://github.com/smcguirekuhn)
