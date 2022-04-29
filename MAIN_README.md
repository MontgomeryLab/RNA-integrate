# Montgomery Lab DESeq2 Pipeline

Welcome to the Montgomery DESeq2 Pipeline in R Markdown and R Shiny, developed by the Montgomery Lab at Colorado State University. Our pipelines use tabulated counts files from a variety of software types to produce data analysis and visualization results for differential gene expression experiments. Our R Shiny app provides a user-friendly and interactive interface for first-time users, while our markdown pipeline provides rapid, malleable analysis for more experienced users. 

The quick start guides below illustrate our basic execution process for each version of the pipeline. 

## Shiny App Pipeline

1. Tabulate Counts Files/[Counts Matrix](Shiny_App/count_matrix.csv)
2. (optional) Prepare [Gene Table CSV](Shiny_App/gene_table.csv)
3. Launch [App](Shiny_App/Montgomery_DESeq2_App.R)
4. Generate or Import [Metadata CSV](Shiny_App/metadata.csv)
5. (optional) Import Gene Table CSV
6. (optional) Generate or Import [Class Parameters CSV](Shiny_App/class_parameters.csv)
7. Run Analysis

See the [Shiny App Directory](Shiny_App/) for more details.

## Markdown Pipeline
1. Tabulate Counts Files/[Counts Matrix](Markdown_Pipeline/count_matrix.csv)
2. Prepare [Metadata CSV](Markdown_Pipeline/metadata.csv)
3. (optional) Prepare [Gene Table CSV](Markdown_Pipeline/gene_table.csv)
4. (optional) Prepare [Class Parameters CSV](Markdown_Pipeline/class_parameters.csv)
5. Prepare [Parameters YAML](Markdown_Pipeline/params.yml)
6. Execute [Markdown Script](Markdown_Pipeline/Montgomery_DESeq2_Pipeline.Rmd)

See the [Markdown Pipeline Directory](Markdown_Pipeline/) for more details.

## Authors

* **Dr. Taiowa Montgomery** - 05/2021-present - Colorado State University - [taimontgomery](https://github.com/taimontgomery)
* **Spencer Kuhn** - 06/2021-present - Colorado State University - [smcguirekuhn](https://github.com/smcguirekuhn)
