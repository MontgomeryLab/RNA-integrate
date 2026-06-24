# RNA-Integrate (RNA-i): Integrative Analysis of Small RNA and mRNA High-Throughput Sequencing Data

RNA-Integrate (RNA-i) is a reproducible workflow for differential expression analysis of **small RNA** and **mRNA sequencing data** using **DESeq2**, with optional integrative analysis to identify **small RNA–mRNA target relationships**.

The workflow supports both **standalone differential expression analysis** (**mRNA** or **small RNA independently**) and **integrative analysis**, where differential expression results are combined with a user-supplied small RNA–mRNA target table to generate integrative visualizations, statistical summaries, and target-pair results tables.

RNA-i can be run through either:

-   a **Shiny application**, providing a user-friendly and interactive interface for first-time users, or
-   a highly configurable **R Markdown workflow**, designed for rapid, reproducible, and customizable execution.

Input files include sequencing count tables, experimental metadata, and (for integrative analyses) a small RNA–mRNA target table. Outputs include differential expression tables, quality-control plots, expression visualizations, and integrative analyses of small RNA–mRNA relationships.

The workflow can be run individually for **mRNA** or **small RNA** differential expression analysis or in **integrative mode**. The quick start guides below illustrate the basic execution process for each version of the workflow.

## Workflow Modes

### Differential Expression Analysis (mRNA)

mRNA differential expression analysis includes:

-   Differential expression results tables
-   PCA plots
-   Intra-condition replicate scatterplots
-   Mean reads scatter plots
-   MA plots
-   Heatmaps

**Supported mRNA Input Formats**

-   HTSeq
-   Counts Matrix
-   RSEM
-   Salmon
-   Kallisto

### Differential Expression Analysis (Small RNA)

Small RNA differential expression analysis includes all standard DESeq2 visualizations and optional class-based customization.

**Supported Small RNA Input Formats**

-   Counts Matrix
-   tinyRNA output tables

Small RNA analyses support automatic handling of annotation columns, feature classes, and common-name extraction for visualization and integration.

### Integrative Analysis

Integrative analysis combines differential expression results from small RNA and mRNA datasets to evaluate small RNA–mRNA target relationships.

**Integrative Outputs**

-   Integrative results tables
-   Cosmic plots
-   Slope plots
-   Statistical analyses

## Shiny App Pipeline *(Under Construction!)*

1.  Prepare count files for **mRNA** or **small RNA** experiments
2.  *(Optional)* Prepare a [Gene Table CSV](Shiny_App/template_gene_table.csv) for common-name conversion and class-based customization
3.  Launch [App](Shiny_App/Montgomery_DESeq2_App.R) from **RStudio** or the **command line**
4.  Generate or import a [Metadata CSV](Shiny_App/template_metadata.csv) to specify experimental design
5.  *(Optional)* Import a Gene Table CSV within the app
6.  *(Optional)* Generate or import a [Class Parameters CSV](Shiny_App/template_class_parameters.csv) to customize mean reads scatter plots by gene class
7.  Run analysis (**5–10 minute runtime**) to produce, save, and interact with tables and plots

See the [Shiny App Directory](Shiny_App/) for more details.

## Markdown Pipeline

1.  Prepare sequencing count data
2.  Prepare a [Metadata CSV](Markdown_Pipeline/template_metadata.csv) to specify experimental design
3.  *(Optional)* Prepare a [Gene Table CSV](Markdown_Pipeline/template_gene_table.csv) for common names conversion and gene-class annotations
4.  *(Optional)* Prepare a [Class Parameters CSV](Markdown_Pipeline/template_class_parameters.csv) to customize mean reads scatter plots by gene class
5.  Prepare a [Parameters YAML](Markdown_Pipeline/template_params.yml) to specify analysis settings and output options
6.  Execute [Markdown Script](Markdown_Pipeline/DESeq2_Report.Rmd) from **RStudio** or the **command line** (**3–5 minute runtime**) to produce and save a variety of tables and plots, including optional interactive HTML files

See the [Markdown Pipeline Directory](Markdown_Pipeline/) for more details.

## Integrative Pipeline

1.  Prepare **mRNA count data** (**HTSeq, Counts Matrix, RSEM, Salmon, or Kallisto**)
2.  Prepare **small RNA count data** (**Counts Matrix or tinyRNA format**)
3.  Prepare a [Gene Table CSV](Integrative_Pipeline/template_gene_table.csv) containing small RNA–mRNA target relationships
4.  Prepare [Metadata CSV](Integrative_Pipeline/template_metadata.csv) files for **mRNA** and **small RNA** experiments to specify experimental design
5.  *(Optional)* Prepare [Class Parameters CSV](Integrative_Pipeline/template_class_parameters.csv) files for **mRNA** and/or **small RNA** experiments to customize mean reads scatter plots by gene class
6.  Prepare Integrative [Parameters YAML](Integrative_Pipeline/template_int_params.yml) to specify comparisons and output settings
7.  Execute [Markdown Script](Integrative_Pipeline/RNA_Integrate.Rmd) from **RStudio** or the **command line** (**10–15 minute runtime**) to produce and save a variety of tables and plots, including interactive HTML files

See the [Integrative Pipeline Directory](Integrative_Pipeline/) for more details.

## Authors

-   **Taiowa Montgomery** — 05/2021–present — Colorado State University — [taimontgomery](https://github.com/taimontgomery)
-   **Spencer Kuhn** — 06/2021–present — Colorado State University — [smcguirekuhn](https://github.com/smcguirekuhn)
-   **Paige Lillibridge** — 08/2025–present — Colorado State University — [paigelily](https://github.com/paigelily)
