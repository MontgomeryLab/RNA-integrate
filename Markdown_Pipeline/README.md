# Montgomery R Markdown Pipeline for mRNA and small RNA DGE Experiments

:warning: **Under Development & Testing** :warning:

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Configuration Files](#configuration-files)
  - [Parameters YAML](#parameters-yaml)
  - [Metadata CSV](#metadata-csv)
  - [Gene Table](#gene-table)
  - [Plot Parameters CSV](#plot-parameters-csv)
- [Output](#output)
  - [Results Tables](#results-tables)
  - [PCA Plot](#pca-plot)
  - [Intra-Condition Scatter Plot](#intra--condition-scatter-plot)
  - [Mean Reads Scatter Plots](#mean-reads-scatter-plots)
  - [MA Plots](#ma-plots)
  - [Heatmap](#heatmap)
- [Authors](#authors)


## Introduction

Welcome to the R Markdown DESeq2 graphics pipeline by the Montgomery Lab at Colorado State University. This pipeline uses tabulated counts files from mRNA or small RNA sequencing experiments to produce differential gene expression (DGE) analysis results and visualizations. In addition to a knitted R Markdown report, the pipeline saves pdf or html files to an output folder within the local directory as well as the parameters YAML used for that experimental run.

The pipeline can be executed within R studio using the following command:

`rmarkdown::render("Montgomery_DESeq2_Pipeline.Rmd",params = yaml::read_yaml("params.yml"))`

Or in the command line as such:

`R -e 'rmarkdown::render("Montgomery_DESeq2_Pipeline.Rmd",params = yaml::read_yaml("params.yml"))'`

## Installation

## Configuration Files

### Parameters YAML

### Metadata CSV

### Gene Table

### Plot Parameters CSV

## Output

### Results Tables

Results tables display the gene-wise counts from each replicate in a given contrast, followed by the fold change value of each gene and the associated p-value (adjusted) of the negative binomial hypothesis test conducted by DESeq2. Lower p-values indicate a lower probability of the null hypothesis that counts between the two conditions are derived from the same distributional parameters. Furthermore, common gene names from any uploaded gene names table will appear in a column alongside standard gene IDs for the genes that aren't associated with a common name. In addition to pdf outputs, these tables can be saved as html widgets with <strong>sorting, searching, and page size customization</strong> features.

### PCA Plot

### Intra-Condition Scatter Plot

### Mean Reads Scatter Plots

### MA Plots

### Heatmap


## Authors

* **Dr. Taiowa Montgomery** - 05/2021-present - Colorado State University - [taimontgomery](https://github.com/taimontgomery)
* **Spencer Kuhn** - 06/2021-present - Colorado State University - [smcguirekuhn](https://github.com/smcguirekuhn)
