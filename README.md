# DESeq2 Results & Graphics Program for Differential Gene Expression

Dr. Taiowa Montgomery and Mr. Spencer Kuhn - Colorado State University RNA Biology

## Table of Contents

[Introduction](#Introduction)

[File Uploads and Samples csv](#File-Uploads-and-Samples-csv)

[Experimental Design and DESeq2](#Experimental-Design-and-DESeq2)

[PCA Plot](#PCA-Plot)

[MA Plots](#MA-Plots)

[Intra-Condition Scatterplots](#Intra--Condition-Scatterplots)

[Mean Reads Scatterplots](#Mean-Reads-Scatterplots)


## Introduction

The programs in this repository run high-throughput mRNA sequencing files through code wrapped around the DESeq2 R programming package for Differential Gene Expression (DGE) experiments. [AllCompare](AllCompareDESeqPipeline.Rmd) and [ControlCompare](ControlCompareDESeqPipeline.Rmd) are R Markdown programs that generate csv formatted results files from the DESeq2 analysis in addition to generating a PCA Plot, standard and Log Fold Change (LFC) shrinkage MA Plots, Intra-Condition Scatterplots for counts between samples, and Mean Reads Scatterplots of average counts across all samples between contrasting conditions. [AllCompare](AllCompareDESeqPipeline.Rmd) constructs all possible contrasts between the defined experimental conditions. For Instance, if four conditions are specified as A, B, C, and D, then results files, MA Plots, and Mean Reads Scatterplots will be created for the following contrasts: A/B, A/C, A/D, B/C, B/D, and C/D. [ControlCompare](ControlCompareDESeqPipeline.Rmd) constructs contrasts between a control condition and the remaining conditions. For Instance, if the same four conditions are specified, A, B, C, and D, then results files, MA Plots, and Mean Reads Scatterplots will be created for the following contrasts: A/B, A/C, and A/D. The first condition by alphabetical order is assumed to be the control condition be default. The outputs of both programs include an overall counts table csv from the general DESeq2 output, a counts table csv for each contrast, and pdf files for the PCA Plot, contrast-specific MA Plots (normal and LFC shrinkage), Intra-Condition Scatterplot, and contrast-specific Mean Reads Scatterplots. 

## File Uploads and Samples csv

## Experimental Design and DESeq2

## PCA Plot

## MA Plots

## Intra-Condition Scatterplots

## Mean Reads Scatterplots

