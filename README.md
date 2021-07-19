# DESeq2 Results & Graphics Wrapping Program for Differential Gene Expression

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

The programs in this repository run high-throughput mRNA sequencing files through code wrapped around the DESeq2 R programming package for Differential Gene Expression (DGE) experiments. [AllCompare](AllCompareDESeqPipeline.Rmd) and [ControlCompare](ControlCompareDESeqPipeline.Rmd) are R Markdown programs that generate csv formatted results files from the DESeq2 analysis in addition to generating a PCA Plot, standard and Log Fold Change (LFC) shrinkage MA Plots, Intra-Condition Scatterplots for counts between samples, and Mean Reads Scatterplots of average counts across all samples between contrasting conditions. [AllCompare](AllCompareDESeqPipeline.Rmd) analyzes all possible contrasts between the defined experimental conditions. For Instance, if four conditions are specified as A, B, C, and D, then results files, MA Plots, and Mean Reads Scatterplots will be created for the following contrasts: A/B, A/C, A/D, B/C, B/D, and C/D. [ControlCompare](ControlCompareDESeqPipeline.Rmd) analyzes contrasts between a control condition and the remaining conditions. For Instance, if the same four conditions are specified, A, B, C, and D, then results files, MA Plots, and Mean Reads Scatterplots will be created for the following contrasts: A/B, A/C, and A/D. The first condition by alphabetical order is assumed to be the control condition be default. The outputs of both programs include an overall counts table csv from the general DESeq2 output, a counts table csv for each contrast, and pdf files for the PCA Plot, contrast-specific MA Plots (normal and LFC shrinkage), Intra-Condition Scatterplot, and contrast-specific Mean Reads Scatterplots. Finally, a knitted report document can be created from the markdown files that compiles all graphics into a single pdf, with raw code included as an appendix. 

## File Uploads and Samples csv

Counts files must be present in the working directory of R or RStudio for the programs to run. 

A file called samples.csv must also be present in the working directory for the programs to run. This file contains the experimental design information necessary for providing structure to the DESeq analysis. The first column includes the names of the files to analyze. The second column includes information on the replicate number and type for each file. For example, the first replicate of a wild type condition might be written "wt_1" in the second column. The third column of samples.csv includes the condition type for each file. For the ControlCompare program, ensure that the control condition is the first condition when listed alphabetically. For the AllCompare program, the names of each condition have no restrictions. An example of a samples.csv file can be seen here

## Experimental Design and DESeq2

Three R Packages are required for this program: DESeq2, tximport, and knitr. Ensure these packages and all their dependencies are installed. 

The first chunk of code in each program constructs several accessory objects to feed into the DESeq2 analysis, including a list of the files named according to their replicate number and type, a factor vector summarizing of the different conditions present, and a data frame associated the conditions with the replicate numbers and type (essentially the last two columns of the samples.csv file). The tximport package is used to read in the counts files for DESeq2 analysis, and then DESeq2 is run, creating the object named "dds". The final line of the chunk outputs a csv file of the general counts table produced by DESeq2. 

## PCA Plot



## MA Plots

## Intra-Condition Scatterplots

## Mean Reads Scatterplots

