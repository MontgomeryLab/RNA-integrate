
# Read Filtering, Mapping, and Feature Counting

## COMING SOON

# DESeq2 Results & Graphics Wrapping Program for Differential Gene Expression

Dr. Taiowa Montgomery and Mr. Spencer Kuhn - Colorado State University RNA Biology

## Table of Contents

[Introduction](#Introduction)

[File Uploads and Samples csv](#File-Uploads-and-Samples-csv)

[Required Packages](#Required-Packages)

[Experimental Design and DESeq2](#Experimental-Design-and-DESeq2)

[PCA Plot](#PCA-Plot)

[MA Plots](#MA-Plots)

[Intra-Condition Scatterplots](#Intra--Condition-Scatterplots)

[Mean Reads Scatterplots](#Mean-Reads-Scatterplots)


## Introduction

The programs in this repository run high-throughput mRNA sequencing files through code wrapped around the DESeq2 R programming package for Differential Gene Expression (DGE) experiments. [AllCompare](AllCompareDESeqPipeline.Rmd) and [ControlCompare](ControlCompareDESeqPipeline.Rmd) are R Markdown programs that generate csv formatted results files from the DESeq2 analysis in addition to generating a PCA Plot, standard and Log Fold Change (LFC) shrinkage MA Plots, Intra-Condition Scatterplots for counts between samples, and Mean Reads Scatterplots of average counts across all samples between contrasting conditions. [AllCompare](AllCompareDESeqPipeline.Rmd) analyzes all possible contrasts between the defined experimental conditions. For Instance, if four conditions are specified as A, B, C, and D, then results files, MA Plots, and Mean Reads Scatterplots will be created for the following contrasts: A/B, A/C, A/D, B/C, B/D, and C/D. [ControlCompare](ControlCompareDESeqPipeline.Rmd) analyzes contrasts between a control condition and the remaining conditions. For Instance, if the same four conditions are specified, A, B, C, and D, then results files, MA Plots, and Mean Reads Scatterplots will be created for the following contrasts: A/B, A/C, and A/D. The first condition by alphabetical order is assumed to be the control condition be default. The outputs of both programs include an overall counts table csv from the general DESeq2 output, a counts table csv for each contrast, and pdf files for the PCA Plot, contrast-specific MA Plots (normal and LFC shrinkage), Intra-Condition Scatterplot, and contrast-specific Mean Reads Scatterplots. Finally, a knitted report document can be created from the markdown files that compiles all graphics into a single pdf, with raw code included as an appendix. To knit the report, press the knit button in RStudio at the top of the code editor window. Run time can be up to several minutes. Individual code chunks within the markdown file can also be run by pressing the green run button in the top corner of each chunk. It is best practice to run chunks in order of their appearance in the script. Several chunks are used only to generate the R Markdown report. These chunks are those with an "include_graphics" command.

## File Uploads and Samples csv

Counts files from Rsem must be present in the working directory of R or RStudio for the programs to run. These files are filtered, mapped to a reference genome, and counted from fastq high-throughput sequencing files. 

A file called samples.csv must also be present in the working directory for the programs to run. This file contains the experimental design information necessary for providing structure to the DESeq analysis. The first column includes the names of the files to analyze. The second column includes information on the replicate number and type for each file. For example, the first replicate of a wild type condition might be written "wt_1" in the second column. The third column of samples.csv includes the condition type for each file. For the ControlCompare program, ensure that the control condition is the first condition when listed alphabetically. For the AllCompare program, the names of each condition have no restrictions. An example of a samples.csv file can be seen here

## Required Packages

Three R Packages are required for this program: DESeq2, tximport, and knitr. Ensure these packages and all their dependencies are installed before proceeding.

## Experimental Design and DESeq2 

In both programs, the first chunk of code following the package-calling chunk constructs several accessory objects to feed into the DESeq2 analysis, including a list of the files named according to their replicate number and type, a factor vector summarizing of the different conditions present, and a data frame associated the conditions with the replicate numbers and type (essentially the last two columns of the samples.csv file). The tximport package is used to read in the counts files for DESeq2 analysis, and then DESeq2 is run, creating the object named "dds". The final line of the chunk outputs a csv file of the general counts table produced by DESeq2. 

## PCA Plot

Both programs use the plotPCA function of the DESeq2 package to generate a PCA plot of the first two principal components for each file, colored by condition and formatted according to default ggplot arrangements. This plot is output as a pdf file with a preceding date stamp. 

## MA Plots

Both AllCompare and ControlCompare use a CompareMA, a function defined within the program that creates two MA Plot pdfs and produces a csv counts table for a given contrast between conditions. With the two contrasting conditions as the function inputs, CompareMA generates a results file from the dds object. This results file is used for the plotMA function from DESeq2. The first pdf consists of a standard MA plot, while the second consists of an MA plot after LFC shrinkage using the lfcShrink function from DESeq2. The LFC Shrinkage for both programs is type 'normal' and takes a contrast argument rather than a coefficient argument. AllCompare creates MA Plots for all possible contrasts given the condition set (CompareMA is called in a for loop iterating over all possible two-condition combinations), whereas ControlCompare contrasts the first condition listed alphabetically with the remaining conditions (CompareMA is called in a for loop iterating over combinations between the first alphabetical condition and the remaining conditions). Output pdfs have a preceding date stamp.

## Intra-Condition Scatterplots

AllCompare and ControlCompare both utilize a program-defined function called scatterplot_by_condition, which takes the condition set factor vector as its only input. The function iterates through every condition and produces a counts scatterplot for each combination of two intra-condition replicates. For instance, if there are two conditions, A and B, and three replicates within each condition, then six total scatterplots will be produced, with the first row of plots comparing replicates 1/2, 1/3, 2/3 for condition A and the second row comparing replicates 1/2, 1/3, and 2/3 for condition B. These scatterplots are compiled into a single pdf with a preceding date stamp. 

## Mean Reads Scatterplots

For each condition, mean reads across all replicates are calculated in both programs by iterating through the general counts table from the dds object. Within the program, the mean_scatterplot function is defined to plot mean reads given two input conditions to contrast. Significant mean counts are separated according to adjusted p-values of < 0.05 and Log Fold Changes of greater than 0.378511623 or less than -0.378511623. These significant points are slightly larger and colored blue, whereas the remaining points are colored gray, although the colors can be changed within the Plot_Colors object. Vertical and horizontal dotted guidelines are added to the plot separating average reads greater than 10. Diagonal guidelines are also added to help delineate genes in which one condition's mean differs from the other condition's mean by a factor of two. Along the plot axes, tick marks are added according to a log2 scale, with major ticks at 1/16, 1, 16, 256, 4096, 65536, and 1,048,596. Minor ticks are added according to regular whole number intervals between major ticks, with four minor ticks within each interval. AllCompare runs this function for all two-condition combinations, whereas ControlCompare runs the function for all combinations between the control condition and the remaining conditions. Each plot is generated as a separate pdf with a preceding date stamp.
