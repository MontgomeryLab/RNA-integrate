# DESeq2 Graphics App - Montgomery Lab at Colorado State University
# Dr. Taiowa Montgomery and Spencer Kuhn

## Environment Setup

# Package Installation (if Necessary)
if (require(shiny) == FALSE) {
  install.packages("shiny")
}
if (require(shinycssloaders) == FALSE) {
  install.packages("shinycssloaders")
}
if (require(DT) == FALSE) {
  install.packages("DT")
}
if (require(heatmaply) == FALSE) {
  install.packages("heatmaply")
}
if (require(ggplot2) == FALSE) {
  install.packages("ggplot2")
}
if (require(BiocManager) == FALSE) {
  install.packages("BiocManager")
}
if (require(tximport) == FALSE) {
  BiocManager::install("tximport",force = TRUE,update = FALSE)
}
if (require(DESeq2) == FALSE) {
  BiocManager::install("DESeq2",force = TRUE,update = FALSE)
}

# Load Installed Packages
library(shiny)
library(shinycssloaders)
library(DESeq2)
library(tximport)
library(DT)
library(heatmaply)
library(ggplot2)

# Source All Functions
for (i in list.files(path = "\\.R$",full.names = TRUE)) {
  source(i)
}
