# Demo Data Usage

This directory contains example datasets and configuration files for running both the standalone **DESeq2 Markdown Pipeline** and the **RNA-Integrate Pipeline**.

The demo files are intended to help users confirm successful installation and learn the expected formatting for metadata, counts files, gene tables, class parameters, and YAML configuration files.

## First-Time Usage

1.  Install the latest versions of **R** (**4.2.1 or newer recommended**) and **RStudio**.

2.  Clone the repository to your local device.

3.  Open one of the following files in **RStudio**:

### DESeq2 Markdown Pipeline

Open:

``` text
Markdown_Pipeline/DESeq2_Report.Rmd
```

This pipeline performs standalone differential expression analysis for **mRNA** or **small RNA** datasets.

### RNA-Integrate Pipeline

Open:

``` text
Integrative_Pipeline/RNA_Integrate.Rmd
```

This pipeline performs integrated analysis between **small RNA** and **mRNA** experiments.

4.  If prompted, install any missing packages required by the pipeline.

In RStudio, missing package installation prompts typically appear in a yellow notification bar below the toolbar.

5.  Run the pipeline by selecting **"Knit"** in the RStudio toolbar.

Approximate runtimes:

| Pipeline                 | Runtime       |
|--------------------------|---------------|
| DESeq2 Markdown Pipeline | 3–5 minutes   |
| RNA-Integrate Pipeline   | 10–15 minutes |

6.  For best results, open the generated HTML report in a web browser after rendering completes.

## Returning Usage

After required packages are installed, the demo workflows can also be run from the command line.

### RNA-Integrate Pipeline

Run from the `Demo/` directory:

``` bash
R -e 'rmarkdown::render("RNA_Integrate.Rmd", params = yaml::read_yaml("int_params.yml"))'
```

### DESeq2 Markdown Pipeline

Run from the `Demo/` directory:

``` bash
R -e 'rmarkdown::render("DESeq2_Report.Rmd", params = yaml::read_yaml("deseq2_params.yml"))'
```

### macOS and Linux

Command-line execution is supported and recommended for reproducible runs.

### Windows

Running the pipelines through **RStudio** is recommended for easier package management and troubleshooting.

## Demo Files Included

The demo directory contains example files for both mRNA and small RNA workflows, including:

-   metadata CSV files
-   gene tables
-   class parameter CSV files
-   mRNA count files
-   small RNA count matrices
-   example YAML parameter files
-   example output reports and result folders

These files can be modified and reused as templates for new analyses.

## Authors

-   **Dr. Taiowa Montgomery** - 05/2021-present - Colorado State University - [taimontgomery](https://github.com/taimontgomery)
-   **Spencer Kuhn** - 06/2021-present - Colorado State University - [smcguirekuhn](https://github.com/smcguirekuhn)
-   **Paige Lillibridge** - 2025-present - Colorado State University - [paigelily](https://github.com/paigelily)
