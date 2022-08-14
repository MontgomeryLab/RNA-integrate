# Demo Data Usage

## First Time Usage

1. Ensure that the latest versions of R (4.2.1 or newer) and Rstudio are installed
2. Clone the repository to your local device
3. Open RNA_Integrate.Rmd File (integrative analysis) or Montgomery_DESeq2_Pipeline.Rmd file (basic markdown analysis) in Rstudio
4. Select the option to install the uninstalled packages found within the script (the option should appear in a yellow bar below the toolbar)
5. Select the toolbar option to "Knit" the Markdown file. This step may take 5-10 minutes for the basic markdown analysis and 10-15 minutes for the integrative analysis
6. For best results, view the output HTML files in a browser window

## Returning Usage:

When the appropriate packages are installed, the command line (macOS) can also be used to analyze demo data with the following commands:

Integrative Pipeline: `R -e 'rmarkdown::render("RNA_Integrate.Rmd",params = yaml::read_yaml("int_params.yml"))'`
Basic Markdown Pipeline: `R -e 'rmarkdown::render("Montgomery_DESeq2_Pipeline.Rmd",params = yaml::read_yaml("params.yml"))'`

For Windows Users, it is advised to run the pipelines from within Rstudio

## Authors

* **Dr. Taiowa Montgomery** - 05/2021-present - Colorado State University - [taimontgomery](https://github.com/taimontgomery)
* **Spencer Kuhn** - 06/2021-present - Colorado State University - [smcguirekuhn](https://github.com/smcguirekuhn)
