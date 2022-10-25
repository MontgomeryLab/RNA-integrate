# Define UI for application that draws a histogram
ui = fluidPage(
  
  # Title
  titlePanel("DESeq2 Graphics App, Montgomery Lab at Colorado State University"),
  
  # Body
  mainPanel(
    tabsetPanel(
      type = "tabs",
      
      # Parameters Tab
      tabPanel(
        "Parameters",
        
        # Header
        h2("Parameters"),br(),
        
        # Introduction
        HTML(
          "<p>Welcome to the Montgomery DESeq2 Graphics and Analysis App, developed by the Montgomery Lab at Colorado State University. 
          This app is designed to convert tabulated counts files from mRNA or small RNA sequencing Differential Gene Expression (DGE) experiments into user-friendly results tables and plots. 
          The app does not assist with upstream processes of the RNA-sequencing pipeline, such as filtering or aligment. 
          Compatible tabulation software includes Salmon, Kallisto, RSEM, and HTSeq. Furthermore, standard count matrices, 
          (such as the kind produced by featureCounts), as well as tiny RNA matrices (developed by the Montgomery Lab) are acceptable input types.</p>"
        ),
        
        br(),
        
        # Experiment ID
        h2("Step One: Experiment ID"),
        
        HTML(
          "<p>Enter a brief identifier (no spaces or special characters) for this experimental run (ex. <strong>elegans1</strong>).</p> 
          <p>The experiment ID will precede file names of saved tables and plots.</p>"
        ),
        
        textInput(inputId = "experiment_id",label = NULL),
        
        # Counts Files Parameters
        h2("Step Two: Counts Files Parameters"),
        
        HTML(
          "<p>First, select the tabulation software used to produce counts files. Options include <strong>RSEM, Salmon, Kallisto, and HTSeq</strong>.
          for FeatureCounts matrices and general counts matrices for which samples are separated by column, select <strong>Counts Matrix</strong>.</p>
          <p>Next, if <strong>Counts Matrix</strong> or <strong>tiny RNA</strong> is selected, choose the appropriate matrix file from the working directory.</p>"
        ),
        
        br(),
        
        fluidRow(
          column(
            width = 6,
            
            # Tabulation Software
            selectInput(
              inputId = "software_method",
              label = "Tabulation Software",
              choices = c("RSEM" = "rsem","Salmon" = "salmon","Kallisto" = "kallisto","HTSeq" = "htseq","Counts Matrix" = "counts_matrix","tiny RNA" = "tiny_rna"),
              selected = "rsem",
              multiple = FALSE
            )
          ),
          
          column(
            width = 6,
            
            # Conditional Counts Matrix Selection from Working Directory
            conditionalPanel(
              condition = "input.software_method == 'counts_matrix' | input.software_method == 'tiny_rna'",
              selectInput(
                inputId = "selected_counts_matrix",
                label = "Select counts matrix for this experiment",
                choices = c("select csv file" = "",dir()),
                selected = character(0),
                multiple = FALSE
              )
            )
          )
        ),
        
        # Metadata and Gene Table Files
        h2("Step Three: Metadata and Gene Table"),
        
        HTML(
          "<p>First, select the metadata file corresponding to the experimental run, or choose to generate a new metadata file within the app.
          To generate a new metadata file, select <strong>generate new</strong> as the metadata method and then follow the prompts that appear.
          Otherwise, choose <strong>working directory</strong> and choose the appropriate csv file.</p>
          <p>Next, choose what type of gene table, if any, should be uploaded for the purposes of customizing tables and plots.
          an uploaded gene table csv can include gene IDs and their <strong>Common Names</strong> (public names), <strong>Class Names</strong> (for small RNA or other classifications), or <strong>Both</strong>.
          To proceed withouth using a gene table, select <strong>No Table</strong>.</p>"
        ),
        
        br(),
        
        fluidRow(
          column(
            width = 6,
            
            # Header Text
            h2("Metadata Parameters"),
            
            # Metadata Method
            selectInput(
              inputId = "metadata_method",
              label = "Use existing csv or generate a new csv",
              choices = c("select method" = "","working directory" = "working_directory","generate new" = "generate"),
              selected = character(0)
            ),
            
            # Conditional Csv Selection from Working Directory
            conditionalPanel(
              condition = "input.metadata_method == 'working_directory'",
              selectInput(
                inputId = "selected_metadata",
                label = "Select csv for this experiment",
                choices = c("select csv file" = "",dir()),
                selected = character(0),
                multiple = FALSE
              )
            ),
            
            # Condition Number Selection
            conditionalPanel(
              condition = "input.metadata_method == 'generate'",
              radioButtons(
                inputId = "condition_number",
                label = "Select number of experimental conditions",
                choices = 2:10
              )
            ),
            
            # Condition Parameter Entries UI
            uiOutput("condition_ui")
          ),
          
          column(
            width = 6,
            
            # Header Text
            h2("Gene Table Parameters"),
            
            # Gene Table Method
            selectInput(
              inputId = "gene_table_method",
              label = "Choose the type of gene table to import",
              choices = c("select method" = "","Full Table" = "full_table","Common Names Only" = "common_names_only","Gene Class Only" = "gene_class_only","No Table" = "no_table"),
              selected = character(0),
              multiple = FALSE
            ),
            
            # Conditional Gene Table Selection from Working Directory
            conditionalPanel(
              condition = "input.gene_table_method != 'no_table' & input.gene_table_method != ''",
              selectInput(
                inputId = "selected_gene_table",
                label = "Select gene table csv for this Experiment",
                choices = c("select csv file" = "",dir()),
                selected = character(0),
                multiple = FALSE
              )
            )
          )
        ),
        
        # Plot Selection
        h2("Step Four: Plot Selection"),
        
        HTML("<p>Choose which plots to render and save for this experiment.</p>"),
        
        checkboxInput(inputId = "results_tables_render",label = "Results Tables",value = TRUE),
        checkboxInput(inputId = "pca_render",label = "PCA Plot",value = TRUE),
        checkboxInput(inputId = "intra_condition_render",label = "Intra-Condition Scatter Plot",value = TRUE),
        checkboxInput(inputId = "mean_reads_render",label = "Mean Reads Scatter Plots",value = TRUE),
        checkboxInput(inputId = "ma_render",label = "MA Plots",value = TRUE),
        checkboxInput(inputId = "heatmap_render",label = "Heatmap",value = TRUE),
        
        # Plot Customization
        h2("Step Five: Plot Customization"),
        
        HTML(
          "<p>Follow the prompts to customize mean reads scatter plots and/or the heatmap. 
          Without importing a gene table, only the p-value, fold change, and transparency thresholds of the mean reads scatter plots will be customizable.
          With a gene table imported, both the mean reads scatter plots and the heatmap can be customized by the classes included in the gene table.</p> 
          <p>A <strong>class parameters</strong> csv can be imported from the working directory or generated within the app to color or size points by class.
          Furthermore, genes with insignificant p-values can be colored grey regardless of their classification color.</p>
          <p>Three heatmap types are available if a gene table is imported. The <strong>complete</strong> heatmap clusters all genes in the experiment regardless of class.
          The <strong>All Classes</strong> heatmap clusters genes within each class and compiles the separate clusters into a single plot.
          The <strong>Selected Classes</strong> heatmap isolates user-chosen classes into a series of individually-clustered subplots.</p>"
        ),
        
        br(),
        
        fluidRow(
          column(
            width = 7,
            
            # Mean Reads Scatter Plot Parameters
            conditionalPanel(
              condition = "input.mean_reads_render",
              
              # Header Text
              h2("Mean Reads Scatter Plot Parameters"),
              
              # P-Value Threshold
              sliderInput(
                inputId = "p_value_threshold",
                label = "select significant gene p-value threshold",
                value = 0.05,min = 0.01,max = 0.1,step = 0.01
              ),
              
              # Fold Change Threshold
              sliderInput(
                inputId = "fold_change_threshold",
                label = "select significant gene fold change threshold",
                value = 1.3,min = 1.1,max = 2.0,step = 0.05
              ),
              
              # Lower Transparency Threshold
              sliderInput(
                inputId = "lower_transparency",
                label = "transparency level for insignificant genes",
                value = 0.2,min = 0.1,max = 1.0,step = 0.05
              ),
              
              # Fold Change Threshold
              sliderInput(
                inputId = "upper_transparency",
                label = "transparency level for significant genes",
                value = 0.9,min = 0.1,max = 1.0,step = 0.05
              )
            ),
            
            # Conditional Class Parameters Method Selection
            conditionalPanel(
              condition = "input.selected_gene_table != '' & input.mean_reads_render",
              checkboxInput(
                inputId = "customize_by_class",
                label = "customize plots by class?",
                value = FALSE
              ),
              
              checkboxInput(
                inputId = "customize_by_significance",
                label = "color insignificant genes grey?",
                value = FALSE
              )
            ),
            
            # Conditional Class Parameters Generation Option
            conditionalPanel(
              condition = "input.customize_by_class",
              checkboxInput(
                inputId = "generate_class_parameters",
                label = "generate new class parameters csv?",
                value = FALSE
              )
            ),
            
            # Conditional Csv Selection from Working Directory
            conditionalPanel(
              condition = "input.customize_by_class & !input.generate_class_parameters",
              selectInput(
                inputId = "select_class_parameters",
                label = "Select class parameters csv or leave blank for default csv",
                choices = c("select csv file" = "",dir()),
                selected = character(0),
                multiple = FALSE
              )
            ),
            
            # Conditional Plot Parameters Class Selection
            uiOutput("class_parameters_ui"),
            
            # Condition Parameter Entries UI
            uiOutput("class_ui")
          ),
          
          column(
            width = 5,
            
            # Heatmap Parameters
            conditionalPanel(
              condition = "input.selected_gene_table != '' & input.heatmap_render",
              
              # Header Text
              h2("Heatmap Parameters"),
              
              # Heatmap Type Selection
              selectInput(
                inputId = "heatmap_type",
                label = "choose heatmap type",
                choices = c("complete" = "complete","class-separated: all classes" = "all_classes","class-separated: selected classes" = "selected_classes"),
                selected = "complete"
              )
            ),
            
            # Conditional Heatmap Class Selection
            uiOutput("heatmap_classes_ui")
          )
        ),
        
        # Run Button
        h2("Step Six: Run Analysis"),
        
        HTML(
          "<p>Running a full analysis may take several minutes depending on the size of the data and the level of customization.</p>
          <p>A <strong>progress bar</strong> will appear in the lower right corner once plot saving has started.</p>"
        ),
        
        br(),
        
        actionButton(
          inputId = "run",
          label = "Analyze",
          icon = icon("flask")
        ),
        
        br()
      ),
      
      tabPanel(
        "Results Tables",
        
        h2("Results Tables"),
        
        HTML(
          "Results tables display the gene-wise counts from each replicate in a given contrast, followed by the fold change value of each gene and the associated p-value (adjusted) of the negative binomial hypothesis test conducted by DESeq2. 
          Lower p-values indicate a lower probability of the null hypothesis that counts between the two conditions are derived from the same distributional parameters. 
          Furthermore, any common gene names and/or classes from an uploaded gene table will be included in columns beside the Gene ID column. 
          In addition to csv outputs, these tables are rendered within the app as html widgets with <strong>sorting, searching, and page size customization</strong> features."
        ),
        
        br(),
        
        uiOutput("results_tables_ui"),
        
        br()
      ),
      
      tabPanel(
        "PCA Plot",
        
        h2("PCA Plot"),
        HTML(
          "<p>The PCA plot displays loadings of the first two principal components for each sample/biological replicate in the experimental design. 
          Colored by experimental condition, the points of the PCA plot provide a visualization of clustering amongst the samples, both within conditions and across conditions.</p>"
        ),
        
        br(),
        
        uiOutput("pca_ui"),
        
        br()
      ),
      
      tabPanel(
        "Intra-Condition Scatter Plot",
        
        h2("Intra-Condition Scatter Plot"),
        
        HTML(
          "<p>The Intra-Condition Scatter Plots display log2 counts between pairs of biological replicates within each condition of the experimental design. 
          Any counts in the data below 1 are replaced with a value of 1 to simplify the log2 transformation. All replicate pairs in each condition are displayed.</p>"
        ),
        
        br(),
        
        uiOutput("intra_condition_ui"),
        
        br()
      ),
      
      tabPanel(
        "Mean Reads Scatter Plots",
        
        h2("Mean Reads Scatter Plots"),
        
        HTML(
          "<p>The Mean Reads Scatter Plots display average log2 counts across biological replicates of experimental condition pairs for a provided contrast. 
          Average counts with a value of 0 are assigned the value -4 following the log2 transformation. Guidelines are added to assist with visualization. 
          Statistically significant genes are less transparent than insignificant genes, and in default plots, 
          they are colored <strong>blue</strong> and sized with a value of <strong>0.5</strong> (insignificant genes are sized with a value of <strong>0.3</strong>). 
          Axes are scaled to more closely resemble the log2 transformation.Unlabeled tick marks, therefore, do not always represent whole number intervals between labeled tick marks. 
          In addition to pdf outputs, Mean Reads Scatter Plots can be saved as html widgets with <strong>draggable zooming, panning, and hover text</strong> features using the download button. 
          To reset the scatter plot axes within the widget, the <strong>home</strong> button in the top right corner of the widget can be pressed.</p>
          <p>The mean reads scatter plots can be customized by changing the <strong>p-value</strong> and <strong>fold change</strong> thresholds for distinguishing statistically significant genes. 
          Furthermore, <strong>upper and lower transparency</strong> thresholds can be set for distinguishing statistical significance. 
          The lower threshold corresponds to insignificant genes. The points of the plot can also be <strong>colored and sized</strong> according to their gene classes as specified by the <strong>gene table</strong>. 
          If no gene table has been selected/uploaded, only the p-value, fold change, and transparency thresholds will be customizable. 
          If a gene table has been selected/uploaded, a class parameters csv can be selected from the working directory or built within the app. If used, only points/genes corresponding to classes in the table will be plotted. 
          If no class parameters csv is selected, every class will be automatically colored from a list of 15 colors, and every point will be sized with a cex value of 0.3. 
          Insignificant genes will be colored grey if the corresponding checkbox option is selected.</p>"
        ),
        
        br(),
        
        uiOutput("mean_reads_ui"),
        
        br()
      ),
      
      tabPanel(
        "MA Plots",
        
        h2("MA Plots"),
        
        HTML(
          "<p>The Standard MA Plots are built from the DESeq2 package and display the fold change of a gene over its mean counts value (normalized) for a provided contrast between experimental conditions. 
          Genes with a statistically significant p-value (p < 0.05) are colored in <strong>blue</strong>.</p>"
        ),
        
        br(),
        
        uiOutput("ma_ui"),
        
        br()
      ),
      
      tabPanel(
        "Heatmap",
        
        h2("Heatmap"),
        
        HTML(
          "<p>The <strong>Complete</strong> heatmap shows log2 counts across all samples for any genes above a certain mean count threshold of three (meaning an average of eight counts across all samples). 
          Using the package <strong>heatmaply</strong>, an interactive html widget is rendered with <strong>draggable zooming, panning, and hover text</strong> features. 
          To reset the heatmap axes, the home button in the top right corner of the widget can be pressed. Darker blue cells indicate lower log2 counts values, while darker red cells indicate higher log2 counts values. 
          Furthermore, rows (genes) are clustered using the complete hierarchical clustering method in R (hclust) with eudclidean distances.</p>
          <p>The <strong>Class-Separated: All Classes</strong> heatmap relies on the gene table to map all genes with an associated class. 
          Such genes are grouped alphabetically by class and clustered hierarchically within their respective group. Hover text indicates the class assigned to each gene.</p>
          <p>The <strong>Class-Separated: Selected Classes</strong> heatmap produces subplots for each class listed in the associated yaml parameter entry. 
          Heatmaps for each listed class are scaled according to the limits of the complete heatmap, although cell sizes are scaled according to the number of genes associated with each class.</p>"
        ),
        
        br(),
        
        uiOutput("heatmap_ui"),
        
        br()
      )
    )
  )
)