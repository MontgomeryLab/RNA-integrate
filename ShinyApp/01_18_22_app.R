# DESeq2 Graphics and Figures Pipeline for mRNA-based DGE Experiments

# Environment Setup
library(shiny)
library(shinycssloaders)
library(DESeq2)
library(tximport)
library(heatmaply)
library(ggplot2)
library(DT)

# User Interface
ui = fluidPage(
    titlePanel("DESeq2 Graphics App, Montgomery Lab at Colorado State University"),
    sidebarLayout(position = "right",
                  sidebarPanel(
                      
                      # Header Text
                      h2("Counts Files Parameters"),
                      
                      # File Selection Method
                      selectInput("file_method", "use working directory or upload counts files",
                                  c("working directory" = "working_directory","upload files" = "upload_files")),
                      
                      # File Tabulation Software Method
                      selectInput("software_method", "what software was used to tabulate counts?",
                                  c("Salmon" = "salmon","Kallisto" = "kallisto","RSEM" = "rsem","HTSeq" = "htseq"),
                                  selected = "rsem",multiple = FALSE),
                      
                      # Header Text
                      h2("Metadata Parameters"),
                      
                      # Samples.csv Method
                      conditionalPanel(condition = "input.file_method == 'working_directory'",
                                       selectInput("samples_method", "use existing csv or generate new csv",
                                                   c("select method" = "",
                                                     "working directory" = "working_directory",
                                                     "upload file" = "upload",
                                                     "generate new" = "generate"),
                                                   selected = character(0))),
                      
                      # Conditional Samples.csv Selection from Working Directory
                      conditionalPanel(condition = "input.samples_method == 'working_directory'",
                                       selectInput("selected_samples","select csv for this experiment",
                                                   c("select csv file" = "",dir()),selected = FALSE,multiple = FALSE)),
                      
                      # Conditional Samples.csv Upload
                      conditionalPanel(condition = "input.samples_method == 'upload' & input.file_method == 'working_directory'",
                                       fileInput("uploaded_samples","upload csv",multiple = FALSE)),
                      
                      # Condition number selection
                      conditionalPanel(condition = "input.file_method == 'upload_files' | input.samples_method == 'generate'",
                                       radioButtons("condition_number","Select number of experimental conditions",2:10)),
                      
                      # Condition Parameter Entries UI
                      uiOutput("condition_ui"),
                      
                      # Header Text
                      h2("Gene Table Parameters (optional)"),
                      
                      # Gene Table Method
                      selectInput("gene_method", "use gene names table to convert common/public names?",
                                  c("select method" = "","working directory" = "working_directory","upload file" = "upload_table","none" = "no_table"),
                                  selected = character(0),multiple = FALSE),
                      
                      # Conditional Gene Table Selection from Working Directory
                      conditionalPanel(condition = "input.gene_method == 'working_directory'",
                                       selectInput("selected_gene_table","select gene table csv for this Experiment",
                                                   c("select csv file" = "",dir()),selected = character(0),multiple = FALSE)),
                      
                      # Conditional Gene Table Upload
                      conditionalPanel(condition = "input.gene_method == 'upload_table'",
                                       fileInput("uploaded_gene_table","upload gene table",multiple = FALSE)),
                      
                      # Desired Plot Selection
                      checkboxGroupInput("plot_selection","Which outputs should be generated?",
                                         choices = c("Results Tables","PCA Plot","Heatmap","Intra-Condition Scatter Plots",
                                                     "Standard MA Plots","LFC Shrinkage MA Plots","Mean Reads Scatter Plots"),
                                         selected = c("Results Tables","PCA Plot","Heatmap","Intra-Condition Scatter Plots",
                                                      "Standard MA Plots","LFC Shrinkage MA Plots","Mean Reads Scatter Plots")),
                      
                      # Header Text
                      h2("Ready to Run DGE Analysis?"),
                      
                      # Run Button and Progress Loader
                      radioButtons("run","",c("Start Analysis","Pause Analysis/Set Parameters"),selected = "Pause Analysis/Set Parameters",inline = TRUE),
                      withSpinner(textOutput("Plots_Save"),type = 7,size = 0.5,proxy.height = "25px"),
                      ),
                  mainPanel(
                      
                      tabsetPanel(type = "tabs",
                                  tabPanel("Getting Started",
                                           h2("Instructions"),
                                           HTML("<p>Welcome to the Montgomery DESeq2 Graphics Output App, developed by the Montgomery Lab at Colorado State University. 
                                             This app is designed to convert tabulated counts files from mRNA sequencing Differential Gene Expression (DGE) experiments into user-friendly results tables and plots. 
                                             The app does not assist with upstream processes of the RNA-sequencing pipeline, such as filtering or aligment. 
                                             Compatible tabulation software includes Salmon, Kallisto, RSEM, and HTSeq.</p>"),
                                           br(),
                                           
                                           h3("Step One: Counts Files Parameters"),
                                           br(),
                                           h4("File Selection Method"),
                                           HTML("<p>If counts files are placed in the same working directory as the app.R file <strong>(best option)</strong>, select <strong>working directory</strong>.</p>
                                                <p>If counts files are in a separate directory, they can be uploaded to a temporary directory that the app can access. To use this option, select <strong>upload files</strong>.</p>"),
                                           br(),
                                           h4("Tabulation Software"),
                                           HTML("<p>Choose the software from which counts files were tabulated. Options include <strong>Salmon, Kallisto, RSEM, and HTSeq</strong>.</p>"),
                                           br(),
                                           
                                           h3("Step Two: Metadata Parameters"),
                                           br(),
                                           HTML("<p>If <strong>working directory</strong> is selected for the counts files method, a metadata csv file can be selected from the working directory, uploaded, or generated from within the app.</p> 
                                                <p>If <strong>upload files</strong> is selected for the counts files method, the metadata csv <em>must</em> be generated from scratch within the app.</p>"), 
                                           br(),
                                           HTML("<p>Any metadata files that are not built within the app interface should be csv files in the format of the following example:</p>"),
                                           tableOutput("demo_table"),
                                           HTML("<p>The first column consists of the file names of each counts file/sample, the second column describes the replicate number attached to each sample of a particular experimental condition, 
                                                the third column provides only the name of the experimental condition for each sample, and the final column describes whether each sample corresponds to a control group condition or not (TRUE/FALSE).</p>"),
                                           br(),
                                           h4("Generating the Metadata csv File (if selected)"),
                                           HTML("<p>If the counts file method is set to <strong>upload files</strong> or the metadata method is set to <strong>generate new</strong>, follow the sidebar input prompts to produce a metadata file from scratch.</p>
                                                <p>First, enter the number of <strong>experimental conditions</strong>. For each condition, the sidebar inputs will ask for <strong>(1)</strong> a name for the condition, <strong>(2)</strong> whether that condition is a control group, 
                                                and <strong>(3)</strong> which files associated with that condition. File uploads and working directory selection methods are both available depending on the file method and metadata methods selected previously 
                                                Once the DESeq2 analysis has been run, the generated metadata file will be <strong>saved in the working directory</strong> with a preceding date stamp in the file name. 
                                                Unless the counts files were uploaded within the app, this generated metadata file can be used in future analysis runs by selecting it from the working directory.</p>"),
                                           br(),
                                           
                                           h3("Step Three: Gene Table Parameters (optional)"),
                                           HTML("<p>If desired, a table can be uploaded describing the <strong>common/public gene names</strong> corresponding to the sequencing results gene names of any list of genes. 
                                             In the <strong>results tables, heatmap, and mean reads scatter plots</strong>, any common names listed in the table will replace their associated gene IDs. 
                                             This table should be a csv file with two columns in the format shown below:</p>"),
                                           tableOutput("demo_genes"),
                                           HTML("<p>Any genes for which common names are preferred should be listed in the table. The table does not need to include genes without a common name, 
                                                although if it does, the <strong>gene_ID</strong> entry for that gene should be repeated in the <strong>public_name</strong> column.</p>"),
                                           br(),
                                           
                                           h3("Step Four: Select Plots of Interest"),
                                           HTML("<p>Using the checkboxes, choose the plot types to be run or excluded from the analysis.</p>"),
                                           br(),
                                           
                                           h3("Step Five: Run DESeq2 Analysis"),
                                           HTML("<p>Once the file method, metadata method, and gene table method parameters have been set, the data is ready for DESeq2 analysis. 
                                                Press the <strong>Start Analysis</strong> button to proceed. This step can take up to <strong>5-10 minutes</strong> to complete depending on the complexity of the experimental design and the size of the counts files. 
                                                All generated plots selected for analysis will be <strong>saved in the working directory</strong>, and a progress bar will indicate the completion of each saved plot.</p>
                                                <p>To change plot selections or run a new analysis, first press the <strong>Freeze/Change Inputs</strong> button, then alter the parameters before pressing <strong>Start Analysis</strong> again.</p>"),
                                           br()
                                  ),
                                  tabPanel("Results Tables",
                                           
                                           # Header Text
                                           h2("Results Tables by Contrast"),
                                           p("Results tables display the gene-wise counts from each replicate in a given contrast, followed by the fold change value of each gene and the associated p-value (adjusted) of the negative binomial hypothesis test conducted by DESeq2. 
                                             Lower p-values indicate a lower probability of the null hypothesis that counts between the two conditions are derived from the same distributional parameters. 
                                             Furthermore, common gene names from any uploaded gene names table will appear in a column alongside standard gene IDs for the genes that aren't associated with a common name. 
                                             Using the 'DT' package, these tables are rendered as html widgets with sorting, searching, and page size customization features."),
                                           
                                           # Results table output, saved tables verification, and contrast selection
                                           withSpinner(uiOutput("results_table_ui"),type = 1),
                                           br(),
                                           fluidRow(column(4,uiOutput("control1")),column(4,uiOutput("treatment1"))),
                                           fluidRow(column(2,uiOutput("download_results_table_ui"))),
                                           br(),
                                           ),
                                  tabPanel("PCA Plot",
                                           
                                           # Header text
                                           h2("PCA Plot"),
                                           
                                           # HTML description
                                           HTML("<p>The PCA plot displays loadings of the first two principal components for each sample/biological replicate in the experimental design. 
                                           Colored by experimental condition, the points of the PCA plot provide a visualization of clustering amongst the samples, 
                                           both within conditions and across conditions.</p>"),
                                           br(),
                                           
                                           # Plot output
                                           withSpinner(uiOutput("pca_ui"),type = 1),
                                           br(),
                                           
                                           # Download button
                                           uiOutput("download_pca_ui"),
                                           br()
                                           
                                           ),
                                  tabPanel("Heatmap",
                                           # Header Text
                                           h2("Comprehensive Heatmap"),
                                           HTML("<p>The comprehensive heatmap shows log2 counts across all samples for any genes above a certain mean count threshold of three 
                                           (meaning an average of eight counts across all samples). Using the package <strong>heatmaply</strong>, 
                                           an interactive html widget is rendered with draggable zooming, panning, and hover text features. 
                                           To reset the heatmap axes, the home button in the top right corner of the widget can be pressed. Darker cells indicate higher log2 counts values.
                                           Furthermore, rows (genes) are clustered using the complete hierarchical clustering method in R (hclust) with eudclidean distances."),
                                           br(),
                                           
                                           # Heatmap output and download button
                                           withSpinner(uiOutput("heatmap_ui"),type = 1),
                                           br(),
                                           
                                           # Download button
                                           uiOutput("download_heatmap_ui"),
                                           br()
                                           ),
                                  tabPanel("Intra-Condition Scatter Plots",
                                           # Header Text
                                           h2("Intra-Condition Scatter Plots"),
                                           p("The Intra-Condition Scatter Plots display log2 counts between pairs of biological replicates within each condition of the experimental design. 
                                             Any counts in the data below 1 are replaced with a value of 1 to simplify the log2 transformation. All replicate pairs in each condition are displayed."),
                                           
                                           # Intra-Condition output and download button
                                           withSpinner(uiOutput("intra_condition_ui"),type = 1),
                                           br(),
                                           uiOutput("download_intra_ui")
                                           ),
                                  tabPanel("MA Plots",
                                           # Header Text
                                           h2("Standard MA Plots"),
                                           p("The Standard MA Plots are built from the DESeq2 package and display the fold change of a gene over its mean counts value (normalized) for a provided contrast between experimental conditions. 
                                             Genes with a statistically significant p-value (p < 0.05) are colored in blue."),
                                           
                                           # Standard MA plot output, saved plots verification, contrast selection, and download button
                                           withSpinner(uiOutput("ma_standard_ui"),type = 1),
                                           br(),
                                           fluidRow(column(4,uiOutput("control2")),column(4,uiOutput("treatment2"))),
                                           fluidRow(column(2,uiOutput("download_ma_standard_ui"))),
                                           br(),
                                           
                                           # Header Text
                                           h2("LFC Shrinkage MA Plots"),
                                           p("LFC Shrinkage MA Plots, also built from the DESeq2 package, are similar to standard MA Plots, only genes with low counts or high dispersion are assigned transformed (lower) fold changes. 
                                             Genes with a statistically significant p-value (p < 0.05) are colored in blue."),
                                           
                                           # LFC Shrinkage MA plot output, contrast selection, and download button
                                           withSpinner(uiOutput("ma_shrunk_ui"),type = 1),
                                           br(),
                                           fluidRow(column(4,uiOutput("control3")),column(4,uiOutput("treatment3"))),
                                           fluidRow(column(2,uiOutput("download_ma_shrunk_ui"))),
                                           br()
                                           ),
                                  tabPanel("Mean Reads Scatter Plots",
                                           # Header Text
                                           h2("Mean Reads Scatter Plots"),
                                           p("The Mean Reads Scatter Plots display average log2 counts across biological replicates of experimental condition pairs for a provided contrast. 
                                             Average counts with a value of 0 are assigned the value -4 following the log2 transformation. Guidelines are added to assist with visualization. 
                                             Statistically significant genes are colored blue. Axes are scaled to more closely resemble the log2 transformation. 
                                             Unlabeled tick marks, therefore, do not always represent whole number intervals between labeled tick marks. 
                                             Using the 'ggplot2' package, Mean Reads Scatter Plots are rendered as an html widget with draggable zooming, panning, and hover text features. 
                                             To reset the scatter plot axes, the home button in the top right corner of the widget can be pressed."),
                                           
                                           # Mean Reads plot output, saved plots verification, contrast selection, and download button
                                           uiOutput("mean_reads_ui"),
                                           br(),
                                           fluidRow(column(4,uiOutput("control4")),column(4,uiOutput("treatment4"))),
                                           fluidRow(column(2,uiOutput("download_mean_reads_ui"))),
                                           br()
                                           )
                                  )
                  )
                  
    )

)

server = function(input,output) {
    
# EXPERIMENTAL SETUP AND METADATA
    
    # demo samples.csv output
    output$demo_table = renderTable({
        demo_files = c("file1.results","file2.results","file3.results","file4.results","file5.results","file6.results")
        demo_replicates = c("wt_1","wt_2","wt_3","treatment_1","treatment_2","treatment_3")
        demo_conditions = c("wt","wt","wt","treatment","treatment","treatment")
        demo_controls = c("TRUE","TRUE","TRUE","FALSE","FALSE","FALSE")
        demo = data.frame(files = demo_files,replicates = demo_replicates,condition = demo_conditions,control_condition = demo_controls)
        return(demo)
    })
    
    # demo gene table output
    output$demo_genes = renderTable({
        demo = data.frame(gene_ID = c("F22F7.5","..."),public_name = c("ckb-4","..."))
        return(demo)
    })
    
    # ui for specifying the parameters of each condition
    output$condition_ui = renderUI({
        lapply(1:input$condition_number, function(i) {
            fluidRow(
                conditionalPanel(condition = "input.file_method == 'upload_files' | input.samples_method == 'generate'",
                                 textInput(paste0("condition_",i),h3(paste0("Enter name of condition ",i)))),
                conditionalPanel(condition = "input.file_method == 'upload_files' | input.samples_method == 'generate'",
                                 radioButtons(paste0("control_",i),"Is this condition a control group?",
                                              choiceNames = c("yes","no"),choiceValues = c(TRUE,FALSE),selected = FALSE)),
                conditionalPanel(condition = "input.samples_method == 'generate' & input.file_method == 'working_directory'",
                                 selectInput(paste0("files_",i,"_w"),"Select files that belong to this condition",
                                             dir(),selected = NULL,multiple = TRUE)),
                conditionalPanel(condition = "input.file_method == 'upload_files'",
                                 fileInput(paste0("files_",i,"_u"),"Upload counts files for this condition",multiple = TRUE))
            )
        })
    })
    
    # samples_generated.csv setup function based on user inputs
    samples_csv_generator = eventReactive(input$run == "Start Analysis",{
        files_list = NULL
        replicates_list = NULL
        conditions_list = NULL
        control_list = NULL
        if (input$file_method == "working_directory" & input$samples_method == "generate") {
            for (i in 1:input$condition_number) {
                files_list = c(files_list,eval(parse(text = paste0("input$files_",i,"_w"))))
                replicates_list = c(replicates_list,
                                    paste0(rep(eval(parse(text = paste0("input$condition_",i))),
                                               length(eval(parse(text = paste0("input$files_",i,"_w"))))),
                                           "_",1:length(eval(parse(text = paste0("input$files_",i,"_w"))))))
                conditions_list = c(conditions_list,
                                    rep(eval(parse(text = paste0("input$condition_",i))),
                                        length(eval(parse(text = paste0("input$files_",i,"_w"))))))
                control_list = c(control_list,
                                 rep(eval(parse(text = paste0("input$control_",i))),
                                     length(eval(parse(text = paste0("input$files_",i,"_w"))))))
            }
            return(cbind(files_list,replicates_list,conditions_list,control_list))
        } else if (input$file_method == "upload_files" & input$samples_method == "generate") {
            for (i in 1:input$condition_number) {
                files_list = c(files_list,eval(parse(text = paste0("input$files_",i,"_u$datapath"))))
                replicates_list = c(replicates_list,
                                    paste0(rep(eval(parse(text = paste0("input$condition_",i))),
                                               length(eval(parse(text = paste0("input$files_",i,"_u$datapath"))))),
                                           "_",1:length(eval(parse(text = paste0("input$files_",i,"_u$datapath"))))))
                conditions_list = c(conditions_list,
                                    rep(eval(parse(text = paste0("input$condition_",i))),
                                        length(eval(parse(text = paste0("input$files_",i,"_u$datapath"))))))
                control_list = c(control_list,
                                 rep(eval(parse(text = paste0("input$control_",i))),
                                     length(eval(parse(text = paste0("input$files_",i,"_u$datapath"))))))
            }
            return(cbind(files_list,replicates_list,conditions_list,control_list))
        } else {
            return(NULL)
        }
    })
    
    # If samples.csv is newly generated, it will be saved in the working directory here
    observeEvent(input$run == "Start Analysis",{
        if (input$file_method == "upload_files" | input$samples_method == "generate") {
            write.csv(samples_csv_generator(),
                      paste0(format(Sys.Date(),format = "%d_%m_%y_"),"samples_generated.csv"),
                      row.names = FALSE,quote = TRUE)
        }
    })
    
    # Samples.csv reading function
    deseq2_samples = eventReactive(input$run == "Start Analysis",{
        
        if (input$samples_method == "working_directory") {
            samples = read.csv(input$selected_samples)
        } else if (input$samples_method == "upload") {
            samples = read.csv(input$uploaded_samples$datapath)
        } else if (input$samples_method == "generate") {
            samples = samples_csv_generator()
        }
        return(samples)
    })
    
# DESEQ2 AND CONTRAST PREPARATION
    
    # Base DESeq object function
    deseq2 = eventReactive(input$run == "Start Analysis",{
        
        samples = deseq2_samples()
        files = samples[,1]
        names(files) = samples[,2]
        conditions_set = as.factor(samples[,3])
        sampleTable = data.frame(condition = conditions_set)
        rownames(sampleTable) = samples[,2]
        
        if (input$software_method == "htseq") {
            samples = data.frame(sampleName = samples[,2],fileName = samples[,1],condition = as.factor(samples[,3]))
            dds = DESeqDataSetFromHTSeqCount(samples,getwd(),~condition)
        } else {
            txi = tximport(files,type = input$software_method,txIn = FALSE,txOut = FALSE)
            txi$length[txi$length == 0] = 1
            dds = DESeqDataSetFromTximport(txi, sampleTable, ~condition)
        }
        
        dds = DESeq(dds)
        return(dds)
    })
    
    # Function generating factor variable of condition types for contrast-based plots
    deseq2_contrasts = eventReactive(input$run == "Start Analysis",{
        samples = deseq2_samples()
        return(levels(as.factor(samples[,3])))
    })
    
    # Combination generator function for experiments with a single control
    combination_generator = function(control,conditions_set){
        rbind(rep(as.character(control),length(levels(conditions_set))-1),
              levels(conditions_set)[-which(levels(conditions_set) == as.character(control))])
    }
    
    # Adaptive combination generator for any set of control/treatment groups
    deseq2_combinations = eventReactive(input$run == "Start Analysis",{
        samples = deseq2_samples()
        conditions_set = as.factor(samples[,3])
        control_condition = as.logical(samples[,4])
        names(control_condition) = conditions_set
        control_condition = control_condition[unique(names(control_condition))]
        
        if (all(!control_condition) | sum(control_condition > 1)) {
            combinations = combn(levels(conditions_set),2)
        } else {
            combinations = combination_generator(names(which(control_condition == TRUE)),conditions_set)
        }
        return(combinations)
    })

# CSV COUNTS TABLES

    # Base table generating function
    results_table = function(cts,group1,group2) {
        if (input$gene_method == "working_directory") {
            GeneIDs = read.csv(input$selected_gene_table)
            GeneUnique = NULL
            for (i in 1:length(GeneIDs[,1])){
                if (as.character(GeneIDs[i,1]) != as.character(GeneIDs[i,2])) {
                    GeneUnique = rbind(GeneUnique,GeneIDs[i,c(1,2)])
                }
            }
        } else if (input$gene_method == "upload_table") {
            GeneIDs = read.csv(input$uploaded_gene_table$datapath)
            GeneUnique = NULL
            for (i in 1:length(GeneIDs[,1])){
                if (as.character(GeneIDs[i,1]) != as.character(GeneIDs[i,2])) {
                    GeneUnique = rbind(GeneUnique,GeneIDs[i,c(1,2)])
                }
            }
        } else {
            GeneUnique = NULL
        }
        conditions_set = as.factor(deseq2_samples()[,3])
        sampleTable = data.frame(condition = conditions_set)
        rownames(sampleTable) = deseq2_samples()[,2]
        res = results(deseq2(),contrast=c("condition",as.character(group2),as.character(group1)))
        resOrdered = res[order(res$pvalue),]
        resOrdered = as.data.frame(resOrdered)
        commonnames = rep(NA,length(rownames(resOrdered)))
        for (i in 1:length(rownames(resOrdered))) {
            if (rownames(resOrdered)[i] %in% GeneUnique[,1]) {
                commonnames[i] = GeneUnique[which(GeneUnique[,1] == rownames(resOrdered)[i]),2]
            } else {
                commonnames[i] = rownames(resOrdered)[i]
            }
        }
        resOrdered = cbind(commonnames,resOrdered)
        resOrdered = resOrdered[,-c(2,4,5,6)]
        resCounts = matrix(rep(NA,dim(cts)[2]*dim(resOrdered)[1]),
                           nrow = dim(resOrdered)[1],ncol = dim(cts)[2])
        for (i in 1:length(resOrdered[,1])) {
            resCounts[i,] = cts[which(rownames(cts) == rownames(resOrdered)[i]),]
        }
        colnames(resCounts) = colnames(cts)
        resCounts = resCounts[,c(rownames(sampleTable)[which(sampleTable == as.character(group1))],
                                 rownames(sampleTable)[which(sampleTable == as.character(group2))])]
        resOrdered = cbind(resOrdered[,1],resCounts,resOrdered[,c(2,3)])
        colnames(resOrdered)[1] = "Common_Name"
        resOrdered$log2FoldChange[which(resOrdered$log2FoldChange > 0)] = 
            2^resOrdered$log2FoldChange[which(resOrdered$log2FoldChange > 0)]
        resOrdered$log2FoldChange[which(resOrdered$log2FoldChange < 0)] = 
            -1/(2^resOrdered$log2FoldChange[which(resOrdered$log2FoldChange < 0)])
        colnames(resOrdered)[length(colnames(resOrdered))-1] = "Fold_Change"
        resOrdered[,-1] = round(resOrdered[,-1],digits = 3)
        return(resOrdered)
    }

    # Control selection
    output$control1 = renderUI({
        req(input$run == "Start Analysis")
        if ("Results Tables" %in% input$plot_selection) {
            selectInput("control1","Select Control Group",deseq2_contrasts(),selected = deseq2_combinations()[1,1])
        }
    })
    
    # Treatment selection
    output$treatment1 = renderUI({
        req(input$run == "Start Analysis")
        if ("Results Tables" %in% input$plot_selection) {
            selectInput("treatment1","Select Treatment Group",deseq2_contrasts(),selected = deseq2_combinations()[2,1])
        }
    })
    
    # Render table for a selected contrast
    output$results_table = renderDataTable({
        req(input$run == "Start Analysis")
        datatable(results_table(counts(deseq2(), normalized=TRUE),
                                as.character(input$control1),
                                as.character(input$treatment1)))
    })
    
    # Render table ui
    output$results_table_ui = renderUI({
        req(input$run == "Start Analysis")
        if ("Results Tables" %in% input$plot_selection) {
            withSpinner(dataTableOutput("results_table"),type = 1)
        } else {
            p("Select 'Results Tables' on the sidebar and re-run the DESeq2 Analysis to view Results Tables")
        }
    })
    
    # Download button for selected contrast
    output$download_results_table = downloadHandler(
        filename = paste0(format(Sys.Date(),format = "%d_%m_%y_"),"results.csv"),
        content = function(file) {
            cts = counts(deseq2(), normalized=TRUE)
            write.csv(results_table(cts,as.character(input$control1),as.character(input$treatment1)),file)
        }
    )
    
    # Download button ui
    output$download_results_table_ui = renderUI({
        req(input$run == "Start Analysis")
        if ("Results Tables" %in% input$plot_selection) {
            conditionalPanel(condition = "input.run == 'Start Analysis'",downloadButton("download_results_table","Download Results Table"))
        }
    })
    
# PCA PLOTTING
    
    # Plotting function
    pca_input = function() {
        vsd = vst(deseq2(),blind = FALSE)
        plotPCA(vsd,intgroup="condition")
    }
    
    # Render plot
    output$pca = renderPlot({
        req(input$run == "Start Analysis")
        pca_input()
    })
    
    # Render plot ui
    output$pca_ui = renderUI({
        req(input$run == "Start Analysis")
        if ("PCA Plot" %in% input$plot_selection) {
            withSpinner(plotOutput("pca",height = "600px",width = "800px"),type = 1)
        } else {
            HTML("<p>Select <strong>PCA Plot</strong> on the sidebar and re-run the DESeq2 analysis to view PCA Plot results</p>")
        }
    })
    
    # Download button
    output$download_pca = downloadHandler(
        filename = paste0(format(Sys.Date(),format = "%d_%m_%y_"),"PCA_plot.pdf"),
        content = function(file) {
            ggsave(file,plot = pca_input(),device = "pdf")
        }
    )
    
    # Download button ui
    output$download_pca_ui = renderUI({
        req(input$run == "Start Analysis")
        if ("PCA Plot" %in% input$plot_selection) {
            downloadButton("download_pca","Download PCA Plot")
        }
    })

# HEATMAP PLOTTING
    
    # Heatmap data function
    deseq2_heat = eventReactive(input$run == "Start Analysis",{
        
        # import counts data
        cts = counts(deseq2(), normalized=TRUE)
        
        # log2 transformation and removal of low-count rows
        cts_heat = cts
        cts_heat = log2(cts_heat)
        cts_heat = cts_heat[which(rowMeans(cts_heat) >= 3),]
        
        # create clustering object
        d = dist(cts_heat,method = "euclidean")
        h = hclust(d,method = "complete")
        
        # order heatmap data by clustering order
        cts_heat = cts_heat[h$order,]
        
        # output
        return(cts_heat)
    })
    
    # Base plotting function
    heatmap_input = function() {
        heatmaply(deseq2_heat(),
                  colors = c("#F7FBFF","#DEEBF7","#C6DBEF","#9ECAE1","#6BAED6","#4292C6",
                             "#2171B5","#08519C","#08306B","#03095B","#020532","#000000"),
                  showticklabels = c(TRUE,FALSE),
                  Colv = FALSE,
                  Rowv = FALSE,
                  main = "Overall Heatmap")
    }
    
    # Render plot
    output$heatmap = renderPlotly({
        req(input$run == "Start Analysis")
        heatmap_input()
    })
    
    # Render plot ui
    output$heatmap_ui = renderUI({
        req(input$run == "Start Analysis")
        if ("Heatmap" %in% input$plot_selection) {
            withSpinner(plotlyOutput("heatmap",height = "600px",width = "800px"),type = 1)
        } else {
            HTML("<p>Select <strong>Heatmap</strong> on the sidebar and re-run the DESeq2 Analysis to view Heatmap results</p>")
        }
    })
    
    # Download button
    output$download_heatmap = downloadHandler(
        filename = paste0(format(Sys.Date(),format = "%d_%m_%y_"),"Heatmap.html"),
        content = function(file) {
            htmlwidgets::saveWidget(heatmap_input(),file = file)
        }
    )
    
    # Download button ui
    output$download_heatmap_ui = renderUI({
        req(input$run == "Start Analysis")
        if ("Heatmap" %in% input$plot_selection) {
            downloadButton("download_heatmap","Download Heatmap")
        }
    })
    
# INTRA-CONDITION SCATTER PLOTTING
    
    # Base plotting function
    intra_condition_input = function() {
        
        # import metadata and counts data
        samples = deseq2_samples()
        conditions_set = as.factor(samples[,3])
        cts = counts(deseq2(), normalized=TRUE)
        
        # produce intra-condition combinations
        combinations = NULL
        for (i in levels(conditions_set)){
            columns = c(samples[which(samples[,3] == as.character(i)),2])
            combinations = cbind(combinations,combn(columns,2))
        }
        
        # produce plots in grid for each combination
        par(mfrow=c(length(levels(conditions_set)),3),mar = c(3,1,1,1), mgp = c(2,0.5,0),pty="s")
        for (i in 1:length(combinations[1,])){
            plot(log2(replace(cts,cts < 1,1))[,c(combinations[1,i],combinations[2,i])],
                 col="lightskyblue4",pch=15,cex=0.2)
        }
    }
    
    # Render plot
    output$intra_condition = renderPlot({
        req(input$run == "Start Analysis")
        intra_condition_input()
    })
    
    # Render plot ui
    output$intra_condition_ui = renderUI({
        req(input$run == "Start Analysis")
        if ("Intra-Condition Scatter Plots" %in% input$plot_selection) {
            withSpinner(plotOutput("intra_condition",height = "800px",width = "800px"),type = 1)
        } else {
            p("Select 'Intra-Condition Scatter Plot' on the sidebar and re-run the DESeq2 Analysis to view Intra-Condition Scatter Plot results")
        }
    })
    
    # Download button
    output$download_intra = downloadHandler(
        filename = paste0(format(Sys.Date(),format = "%d_%m_%y_"),"Intra_Condition_Scatterplot.pdf"),
        content = function(file) {
            ggsave(file,plot = intra_condition_input(),device = "pdf")
        }
    )
    
    # Download button ui
    output$download_intra_ui = renderUI({
        req(input$run == "Start Analysis")
        if ("Intra-Condition Scatter Plots" %in% input$plot_selection) {
            downloadButton("download_intra","Download Intra-Condition Scatterplot")
        }
    })

# STANDARD MA PLOTTING    

    # Base plotting function
    ma_standard_input = function(group1,group2) {
        res = results(deseq2(),contrast=c("condition",group2,group1))
        plotMA(res, ylim=c(-4,4),main = paste0(group2,' vs ',group1))
    }
    
    # UI contrast selection
    output$control2 = renderUI({
        req(input$run == "Start Analysis")
        if ("Standard MA Plots" %in% input$plot_selection) {
            selectInput("control2","Select Control Group",deseq2_contrasts(),selected = deseq2_combinations()[1,1])
        }
    })
    
    output$treatment2 = renderUI({
        req(input$run == "Start Analysis")
        if ("Standard MA Plots" %in% input$plot_selection) {
            selectInput("treatment2","Select Treatment Group",deseq2_contrasts(),selected = deseq2_combinations()[2,1])
        }
    })
    
    # Render plot for a selected contrast in ui
    output$ma_standard = renderPlot({
        req(input$run == "Start Analysis")
        ma_standard_input(as.character(input$control2),as.character(input$treatment2))
    })
    
    # UI object for rendered plot
    output$ma_standard_ui = renderUI({
        req(input$run == "Start Analysis")
        if ("Standard MA Plots" %in% input$plot_selection) {
            withSpinner(plotOutput("ma_standard",height = "800px",width = "800px"),type = 1)
        } else {
            p("Select 'Standard MA Plots' on the sidebar and re-run the DESeq2 Analysis to view Standard MA Plot results")
        }
    })
    
    # Download button for selected plot in ui
    output$download_ma_standard = downloadHandler(
        filename = paste0(format(Sys.Date(),format = "%d_%m_%y_"),"MA_plot.pdf"),
        content = function(file) {
            ggsave(file,plot = ma_standard_input(as.character(input$control2),as.character(input$treatment2)),
                   device = "pdf")
        }
    )
    
    # UI object for download button
    output$download_ma_standard_ui = renderUI({
        req(input$run == "Start Analysis")
        if ("Standard MA Plots" %in% input$plot_selection) {
            conditionalPanel(condition = "input.run == 'Start Analysis'",downloadButton("download_ma_standard","Download Standard MA Plot"))
        }
    })
    
# LFC SHRINKAGE MA PLOTTING
    
    # Base plotting function
    ma_shrunk_input = function(group1,group2) {
        resLFC = lfcShrink(deseq2(),contrast=c("condition",group2,group1),type="normal")
        plotMA(resLFC, ylim=c(-4,4),main = paste0(group2,' vs ',group1))
    }
    
    # UI contrast selection
    output$control3 = renderUI({
        req(input$run == "Start Analysis")
        if ("LFC Shrinkage MA Plots" %in% input$plot_selection) {
            selectInput("control3","Select Control Group",deseq2_contrasts(),selected = deseq2_combinations()[1,1])
        }
    })
    
    output$treatment3 = renderUI({
        req(input$run == "Start Analysis")
        if ("LFC Shrinkage MA Plots" %in% input$plot_selection) {
            selectInput("treatment3","Select Treatment Group",deseq2_contrasts(),selected = deseq2_combinations()[2,1])
        }
    })
    
    # Render plot for a selected contrast in ui
    output$ma_shrunk = renderPlot({
        req(input$run == "Start Analysis")
        ma_shrunk_input(as.character(input$control3),as.character(input$treatment3))
    })
    
    # UI object for rendered plot
    output$ma_shrunk_ui = renderUI({
        req(input$run == "Start Analysis")
        if ("LFC Shrinkage MA Plots" %in% input$plot_selection) {
            withSpinner(plotOutput("ma_shrunk",height = "800px",width = "800px"),type = 1)
        } else {
            p("Select 'LFC Shrinkage MA Plots' on the sidebar and re-run the DESeq2 Analysis to view LFC Shrinkage MA Plot results")
        }
    })
    
    # Download button for selected plot in ui
    output$download_ma_shrunk = downloadHandler(
        filename = paste0(format(Sys.Date(),format = "%d_%m_%y_"),"MA_Shrunk_plot.pdf"),
        content = function(file) {
            ggsave(file,plot = ma_shrunk_input(as.character(input$control3),as.character(input$treatment3)),
                   device = "pdf")
        }
    )
    
    # UI object for download button
    output$download_ma_shrunk_ui = renderUI({
        req(input$run == "Start Analysis")
        if ("LFC Shrinkage MA Plots" %in% input$plot_selection) {
            conditionalPanel(condition = "input.run == 'Start Analysis'",downloadButton("download_ma_shrunk","Download LFC Shrinkage MA Plot"))
        }
    })
    
# MEAN READS SCATTER PLOTTING
    
    # Base plotting function
    mean_reads_input = function(group1,group2) {
        cts = counts(deseq2(), normalized=TRUE)
        samples = deseq2_samples()
        conditions_set = levels(as.factor(samples[,3]))
        cts_avg = NULL
        for (i in conditions_set) {
            cts_avg = cbind(cts_avg,
                            log2(rowMeans(cts[,samples[,2][samples[,3] == i]])))
        }
        cts_avg = replace(cts_avg,cts_avg=='-Inf',-4)
        colnames(cts_avg) = conditions_set
        plot_colors = c("lightgray","deepskyblue3")
        
        res = results(deseq2(),contrast=c("condition",group2,group1))
        res_sig = as.data.frame(subset(res,padj < 0.05 & abs(log2FoldChange) > 0.378511623))
        sig_logical = rownames(cts_avg) %in% rownames(res_sig)
        
        plot_colors = c("lightgray","deepskyblue3")
        
        tick_sep = seq(-4,20,0.8)
        tick_sep = 2^tick_sep
        tick_sep[1:6] = seq(tick_sep[1],tick_sep[6],length.out = 6)
        tick_sep[6:11] = seq(tick_sep[6],tick_sep[11],length.out = 6)
        tick_sep[11:16] = seq(tick_sep[11],tick_sep[16],length.out = 6)
        tick_sep[16:21] = seq(tick_sep[16],tick_sep[21],length.out = 6)
        tick_sep[21:26] = seq(tick_sep[21],tick_sep[26],length.out = 6)
        tick_sep[26:31] = seq(tick_sep[26],tick_sep[31],length.out = 6)
        tick_labs = rep("",31)
        tick_labs[seq(1,31,5)] = prettyNum(log2(tick_sep[seq(1,31,5)]),digits=2,format="f")
        
        x = cts_avg[,group1]
        y = cts_avg[,group2]
        n = rownames(cts_avg)
        for (i in 1:length(n)) {
            n[i] = paste0("Gene: ",n[i])
        }
        
        mean_reads = ggplot() + 
            geom_point(aes(x,y,text = n),size = 0.5+(0.2*sig_logical),
                       color = plot_colors[1+sig_logical]) +
            geom_abline(intercept = 0,slope = 1,color = "grey60") + 
            geom_abline(intercept = 1,slope = 1,color = "grey60") +
            geom_abline(intercept = -1,slope = 1,color = "grey60") + 
            geom_hline(yintercept = log2(10),color = "grey60",linetype = "dashed") + 
            geom_vline(xintercept = log2(10),color = "grey60",linetype = "dashed") + 
            scale_y_continuous(breaks = log2(tick_sep),labels = tick_labs,
                               limits = c(1/16,max(cts_avg[,c(group1,group2)]))) + 
            scale_x_continuous(breaks = log2(tick_sep),labels = tick_labs,
                               limits = c(1/16,max(cts_avg[,c(group1,group2)]))) + 
            ggtitle(paste0(group2,' vs ',group1)) +
            xlab(paste0('Log2 Average reads in ',group1)) +
            ylab(paste0('Log2 Average reads in ',group2)) +
            theme_classic() + 
            theme(plot.title = element_text(hjust = 0.5))
    }
    
    # UI contrast selection
    output$control4 = renderUI({
        req(input$run == "Start Analysis")
        if ("Mean Reads Scatter Plots" %in% input$plot_selection) {
            selectInput("control4","Select Control Group",deseq2_contrasts(),selected = deseq2_combinations()[1,1])
        }
    })
    
    output$treatment4 = renderUI({
        req(input$run == "Start Analysis")
        if ("Mean Reads Scatter Plots" %in% input$plot_selection) {
            selectInput("treatment4","Select Treatment Group",deseq2_contrasts(),selected = deseq2_combinations()[2,1])
        }
    })

    # Render plot for a selected contrast in ui
    output$mean_reads = renderPlotly({
        req(input$run == "Start Analysis")
        ggplotly(mean_reads_input(input$control4,input$treatment4))
    })
    
    # UI object for rendered plot
    output$mean_reads_ui = renderUI({
        req(input$run == "Start Analysis")
        if ("Mean Reads Scatter Plots" %in% input$plot_selection) {
            withSpinner(plotlyOutput("mean_reads",height = "600px",width = "800px"),type = 1)
        } else {
            p("Select 'Mean Reads Scatter Plots' on the sidebar and re-run the DESeq2 Analysis to view Mean Reads Scatter Plot results")
        }
    })
    
    # Download button for selected plot in ui
    output$download_mean_reads = downloadHandler(
        filename = paste0(format(Sys.Date(),format = "%d_%m_%y_"),"Mean_Reads_plot.html"),
        content = function(file) {
            htmlwidgets::saveWidget(ggplotly(mean_reads_input(as.character(input$control4),as.character(input$treatment4))),file = file)
        }
    )
    
    # UI object for download button
    output$download_mean_reads_ui = renderUI({
        req(input$run == "Start Analysis")
        if ("Mean Reads Scatter Plots" %in% input$plot_selection) {
            conditionalPanel(condition = "input.run == 'Start Analysis'",downloadButton("download_mean_reads","Download Mean Reads Scatterplot"))
        }
    })
    
    # Save all plots in working directory
    
    output$Plots_Save = renderText({
        
        # Setup
        req(input$run == "Start Analysis")
        combinations = deseq2_combinations()
        cts = counts(deseq2(), normalized=TRUE)
        
        # Progress Bar Initialization
        withProgress(message = "Saving PCA Plot",{
            
            # Results Tables Saving
            if ("Results Tables" %in% input$plot_selection) {
                incProgress(0.1,message = "Saving Results Tables")
                for (i in 1:length(combinations[1,])){
                    res = results_table(cts,combinations[1,i],combinations[2,i])
                    write.csv(res,paste0(format(Sys.Date(),format = "%d_%m_%y_"),
                                         as.character(combinations[2,i]),'_vs_',
                                         as.character(combinations[1,i]),'.csv'))
                }
            } else {
                incProgress(0.1,message = "Results Tables Not Selected")
            }
            
            Sys.sleep(0.5)

            # PCA Plot Saving
            if ("PCA Plot" %in% input$plot_selection) {
                incProgress(0.1,message = "Saving PCA Plot")
                ggsave(paste0(format(Sys.Date(),format = "%d_%m_%y_"),"PCA_plot.pdf"),
                       plot = pca_input(),device = "pdf",
                       width = 6.5,height = 6.5,units = "in")
            } else {
                incProgress(0.1,message = "PCA Plot Not Selected")
            }
            
            Sys.sleep(0.5)
            
            # Heatmap Saving
            if ("Heatmap" %in% input$plot_selection) {
                incProgress(0.1,message = "Saving Heatmap")
                htmlwidgets::saveWidget(heatmap_input(),
                                        file = paste0(format(Sys.Date(),format = "%d_%m_%y_"),"Heatmap.html"))
            } else {
                incProgress(0.1,message = "Heatmap Not Selected")
            }
            
            Sys.sleep(0.5)
            
            # Intra-Condition Scatter Plots Saving
            if ("Intra-Condition Scatter Plots" %in% input$plot_selection) {
                incProgress(0.1,message = "Saving Intra-Condition Plot")
                ggsave(paste0(format(Sys.Date(),format = "%d_%m_%y_"),"Intra_Condition_Scatterplot.pdf"),
                       plot = intra_condition_input(),device = "pdf",
                       width = 6.5,height = 6.5,units = "in")
            } else {
                incProgress(0.1,message = "Intra-Condition Plots Not Selected")
            }
            
            Sys.sleep(0.5)
            
            # Standard MA Plots Saving
            if ("Standard MA Plots" %in% input$plot_selection) {
                incProgress(0.1,message = "Saving Standard MA Plots")
                for (i in 1:length(combinations[1,])){
                    ggsave(paste0(format(Sys.Date(),format = "%d_%m_%y_"),
                                  combinations[2,i],'_vs_',combinations[1,i],'_MA.pdf'),
                           plot = ma_standard_input(combinations[1,i],combinations[2,i]),
                           device = "pdf",width = 6.5,height = 6.5,units = "in")
                }
            } else {
                incProgress(0.1,message = "Standard MA Plots Not Selected")
            }
            
            Sys.sleep(0.5)
            
            # LFC Shrinkage MA Plots Saving
            if ("LFC Shrinkage MA Plots" %in% input$plot_selection) {
                incProgress(0.1,message = "Saving LFC Shrinkage MA Plots")
                for (i in 1:length(combinations[1,])){
                    ggsave(paste0(format(Sys.Date(),format = "%d_%m_%y_"),
                                  combinations[2,i],'_vs_',combinations[1,i],'_MA_Shrunk.pdf'),
                           plot = ma_shrunk_input(combinations[1,i],combinations[2,i]),
                           device = "pdf",width = 6.5,height = 6.5,units = "in")
                }
            } else {
                incProgress(0.1,message = "LFC Shrinkage Plots Not Selected")
            }
            
            Sys.sleep(0.5)

            # Mean Reads Scatter Plots Saving
            if ("Mean Reads Scatter Plots" %in% input$plot_selection) {
                incProgress(0.1,message = "Saving Mean Reads Scatter Plots")
                for (i in 1:length(combinations[1,])){
                    ggsave(paste0(format(Sys.Date(),format = "%d_%m_%y_"),
                                  combinations[2,i],'_vs_',combinations[1,i],'_Mean_Reads.pdf'),
                           plot = mean_reads_input(combinations[1,i],combinations[2,i]),
                           device = "pdf",width = 6.5,height = 6.5,units = "in")
                }
            } else {
                incProgress(0.1,message = "Mean Reads Scatter Plots Not Selected")
            }
            
            Sys.sleep(0.5)
            
            # Completion Message
            incProgress(0.2,message = "Complete")
        })
        return("All Selected Tables and Plots Saved!")
    })   
}

# Run the application
shinyApp(ui = ui,server = server)
