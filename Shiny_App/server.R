# Define server logic required to draw a histogram
server = function(input,output) {
  
  rv = reactiveValues()
  rv$metadata = NULL
  rv$gene_table = NULL
  rv$common_names = NULL
  rv$class_names = NULL
  rv$class_parameters = NULL
  rv$point_color = NULL
  rv$point_size = NULL
  rv$sample_table = NULL
  rv$dds = NULL
  rv$combinations = NULL
  rv$res = NULL
  rv$cts = NULL
  rv$cts_names = NULL
  rv$cts_avg = NULL
  rv$cts_heat = NULL
  rv$tables = NULL
  rv$mean_reads_plots = NULL
  
  ## Experimental Setup and Metadata
  
  observeEvent(input$run,{
    
    # Table and Plot Saving
    withProgress(message = "Initializing",{
      
      # Create Output Directory
      dir.create(paste0(getwd(),"/",format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/"))
      
      # Compile List of Parameters
      
      # Conditional Metadata Entry
      if (input$metadata_method == "generate") {
        metadata_param = paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_metadata_generated.csv")
      } else {
        metadata_param = input$selected_metadata
      }
      
      # Conditional Class Parameters Entry
      if (input$generate_class_parameters) {
        class_parameters_param = paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_class_parameters_generated.csv")
      } else {
        class_parameters_param = input$select_class_parameters
      }
      
      # Parameter List
      params = list(
        "experiment_id" = as.character(input$experiment_id),
        "software_method" = input$software_method,
        "counts_matrix" = input$selected_counts_matrix,
        "metadata" = metadata_param,
        "gene_table_method" = input$gene_table_method,
        "gene_table" = input$selected_gene_table,
        "generate_results_tables" = input$results_tables_render,
        "generate_pca" = input$pca_render,
        "generate_intra_condition" = input$intra_condition_render,
        "generate_mean_reads" = input$mean_reads_render,
        "save_mean_reads_interactive" = FALSE,
        "p_value_threshold" = input$p_value_threshold,
        "fold_change_threshold" = input$fold_change_threshold,
        "lower_transparency" = input$lower_transparency,
        "upper_transparency" = input$upper_transparency,
        "customize_by_class" = input$customize_by_class,
        "customize_by_significance" = input$customize_by_significance,
        "class_parameters" = class_parameters_param,
        "generate_ma" = input$ma_render,
        "generate_heatmap" = input$heatmap_render,
        "heatmap_type" = input$heatmap_type,
        "heatmap_selected_classes" = paste0(input$heatmap_selected_classes,collapse = ",")
      )
      
      # Save Parameters YAML
      yaml::write_yaml(
        params,
        paste0(
          format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
          format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_params.yml"
        )
      )
      
      # Generate Metadata if Necessary
      incProgress(0.05,message = "compiling metadata")
      if (input$metadata_method == "generate") {
        
        # Column Setup
        files_list = NULL
        replicates_list = NULL
        conditions_list = NULL
        control_list = NULL
        
        for (i in 1:input$condition_number) {
          
          # Files Column
          files_list = c(files_list,eval(parse(text = paste0("input$files_",i))))
          
          # Replicates Column
          replicates_list = c(
            replicates_list,
            paste0(
              rep(eval(parse(text = paste0("input$condition_",i))),length(eval(parse(text = paste0("input$files_",i))))),
              "_",1:length(eval(parse(text = paste0("input$files_",i))))
            )
          )
          
          # Conditions Column
          conditions_list = c(
            conditions_list,
            rep(eval(parse(text = paste0("input$condition_",i))),length(eval(parse(text = paste0("input$files_",i)))))
          )
          
          # Controls Column
          control_list = c(
            control_list,
            rep(eval(parse(text = paste0("input$control_",i))),length(eval(parse(text = paste0("input$files_",i)))))
          )
        }
        
        # Object Return
        rv$metadata = cbind(files_list,replicates_list,conditions_list,control_list)
        
        # Save Metadata
        write.csv(
          rv$metadata,
          paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_metadata_generated.csv"),
          row.names = FALSE,quote = TRUE
        )
      } else if (input$metadata_method == "working_directory") {
        rv$metadata = read.csv(input$selected_metadata)
      }
      
      # Conditional Gene Table Import
      incProgress(0.05,message = "compiling gene table")
      if (input$selected_gene_table == "") {
        rv$gene_table = NULL
      } else {
        rv$gene_table = read.csv(input$selected_gene_table)
      }
      
      # Gene Table Values
      if (is.null(rv$gene_table)) {
        
        # No Gene Table Option
        rv$common_names = NULL
        rv$class_names = NULL
        
      } else {
        if (input$gene_table_method == "full_table") {
          
          # Full Gene Table Option
          rv$common_names = repeated_values_solver(rv$gene_table[,2])
          names(rv$common_names) = repeated_values_solver(rv$gene_table[,1])
          rv$class_names = rv$gene_table[,3]
          names(rv$class_names) = repeated_values_solver(rv$gene_table[,1])
          
        } else if (input$gene_table_method == "common_names_only") {
          
          # Common Names Only Option
          rv$common_names = repeated_values_solver(rv$gene_table[,2])
          names(rv$common_names) = repeated_values_solver(rv$gene_table[,1])
          rv$class_names = NULL
          
        } else if (input$gene_table_method == "gene_class_only") {
          
          # Gene Class Only Option
          rv$common_names = NULL
          rv$class_names = rv$gene_table[,2]
          names(rv$class_names) = repeated_values_solver(rv$gene_table[,1])
        }
      }
      
      # Class Parameters Metadata
      incProgress(0.05,message = "compiling class parameters")
      if (input$generate_class_parameters & !is.null(input$mean_reads_selected_classes[1])) {
        
        # Generated Class Parameters Option
        
        color_list = NULL
        size_list = NULL
        
        for (i in 1:length(input$mean_reads_selected_classes)) {
          color_list = c(color_list,eval(parse(text = paste0("input$color_",i))))
          size_list = c(size_list,eval(parse(text = paste0("input$size_",i))))
        }
        
        # Object Return
        rv$class_parameters = data.frame(
          "point_class" = input$mean_reads_selected_classes,
          "point_colors" = color_list,
          "point_sizes" = size_list
        )
        
        # Save Class Parameters
        write.csv(
          rv$class_parameters,
          paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_class_parameters_generated.csv"),
          row.names = FALSE,quote = TRUE
        )
        
      } else if (input$select_class_parameters != "") {
        
        # Imported Class Parameters Option
        rv$class_parameters = read.csv(input$select_class_parameters)
        
      } else if (input$select_class_parameters == "" & input$customize_by_class) {
        
        # Default Class Parameters Option
        
        # Default Colors List
        class_colors = c(
          "#1F78B4","#1B9E77","#33A02C","#014421","#E6AB02","#FF7F00","#FF9E84","#8B0000",
          "#E31A1C","#E7298A","#6A3D9A","#907567","#A6761D","#525252","#000000"
        )
        
        # Default Table Generation
        rv$class_parameters = data.frame(
          point_class = levels(as.factor(gene_table()[[2]])),
          point_colors = rep(class_colors,ceiling(length(levels(as.factor(gene_table()[[2]])))/length(class_colors)))[1:length(levels(as.factor(gene_table()[[2]])))],
          point_sizes = rep(0.5,length(levels(as.factor(gene_table()[[2]]))))
        )
        
        # Default Table Saving
        write.csv(
          rv$class_parameters,
          paste0(format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_class_parameters_generated.csv"),
          row.names = FALSE
        )
        
      } else {
        
        # Null Class Parameters
        rv$class_parameters = NULL
        
      }
      
      # Plot Parameter Values Preparation
      if (is.null(rv$class_parameters)) {
        
        # No Plot Parameters Option
        rv$point_color = NULL
        rv$point_size = NULL
        
      } else {
        
        # Full Plot Parameters Option
        rv$point_color = rv$class_parameters[,2]
        names(rv$point_color) = rv$class_parameters[,1]
        rv$point_size = rv$class_parameters[,3]
        names(rv$point_size) = rv$class_parameters[,1]
      }
      
      # Sample Table Preparation
      rv$sample_table = data.frame(condition = as.factor(rv$metadata[,3]))
      rownames(rv$sample_table) = rv$metadata[,2]
      
      # Conditional DESeq Data Import Calls
      incProgress(0.05,message = "assembling counts data")
      if (input$software_method == "htseq") {
        
        # HTSeq Option
        metadata_htseq = data.frame(sampleName = rv$metadata[,2],fileName = rv$metadata[,1],condition = factor(rv$metadata[,3]))
        rv$dds = DESeqDataSetFromHTSeqCount(rv$metadata_htseq,getwd(),~condition)
        
      } else if (input$software_method == "counts_matrix") {
        
        # Counts Matrix/FeatureCounts Option
        counts_matrix = counts_matrix_prep(read.csv(input$selected_counts_matrix))
        colnames(counts_matrix) = rownames(rv$sample_table)
        rv$dds = DESeqDataSetFromMatrix(counts_matrix,rv$sample_table,~condition)
        
      } else if (input$software_method == "tiny_rna") {
        
        # tiny RNA Option
        tiny_rna_matrix = tiny_rna_matrix_prep(read.csv(input$selected_counts_matrix))[,-c(1,2)]
        colnames(tiny_rna_matrix) = rownames(rv$sample_table)
        rv$dds = DESeqDataSetFromMatrix(tiny_rna_matrix,rv$sample_table,~condition)
        
      } else {
        
        # RSEM, Salmon, and Kallisto Options
        txi = tximport(rv$metadata[,1],type = input$software_method,txIn = FALSE,txOut = FALSE)
        txi$length[txi$length == 0] = 1
        rv$dds = DESeqDataSetFromTximport(txi,rv$sample_table,~condition)
        
      }
      
      # DESeq2 Analysis
      incProgress(0.05,message = "running DESeq2")
      rv$dds = DESeq(rv$dds)
      
      # Control Condition Status Preparation
      conditions_set = as.factor(rv$metadata[,3])
      control_status = as.logical(rv$metadata[,4])
      names(control_status) = conditions_set
      control_status = control_status[unique(names(control_status))]
      control_group = names(which(control_status == TRUE))
      
      # Combination Generation
      if (all(!control_status) | sum(control_status > 1)) {
        rv$combinations = combn(levels(conditions_set),2)
      } else {
        rv$combinations = data.frame(
          rbind(
            rep(control_group,length(levels(conditions_set))-1),
            levels(conditions_set)[-which(levels(conditions_set) == control_group)]
          )
        )
      }
      
      # Column Naming
      for (i in 1:length(rv$combinations[1,])) {
        colnames(rv$combinations)[i] = paste0(rv$combinations[2,i]," vs ",rv$combinations[1,i])
      }
      
      # List Generation
      rv$res = list()
      for (i in 1:length(rv$combinations[1,])) {
        rv$res = append(rv$res,list(results(rv$dds,contrast=c("condition",rv$combinations[2,i],rv$combinations[1,i]))))
      }
      
      rv$cts = counts(rv$dds,normalized = TRUE)
      
      # Replace Gene IDs with Common Names (if Applicable)
      incProgress(0.05,message = "compiling common names")
      rv$cts_names = rownames(rv$cts)
      if (!is.null(rv$common_names)) {
        for (i in 1:length(rv$cts_names)) {
          if (rv$cts_names[i] %in% names(rv$common_names)) {
            rv$cts_names[i] = rv$common_names[rv$cts_names[i]]
          }
        }
      }
      
      # Reconcile Repeated Values
      rv$cts_names = repeated_values_solver(rv$cts_names)
      
      # Assign Gene IDs to Common Names
      names(rv$cts_names) = rownames(rv$cts)
      
      # Average Counts Matrix Setup
      rv$cts_avg = NULL
      for (i in levels(as.factor(rv$metadata[,3]))) {
        rv$cts_avg = cbind(rv$cts_avg,log2(rowMeans(rv$cts[,rv$metadata[,2][rv$metadata[,3] == i]])))
      }
      
      # Replace Infinite Values With -4
      rv$cts_avg = replace(rv$cts_avg,rv$cts_avg=='-Inf',-4)
      
      # Assign Appropriate Column Names
      colnames(rv$cts_avg) = levels(as.factor(rv$metadata[,3]))
      
      # Heatmap Counts Matrix
      rv$cts_heat = log2(counts(rv$dds,normalized = TRUE))[which(rowMeans(log2(rv$cts)) >= 3),]
      
      # Tables List Setup
      incProgress(0.05,message = "generating results tables")
      rv$tables = list()
      
      # Conditional Rendering with Progress Bar
      if (input$results_tables_render) {
        for (i in 1:length(rv$combinations[1,])) {
          rv$tables = append(rv$tables,list(results_table(rv$sample_table,rv$cts,rv$common_names,rv$class_names,rv$res[[i]],rv$combinations[1,i],rv$combinations[2,i])))
        }
      }
      
      # Isolate Specific Classes
      if (!is.null(rv$point_color)) {
        rv$cts_avg = rv$cts_avg[intersect(names(rv$class_names[rv$class_names %in% names(rv$point_color)]),names(rv$cts_names)),]
      }
      
      # Plots List Setup
      incProgress(0.05,message = "generating mean reads plots")
      rv$mean_reads_plots = list()
      
      # Conditional Plot Saving
      if (input$mean_reads_render) {
        for (i in 1:length(rv$combinations[1,])) {
          rv$mean_reads_plots = append(
            rv$mean_reads_plots,
            list(
              mean_reads_input(
                input$p_value_threshold,
                input$fold_change_threshold,
                input$lower_transparency,
                input$upper_transparency,
                rv$cts_avg,
                rv$cts_names,
                rv$class_names,
                input$customize_by_significance,
                rv$point_color,
                rv$point_size,
                rv$res[[i]],
                rv$combinations[1,i],
                rv$combinations[2,i]
              )
            )
          )
        }
      }
      
      if (input$results_tables_render) {
        for (i in 1:length(rv$combinations[1,])){
          incProgress(0.10/length(rv$combinations[1,]),message = "saving results tables")
          write.csv(
            rv$tables[[i]],
            paste0(
              format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
              format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_",
              rv$combinations[2,i],"_vs_",rv$combinations[1,i],"_Results.csv"
            )
          )
        }
      } else {
        incProgress(0.10,message = "results tables not selected")
      }
      
      if (input$pca_render) {
        incProgress(0.05,message = "saving pca plot")
        pdf(
          paste0(
            format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
            format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_PCA_Plot.pdf"
          ),
          title = paste0(input$experiment_id," PCA Plot"),
          height = 6.5,width = 8)
        
        print(pca_input(rv$dds))
        
        dev.off()
      } else {
        incProgress(0.05,message = "pca plot not selected")
      }
      
      if (input$intra_condition_render) {
        incProgress(0.05,message = "saving intra-condition scatter plot")
        pdf(
          paste0(
            format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
            format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_Intra_Condition.pdf"
          ),
          title = paste0(input$experiment_id," Intra Condition Scatter Plot"),
          height = 6.5,width = 8
        )
        
        print(intra_condition_input(rv$sample_table,rv$cts))
        
        dev.off()
      } else {
        incProgress(0.05,message = "intra-condition scatter plots not selected")
      }
      
      if (input$mean_reads_render) {
        for (i in 1:length(rv$combinations[1,])){
          incProgress(0.10/length(rv$combinations[1,]),message = "saving mean reads scatter plots")
          pdf(
            paste0(
              format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
              format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_",
              rv$combinations[2,i],"_vs_",rv$combinations[1,i],"_Mean_Reads.pdf"
            ),
            title = paste0(input$experiment_id," ",rv$combinations[2,i]," vs ",rv$combinations[1,i]," Mean Reads Scatter Plot"),
            height = 6.5,width = 8
          )
          
          print(rv$mean_reads_plots[[i]])
          
          dev.off()
        }
      } else {
        incProgress(0.10,message = "mean reads scatter plots not selected")
      }
      
      if (input$ma_render) {
        for (i in 1:length(rv$combinations[1,])){
          incProgress(0.10/length(rv$combinations[1,]),message = "saving ma plots")
          pdf(
            paste0(
              format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
              format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_",
              rv$combinations[2,i],"_vs_",rv$combinations[1,i],"_MA.pdf"
            ),
            
            title = paste0(input$experiment_id," ",rv$combinations[2,i]," vs ",rv$combinations[1,i]," MA Plot"),
            height = 6.5,width = 8
          )
          
          print(ma_input(rv$res[[i]],rv$combinations[1,i],rv$combinations[2,i]))
          
          dev.off()
        }
      } else {
        incProgress(0.10,message = "ma plots not selected")
      }
      
      if (input$heatmap_render) {
        incProgress(0.10,message = "saving heatmap")
        if (input$heatmap_type == "complete") {
          
          # Complete Heatmap Plot Saving
          htmlwidgets::saveWidget(
            complete_heatmap(rv$cts_heat,rv$cts_names),
            file = paste0(
              format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
              format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_Heatmap.html"
            )
          )
          
        } else if (input$heatmap_type == "all_classes") {
          
          # All Classes Heatmap Plot Saving
          htmlwidgets::saveWidget(
            all_classes_heatmap(rv$cts_heat,rv$cts_names,rv$class_names),
            file = paste0(
              format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
              format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_Heatmap.html"
            )
          )
          
        } else if (input$heatmap_type == "selected_classes") {
          
          # Selected Classes Heatmap Plot Saving
          htmlwidgets::saveWidget(
            selected_classes_heatmap(rv$cts_heat,rv$cts_names,rv$class_names,input$heatmap_selected_classes),
            file = paste0(
              format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
              format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_Heatmap.html"
            )
          )
        }
      } else {
        incProgress(0.10,message = "heatmap not selected")
      }
      
      # Delete Accessory Files
      unlink(
        paste0(
          format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_output/",
          format(Sys.Date(),format = "%d_%m_%y_"),input$experiment_id,"_Heatmap_files/"
        ),
        recursive = TRUE
      )
    })
  })
  
  # UI for Specifying the Parameters of Each Condition
  output$condition_ui = renderUI({
    
    # Looping Condition Parameters UI
    lapply(1:input$condition_number,function(i) {
      fluidRow(
        
        # Condition Name
        conditionalPanel(
          condition = "input.metadata_method == 'generate'",
          textInput(
            inputId = paste0("condition_",i),
            label = h3(paste0("Enter name of condition ",i))
          )
        ),
        
        # Control Group Specification
        conditionalPanel(
          condition = "input.metadata_method == 'generate'",
          radioButtons(
            inputId = paste0("control_",i),
            label = "Is this condition a control group?",
            choiceNames = c("yes","no"),
            choiceValues = c(TRUE,FALSE),
            selected = FALSE
          )
        ),
        
        # Working Directory File Selections
        conditionalPanel(
          condition = "input.metadata_method == 'generate'",
          selectInput(
            inputId = paste0("files_",i),
            label = "Select files that belong to this condition",
            choices = c("select files" = "",dir()),
            selected = character(0),
            multiple = TRUE
          )
        )
      )
    })
  })
  
  output$class_parameters_ui = renderUI({
    conditionalPanel(
      condition = "input.generate_class_parameters & input.selected_gene_table != ''",
      checkboxGroupInput(
        inputId = "mean_reads_selected_classes",
        label = "select classes to customize",
        choices = levels(as.factor(rv$class_names))
      )
    )
  })
  
  # UI for Specifying the Mean Reads Scatter Plot Parameters
  output$class_ui = renderUI({
    
    if (!is.null(input$mean_reads_selected_classes[1])) {
      
      # Looping Condition Parameters UI
      lapply(1:length(input$mean_reads_selected_classes),function(i) {
        
        fluidRow(
          
          # Color Specification
          conditionalPanel(
            condition = "input.generate_class_parameters",
            selectInput(
              inputId = paste0("color_",i),
              label = paste0("Color of genes/points for ",input$mean_reads_selected_classes[i]),
              choices = c("select a color" = "","light blue" = "#A6CEE3","blue" = "#1F78B4","teal" = "#1B9E77",
                          "light green" = "#B2DF8A","green" = "#33A02C","lime green" = "#66A61E",
                          "light red" = "#FB9A99","red" = "#E31A1C","bright pink" = "#E7298A",
                          "light orange" = "#FDBF6F","orange" = "#FF7F00","burnt orange" = "#D95F02",
                          "light purple" = "#CAB2D6","purple" = "#6A3D9A","plum" = "#7570B3",
                          "gold" = "#E6AB02","bronze" = "#A6761D","dark grey" = "#525252","black" = "#000000"),
              selected = character(0)
            )
          ),
          
          # Size Specification
          conditionalPanel(
            condition = "input.generate_class_parameters",
            sliderInput(
              inputId = paste0("size_",i),
              label = paste0("Size of genes/points for class ",input$mean_reads_selected_classes[i]),
              value = 0.5,min = 0.1,max = 0.9,step = 0.05
            )
          )
        )
      })
    }
  })
  
  # Heatmap Class Selection UI
  output$heatmap_classes_ui = renderUI({
    if (input$heatmap_type == "selected_classes") {
      if (input$gene_table_method == "full_table") {
        selectInput(
          inputId = "heatmap_selected_classes",
          label = "select classes to isolate",
          choices = c("select classes" = "",levels(as.factor(read.csv(input$selected_gene_table)[,3]))),
          selected = character(0),
          multiple = TRUE
        )
      } else if (input$gene_table_method == "gene_class_only") {
        selectInput(
          inputId = "heatmap_selected_classes",
          label = "select classes to isolate",
          choices = c("select classes" = "",levels(as.factor(read.csv(input$selected_gene_table)[,2]))),
          selected = character(0),
          multiple = TRUE
        )
      } else {
        HTML("Please upload valid gene table before selecting classes ('full table' or 'gene class only' table)")
      }
    }
  })
  
  # Results Tables
  
  # Results Tables UI
  output$results_tables_ui = renderUI({
    req(input$run)
    if (input$results_tables_render) {
      fluidRow(selectInput("results_tables_contrast","Select contrast",colnames(rv$combinations),selected = 1),dataTableOutput("results_tables_render"))
    } else {
      HTML("To view results tables, check the corresponding box under <strong>Plot Selection</strong> and run the analysis again.")
    }
  })
  
  # Results Tables Rendering by Contrast
  output$results_tables_render = renderDataTable({
    req(input$run)
    return(datatable(rv$tables[[which(colnames(rv$combinations) == input$results_tables_contrast)]]))
  })
  
  # PCA Plotting
  
  # PCA Plot UI
  output$pca_ui = renderUI({
    req(input$run)
    if (input$pca_render) {
      withSpinner(plotOutput("pca_plot_render",height = "600px",width = "800px"),type = 1)
    } else {
      HTML("To view the pca plot, check the corresponding box under <strong>Plot Selection</strong> and run the analysis again.")
    }
  })
  
  # PCA Plot Rendering
  output$pca_plot_render = renderPlot({
    req(input$run)
    return(pca_input(rv$dds))
  })
  
  # Intra-Condition Scatter Plotting
  
  # Intra-Condition UI
  output$intra_condition_ui = renderUI({
    req(input$run)
    if (input$intra_condition_render) {
      withSpinner(plotOutput("intra_condition_plot_render",height = "800px",width = "800px"),type = 1)
    } else {
      HTML("To view the intra-condition scatter plot, check the corresponding box under <strong>Plot Selection</strong> and run the analysis again.")
    }
  })
  
  # Intra-Condition Rendering
  output$intra_condition_plot_render = renderPlot({
    req(input$run)
    return(intra_condition_input(rv$sample_table,counts(rv$dds,normalized = TRUE)))
  })
  
  # Mean Reads Scatter Plotting
  
  # Mean Reads Scatter Plots UI
  output$mean_reads_ui = renderUI({
    req(input$run)
    if (input$mean_reads_render) {
      fluidRow(
        selectInput(
          inputId = "mean_reads_contrast",
          label = "Select contrast",
          choices = colnames(rv$combinations),
          selected = 1
        ),
        
        withSpinner(plotOutput("mean_reads_render",height = "650px",width = "800px"),type = 1),
        
        downloadButton("download_mean_reads","Download Interactive Version")
      )
    } else {
      HTML("To view mean reads scatter plots, check the corresponding box under <strong>Plot Selection</strong> and run the analysis again.")
    }
  })
  
  # Plot Rendering by Contrast
  output$mean_reads_render = renderPlot({
    req(input$run)
    return(rv$mean_reads_plots[[which(colnames(rv$combinations) == input$mean_reads_contrast)]])
  })
  
  # Interactive Plot Download
  output$download_mean_reads = downloadHandler(
    filename = paste0(format(Sys.Date(),format = "%d_%m_%y_"),"Mean_Reads.html"),
    content = function(file) {
      htmlwidgets::saveWidget(
        ggplotly(rv$mean_reads_plots[[which(colnames(rv$combinations) == input$mean_reads_contrast)]],tooltip = "text"),
        file = file
      )
    }
  )
  
  # MA Plotting
  
  # Mean Reads Scatter Plots UI
  output$ma_ui = renderUI({
    req(input$run)
    if (input$ma_render) {
      fluidRow(
        selectInput(
          inputId = "ma_contrast",
          label = "Select contrast",
          choices = colnames(rv$combinations),
          selected = 1
        ),
        
        withSpinner(plotOutput("ma_render",height = "650px",width = "800px"),type = 1)
      )
    } else {
      HTML("To view ma plots, check the corresponding box under <strong>Plot Selection</strong> and run the analysis again.")
    }
  })
  
  # Plot Rendering by Contrast
  output$ma_render = renderPlot({
    req(input$run)
    
    # Plotting Call
    return(
      ma_input(
        rv$res[[which(colnames(rv$combinations) == input$ma_contrast)]],
        rv$combinations[1,which(colnames(rv$combinations) == input$ma_contrast)],
        rv$combinations[2,which(colnames(rv$combinations) == input$ma_contrast)]
      )
    )
  })
  
  # Heatmap Plotting
  
  # Heatmap UI
  output$heatmap_ui = renderUI({
    req(input$run)
    if (input$heatmap_render) {
      withSpinner(plotlyOutput("heatmap_plot",height = "800px",width = "1000px"),type = 1)
    } else {
      HTML("To view the heatmap, check the corresponding box under <strong>Plot Selection</strong> and run the analysis again.")
    }
  })
  
  # Heatmap Rendering
  output$heatmap_plot = renderPlotly({
    req(input$run)
    
    # Conditional Plotting Calls
    if (input$heatmap_type == "complete") {
      complete_heatmap(rv$cts_heat,rv$cts_names)
    } else if (input$heatmap_type == "all_classes") {
      all_classes_heatmap(rv$cts_heat,rv$cts_names,rv$class_names)
    } else if (input$heatmap_type == "selected_classes") {
      selected_classes_heatmap(rv$cts_heat,rv$cts_names,rv$class_names,input$heatmap_selected_classes)
    }
  })
}