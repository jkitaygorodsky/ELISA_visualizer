#ELISA Participant Visualization - Shiny app
#Julia Kitaygorodsky
#September 15, 2022

# ipak <- function(pkg){
#   new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#   if (length(new.pkg))
#     install.packages(new.pkg, dependencies = TRUE)
#   sapply(pkg, require, character.only = TRUE)
# }
# 
# packages <- c("rmarkdown", "rvest", "tidyverse", "stringr", "gplots", "ggplot2", "shiny", "shinymanager", "bslib", "rsconnect", "grid", "gtable", "gridExtra", "gridtext", "mdthemes")
# ipak(packages)



library("shiny")
library("shinymanager")
library("ggplot2")
library("tidyverse")
library('grid')
library('gtable')
library('gridExtra')
library('gridtext')
# library('mdthemes')

seroconversionThresholds <- c(SmT1 = 0.19, RBD = 0.186, NP = 0.396)
convalescentMedians <- c(SmT1 = 1.38, RBD = 1.25, NP = 1.13)

DBS_Thresholds <- c(SmT1 = 0.482, RBD = 0.324, NP = 0.642)

colours <- c("#F8766D" = "Your antibodies fall in the 'red zone'. This means we could not detect any antibodies to the virus. This probably means that you are susceptible to future infections with SARS-CoV-2/COVID-19; however, we and others are trying to determine if you might have other immune cells that can help protect you.",
             "#00BFC4" = "Your antibodies fall in the 'blue zone'. This means you did not make as many antibodies as a person who had COVID-19. We do not yet know how many antibodies you need to be protected from infections so we cannot advise you as to whether you are protected from future infections.",
             "#7CAE00" = "Your antibody levels fall in the 'green zone'. This means you made about the same level of antibodies as a person who had COVID-19.")

zones <- c("F8766D" = "red", "#00BFC4" = "blue", "#7CAE00" = "green")

get_participantID <- function(ID, subjects, externalID){
  if(ID == ""){
    return(ID)
  } else if(ID %in% unique(subjects)){
    return(ID)
  } else if(ID %in% externalID){
    if(length(grep("[\\(\\)]", ID) > 0)){
      ID_tmp <- gsub("\\(", "\\\\(", ID)
      ID_tmp <- gsub("\\)", "\\\\)", ID_tmp)
      return(subjects[grep(ID_tmp, externalID)[1]])
    } else{
      return(subjects[grep(ID, externalID)[1]])
    }
  }
  strsplit(ID, "-")[[1]] %>% .[1:(length(.)-1)] %>% paste0(., collapse = "-")
}


dilution_factors <- data.frame(amount = c(2.5, 0.625, 0.15625, 0.0625, 0.0156, 0.00390625), dilution = c(4, 16, 64, 160, 640, 2560), sourceType = c("DBS", "DBS", "DBS", "plasma_serum", "plasma_serum", "plasma_serum"))
dilution_factors$amount <- round(dilution_factors$amount, digits = 4)
BAUmLLinearRange <- list(upper = c(SmT1 = 1, RBD = 1, NP = 2), lower = c(SmT1 = 0.03125, RBD = 0.03125, NP = 0.0625))

BAUmL_conversion <- function(RR, amount, antigen){
  if(!(amount %in% dilution_factors$amount)){stop("The amount does not match a standard dilution.")}
  if(!(antigen %in% names(seroconversionThresholds))){stop("The antigen is not in the study.")}
  if((RR >= BAUmLLinearRange$upper[[antigen]]) | (RR <= BAUmLLinearRange$lower[[antigen]])){ ##cap values at upper limit or lower limit
    RR = min(max(RR, BAUmLLinearRange$lower[[antigen]]), BAUmLLinearRange$upper[[antigen]])
    # return("not_linear")
  }
  
  # else if(antigen == "NP"){
  if(antigen == "NP"){
    
    2^((log2(RR) - 0.243) / 0.713 + log2(dilution_factors$dilution[which(dilution_factors$amount == amount)]))
  } else if(antigen == "RBD"){
    2^((log2(RR) + 0.612) / 0.766 + log2(dilution_factors$dilution[which(dilution_factors$amount == amount)]))
  } else if(antigen == "SmT1"){
    2^((log2(RR) - 0.604) / 0.784 + log2(dilution_factors$dilution[which(dilution_factors$amount == amount)]))
  }
}

BAUmLThresholds <- c(BAUmL_conversion(seroconversionThresholds["SmT1"], 0.0625, "SmT1"),
                     BAUmL_conversion(seroconversionThresholds["RBD"], 0.0625, "RBD"),
                     BAUmL_conversion(seroconversionThresholds["NP"], 0.0625, "NP"))
BAUmLConvalescentMedians <- c(BAUmL_conversion(convalescentMedians["SmT1"], 0.0625, "SmT1"),
                              BAUmL_conversion(convalescentMedians["RBD"], 0.0625, "RBD"),
                              BAUmL_conversion(convalescentMedians["NP"], 0.0625, "NP"))
names(BAUmLConvalescentMedians) = names(convalescentMedians)


ui <- fluidPage(
  tags$head(
    HTML(
      "
          <script>
          var socket_timeout_interval
          var n = 0
          $(document).on('shiny:connected', function(event) {
          socket_timeout_interval = setInterval(function(){
          Shiny.onInputChange('count', n++)
          }, 15000)
          });
          $(document).on('shiny:disconnected', function(event) {
          clearInterval(socket_timeout_interval)
          });
          </script>
          "
    )
  ),
  textOutput("keepAlive"),
  titlePanel("Serology ELISA Data Visualization"),
  fluidRow(column(4,
                  fileInput("file1", "Choose CSV File",
                            multiple = FALSE,
                            accept = c("text/csv",
                                       "text/comma-separated-values,text/plain",
                                       ".csv")),
                  br(),
                  textInput("participantID", label = "Participant ID:"),
                  br(),
                  uiOutput("select_timepoint"),
                  br(),
                  selectInput("dataType", label = "Select dilution or BAU/mL:", choices = c("0.00390625", "0.0625", "BAU/ml"), selected = "0.0625"),
                  br(),
                  checkboxInput("log2transf", "Log2 Data Transformation", value = FALSE),
                  br(),
                  fluidRow(column(4, tableOutput('table')))),
           column(8, plotOutput("studyPlot"),
                  br(),
                  textOutput('result_description'))),
  br(),
  br(),
  fluidRow(column(1, downloadLink('export', 'Download')))
)



server <- function(input, output, session) {
  output$keepAlive <- renderText({
    req(input$count)
    paste("keep alive ", input$count)
  })
  ## load the appropriate data set (user-reactive)
  filtered_data <- reactive({
    req(input$file1)
    
    tryCatch(
      {
        df <- read.csv(input$file1$datapath,
                       header = TRUE,
                       sep = ",",
                       quote = '"')
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    # print(colnames(df))

    # print(as.character(input$dataType))
    columns_to_use <<- c()
    if(input$dataType == "BAU/ml"){
      thresholds <<- BAUmLThresholds
      medians <<- BAUmLConvalescentMedians
      for(antigen in c("RBD", "SmT1")){
        columns_to_use <<- c(columns_to_use, paste0(antigen, ".BAU.ml.aggregated"))
        
      }
    } else{
      thresholds <<- seroconversionThresholds
      medians <<- convalescentMedians
      for(antigen in c("RBD", "SmT1")){
        columns_to_use <<- c(columns_to_use, paste0(antigen, ".IgG.", as.character(input$dataType)))
      }
    }
    
    if(!("Timepoint" %in% colnames(df))){
      df$Timepoint = NA
    }
    if(!("Participant_Id" %in% colnames(df))){
      df$Participant_Id = df$External_Id
    }
    
    to_plot <<- df[, c("External_Id", "Participant_Id", columns_to_use, "Timepoint")]
    colnames(to_plot) <<- c("External_Id", "Participant", "RBD Level", "Spike Level", "Timepoint")
    
    if(input$dataType == "BAU/ml"){
      to_plot$`Spike Level` <<- lapply(to_plot$`Spike Level`, function(x){as.character(x) %>% strsplit("[ ><]") %>% .[[1]] %>% .[. != ""] %>% as.numeric()}) %>% unlist()
      to_plot$`RBD Level` <<- lapply(to_plot$`RBD Level`, function(x){as.character(x) %>% strsplit("[ ><]") %>% .[[1]] %>% .[. != ""] %>% as.numeric()}) %>% unlist()
    }
    
    if(input$dataType == "0.00390625"){
      to_plot$colour <<- "grey"
    } else{
      to_plot$colour <<- apply(to_plot, MARGIN = 1, function(x){
          if(any(as.numeric(x[c("Spike Level", "RBD Level")]) < thresholds[c("SmT1", "RBD")]))
                 {names(colours)[1]}
          else if((input$dataType == "0.0625") &&
            (all(as.numeric(x[c("Spike Level", "RBD Level")]) > thresholds[c("SmT1", "RBD")]) &
                  any(as.numeric(x[c("Spike Level", "RBD Level")]) < medians[c("SmT1", "RBD")])))
            {names(colours)[2]}
          else{
            names(colours)[3]}
        })
    }
    
    if(input$log2transf == TRUE){
      to_plot$`Spike Level` <<- log2(to_plot$`Spike Level`)
      to_plot$`RBD Level` <<- log2(to_plot$`RBD Level`)
    }
    
    to_plot$pointType <<- 19
    
    to_plot$size <<- 1.5
    
    # print(head(to_plot))
  
    if((length(input$Timepoint) > 0) && (input$Timepoint != "all")){
      dplyr::filter(to_plot, Timepoint == input$Timepoint)
    } else{
      to_plot
    }
})
      
      
      
  output$select_timepoint <- renderUI({
    wake_filtered_data <- filtered_data()[1, ]
    selectInput('Timepoint', 'Please Select Desired Time Point:', choices = c("all", as.character(sort(unique(to_plot$Timepoint)))), selected = input$Timepoint)
  })



  plotFunc <- reactive({
    # req(input$file1)
    prepared_plot <- ggplot(filtered_data(), aes(x=`RBD Level`, y=`Spike Level`, color=colour)) +
      geom_point(color=filtered_data()$colour, shape=filtered_data()$pointType, size=filtered_data()$size) +
      theme_classic()
    
    if(input$dataType != "0.00390625"){
      prepared_plot <- prepared_plot +
        geom_hline(yintercept = if(input$log2transf == TRUE){log2(thresholds["SmT1"])}else{thresholds["SmT1"]}, color = "grey") +
        geom_vline(xintercept = if(input$log2transf == TRUE){log2(thresholds["RBD"])}else{thresholds["RBD"]}, color = "grey")
    }
    if(input$dataType == "0.0625"){
      prepared_plot <- prepared_plot +
        geom_hline(yintercept = if(input$log2transf == TRUE){log2(medians["SmT1"])}else{medians["SmT1"]}, color = "grey", linetype = "dashed") +
        geom_vline(xintercept = if(input$log2transf == TRUE){log2(medians["RBD"])}else{medians["RBD"]}, color = "grey", linetype = "dashed")
    }

    if(input$participantID == ""){
      show_all = FALSE
    } else{
      #check if this is a subject or a specific timepoint ID
      if(get_participantID(input$participantID, to_plot$Participant, to_plot$External_Id) == input$participantID){show_all = TRUE} else{show_all = FALSE}
    }

    if(show_all == TRUE){
      prepared_plot +
        geom_point(data = filtered_data()[filtered_data()$Participant == get_participantID(input$participantID, to_plot$Participant, to_plot$External_Id), ],
                   color = "black", size = 5, shape=8, stroke = 3)
    } else{
      prepared_plot +
        geom_point(data = filtered_data()[filtered_data()$External_Id == input$participantID, ],
                   color = "black", size = 5, shape=8, stroke = 3)
    }
  })


  # output the filtered plot
   output$studyPlot <- renderPlot({
     # req(input$file1)
     plotFunc()
     })

   values <- reactive({
     # req(input$file1)
     wake_filtered_data <- filtered_data()[1, ]
     
       if(input$participantID == ""){
         participant = ""
       } else if(input$participantID %in% to_plot$External_Id){
         participant = to_plot$Participant[to_plot$External_Id == input$participantID]
       } else{
         participant = get_participantID(input$participantID, to_plot$Participant, to_plot$External_Id)
       }
      valuesTable <<- to_plot[to_plot$Participant == participant, c("External_Id", "Participant", "Timepoint", "RBD Level", "Spike Level", "colour")] %>% mutate(across(where(is.numeric), round, 2))

      valuesTable[, c("External_Id", "Timepoint", "RBD Level", "Spike Level")]
     # }
   })
   # show the values of the participant's Spike and RBD levels
   output$table <- renderTable(
     values(), width = "200px"
   )

   # explain the meaning of the red/blue/green categories
   output$result_description <- renderText({
     # req(input$file1)
     if(input$participantID != ""){
       if(input$participantID %in% filtered_data()$External_Id){
         result_colour <<- filtered_data()$colour[filtered_data()$External_Id == input$participantID]
       } else if((input$Timepoint != "all") && (input$participantID %in% filtered_data()$Participant)){
         result_colour <<- filtered_data() %>% filter(Participant == input$participantID) %>% filter(Timepoint == input$Timepoint) %>% select(colour)
       } else{
         result_colour <<- 0
       }
     } else{
       result_colour <<- 0
     }

     if(length(result_colour) < 1){
       result_colour <<- 0
     }

     description <<- if(unlist(result_colour) %in% names(colours)){colours[[unlist(result_colour)]]} else{""}

     description
     })
   
   output$export = downloadHandler(
     filename = function() {"Plot output.pdf"},
     content = function(file) {
       pdf(file, paper = "default")
       # grid.text(paste0("Serology Results"),  x=0.5, y=.9, gp=gpar(fontsize=18), check=TRUE)
       #
       # grid.newpage()
       
       sample_vp_1 <- viewport(x = 0, y = 1,
                               width = 1, height = 0.5,
                               just = c("left", "top"))
       pushViewport(sample_vp_1)
       
       grid.draw(ggplotGrob(plotFunc()))
       
       popViewport()
       
       
       if(input$participantID != ""){
         sample_vp_2 <- viewport(x = 0, y = 0.5,
                                 width = 1, height = 0.3,
                                 just = c("left", "top"))
         pushViewport(sample_vp_2)
         
         tableForPDF <<- valuesTable
         tableForPDF$Zone <<- unlist(lapply(tableForPDF$colour, function(x){zones[[x]]}))
         
         tableForPDF <<- tableForPDF[, c("Participant", "Timepoint", "RBD Level", "Spike Level", "Zone")]
         
         g_values <- tableGrob(tableForPDF, rows = NULL, theme=ttheme_default(base_size = 15))
         title <- textGrob("Antibody Levels:", gp = gpar(fontsize = 16))
         padding <- unit(0.5,"line")
         g_values <- gtable_add_rows(g_values, heights = grobHeight(title) + padding, pos = 0)
         g_values <- gtable_add_grob(g_values, list(title), t = 1, l = 1, r = ncol(g_values))
         g_values <- gtable_add_grob(g_values,
                                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                     t = 2, b = nrow(g_values), l = 1, r = ncol(g_values))
         g_values <- gtable_add_grob(g_values,
                                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                                     t = 1, l = 1, r = ncol(g_values))
         grid.draw(g_values)
         
         popViewport()
         
         
         sample_vp_3 <- viewport(x = 0, y = 0.25,
                                 width = 1, height = 0.5,
                                 just = c("left", "top"))
         pushViewport(sample_vp_3)
         
         
         g_descrip <- textbox_grob(
           description,
           x = unit(0.5, "npc"), y = unit(0.7, "npc"),
           gp = gpar(fontsize = 15),
           r = unit(5, "pt"),
           padding = unit(c(10, 10, 10, 10), "pt"),
           margin = unit(c(0, 10, 0, 10), "pt")
         )
         
         grid.draw(g_descrip)
         
         
         popViewport()
       }
       dev.off()
     })
}

shinyApp(ui = ui, server = server, options = list(height = 600))


