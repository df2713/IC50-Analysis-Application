#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(data.table)
library(tidyverse)
library(shinyFiles)
library(shinythemes)
library(gridExtra)
library(shinydashboard)
library(dplyr)
library(drc)
library(DescTools)

file_upload <- fileInput("data",
                       label = "Please Upload the Data as a CSV File",
                       multiple = FALSE,
                       accept = c(".csv"),
                       buttonLabel = "Upload")

# ---------------------------------------------------------------------
# - FUNCTION THAT ANALYZES THE DATA - FUNCTION THAT ANALYZES THE DATA -  
# ---------------------------------------------------------------------

# list of plots to output in the Dose-Response Tab

analyze <- function(csv_file) {
  raw_data <- read.csv(csv_file)
  plots_list <- list()
  auc_list <- list()
  auc_labels <- list()
  parameter_list <- data.frame(Plot=character(),Hill_Slope=double(), Min = double(), Max = double(), IC50 = double(), AUC = double(), R1 = double(), R2 = double())
  parameter_names <- c("Plot", "Hill Slope", "Max", "Min", "IC50", "AUC", "IC50 over AUC", "AUC over IC50")
  colnames(parameter_list) <- parameter_names
  
  # ------------------------
  # - DATA SETUP/CLEANING  -
  # ------------------------
  
  # -------------------------------
  # ---- IF one compound plate ----
  # -------------------------------
  
  if (ncol(raw_data) > 15) {
    # Makes data frame so that we are just one long column of tables if there are two columns of tables
    data1 <- raw_data[,1:16]
    data2 <- raw_data[,18:33]
    # empty rows that maintain the format, plus adding 1 since the column labels get cut off
    # Feb 15th, 2024 - If the last place for data1 is a plate with two concetrations, then our spacer should be 3 not 4
    # we do this by checking for zeros in the column named 2 but the third actual column
    
    if (as.numeric(data1[nrow(data1) ,3]) == 0  & as.numeric(data1[nrow(data1)-1 ,3]) == 0) {
      spacer <- as.data.table(matrix(nrow = 3, ncol = 16))
    } else {
      spacer <- as.data.table(matrix(nrow = 4, ncol = 16))
    }
    
    extra <- data.frame(t(c("Letter", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "uM", "Comp", "Line")))
    colnames(extra) <- c("Letter", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "uM", "Comp", "Line")
    colnames(data1) <- c("Letter", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "uM", "Comp", "Line")
    colnames(data2) <- c("Letter", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "uM", "Comp", "Line")
    colnames(spacer) <- c("Letter", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "uM", "Comp", "Line")
    fulldata <- rbind(data1, spacer, extra, data2, spacer)
    # adding a spacer at the end to account for if the last plate is a two concentration plate
    
  } else {
    # IF there is just one column of tables...
    # then we shouldn't need to do anything really
    fulldata <- raw_data[,1:16]
  }
  
  # ------------------------
  # - ANALYSIS - FOR PRISM -
  # ------------------------
  
  # setting variable for first row of current table
  curr_row <- 1
  
  while (curr_row < nrow(fulldata)) {
    # there are two cases to consider, one where the plate is 3 of 1 substance and 3 of another substance
    # the other case is when there is just one compound per plate
    
    # -------------------------------
    # ---- IF one compound plate ----
    # -------------------------------
    
    # We check to see that the labels in the first three match the labels in the last three
    if (fulldata[curr_row + 1, 15] == fulldata[curr_row + 4,15]) {
      
      # first compute the control mean value for the plate
      control <- mean(c(as.numeric(fulldata[curr_row + 1, 3]), as.numeric(fulldata[curr_row + 2, 3]), as.numeric(fulldata[curr_row + 3, 3]), 
                        as.numeric(fulldata[curr_row + 4, 3]), as.numeric(fulldata[curr_row + 5, 3]), as.numeric(fulldata[curr_row + 6, 3])))
      
      # divide everything in the table by the control value
      for (i in curr_row : (curr_row + 7)) {
        for (j in 2:13) {
          fulldata[i,j] <- as.numeric(fulldata[i,j])/control
          
          # multiply the value by 100
          fulldata[i,j] <- as.numeric(fulldata[i,j]) * 100
        }
      }
      # clearing out the top and bottom columns in with NAs
      for (k in 2:13) {
        fulldata[curr_row, k] <- NA
        fulldata[curr_row + 7,k] <- NA
      } 
    } else {
      # ---------------------------------
      # ---- ELSE two compound plate ----
      # ---------------------------------
      # compute control for both compounds
      control_top <- mean(c(as.numeric(fulldata[curr_row + 1, 3]), as.numeric(fulldata[curr_row + 2, 3]), as.numeric(fulldata[curr_row + 3, 3])))
      control_bottom <- mean(c(as.numeric(fulldata[curr_row + 4, 3]), as.numeric(fulldata[curr_row + 5, 3]), as.numeric(fulldata[curr_row + 6, 3])))
      
      # start dividing the top rows by the control
      for (i in curr_row : (curr_row + 3)) {
        for (j in 2:13) {
          fulldata[i,j] <- as.numeric(fulldata[i,j])/control_top
          
          # multiply the value by 100
          fulldata[i,j] <- as.numeric(fulldata[i,j]) * 100
        }
      }
      
      # start dividing the bottom rows by the control 
      for (i in (curr_row + 4) : (curr_row + 7)) {
        for (j in 2:13) {
          fulldata[i,j] <- as.numeric(fulldata[i,j])/control_bottom
          
          # multiply the value by 100
          fulldata[i,j] <- as.numeric(fulldata[i,j]) * 100
        }
      }
      # clearing out the top and bottom columns in with NAs
      for (k in 2:13) {
        fulldata[curr_row, k] <- NA
        fulldata[curr_row + 7,k] <- NA
      } 
    }
    
    # increment curr_row accordingly
    curr_row <- curr_row + 14
  }
  
  # outputs data normlaized to 100
  write.csv(fulldata, "C:/Users/david/OneDrive/Desktop/ic50/SKOV3_normalized.csv", row.names=FALSE)
  
  # -----------------------------------------------------------------------------------------
  # - ANALYSIS - FOR IC50 - ANALYSIS - FOR IC50 - ANALYSIS - FOR IC50 - ANALYSIS - FOR IC50 -
  # -----------------------------------------------------------------------------------------
  
  # reset row counter
  curr_row <- 1
  num_plots <- 1
  
  while (curr_row < nrow(fulldata)) {
    
    # -------------------------------------------------------------------------------------------------------------
    # IF it is a 1 compound plate || 1 compound plate || 1 compound plate || 1 compound plate || 1 compound plate -
    # -------------------------------------------------------------------------------------------------------------
    
    if (fulldata[curr_row + 1, 15] == fulldata[curr_row + 4, 15]) {
      # isolate the data we need
      drc <- fulldata[(curr_row + 1):(curr_row+6),3:12]
      
      # appends dosages to the data frame
      dosages_val <- t(c(as.numeric(fulldata[curr_row + 8,3:12])))
      colnames(dosages_val) <- colnames(drc)
      drc <- rbind(dosages_val, drc)
      
      # Transpose Dataframe and add labels for sample
      drc_t <- data.frame(t(drc))
      names <- c("dose", "1", "2", "3", "4", "5", "6")
      colnames(drc_t) <- names
      
      # Put into long format
      drc_t %>%
        pivot_longer(
          col = -dose,
          names_to = "sample",
          values_to = "response"
        ) -> drc_tidy
      
      # ensures that data is in proper type for use in model
      resp_numeric <- sapply(drc_tidy[,3], as.numeric)
      drc_tidy[,3] <- resp_numeric
      
      # Computes averages and standard error for each dosage
      drc_tidy %>%
        group_by(dose) %>%
        summarize(
          sem = sd(response)/sqrt(n()),
          response = mean(response)
        ) -> drc_per_dose
      
      # ensures dosages are in numeric form
      drc_per_dose$dose <- as.numeric(drc_per_dose$dose)
      
      # Fits the model
      drm(
        data = drc_per_dose,
        formula = response ~ dose,
        fct = LL.4(names = c("Hill Slope", "Min", "Max", "IC50"))
      ) -> model
      
      # saves the paramters and gets the IC50
      coeff <- data.frame(model[28])
      IC50 <- coeff[4,1]
      
      # puts IC50 in dataframe
      fulldata[curr_row + 1, ncol(fulldata)] <- paste("IC50=", IC50)
      
      # -----------------------------------------------------------------------
      # - PLOT - PLOT - PLOT - PLOT - PLOT - PLOT - PLOT - PLOT - PLOT - PLOT -
      # -----------------------------------------------------------------------
      
      data.frame(
        dose = seq(from = 0, to = max(drc_per_dose$dose), by = 0.1)
      ) -> predicted_data
      predict(model, newdata=predicted_data) -> predicted_data$prediction
      
      # ensures that things are the correct data type to aboiv non-numeric argument errors 2/7/2024
      drc_tidy$dose <- as.numeric(drc_tidy$dose)
      
      # -------------------------------------------
      # - AUC - AUC - AUC - AUC - AUC - AUC - AUC -
      # -------------------------------------------
      auc <- AUC(predicted_data$dose, predicted_data$prediction, method = "trapezoid")
      auc_list <- append(auc_list, auc)
      auc_labels <- append(auc_labels, paste("Plot", num_plots, "AUC"))
      
      # plotting
      drc_tidy %>%
        ggplot(aes(x=dose, y=response)) +
        geom_point() +
        scale_x_log10() + 
        geom_line(
          data = predicted_data,
          aes(x = dose, y = prediction),
          size = 1, color = "#339999"
        ) +
        labs(title = paste("Plot", num_plots, fulldata[curr_row + 1, 15]), fulldata[curr_row + 1, 16]) -> plot
      
      # add to plot list
      plots_list <- append(plots_list, list(plot))
      
      # annotates fulldata for reference, both plot number and AUC
      fulldata[curr_row + 2, 16] <- paste("Plot", num_plots)
      fulldata[curr_row + 3, 16] <- paste("AUC", auc)
      
      # PARAMETERS
      parameters <- as.data.frame(t(c(paste("Plot", num_plots, fulldata[curr_row + 1, 15], fulldata[curr_row, 16]), 
                                      round(coeff[1,1],3), round(coeff[2,1],3), round(coeff[3,1],3), round(coeff[4,1],3), round(auc,3),
                                      round(coeff[4,1]/auc,6), round(auc/coeff[4,1],6)
                                      )))
      colnames(parameters) <- parameter_names
      parameter_list <- rbind(parameter_list, parameters)
      
      # increment num_plots
      num_plots <- num_plots + 1
      
    } else {
      # -------------------------------------------------------------------------------------------
      # ELSE if is a 2 compound plate || 2 compound plate || 2 compound plate || 2 compound plate -
      # -------------------------------------------------------------------------------------------  
      
      # split data into one for the top compound and one for the bottom compound
      drc_top <- fulldata[(curr_row + 1):(curr_row+3),3:12]
      drc_bottom <- fulldata[(curr_row + 4):(curr_row+6),3:12]
      
      # -------------------------------------------------
      # - ADDING PROPER DOSAGES - ADDING PROPER DOSAGES -
      # -------------------------------------------------
      # IF there are two concentrations on the plate -
      # ----------------------------------------------
      # This part is a bit tricky. First, I need to check if there are two concentrations. This is done by checking whether or not the 
      # first after the concentration is NA
      if (!is.na(fulldata[curr_row + 9, 3])) {
        if (as.numeric(fulldata[curr_row + 8,3]) == 0 & as.numeric(fulldata[curr_row + 9,3]) == 0)  {
          dosages_val_top <- t(c(as.numeric(fulldata[curr_row + 8,3:12])))
          dosages_val_bottom <- t(c(as.numeric(fulldata[curr_row + 9,3:12])))
          
          # ensures column names are consistent
          colnames(dosages_val_top) <- colnames(drc_top)
          colnames(dosages_val_bottom) <- colnames(drc_top)
          
          # creates data sets for each
          drc_top <- rbind(dosages_val_top, drc_top)
          drc_bottom <- rbind(dosages_val_bottom, drc_bottom)
        }
      } else {
        dosages_val <- t(c(as.numeric(fulldata[curr_row + 8,3:12])))
        colnames(dosages_val) <- colnames(drc_top)
        
        drc_top <- rbind(dosages_val, drc_top)
        drc_bottom <- rbind(dosages_val, drc_bottom)
      }
      
      # Transpose Dataframe and add labels for samples
      drc_top_t <- data.frame(t(drc_top))
      drc_bottom_t <- data.frame(t(drc_bottom))
      names <- c("dose", "1", "2", "3")
      colnames(drc_top_t) <- names
      colnames(drc_bottom_t) <- names
      
      # Put into long format
      drc_top_t %>%
        pivot_longer(
          col = -dose,
          names_to = "sample",
          values_to = "response"
        ) -> drc_top_tidy
      
      drc_bottom_t %>%
        pivot_longer(
          col = -dose,
          names_to = "sample",
          values_to = "response"
        ) -> drc_bottom_tidy
      
      # Summarizes data and prepares to graph/model
      drc_top_tidy %>%
        group_by(dose) %>%
        summarize(
          sem = sd(response)/sqrt(n()),
          response = mean(as.numeric(response))
        ) -> drc_top_per_dose
      
      drc_bottom_tidy %>%
        group_by(dose) %>%
        summarize(
          sem = sd(response)/sqrt(n()),
          response = mean(as.numeric(response))
        ) -> drc_bottom_per_dose
      
      # coerce dosages into numeric
      drc_bottom_per_dose$dose <- as.numeric(drc_bottom_per_dose$dose)
      drc_top_per_dose$dose <- as.numeric(drc_top_per_dose$dose)
      
      # generates models and coefficients
      drm(
        data = drc_top_per_dose,
        formula = response ~ dose,
        fct = LL.4(names = c("Hill Slope", "Min", "Max", "IC50"))
      ) -> model_top
      
      drm(
        data = drc_bottom_per_dose,
        formula = response ~ dose,
        fct = LL.4(names = c("Hill Slope", "Min", "Max", "IC50"))
      ) -> model_bottom
      
      # storing the calculated parameters
      coeff_top <- data.frame(model_top[28])
      IC50_top <- coeff_top[4,1]
      coeff_bottom <- data.frame(model_bottom[28])
      IC50_bottom <- coeff_bottom[4,1]
      
      # puts IC50 values in dataframe
      fulldata[curr_row + 1, ncol(fulldata)] <- paste("IC50_",fulldata[curr_row+1,15],"=", IC50_top)
      fulldata[curr_row + 4, ncol(fulldata)] <- paste("IC50_",fulldata[curr_row+4,15],"=", IC50_bottom)
      
      # -----------------------------------------------------------------------
      # - PLOT - PLOT - PLOT - PLOT - PLOT - PLOT - PLOT - PLOT - PLOT - PLOT -
      # -----------------------------------------------------------------------
      # need to create two plots for each of the compounds
      
      # generates the data to predict - make sure that we are using the correct model for each
      data.frame(
        dose = seq(from = 0, to = max(drc_top_per_dose$dose), by = 0.1)
      ) -> predicted_data_top
      predict(model_top, newdata=predicted_data_top) -> predicted_data_top$prediction
      
      data.frame(
        dose = seq(from = 0, to = max(drc_bottom_per_dose$dose), by = 0.1)
      ) -> predicted_data_bottom
      predict(model_bottom, newdata=predicted_data_bottom) -> predicted_data_bottom$prediction
      
      # turn tidy dataframe response/dpse values into numeric 2/7/2024
      drc_top_tidy$dose <- as.numeric(drc_top_tidy$dose)
      drc_bottom_tidy$dose <- as.numeric(drc_bottom_tidy$dose)
      drc_top_tidy$response <- as.numeric(drc_top_tidy$response)
      drc_bottom_tidy$response <- as.numeric(drc_bottom_tidy$response)
      
      # -------------------------------------------
      # - AUC - AUC - AUC - AUC - AUC - AUC - AUC -
      # -------------------------------------------
      auc_top <- AUC(predicted_data_top$dose, predicted_data_top$prediction, method = "trapezoid")
      auc_bottom <- AUC(predicted_data_bottom$dose, predicted_data_bottom$prediction, method = "trapezoid")
      auc_list <- append(auc_list, auc_top)
      auc_list <- append(auc_list, auc_bottom)
      auc_labels <- append(auc_labels, paste("Plot", num_plots, "AUC"))
      auc_labels <- append(auc_labels, paste("Plot", num_plots + 1, "AUC"))
      
      # Using ggplot here
      drc_top_tidy %>%
        ggplot(aes(x=dose, y=response)) +
        geom_point() +
        scale_x_log10() + 
        geom_line(
          data = predicted_data_top,
          aes(x = dose, y = prediction), size = 1, color = "#339999"
        ) +
        labs(title = paste("Plot", num_plots, fulldata[curr_row + 1, 15]), fulldata[curr_row, 16]) -> plot_top
      
      drc_bottom_tidy %>%
        ggplot(aes(x=dose, y=response)) +
        geom_point() +
        scale_x_log10() + 
        geom_line(
          data = predicted_data_bottom,
          aes(x = dose, y = prediction), size = 1, color = "#339999"
        ) +
        labs(title = paste("Plot", num_plots + 1, fulldata[curr_row + 4, 15]), fulldata[curr_row, 16]) -> plot_bottom
      
      # start appending to plots_list
      # add to plot list
      plots_list <- append(plots_list, list(plot_top))
      plots_list <- append(plots_list, list(plot_bottom))
      
      # annotates fulldata for reference, auc included
      fulldata[curr_row + 2, 16] <- paste("Plot", num_plots)
      fulldata[curr_row + 5, 16] <- paste("Plot", num_plots + 1)
      fulldata[curr_row + 3, 16] <- paste("AUC", auc_top)
      fulldata[curr_row + 6, 16] <- paste("AUC", auc_bottom)
      
      # PARAMETERS
      parameters_top <- as.data.frame(t(c(paste("Plot", num_plots, fulldata[curr_row + 1, 15], fulldata[curr_row + 1, 16]), 
                                          round(coeff_top[1,1],3), round(coeff_top[2,1],3), round(coeff_top[3,1],3), round(coeff_top[4,1],3), round(auc_top,3),
                                          round(coeff_top[4,1]/auc_top,6), round(auc_top/coeff_top[4,1],6)
                                        )))
      colnames(parameters_top) <- parameter_names
      parameters_bottom <- as.data.frame(t(c(paste("Plot", num_plots, fulldata[curr_row + 4, 15], fulldata[curr_row + 1, 16]), 
                                          round(coeff_bottom[1,1],3), round(coeff_bottom[2,1],3), round(coeff_bottom[3,1],3), round(coeff_bottom[4,1],3), round(auc_bottom,3),
                                          round(coeff_bottom[4,1]/auc_bottom,6), round(auc_bottom/coeff_bottom[4,1],6)
                                          )))
      colnames(parameters_bottom) <- parameter_names
      parameter_list <- rbind(parameter_list, parameters_top, parameters_bottom)
      
      # since there are two plates, we must increment by 2
      num_plots <- num_plots + 2
    }
    # increment curr_row accordingly
    curr_row <- curr_row + 14
    
  }
  write.csv(parameter_list, "C:/Users/david/Desktop/ic50/APP_PARAMS.csv", row.names=FALSE)
  write.csv(fulldata, "C:/Users/david/Desktop/ic50/APP_OUTPUT.csv", row.names=FALSE)
  AUC_FINAL <- rbind(auc_labels, auc_list)
  AUC_FINAL <- t(as.data.frame(AUC_FINAL))
  
  # returns data frame
  return(list(fulldata, plots_list, AUC_FINAL, parameter_list))
}

ui <- fluidPage(theme = shinytheme("darkly"),

    # Application title
    titlePanel("David Fang's IC50 Analysis App"),
    
    tabsetPanel(
      tabPanel(
        "Upload and Instructions",
        sidebarLayout(
          sidebarPanel(
            file_upload
          ),
          mainPanel(
            # INSTRUCTIONS ON FORMATTING EACH INDIVIDUAL TABLE
            h2("Single Table Formatting"),
            p("Each table must have all wells accounted for. This means the table should span 1-12 across, and A-H down. In addition, there also be an extra row below
            the last row H for the dosages."),
            p("Finally, we should allocate three extra columns for information on the wavelength, compound, and cell line."),
            p(strong("IMPORTANT:"), "cell line information strong ", strong("must"), "be aligned with Row A for each table. Else, the program will cut off the first row's cell line."),
            p(strong("Example:"), "below is a single table in the correct format."), 
            img(src = "correcttable.png", width = "75%", style="display: block; margin-left: auto; margin-right: auto;"),
            p("Notice how the cell line, WT, is aligned with Row A"),
            
            # INSTRUCTIONS ON PROPER DATA FORMAT FOR MULTIPLE TABLES
            h2("Multiple Tables Formatting"),
            p("When uploading data ensure your spreadsheet or CSV is the following format:"),
            p("Tables for each plate should be separated by 4 rows, for example..."),
            img(src = "fourrows.png", width = "50%", style="display: block; margin-left: auto; margin-right: auto;"),
            p("If there are two columns worth of tables, tables for each plate should be separated by 4 columns, including any notes or annotations, for example..."),
            img(src = "fourcolumns.png", width = "50%", style="display: block; margin-left: auto; margin-right: auto;"),
            p("Lastly, ensure that there are no extra rows or gaps at the top of the spreadsheet. The first row should be the labeled column numbers as shown below:"),
            img(src = "firstrow.png", width = "50%", style="display: block; margin-left: auto; margin-right: auto;"),
            
            # FOR WHEN THERE ARE TWO COMPOUNDS ON A SINGLE PLATE
            h2("Plates with two compounds with different concentrations"),
            p("If we are running a plate with two compounds, but the concentrations we use for each compounds differ, please list the concentrations
              directly below the table, with the concentrations for compound run on the top half of the plate directly above the conctrations for the bottom
              half of the plate."),
            p("For example, in the plate below, we run the compounds ATO and AP7 on the same plate, but the concentrations used are
              different for both. SInce we ran ATO on top and AP7 on bottom, we have the concentrations listed in that order as well.
              "),
            img(src = "twoconcentrations.png", width = "75%", style="display: block; margin-left: auto; margin-right: auto;"),
            # properly formatted spreadsheet
            h2("Properly Formatted Sheet"),
            img(src = "propertable.png", width = "90%", style="display: block; margin-left: auto; margin-right: auto;")
          )
        )
      ),
      tabPanel(
        "Data",
        tableOutput("table")
      )
      ,
      tabPanel(
        "Dose-Response Plots",
        plotOutput("plot", height = "1080px")
      )
      ,
      tabPanel(
        "Dose-Response Parameters",
        tableOutput("param")
      )
      ,
      tabPanel(
        "AUC",
        tableOutput("area")
      )
    )
)

server <- function(input, output) {
  
  # file input
  finished_data <- reactive({
    req(input$data)
    file1 <- input$data
    file1$datapath
  })
  
  # fulldata output
  output$table <- renderTable({
    req(input$data)
    if (is.null(finished_data())) {
      return ()
      }
    analyze(finished_data())[1]
  })
  
  # AUC output
  output$area <- renderTable({
    req(input$data)
    if (is.null(finished_data())) {
      return ()
    }
    analyze(finished_data())[3]
  })
  
  # Parameter output
  output$param <- renderTable({
    req(input$data)
    if (is.null(finished_data())) {
      return ()
    }
    analyze(finished_data())[4]
  })
  
  # renderPlot
  output$plot <- renderPlot({
    if (is.null(finished_data())) {
      return ()
    }
    n <- length(analyze(finished_data())[[2]])
    print(grid.arrange(grobs = analyze(finished_data())[[2]], ncol = floor(sqrt(n))))})
}

# Run the application 
shinyApp(ui = ui, server = server)
