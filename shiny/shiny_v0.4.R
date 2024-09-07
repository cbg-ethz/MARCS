

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyWidgets)
  library(parallel)
  library(BiDAG)
  # I had trouble installing this package with R 4.4.0
  # don't use R package manager, use Bioconductor.
  library(clustNet)
  # I also had trouble installing this one with R 4.4.0
  # only worked with Bioconductor using force=TRUE

  # Define UI with shinydashboard
  library(dplyr)
  library(htmlTable)
})

cluster_results <- readRDS("euler_memberships_8k_9clusters.rds")
source("get_classification.R")


ui <- dashboardPage(
  dashboardHeader(title = "MARCS: MDS/AML Aggregative Risk Classification System", titleWidth = 1000),
  dashboardSidebar(
    sidebarMenu(
      menuItem("FOR RESEARCH USE ONLY", tabName = "paper", icon = icon("file-text")),
      tags$li(class = "treeview",
              tags$ul(class = "treeview-menu menu-open", style = "display: block;",
              ),
      tags$li(style = "margin-left: 20px; margin-top: 10px; white-space: normal;",
                              HTML("Roncador, M. <em> et al. </em> (2024). Beyond AML and MDS separation: dissecting the complexity of MDS-AML continuum with the MDS-AML Aggregative Risk Classification System (MARCS) in 7480 patients. <em> Submitted. </em>, 45(2), 123-134.")

              )
      ),
      uiOutput("patientRecap")
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "paper",
              uiOutput("dynamicUI"),
              uiOutput("images")
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  results <- reactiveValues(data = NULL, myindex = NULL, display = FALSE, validata=list(), clinical_data=data.frame())
  displayState <- reactiveVal(TRUE)
  myindex <- reactiveVal(NULL) # Store myindex as a reactive value

  # Define a list of genes and aberrations
  genes <- c("ASXL1", "BCOR", "CBL", "CUX1", "DNMT3A", "ETV6", "EZH2", "GATA2",
             "GNAS", "IDH1", "JAK2", "KIT", "KRAS", "MPL", "NF1", "NPM1", "PHF6",
             "PTPN11", "RAD21", "RUNX1", "SF3B1", "SRSF2", "STAG2", "TET2", "TP53",
             "U2AF1", "WT1", "ZRSR2", "IDH2", "NRAS", "CEBPA", "FLT3")

  aberrations <- c("+8", "complex", "-17", "-5", "-7", "inv(16)", "inv(3)", "t(15;17)",
                   "t(6;9)", "t(8;21)", "t(9;11)", "t(9;22)", "t(v;11)")

  resultsData <- reactiveValues(myindex = NULL, error = FALSE)

  observeEvent(input$submit, {

    errorMessages <- ""

    # Check each condition and append an error message if the condition is met
    if (input$age <= 0 || is.na(input$age) ) {
      errorMessages <- paste(errorMessages, "Age cannot be negative or empty.", sep = "\n")
    }
    if (is.na(input$HB) || input$HB <= 0) {
      errorMessages <- paste(errorMessages, "Hemoglobin must be > 0 and cannot be empty.", sep = "\n")
    }
    if (input$HB > 20) {
      errorMessages <- paste(errorMessages, "Enter Hemoglobin concentration in g/dL, not g/L", sep = "\n")
    }
    if (input$PLT <= 0 || is.na(input$PLT) ) {
      errorMessages <- paste(errorMessages, "Platelets must be > 0 and cannot be empty.", sep = "\n")
    }
    if (input$PLT > 2000) {
      errorMessages <- paste(errorMessages, "Implausibly high PLT value.", sep = "\n")
    }
    if (input$WBC <= 0 || is.na(input$WBC) ) {
      errorMessages <- paste(errorMessages, "WBC must be > 0 and cannot be empty.", sep = "\n")
    }
    if (input$WBC > 200) {
      errorMessages <- paste(errorMessages, "Implausibly high WBC value:", sep="\n")
    }

    # Check if there were any error messages accumulated
    if (nchar(errorMessages) > 0) {
      # Show all error messages in a single notification
      showNotification(errorMessages, type = "error", duration = NULL) # duration = NULL keeps the notification visible until the user closes it
      return() # Prevent further execution of this observer
    }
      # Initialize validata with NA for all specified columns
      validata_temp <- data.frame(ASXL1 = NA, BCOR = NA, CBL = NA, CUX1 = NA, DNMT3A = NA, ETV6 = NA, EZH2 = NA, GATA2 = NA,
                                  GNAS = NA, IDH1 = NA, JAK2 = NA, KIT = NA, KRAS = NA, MPL = NA, NF1 = NA, NPM1 = NA,
                                  PHF6 = NA, PTPN11 = NA, RAD21 = NA, RUNX1 = NA, SF3B1 = NA, SRSF2 = NA, STAG2 = NA,
                                  TET2 = NA, TP53 = NA, U2AF1 = NA, WT1 = NA, ZRSR2 = NA, IDH2 = NA, NRAS = NA, CEBPA = NA,
                                  FLT3 = NA, BM_BLASTS = NA, HB = NA, PLT = NA, WBC = NA, AGE = NA, SEX = NA, stringsAsFactors = FALSE)



      # Replace NAs with entered values for basic inputs
      validata_temp$AGE <- as.numeric(input$age)
      validata_temp$SEX <- input$sex
      validata_temp$BM_BLASTS <- input$BM_BLASTS
      validata_temp$HB <- as.numeric(input$HB)
      validata_temp$PLT <- as.numeric(input$PLT)
      validata_temp$WBC <- as.numeric(input$WBC)

      resultsData$clinical_data = data.frame(
        Variable = c("Age","Sex","BM Blasts","HB","PLT","WBC"),
        Value = c(
          input$age |> as.numeric(),
          ifelse(input$sex==1,"male","female"),
          input$BM_BLASTS,
          input$HB,
          input$PLT,
          input$WBC
        )
      )



      # Dynamically update gene mutations and chromosomal aberrations
      gene_mutations <- c("ASXL1", "BCOR", "CBL", "CUX1", "DNMT3A", "ETV6", "EZH2", "GATA2",
                          "GNAS", "IDH1", "JAK2", "KIT", "KRAS", "MPL", "NF1", "NPM1", "PHF6",
                          "PTPN11", "RAD21", "RUNX1", "SF3B1", "SRSF2", "STAG2", "TET2", "TP53",
                          "U2AF1", "WT1", "ZRSR2", "IDH2", "NRAS", "CEBPA", "FLT3")

      chrom_aberrations <- c("+8", "complex", "-17", "-5", "-7", "inv(16)", "inv(3)", "t(15;17)",
                             "t(6;9)", "t(8;21)", "t(9;11)", "t(9;22)", "t(v;11)")

      # Update gene mutations
      for(gene in gene_mutations) {
        input_id <- paste0("gene_", gene)
        validata_temp[[gene]] <- as.numeric(input[[input_id]])
      }

      # Update chromosomal aberrations, adjusting for input IDs
      for(aberration in chrom_aberrations) {
        input_id <- paste0("aberration_", gsub("[.+]", "", gsub(" ", "_", aberration)))
        validata_temp[[gsub("[.+]", "", gsub(" ", "_", aberration))]] <- as.numeric(input[[input_id]])
      }
      colnames(validata_temp)[39:51] <- c("X.8", "complex", "X.17", "X.5", "X.7", "inv.16.", "inv.3.", "t.15.17.", "t.6.9.", "t.8.21.", "t.9.11.", "t.9.22.", "t.v.11.")

      if (!is.null(validata_temp)) {

        validata <- validata_temp
        validata$SEX <- as.numeric(validata$SEX)
        validata$AGE <- as.numeric(ifelse(validata$AGE >65, 1, 0)) # median from read.table("../data/aml_mds_matrix_8k.csv")
        validata$HB <- as.numeric(ifelse(validata$HB > 9.5, 1, 0))
        validata$PLT <- as.numeric(ifelse(validata$PLT > 85, 1, 0))
        validata$WBC <- as.numeric(ifelse(validata$WBC > 5.9 , 1, 0))
        validata$BM_BLASTS <- ifelse(validata$BM_BLASTS >9, ifelse(validata$BM_BLASTS >19, 2, 1), 0)

        rownames(validata) <- "ShinyPatient"


        validata <- validata %>%
          select(ASXL1, BCOR, CBL, CUX1, DNMT3A, ETV6, EZH2, GATA2, GNAS, IDH1, JAK2, KIT, KRAS, MPL, NF1, NPM1, PHF6, PTPN11, RAD21, RUNX1, SF3B1, SRSF2, STAG2, TET2, TP53, U2AF1, WT1, ZRSR2, IDH2, NRAS, CEBPA, FLT3, X.8, complex, X.17, X.5, X.7, inv.16., inv.3., t.15.17., t.6.9., t.8.21., t.9.11., t.9.22., t.v.11., BM_BLASTS, HB, PLT, WBC, AGE, SEX)
        resultsData$validata <- validata

        validata <- validata %>%
          select_if(~ !any(is.na(.)))



        mds_reclass <- get_classification(cluster_results, data_classify =  as.data.frame(validata))
        print(mds_reclass$clustermembership)

        lab <- c("UHR","HR1","NPM1","HR2","HR3","INT1","HR3","INT2","LR1","LR2")
        myindex <- mds_reclass$clustermembership
        resultMessage <- lab[mds_reclass$clustermembership]

        resultsData$myclass <- lab[mds_reclass$clustermembership]
        resultsData$myindex <- myindex
        results$display <- TRUE

        output$showResults <- reactive({ results$display })
        outputOptions(output, "showResults", suspendWhenHidden = FALSE)

        output$analysisResults <- renderText({
          if (is.null(results$myindex) || identical(results$myindex, integer(0))){
            "Please enter more data."
          } else {paste("Cluster: ", resultMessage)}
        })


      }

    displayState(FALSE) # Switch state to show results or error
  })

  observeEvent(input$new_entry, {
    displayState(TRUE) # Reset to the data entry state
  })

  output$dynamicUI <- renderUI({
    if (displayState()) {
      # Data entry UI with the specific structure provided
      tabItem(tabName = "data_entry",
              fluidRow(
                box(title = "Patient Info", status = "primary", solidHeader = TRUE,
                    numericInput("age", "Age", value = NA, min = 1),
                    radioButtons("sex", "Sex", choices = c("Male" = 1, "Female" = 0), inline=T)
                ),
                box(title = "Clinical Information", status = "warning", solidHeader = TRUE, collapsible = TRUE,
                    sliderInput("BM_BLASTS", "BM_BLASTS", min = 0, max = 100, value = 10),
                    numericInput("HB", "HB [g/dL]", value = 0 , min = 0.01, max=20),
                    numericInput("PLT", HTML(paste0("PLT [10",tags$sup("9"),"/L]")), value = 0, min = 0, max=2000),
                    numericInput("WBC", HTML(paste0("WBC [10",tags$sup("9"),"/L]")), value = 0, min = 0.01,  max=200)
                ),
                box(title = "Gene Mutations", status = "info", solidHeader = TRUE, width = 12,
                    div(style = "height:200px; overflow-y: scroll;",
                        # Generate gene mutation inputs dynamically
                        lapply(c("ASXL1", "BCOR", "CBL", "CUX1", "DNMT3A", "ETV6", "EZH2", "GATA2",
                                 "GNAS", "IDH1", "JAK2", "KIT", "KRAS", "MPL", "NF1", "NPM1", "PHF6",
                                 "PTPN11", "RAD21", "RUNX1", "SF3B1", "SRSF2", "STAG2", "TET2", "TP53",
                                 "U2AF1", "WT1", "ZRSR2", "IDH2", "NRAS", "CEBPA", "FLT3"), function(gene) {
                                   radioButtons(paste0("gene_", gene), label = tags$em(gene), choices = c("NA" = "NA", "Mutated" = 1, "Unmutated" = 0), selected = "NA", inline=T)
                                 })
                    )
                ),
                box(title = "Chromosomal Aberrations", status = "success", solidHeader = TRUE, collapsible = TRUE, width = 12,
                    div(style = "height:200px; overflow-y: scroll;",
                        # Generate chromosomal aberration inputs dynamically
                        lapply(c("+8", "complex", "-17", "-5", "-7", "inv(16)", "inv(3)", "t(15;17)",
                                 "t(6;9)", "t(8;21)", "t(9;11)", "t(9;22)", "t(v;11)"), function(aberration) {
                                   radioButtons(paste0("aberration_", gsub("[+.]", "", aberration)), label = aberration, choices = c("NA" = "NA", "Present" = 1, "Absent" = 0), selected = "NA", inline=T)
                                 })
                    )
                ),
                actionButton("submit", "Submit", class = "btn-primary")
              )
      )
    } else {
      # Results or error message UI
      if (resultsData$error) {
        list(
          fluidRow(
            box(title = "Error", status = "danger", solidHeader = TRUE,
                "Please enter more data.",
                actionButton("new_entry", "New Entry", class = "btn-primary")
            )
          )
        )
      } else {
        list(
          fluidRow(
            box(title = "Results", status = "info", solidHeader = TRUE,
                paste("Cluster: ", resultsData$myclass),
                tags$p(img(src = paste0("risk_icons/",resultsData$myclass,".png"),style = "width:100%; height: auto; display: block;")),
                actionButton("new_entry", "New Entry", class = "btn-primary")
            )
          )
        )
      }
    }
  })



  output$images <- renderUI({
    # Ensure myindex is not NULL or an empty value before trying to display images
      tagList(
        # Display images one under the other, adjusting the height to fit the screen if necessary
        # Assuming images are stored in the www directory, adjust paths if different
        img(src = paste0("kaplan_maier/KM_",resultsData$myclass,".png"), style = "width:100%; height: auto; display: block;"), #, style = "width:50%; height:auto; display: block;")
      )

  })


  # Dynamic sidebar for patient details recap
  output$patientRecap <- renderUI({
    if (!displayState()) {
      geneRecap <- lapply(genes, function(gene) {
        paste0(gene, ":", match(input[[paste0("gene_",gene)]], c("NA" = "NA", "0" = "Unmutated", "1" = "Mutated")))
      })
      aberrationsRecap <- paste("Chromosomal Aberrations:", ifelse(length(input$aberrations) > 0, paste(input$aberrations, collapse=", "), "None"))
      #print(geneRecap)
      patientData <- t(resultsData$validata) |> as_tibble() |> rename(Value="ShinyPatient")
      patientData$Variable <- colnames(resultsData$validata)
      patientData <- patientData |> select(Variable, Value)
      patientData$Variable[33:45] <- aberrations
      patientData <- patientData |> mutate(Value=ifelse(is.na(Value),"n/a",ifelse(Value==1,"yes","no")))
      return(tags$div(list(
        tags$h4("Submitted patient data"),
        tags$h5("Clinical data"),
        htmlTable(resultsData$clinical_data,rnames=F),
        tags$h5("Mutations"),
        htmlTable(patientData[1:32,],rnames=F),
        tags$h5("Cytogenetics"),
        htmlTable(patientData[33:45,],rnames=F)), class="container"))
      list(
        tags$h4("Patient Details Recap:"),
        tags$p(paste("Age:", input$age)),
        tags$p(paste("Sex:", input$sex)),
        tags$p(paste("BM Blasts:", input$bm_blasts)),
        tags$p(paste("HB:", input$hb)),
        tags$p(paste("PLT:", input$plt)),
        tags$p(paste("WBC:", input$wbc)),
        tags$h5("Gene Mutations Recap:"),
        tags$ul(lapply(geneRecap, tags$li)),
        tags$h5(aberrationsRecap)
      )
    }
  })
}

shinyApp(ui, server)
