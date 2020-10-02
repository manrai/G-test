library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(ggplot2)
library(reshape2)
library(DT)
library(BiocManager)
options(repos = BiocManager::repositories())
library(RGGV)
library(dplyr)
library(DBI)
library(RSQLite)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install(version="3.10")


ui <- dashboardPage(skin = "red", dashboardHeader(title = h5(HTML("NHLBI BioData Catalyst Penetrance API"))), 
    dashboardSidebar(
    sidebarMenu(
    menuItem("Gene and variant selection", 
        tabName = "geneselection", startExpanded = TRUE, icon = icon("dna")), 
        selectizeInput("checkGroup", "Select gene(s) of interest", c(MYBPC3 = "MYBPC3", 
            MYH7 = "MYH7", ACTC1 = "ACTC1", MYL2 = "MYL2", MYL3 = "MYL3", 
            TNNI3 = "TNNI3", TNNT2 = "TNNT2", TPM1 = "TPM1"), selected = c("MYBPC3"), 
            multiple = TRUE), 
        checkboxGroupInput("checkbox", "Variant class", 
            c(Pathogenic = "Pathogenic", `Uncertain significance` = "Uncertain significance", 
            `Likely pathogenic` = "Likely pathogenic"), selected = c("Pathogenic", 
             "Uncertain significance")), 
    menuItem("Phenotype Criteria", 
            tabName = "patienthx", startExpanded = TRUE, icon = icon("user")), 
        selectizeInput("checkGroupCohorts", "Select cohort(s)", c(CARDIA = "CARDIA", 
            `Jackson Heart Study` = "Jackson Heart Study", `Heart and Vascular Health Study` = "Heart and Vascular Health Study"), 
            selected = c("CARDIA", "Jackson Heart Study", "Heart and Vascular Health Study"), 
            multiple = TRUE), 
         sliderInput("age", "Age", value = 35, min = 18, 
            max = 80), selectInput("gender", "Gender", choices = c("Male", 
            "Female"), multiple = FALSE), 
         selectInput("race", "Race/Ethnicity", 
            choices = c("African/African-American", "Latino/Admixed American", 
                "Ashkenazi Jewish", "East Asian", "Finnish", "Non-Finnish European", 
                "South Asian", "Other"), multiple = FALSE), 
         switchInput(inputId = "Id015", 
            label = "Hypertension", labelWidth = "80px", onLabel = "Yes", 
            offLabel = "No"), sliderInput("n", "Left ventricular wall thickness (mm)", 
            value = 10, min = 5, max = 20)), 
        submitButton("Submit", width = "230")), 
        
    dashboardBody(fluidRow(
        tags$style(HTML(".skin-blue .sidebar-menu > li.active > a,\n .skin-blue .sidebar-menu > li:hover > a {\n background-color: #000;\n border-left-color: #ff0000;}")), 
        tags$head(tags$style(HTML("/* body */.content-wrapper, .right-side {background-color: white;}"))), 
        tags$head(tags$style(HTML(".span12 {background-color: black;}"))), 
        tags$head(tags$style(HTML(".content-wrapper {float:top; margin-top:0px; padding-left:15px; padding-right:15px}"))), 
        tabsetPanel(
        
        tabPanel("Calculate Penetrance", br(), fluidRow(valueBoxOutput("value1"), valueBoxOutput("value2"), valueBoxOutput("value3")), downloadButton("downloadData", 
            "Download report"), br(), br(), br(), DT::dataTableOutput("mytable1"), 
            helpText("Table for top variants by maximum population allele frequency (AF) for the following ancestries: 
            African/African-American (Afr), Latino/Admixed American (Lat), Ashkenazi Jewish (AJ), East Asian (EA), 
            European Finnish (EF), European Non-Finnish (ENF), and South Asian (SA)."), plotOutput("plot")), 
        tabPanel(title = "Explore Penetrance"), 
        tabPanel(title = "About Method", 
                withMathJax(), helpText("Penetrance, the probability that an individual who carries a genetic variant is affected with the disease, 
                                        can be represented in a Bayesian framework by:"), uiOutput("formula4"), helpText("where P(G|D) represents allele 
                                        frequency in cases, P(D) the lifetime risk of the disease in the general population, and P(G)the allele 
                                        frequency in controls. We modeled each term with a beta-binomial distribution, which is the binomial 
                                        distribution such that the probability of success at each of n non-fixed trials is randomly drawn from a 
                                        beta distribution. For p with a Beta distribution and n trials, the probability density function, mean, and
                                                                                                                         variance can be modeled as below."), 
                h4(textOutput("diagTitle")), 
                uiOutput("formula1"),uiOutput("formula2"), uiOutput("formula3"), 
                helpText("Visualize how your distribution changes with your selection of alpha, beta, and cohort:"), 
                selectizeInput("checkGroupCohorts", "Select cohort(s)", 
                     c(CARDIA = "CARDIA", `Jackson Heart Study` = "Jackson Heart Study", 
                    `Heart and Vascular Health Study` = "Heart and Vascular Health Study"), 
                     selected = c("CARDIA", "Jackson Heart Study", "Heart and Vascular Health Study"), 
                     multiple = TRUE), 
               selectizeInput("checkGroupCohorts2", "Select phenotypes(s)", c(HCM = "HCM"), selected = c("HCM")), 
                textInput(inputId = "Input1", label = "Alpha:", value = "1.0"), 
                textInput(inputId = "Input1", label = "Beta:", value = "1.0"), 
                helpText("Calculated mean is: 0.5881096"), helpText("Calculated variance is:"), 
                plotOutput("distPlot")
                )
             )
          )
       )
    )
                
                
                
                
                