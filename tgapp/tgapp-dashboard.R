library(shiny)
library(shinydashboard)
library(igvShiny)
library(DT)

library(TrenaGene)
library(TrenaGenePlacentaData)
geneSymbol <- "PSG1"
dataDirectory <- system.file(package="TrenaGenePlacentaData", "extdata")

psg1 <- TrenaGene(geneSymbol, "hg38", expressionDataDirectory=dataDirectory)

expression.data.set.names <- getExpressionMatrixNames(psg1)

ui <- dashboardPage(
  dashboardHeader(title = geneSymbol),
  dashboardSidebar(
     sidebarMenu(id="sidebarTabMenu",
        menuItem("Dashboard", tabName = NULL, icon = icon("dashboard")),
        menuItem("Genomic Region", tabName = NULL, icon = icon("th")),
          menuSubItem("Sub-item 1", tabName = NULL),
          menuSubItem("Sub-item 2", tabName = NULL)
          )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "dashboard",
        fluidRow(
           box(igvShinyOutput('igv'), width=8, height=800),
           box(selectInput("chooseExpression", "Choose Expression Data Set",
                            expression.data.set.names),
               DTOutput("condensedGeneModelTable"),width=4)
           )
        ),

      tabItem(tabName = "genomicRegion",
        h2("Widgets tab content"),
        DTOutput("fullGeneModelTable"), width=10)
        )
    ) # dashboardBody

) # dashboardPage
#------------------------------------------------------------------------------------------------------------------------
server <- function(input, output)
{

  observeEvent(input$sidebarTabMenu, {
     newTab <- isolate(input$sidebarTabMenu)
     printf("sidebarTabMenu event: %s", newTab)
     if(newTab == "fullModels"){
        printf("display full table: %s", isolate(input$chooseModel))
        }
     })

  output$igv <- renderIgvShiny(
    igvShiny(options <- list(genomeName="hg38", initialLocus=geneSymbol), width=900, height=800))

  output$condensedGeneModelTable <- renderDT({
      tbl <- prepareModelTable(modelName, "fullTable")
      tbl <- as.data.frame(mtx.summary)
      selection="single"
      tbl$TF <- rownames(tbl)
      tbl <- tbl[, c("TF", colnames(tbl)[-1])]
      rownames(tbl) <- NULL
      return(tbl[, c("TF", "observed", "mean")])
      }, width="200px", selection="single", options=list(dom='t'))

  output$fullGeneModelTable <- renderDT({
      modelName <- isolate(input$chooseModel)
      tbl <- prepareModelTable(modelName, "fullTable")
      tbl <- as.data.frame(mtx.summary)
      selection="single"
      tbl$TF <- rownames(tbl)
      tbl <- tbl[, c("TF", colnames(tbl)[-1])]
      rownames(tbl) <- NULL
      return(tbl)# [, c("TF", "observed", "mean")])
      }, width="800px", selection="single", options=list(dom='t', scrollX=TRUE))


} # server
#------------------------------------------------------------------------------------------------------------------------
prepareModelTable <- function(modelName, displayMode)
{
   tbl <- tables[[modelName]];
   ## if(displayMode ==


} # prepareModelTable
#------------------------------------------------------------------------------------------------------------------------

app <- shinyApp(ui, server)

