library(shiny)
library(shinydashboard)
library(shinyjs)
library(later)
library(igvShiny)
library(DT)
#------------------------------------------------------------------------------------------------------------------------
source("serverFunctions.R")
source("addTracks.R")
source("buildModel.R")
#------------------------------------------------------------------------------------------------------------------------
library(trenaSGM)
library(TrenaGene)
library(TrenaGeneDataPlacenta)
#------------------------------------------------------------------------------------------------------------------------
state <- new.env(parent=emptyenv())
state$models <- list()
#------------------------------------------------------------------------------------------------------------------------
targetGene <- list(hugo="PSG1", ensembl="ENSG00000231924", combined="ENSG00000231924|PSG1")
genomeName <- "hg38"
dataDirectory <- system.file(package="TrenaGenePlacentaData", "extdata")
psg1.data <- TrenaGeneDataPlacenta("PSG1")

trenaGene <- TrenaGene(psg1.data)
expressionMatrixNames <- getExpressionMatrixNames(trenaGene)
tbl.enhancers <- getEnhancers(trenaGene)
tbl.dhs <- getEncodeDHS(trenaGene)
footprintDatabases <- getFootprintDatabaseNames(getGeneData(trenaGene))
tbl.transcripts <- getTranscriptsTable(trenaGene)

#------------------------------------------------------------------------------------------------------------------------
# calculate the initial igv region, using the (typically? always? single transcripted designated as "gene"
#------------------------------------------------------------------------------------------------------------------------
tbl.gene <- subset(tbl.transcripts, moleculetype=="gene")[1,]
start <- tbl.gene$start
end   <- tbl.gene$endpos
chrom <- tbl.gene$chr
padding <- round(0.25 * (end-start))
initial.chromLoc <- sprintf("%s:%d-%d", chrom, start-padding, end+padding)
#------------------------------------------------------------------------------------------------------------------------
createSidebar <- function()
{
  dashboardSidebar(
    sidebarMenu(id="sidebarMenu",
       menuItem("IGV and Current Models", tabName = "igvAndTable"),
       menuItem("Build a new Model",      tabName = "buildModels")
      ),
    conditionalPanel(id="igvTableWidgets",
        condition = "input.sidebarMenu == 'igvAndTable'",
        actionButton(inputId = "igvHideButton", label = "Toggle IGV"),
        actionButton(inputId = "tableHideButton", label = "Toggle table"),
        selectInput("addTrack", "Add Track:",
                    c("",
                     "Enhancers"="enhancers",
                     "DHS open chromatin (encode, clustered)"="dhs"
                     )),
        selectInput("displayGenomicRegion", "Display Genomic Region:",
                    c("",
                      "Tradtional Promoter" = "traditionalPromoter",
                      "Enhancers Region" = "enhancersRegion"))
        ) # conditionalPanel
     ) # dashboardSidebar

} # createSidebar
#------------------------------------------------------------------------------------------------------------------------
createBody <- function()
{
   body <- #fluidRow(
      tabItems(
         tabItem(tabName="igvAndTable",
            fluidRow(
               column(width=3, offset=9, id="modelSelectorColumn",
                      selectInput("modelSelector", NULL,  "- no models yet -"))),
            fluidRow(
               column(width=9, id="igvColumn", igvShinyOutput('igvShiny', height="800px")),
               column(width=3, id="dataTableColumn", DTOutput("table"), type=1)
               ) # fluidRow
            ), # tabItem 1
         createBuildModelTab()
         ) # tabItems

   return(body)

} # createBody
#------------------------------------------------------------------------------------------------------------------------
ui <- dashboardPage(

  dashboardHeader(title=sprintf("trena %s", targetGene$hugo)),
  createSidebar(),
  dashboardBody(
      tags$head(tags$style(HTML('
        .main-header .logo {
          font-family: "Georgia", Times, "Times New Roman", serif;
          font-weight: bold;
          font-size: 24px;
          }
        #igvShiny{
          background: white;
          border: 1px solid black;
          border-radius: 5px;
          margin: 5px;
          margin-right: 15px;
          }
        #table{
          border: 1px solid black;
          border-radius: 5px;
          }
       '))),
     createBody()),
  useShinyjs()
  ) # dashboardPage

#------------------------------------------------------------------------------------------------------------------------
server <- function(session, input, output){

   output$igvShiny <- renderIgvShiny({
      options <- list(genomeName="hg38",
                      initialLocus=initial.chromLoc,
                      displayMode="EXPANDED",
                      trackHeight=300)
      igvShiny(options) # , height=800)
      })

   output$table = DT::renderDataTable(data.frame(),
                                      width="800px",
                                      class='nowrap display',
                                      options=list(scrollX=TRUE,
                                                   dom='t',
                                                   autowWdth=FALSE,
                                                   pageLength=100))

   observeEvent(input$currentGenomicRegion, {
      new.region <- isolate(input$currentGenomicRegion)
      #printf("new region: %s", new.region)
      state[["chromLocRegion"]] <- new.region
      })

   setupIgvAndTableToggling(session, input);
   setupAddTrack(session, input, output)
   setupBuildModel(session, input, output)
   later(function(){
            state$chromLocRegion <- initial.chromLoc
            showGenomicRegion(session, initial.chromLoc)
            }, 2)

} # server
#------------------------------------------------------------------------------------------------------------------------
app <- shinyApp(ui, server)
