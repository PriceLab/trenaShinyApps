library(shiny)
library(shinyjs)
library(htmlwidgets)
library(igvShiny)

targetGene <- list(hugo="PSG1", ensembl="ENSG00000231924", combined="ENSG00000231924|PSG1")
footprintDatabases <- c("db one", "db two")

createUI <- function()
{
  ui <- fluidPage(
     useShinyjs(),
     titlePanel(sprintf("%s Trena Model & Disruptions", targetGene$hugo)),

     sidebarLayout(
        sidebarPanel(
           width=2,
           selectInput("displayGenomicRegion", "Display Genomic Region:",
                       c("",
                         "Tradtional Promoter" = "traditionalPromoter",
                         "Enhancers Region" = "enhancersRegion")),
           selectInput("addTrack", "Add Track:",
                      c("",
                        "Enhancers"="enhancers",
                        "DHS open chromatin - encode clustered"="dhs"
                        )),
           HTML("<br>"),
           actionButton("removeOptionalTracks", "Remove tracks")
           ),  # sidebarPanel

        mainPanel(width=10,
            tabsetPanel(type="tabs",
                        id="trenaTabs",
                        igvTabPanel(),
                        createGeneModelPanel()
                        ) # tabsetPanel
             ) # mainPanel
         ) # sidebarLayout
       ) # fluidPage

    return(ui)

} # createUI
#------------------------------------------------------------------------------------------------------------------------
igvTabPanel <- function()
{
   panel <- tabPanel(title="IGV", value="igvTab",
                     fluidRow(
                        column(8, igvShinyOutput('igvShiny')),
                        column(4, DT::dataTableOutput('dt'))
                        #column(4, sliderInput("obs", "Number of observations:", min = 1, max = 1000, value = 500))
                     ))
   return(panel)

} # igvTabPanel
#------------------------------------------------------------------------------------------------------------------------
createGeneModelPanel <- function()
{

   panel <- tabPanel(title="Create new trena model",  value="geneModelTab",
                  mainPanel(width=10,
                     br(),
                     h4("(From binding sites In the genomic region currently displayed in the IGV view)"),
                     br(),
                     fluidRow(
                       column(4,checkboxGroupInput("footprintDatabases", "Footprint Databases",
                                                   footprintDatabases, selected=footprintDatabases[1])),
                       column(4, checkboxGroupInput("intersectWithRegions", "Intersect footprints with:",
                                                    c("GeneHancer" = "genehancer",
                                                      "Encode DHS" = "encodeDHS"))),
                       column(4, selectInput("expressionSet", "Choose gene expression dataset:",
                                             c("mtx one", "mtx two")))),
                     fluidRow(column(width=1, offset=4, actionButton("buildModelButton", "Build model")))
                     ))
   return(panel)

} # createGeneModelPanel
#------------------------------------------------------------------------------------------------------------------------
setupBuildModelUI <- function(input, output, session)
{
   observeEvent(input$footprintDatabases, {
                   printf("footprintDatabases checkbox group event")
                })
   observeEvent(input$buildModelButton, ignoreInit = TRUE, {
      buildFootprintModel(2000, 500)
      })


} # setupBuildModelUI
#------------------------------------------------------------------------------------------------------------------------
server = function(input, output){
   output$igvShiny <- renderIgvShiny({
     igvShiny(list(genomeName="hg38",
                   initialLocus="PSG1",
                   displayMode="SQUISHED",
                   trackHeight=150))
     })

   output$dt = DT::renderDataTable(mtcars, server = FALSE, width="800px", options = list(scrollX = TRUE, dom='t'))


} # server

app <- shinyApp(ui=createUI(), server=server)

