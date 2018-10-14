library(shiny)
library(shinydashboard)
library(shinyjs)
library(igvShiny)
#------------------------------------------------------------------------------------------------------------------------
targetGene <- list(hugo="PSG1", ensembl="ENSG00000231924", combined="ENSG00000231924|PSG1")
footprintDatabases <- c("db one", "db two")
#------------------------------------------------------------------------------------------------------------------------
createSidebar <- function()
{
  dashboardSidebar(
     actionButton(inputId = "igvHideButton", label = "Toggle IGV"),
     actionButton(inputId = "tableHideButton", label = "Toggle table")
     )

} # createSidebar
#------------------------------------------------------------------------------------------------------------------------
createBody <- function()
{
   body <- fluidRow(
      tags$head(tags$style(HTML('
        .main-header .logo {
           font-family: "Georgia", Times, "Times New Roman", serif;
          font-weight: bold;
          font-size: 24px;
          }
        #igvShiny{
          background: lightYellow;
          border: 1px solid black;
          border-radius: 5px;
          margin: 5px;
          margin-right: 15px;
          }
       '))),
      column(width=8, id="igvColumn",
             igvShinyOutput('igvShiny', height="800px")
             ),
      column(width=4, id="dataTableColumn", title="mtcars",
             actionButton("button2", "mainPanel button")
             )
      )

   return(body)

} # createBody
#------------------------------------------------------------------------------------------------------------------------
ui <- dashboardPage(
  dashboardHeader(title=sprintf("trena %s", targetGene$hugo)),
  createSidebar(),
  #dashboardSidebar(),
  dashboardBody(createBody()),
  useShinyjs()
  )
#------------------------------------------------------------------------------------------------------------------------
server <- function(input, output){

   output$igvShiny <- renderIgvShiny({
      options <- list(genomeName="hg38",
                      initialLocus="PSG1",
                      displayMode="EXPANDED",
                      trackHeight=300)
      igvShiny(options) # , height=800)
      })

    observeEvent(input$igvHideButton, {
      if(input$igvHideButton %% 2 == 1){
         shinyjs::hide(id = "igvColumn")
       } else {
         shinyjs::show(id = "igvColumn")
         }
       })

   observeEvent(input$tableHideButton, {
      if(input$tableHideButton %% 2 == 1){
         shinyjs::hide(id = "dataTableColumn")
      } else {
        shinyjs::show(id = "dataTableColumn")
         }
      })


} # server
#------------------------------------------------------------------------------------------------------------------------
app <- shinyApp(ui, server)
