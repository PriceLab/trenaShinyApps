library(shiny)
library(shinydashboard)
library(shinyjs)
library(DT)
#------------------------------------------------------------------------------------------------------------------------
targetGene <- list(hugo="PSG1", ensembl="ENSG00000231924", combined="ENSG00000231924|PSG1")
footprintDatabases <- c("db one", "db two")
#------------------------------------------------------------------------------------------------------------------------
createBody <- function()
{
   body <- fluidRow(
      column(8,
             box(plotOutput("plot1", height = 250), width=12)),
      column(4,
             box(selectInput("modelSelector", label=NULL, c("mtcars", "airquality")), height=60, width=12),
             box(DTOutput("table"), width=12))
      )

   body

} # createBody
#------------------------------------------------------------------------------------------------------------------------
ui <- dashboardPage(
  dashboardHeader(title=sprintf("trena %s", targetGene$hugo)),
  dashboardSidebar(),
  dashboardBody(createBody()),
  useShinyjs()
  )
#------------------------------------------------------------------------------------------------------------------------
server <- function(input, output){

   output$table = DT::renderDataTable(mtcars,
                                      width="800px",
                                      class='nowrap display',
                                      options=list(scrollX=TRUE, dom='t'))
} # server
#------------------------------------------------------------------------------------------------------------------------
app <- shinyApp(ui, server)
