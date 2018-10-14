library(shiny)
library(shinydashboard)
library(shinyjs)
library(igvShiny)
#------------------------------------------------------------------------------------------------------------------------
genomeName <- "hg38"
targetGene <- list(hugo="PSG1", ensembl="ENSG00000231924", combined="ENSG00000231924|PSG1")
#------------------------------------------------------------------------------------------------------------------------
createBody <- function(input, output)
{
   body <- fluidRow(
      column(width=8, id="igvColumn", height="800px", igvShinyOutput('igvShiny', height="800px")),
      column(width=4, id="dataTableColumn",
         box(id = "myBox", title = "mtcars", width = '800px',
             DTOutput("table")
             ))
     )

  body

} # createBody
#------------------------------------------------------------------------------------------------------------------------
createSidebar <- function()
{
  dashboardSidebar(
     actionButton(inputId = "igvHideButton", label = "Toggle IGV"),
     actionButton(inputId = "tableHideButton", label = "Toggle table")
     )

} # createSidebar
#------------------------------------------------------------------------------------------------------------------------
ui <- dashboardPage(
    dashboardHeader(title="hide box"),
    createSidebar(),
    dashboardBody(createBody()),
    useShinyjs()
    )


#------------------------------------------------------------------------------------------------------------------------
server <- function(input, output){

   output$table = DT::renderDataTable(mtcars,
                                      width="800px",
                                      class='nowrap display',
                                      options=list(scrollX=TRUE, dom='t', pageLength=100))

   output$igvShiny <- renderIgvShiny({
         igvShiny(options=list(genomeName=genomeName,
                               initialLocus=targetGene$hugo),
                  height=800)
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
