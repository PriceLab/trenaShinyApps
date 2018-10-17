library(shiny)
library(shinydashboard)
library(shinyjs)
library(igvShiny)
library(DT)
#------------------------------------------------------------------------------------------------------------------------
targetGene <- list(hugo="PSG1", ensembl="ENSG00000231924", combined="ENSG00000231924|PSG1")
footprintDatabases <- c("db one", "db two")
expressionMatrixNames <- c("mtx1", "mtx2", "mtx3")
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
        actionButton(inputId = "tableHideButton", label = "Toggle table")
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
                      selectInput("modelSelector", NULL,  c("mtcars", "something else")))),
            fluidRow(
               column(width=9, id="igvColumn", igvShinyOutput('igvShiny', height="800px")),
               column(width=3, id="dataTableColumn", title="mtcars", DTOutput("table")))
            ), # tabItem 1
         tabItem(tabName="buildModels",
             h4("Build trena regulatory model from DNase footprints in the genomic region currently displayed in the IGV view, with these constraints:"),
             br(),
             fluidRow(
                column(2, offset=1, checkboxGroupInput("footprintDatabases", "Footprint Databases",
                                            footprintDatabases, selected=footprintDatabases[1])),
                column(2, checkboxGroupInput("intersectWithRegions", "Intersect footprints with:",
                                             c("GeneHancer" = "genehancer",
                                               "Encode DHS" = "encodeDHS",
                                               "Use all footprints in region" = "allDNAForFootprints"))),
                column(3, selectInput("expressionSet", "Choose gene expression dataset:",
                                      expressionMatrixNames))),
             br(),
             fluidRow(column(width=3, offset=3, textInput("modelNameTextInput", "Model name", width=240))),
             fluidRow(column(width=1, offset=3, actionButton("buildModelButton", "Build model")))
            ) # tabItem 2
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
server <- function(input, output){

   output$igvShiny <- renderIgvShiny({
      options <- list(genomeName="hg38",
                      initialLocus="PSG1",
                      displayMode="EXPANDED",
                      trackHeight=300)
      igvShiny(options) # , height=800)
      })

   observeEvent(input$modelNameTextInput, {
       currentName <- isolate(input$modelNameTextInput)
       printf("currentName size: %d", nchar(currentName))
       if(nchar(currentName) > 5)
          shinyjs::enable("buildModelButton")
       else
          shinyjs::disable("buildModelButton")
       })

   output$table = DT::renderDataTable(mtcars[1:20,],
                                      width="800px",
                                      class='nowrap display',
                                      extensions='FixedColumns',
                                      options=list(scrollX=TRUE, dom='t', pageLength=100,
                                                   fixedColumns=list(leftColumns=2)))



    observeEvent(input$igvHideButton, {
      printf("igvHideButton: %s", input$igvHideButton)
      if(input$igvHideButton %% 2 == 1){
         printf("  --- hiding igv, widening dataTable")
         shinyjs::hide(id = "igvColumn")
         shinyjs::toggleClass("dataTableColumn", "col-sm-3")
         shinyjs::toggleClass("dataTableColumn", "col-sm-12")
       } else {
         printf("  --- showing igv, narrowing dataTable")
         shinyjs::toggleClass("dataTableColumn", "col-sm-12")
         shinyjs::toggleClass("dataTableColumn", "col-sm-3")
         shinyjs::show(id = "igvColumn")
         }
       })

   observeEvent(input$tableHideButton, {
      printf("tableHideButton: %s", input$tableHideButton)
      if(input$tableHideButton %% 2 == 1){
         shinyjs::hide(id = "dataTableColumn")
         shinyjs::hide(id = "modelSelectorColumn")
         shinyjs::toggleClass("igvColumn", "col-sm-9")
         shinyjs::toggleClass("igvColumn", "col-sm-12")
      } else {
         shinyjs::toggleClass("igvColumn", "col-sm-12")
         shinyjs::toggleClass("igvColumn", "col-sm-9")
         shinyjs::show(id = "dataTableColumn")
         shinyjs::show(id = "modelSelectorColumn")
         }
      })


} # server
#------------------------------------------------------------------------------------------------------------------------
app <- shinyApp(ui, server)
