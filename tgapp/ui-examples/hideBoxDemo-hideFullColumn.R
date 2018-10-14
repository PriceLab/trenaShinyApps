library(shiny)
library(shinydashboard)
library(shinyjs)

#------------------------------------------------------------------------------------------------------------------------
createBody <- function(input, output)
{
   body <- fluidRow(
      column(width=10, id="myBoxColumn",
         #actionButton(inputId = "button", label = "Toggle mtcars"),
         box(id = "myBox", title = "mtcars", width = '800px',
             DTOutput("table")
             )
         )
     )

  body

} # createBody
#------------------------------------------------------------------------------------------------------------------------
createSidebar <- function()
{
  actionButton(inputId = "sidebarHideButton", label = "Toggle table")

} # createSidebar
#------------------------------------------------------------------------------------------------------------------------
ui <- dashboardPage(
    dashboardHeader(title="hide box"),
    dashboardSidebar(createSidebar()),
    dashboardBody(createBody()),
    useShinyjs()
    )


server <- function(input, output){

   output$table = DT::renderDataTable(mtcars,
                                      width="800px",
                                      class='nowrap display',
                                      options=list(scrollX=TRUE, dom='t', pageLength=100))

   observeEvent(input$sidebarHideButton, {
   if(input$sidebarHideButton %% 2 == 1){
     #shinyjs::hide(selector='$("#myBox").parent()')
     shinyjs::hide(id = "myBoxColumn")
     } else {
     shinyjs::show(id = "myBoxColumn")
     }
    })
}

app <- shinyApp(ui, server)
