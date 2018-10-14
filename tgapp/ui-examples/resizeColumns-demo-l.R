library(shiny)
library(shinyjs)
ui <- fluidPage(useShinyjs(), br(), wellPanel(fluidRow(
  column(width=4, id="spcol", actionButton("dummy", "dummy")),
  column(width=8, id="main", wellPanel(actionButton("sideBarControl", label = "Show/Hide")))
)))

server <- function(input, output) {
  observeEvent(input$sideBarControl, {
    shinyjs::toggle(id = "spcol")
    shinyjs::toggleClass("main", "col-sm-8")
    shinyjs::toggleClass("main", "col-sm-12")
  })
}

app <- shinyApp(ui = ui, server = server)
