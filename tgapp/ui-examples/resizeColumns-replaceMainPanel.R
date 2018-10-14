library(shiny)
library(shinyjs)
ui <- fluidPage(useShinyjs(), br(), wellPanel(fluidRow(
  column(width=4, id="spcol", actionButton("dummy", "dummy")),
  column(width=8, id="main1", wellPanel(actionButton("sideBarControl", label = "Show/Hide"))),
  column(width=12, id="main2", wellPanel(actionButton("sideBarControl2",  label = "Show/Hide")))
)))

server <- function(input, output) {
  shinyjs::toggle(id = "main2")

  observeEvent(input$sideBarControl+input$sideBarControl2, {
    shinyjs::toggle(id = "spcol")
    shinyjs::toggle(id = "main1")
    shinyjs::toggle(id = "main2")
  })
}
app <- shinyApp(ui = ui, server = server)
