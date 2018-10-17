library(shiny)
library(shinycssloaders)

ui <- fluidPage(
   fluidRow(
      column(width=6,
             wellPanel(
                withSpinner(plotOutput("plot"),type=1)
                )
             ) # column
      ) # fluidRow
   ) # fluidPage

server <- function(input, output,session) {
    output[["plot"]] <- renderPlot({
      Sys.sleep(2) # just for demo so you can enjoy the animation
      plot(x = runif(1e4), y = runif(1e4))
      })
    } # server

app <- shinyApp(ui = ui, server = server)
