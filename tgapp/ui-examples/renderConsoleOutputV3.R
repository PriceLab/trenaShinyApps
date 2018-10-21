ui <- shinyUI(fluidPage(
  sidebarLayout(
    sidebarPanel(
    ),
    mainPanel(verbatimTextOutput("console"))
  )
))


server <- function(input, output)
{
  values <- reactiveValues()

  queryMagic <- function() {
    print("Warning")
    return("Data")
    }

  output$console <- renderPrint({
    logText()
    return(print(values[["log"]]))
    # You could also use grep("Warning", values[["log"]]) to get warning messages and use shinyBS package
    # to create alert message
    })

  logText <- reactive({
    values[["log"]] <- capture.output(data <- queryMagic())

  })
}


app <- shinyApp(ui = ui, server = server)
