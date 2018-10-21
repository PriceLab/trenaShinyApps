library(shiny)

longfunc <- function() {
  message("Thinking...")
  #print("Thinking...")
  Sys.sleep(1)
  message("Still thinking...")
  #print("Still thinking...")
  Sys.sleep(1)
  message("Done")
  #print("Done")
}

withConsoleRedirect <- function(containerId, expr) {
  # Change type="output" to type="message" to catch stderr
  # (messages, warnings, and errors) instead of stdout.
  txt <- capture.output(results <- expr, type = "message")
  if (length(txt) > 0) {
    insertUI(paste0("#", containerId), where = "beforeEnd",
             ui = paste0(txt, "\n", collapse = "")
    )
  }
  results
}

ui <- fluidPage(
  shinyjs::useShinyjs(),
  actionButton("joe", "Joe"),
  actionButton("dean", "Dean"),
  pre(id = "console")
)

server <- function(input, output, session) {

  observeEvent(input$joe, {
     withConsoleRedirect("console", {
      longfunc()
    })
  })

  observeEvent(input$dean, {
    withCallingHandlers(
      longfunc(),
      # can use "warning" instead/on top of "message" to catch warnings too
      message = function(m) {
        shinyjs::html("console", m$message, TRUE)
      }
    )
  })
}

app <- shinyApp(ui, server)
