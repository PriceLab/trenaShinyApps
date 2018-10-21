library(shiny)
library(trenaSGM)
stopifnot(packageVersion("trenaSGM") >= "0.99.56")
genomeName <- "hg38"
targetGene <- "COL1A1"
f <- "../trena.model.1kdemo.buildSpec.RData"
#------------------------------------------------------------------------------------------------------------------------
runTrena <- function()
{
  message(load(f))
  fpBuilder <- FootprintDatabaseModelBuilder(genomeName, targetGene, build.spec, quiet=FALSE)
  x <- build(fpBuilder)
  message(paste(head(x$model$gene, n=10),  collapse=", "))

} # runTrena
#------------------------------------------------------------------------------------------------------------------------
ui <- fluidPage(
  shinyjs::useShinyjs(),
  actionButton("go", "go"),
  wellPanel(style="overflow-y:scroll; height:300px", pre(id = "console"))
  ) # fluidPage

server <- function(input, output, session) {

  observeEvent(input$go, {
    withCallingHandlers(
      runTrena(),
      message = function(m) {   # can use "warning" instead/on top of "message" to catch warnings too
        shinyjs::html(id="console", html=m$message, add=TRUE)
      }
    )
  })
}

app <- shinyApp(ui, server)
