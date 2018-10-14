library(shiny)
library(igvShiny)

#------------------------------------------------------------------------------------------------------------------------
createUI <- function()
{
   ui <- fluidPage(
      includeCSS("custom.css"),
      fluidRow(width=12,
            igvShinyOutput('igvShiny', height="800px"),
            actionButton("dummyButton", "Dummy"))
           )

  return(ui)

} # createUI
#------------------------------------------------------------------------------------------------------------------------
server <- function(input, output)
{
   output$igvShiny <- renderIgvShiny({
      options <- list(genomeName="hg38",
                      initialLocus="PSG1",
                      displayMode="EXPANDED",
                      trackHeight=300)
      igvShiny(options) # , height=800)
      })

} # server
#------------------------------------------------------------------------------------------------------------------------
app <- shinyApp(ui=createUI(), server=server)

