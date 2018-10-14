library(shiny)
library(igvShiny)

#------------------------------------------------------------------------------------------------------------------------
createUI <- function()
{
   ui <- fluidPage(
     includeCSS("custom.css"),
     sidebarLayout(
        sidebarPanel(width=2,
           actionButton("dummyButton", "sidebar button")
           ),
        mainPanel(width=10,
           igvShinyOutput('igvShiny', height="800px"),
           actionButton("button2", "mainPanel button")
           )
         ) # sidebarLayout
      ) # fluidPage

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

