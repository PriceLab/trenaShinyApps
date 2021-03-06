library(shiny)
library(Rmpfr)

shinyUI(
    fluidPage(
        titlePanel("Euler's e in arbitrary precision", "Euler's e"),
        tags$br(),

        fluidRow(
            column(2, sliderInput("precision", "Number of Precision Bits", min = 2, max = 256,
                    value = 10))
        ),
        fluidRow(
            column(12, tags$h1(textOutput("result")))

        ),
        tags$br(),
        tags$p("We wish  to acknowledge Leonhard Euler for his number and Martin Maechler for his Rmpfr package.")
    )
)


shinyServer(function(input, output){


      output$result <- renderText({

            precisionBits <- input$precision
            one <- mpfr(1, precBits = precisionBits)
            e <- exp(one)
            # TODO fix printing...
            x <- capture.output(print(e, ndigits = precisionBits))[2]
            gsub("^\\[1\\] (.+)$", "\\1", x)


          })

      })

#----------------------------------------------------------------------------------------------------
port <- 60030
launchBrowser=TRUE

if(Sys.info()[["nodename"]] == "trena.systemsbiology.net"){
   launchBrowser=FALSE
   printf("running on trena, using port %d", port)
   }

shinyOptions <- list(host="0.0.0.0", port=port, launch.browser=launchBrowser)
# app <- shinyApp(ui=ui, server=server, options=shinyOptions)

