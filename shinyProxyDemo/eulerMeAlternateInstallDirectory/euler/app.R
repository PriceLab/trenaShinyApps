library(shiny)
library(Rmpfr)

ui <- fluidPage(
        titlePanel("Euler's e in arbitrary precision (pshannon 456p 9 jul)", "Euler's e"),
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


server <- function(input, output, session){


      output$result <- renderText({

            precisionBits <- input$precision
            one <- mpfr(1, precBits = precisionBits)
            e <- exp(one)
            # TODO fix printing...
            x <- capture.output(print(e, ndigits = precisionBits))[2]
            gsub("^\\[1\\] (.+)$", "\\1", x)


          })

      } # server

#----------------------------------------------------------------------------------------------------
launchBrowser <- TRUE
port <- 60030

if(Sys.info()[["nodename"]] == "trena.systemsbiology.net"){
   printf("running on trena, using port %d", port)
   launchBrowser <- FALSE
   }
shinyOptions <- list(host="0.0.0.0", port=port, launch.browser=launchBrowser)
app <- shinyApp(ui=ui,server=server, options=shinyOptions)
