library(shiny)
library(DT)
library(igvShiny)
library(htmlwidgets)
library(GenomicRanges)
library(later)
library(MotifDb)
library(RUnit)
#----------------------------------------------------------------------------------------------------
targetGene <- "FRD3"
addResourcePath("tmp", "tmp") # so the shiny webserver can see, and return track files for igv
#----------------------------------------------------------------------------------------------------
chromosome <- "chr3"
loc.min <- 2557465
loc.max <- 2580964
#----------------------------------------------------------------------------------------------------
ui <- fluidPage(

  # includeScript("message-handler.js"),
  tags$script(JS('setInterval(function(){console.log("keepAlive"); Shiny.setInputValue("foo", "var");}, 1 * 60 * 1000)')),

  titlePanel(sprintf("%s Trena Models ", targetGene)),

  sidebarLayout(
     sidebarPanel(
        width=3,
        selectInput("displayGenomicRegion", "Display Genomic Region:",
                    c("All INPP5D enhancers" = "enhancers",
                      "INPP5D +/- 40kb" = "tssPlusMinus40kb",
                      "TSS-proximal interesting snps"="interestingSnpsNearTSS")),

        selectInput("addTrack", "Add Track:",
                   c(" - ",
                     "Enhancers"="enhancers",
                     "IGAP GWAS SNPs"="snps",
                     "SNPs in Enhancers"="enhancer.snps",
                     "Motif-breaking SNPs"="breaking.snps"
                     ))
     ),

     mainPanel(
        tabsetPanel(type="tabs",
                    id="trenaTabs",
                    tabPanel(title="IGV",             value="igvTab",          igvShinyOutput('igvShiny'))
                    ) # tabsetPanel
      )  # mainPanel
    ) # sidebarLayout
  ) # fluidPage

#--------------------------------------------------------------------------------
server <- function(input, output, session) {

   session$sendCustomMessage(type="showGenomicRegion", message=(list(region="FLT1")))
   dt.proxy <- dataTableProxy("geneModelTable")

   observeEvent(input$displayGenomicRegion, {
      regionName <- input$displayGenomicRegion
      #printf("display event: %s", regionName)
      regions <- list(enhancers="chr2:232850784-233,458,236",
                      tssPlusMinus40kb="chr2:233,058,937-233,248,903",
                      interestingSnpsNearTSS="chr2:233,092,624-233,095,615"
                      )
      new.region.name <- input$displayGenomicRegion
      if(!new.region.name %in% names(regions))
         return()
      roi <- regions[[new.region.name]]
      updateTabsetPanel(session, "trenaTabs", select="igvTab");
      myFunc <- function(){
         session$sendCustomMessage(type="showGenomicRegion", message=(list(region=roi)))
         } # myFunc
      later(myFunc, 2)
      later(function() {updateSelectInput(session, "displayGenomicRegion", selected=character(0))}, 1)
      })

   output$igvShiny <- renderIgvShiny({
     igvShiny(list(genomeName="tair10", initialLocus=sprintf("%s:%d-%d", chromosome, loc.min, loc.max)))
     })

  } # server

#----------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
shinyOptions <- list()

if(Sys.info()[["nodename"]] == "trena.systemsbiology.net"){
   port <- 60031
   printf("running on trena, using port %d", port)
   shinyOptions <- list(host="0.0.0.0", port=port, launch.browser=FALSE)
   }
app <- shinyApp(ui=ui,server=server, options=shinyOptions)

