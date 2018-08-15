library(shiny)
library(DT)
library(igvShiny)
library(htmlwidgets)
library(GenomicRanges)
library(later)
library(MotifDb)
library(RUnit)
#----------------------------------------------------------------------------------------------------
load("data/tbls.frd3.dhs.budAndLeaf.RData")   # "tbl.frd3buds" "tbl.frd3leaf"

tbl.buds <- tbl.frd3buds
tbl.buds$sampleID <- paste("seg", seq_len(nrow(tbl.buds)), sep=".")
tbl.buds$chrom <- rep("3", nrow(tbl.buds))
colnames(tbl.buds) <- c("chr", "start", "end", "value", "sampleID")

tbl.leaf <- tbl.frd3leaf
tbl.leaf$sampleID <- paste("seg", seq_len(nrow(tbl.leaf)), sep=".")
tbl.leaf$chrom <- rep("3", nrow(tbl.leaf))
colnames(tbl.leaf) <- c("chr", "start", "end", "value", "sampleID")
#----------------------------------------------------------------------------------------------------
load("data/tbl.model.frd3-transcript2.RData")  # tbl.model
#----------------------------------------------------------------------------------------------------
geneSymbol <- "FRD3"
orf <- "AT3G08040"
tss.orf.1 <- 2569396
tss.orf.2 <- 2572151
targetGene <- orf
chromosome <- "chr3"
loc.start <- 2566277
loc.end <- 2572151
chrom.loc.string <- sprintf("%s:%d-%d", chromosome, loc.start - 1000, loc.end + 1000)
#----------------------------------------------------------------------------------------------------
ui <- fluidPage(

  tags$script(JS('setInterval(function(){console.log("keepAlive"); Shiny.setInputValue("foo", "var");}, 1 * 60 * 1000)')),
  titlePanel(sprintf("%s Trena Models ", geneSymbol)),

  sidebarLayout(
     sidebarPanel(
        width=3,
        selectInput("displayGenomicRegion", "Display Genomic Region:",
                    c("FRD3" = "fullgene",
                      "AT3G08040.1 promoter" = "transcript.1.promoter",
                      "AT3G08040.2 promoter" = "transcript.2.promoter")),

        selectInput("geneModel", "Choose model:",
                    c("transcript 2", "transcript.2.model")),

        selectInput("addTrack", "Add Track:",
                   c(" - ",
                     "bud DHS"="bud.dhs",
                     "leaf DHS"="leaf.dhs"
                     ))
       ),  # sidebarPanel

     mainPanel(
        tabsetPanel(type="tabs",
                    id="trenaTabs",
                    tabPanel(title="IGV",             value="igvTab",          igvShinyOutput('igvShiny')),
                    tabPanel(title="trena model",     value="geneModelTab",
                             mainPanel(radioButtons("tfSelectionChoice", "Row (transcription factor) selection will display:",
                                                    c("XY plot" = "xyPlot",
                                                      "Binding sites" = "displayBindingSites"),
                                                    inline=TRUE),
                                       DTOutput("geneModelTable")
                                      )),
                    tabPanel(title="TF/FRD3 xy plot",  value="plotTab", plotOutput("xyPlot", height=800))
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
      printf("display event: %s", regionName)
      regions <- list(fullgene="3:2,565,078-2,572,981",
                      transcript.1.promoter=sprintf("%s:%d-%d", chromosome, tss.orf.1-499, tss.orf.1+2000),
                      transcript.2.promoter=sprintf("%s:%d-%d", chromosome, tss.orf.2-499, tss.orf.2+2000)
                      )
      new.region.name <- input$displayGenomicRegion
      if(!new.region.name %in% names(regions))
         return()
      roi <- regions[[new.region.name]]
      printf("new region: %s", roi)
      updateTabsetPanel(session, "trenaTabs", select="igvTab");
      myFunc <- function(){
         session$sendCustomMessage(type="showGenomicRegion", message=(list(region=roi)))
         } # myFunc
      later(myFunc, 0.5)
      later(function() {updateSelectInput(session, "displayGenomicRegion", selected=character(0))}, 1)
      })

   observeEvent(input$addTrack, {
      trackName <- input$addTrack
      printf("addTrack: %s", trackName)
      tbl.dhs <- switch(trackName, "bud.dhs"=tbl.buds, "leaf.dhs"=tbl.leaf)
      loadBedGraphTrack(trackName, tbl.dhs)
      })

   output$igvShiny <- renderIgvShiny({
     igvShiny(list(genomeName="tair10", initialLocus=sprintf("%s:%d-%d", chromosome, loc.start, loc.end)))
     })

   output$geneModelTable <- renderDT({
      model.name <- input$geneModel;
      updateTabsetPanel(session, "trenaTabs", select="geneModelTab");
      printf("    rendering DT, presumably because input$geneModel changes, model.name: %s", model.name);
      selection="single"
      tbl.model <- tbl.model
      printf("data.frame dimensions: %d, %d", nrow(tbl.model), ncol(tbl.model))
      if(nrow(tbl.model) > 20)
         tbl.model <- tbl.model[1:20,]
      rownames(tbl.model) <- NULL
      #state[["currentModel"]] <- tbl.model
      #state[["currentModelName"]] <- model.name
      return(tbl.model)
      },
    selection = 'single',
    options=list(pageLength=25, dom="t"))

  } # server

#----------------------------------------------------------------------------------------------------
shinyOptions <- list()

if(Sys.info()[["nodename"]] == "trena.systemsbiology.net"){
   port <- 60031
   printf("running on trena, using port %d", port)
   shinyOptions <- list(host="0.0.0.0", port=port, launch.browser=FALSE)
   }
app <- shinyApp(ui=ui,server=server, options=shinyOptions)

