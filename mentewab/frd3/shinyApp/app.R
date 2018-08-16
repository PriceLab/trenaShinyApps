library(shiny)
library(DT)
library(igvShiny)
library(htmlwidgets)
library(GenomicRanges)
library(later)
library(MotifDb)
library(RUnit)
library(RColorBrewer)
#----------------------------------------------------------------------------------------------------
totalColorCount <- 12
colors <- brewer.pal(totalColorCount, "Paired")
colorNumber <- 0
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

load("data/tbl.bindingSites.frd3-transcript1.RData") #  "tbl.bindingSites.1"
load("data/tbl.bindingSites.frd3-transcript2.RData") #  "tbl.bindingSites.2
#----------------------------------------------------------------------------------------------------
load("data/tbl.model.frd3-transcript1.RData")  # tbl.model.1
load("data/tbl.model.frd3-transcript2.RData")  # tbl.model.2
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

        selectInput("geneModel", "Choose regulatory model:",
                    c(" - " = "no.model",
                      "transcript 1"="transcript.1.model",
                      "transcript 2"="transcript.2.model")),

        selectInput("addTrack", "Add Track:",
                   c(" - ",
                     "bud DHS"="bud.dhs",
                     "leaf DHS"="leaf.dhs"
                     )),

        selectInput("addBindingSites", "Add TF binding sites for regulatory model:",
                   c(" - ",
                     "transcript 1"="transcript.1.model",
                     "transcript 2"="transcript.2.model"
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

   dt.proxy <- dataTableProxy("geneModelTable")

   observeEvent(input$displayGenomicRegion, {
      regionName <- input$displayGenomicRegion
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
      updateTabsetPanel(session, "trenaTabs", select="igvTab");
      if(trackName == " - ") return();
      tbl.dhs <- switch(trackName, "bud.dhs"=tbl.buds, "leaf.dhs"=tbl.leaf)
      myFunc <- function() loadBedGraphTrack(session, trackName, tbl.dhs)
      later(myFunc, delay=1);
      })

   observeEvent(input$addBindingSites, ignoreInit= TRUE, {
      modelName <- input$addBindingSites;
      printf("addBindingSites: %s", modelName);
      tbl.bs <- switch(modelName, "transcript.1.model" = tbl.bindingSites.1, "transcript.2.model" = tbl.bindingSites.2)
      transcript.number <- switch(modelName, "transcript.1.model" = 1, "transcript.2.model" = 2)
      tbl.bs <- tbl.bs[, c("chrom", "start", "stop", "score", "motif", "geneSymbol")]
      colnames(tbl.bs) <- c("chr", "start", "end", "value", "sampleID", "geneSymbol")
      printf("loadingBedGraphTrack, dim: %d %d", nrow(tbl.bs), ncol(tbl.bs))
      myFunc <- function(){
         printf("entering myFunc, rows in tbl: %d", nrow(tbl.bs));
         tfs <- sort(unique(tbl.bs$geneSymbol))
         printf(" %d tfs", length(tfs))
         for(tf in tfs){
            colorNumber <<- (colorNumber + 1) %% totalColorCount
            tbl.bs.sub <- subset(tbl.bs, geneSymbol==tf)
            track.name <- sprintf("%s.%d", tf, transcript.number)
            loadBedTrack(session, track.name, tbl.bs.sub, color=colors[colorNumber], trackHeight=25)
            } # for tf
         }
        # myFunc <- function() loadBedGraphTrack(session, modelName, tbl.bs, color=colors[colorNumber], trackHeight=30)
      updateTabsetPanel(session, "trenaTabs", select="igvTab");
        # allow time for the igv tab to be displayed
      later(myFunc, delay=1);
      })

   output$igvShiny <- renderIgvShiny({
     printf("rendering igvShiny");
     igvShiny(list(genomeName="tair10",
                   initialLocus=sprintf("%s:%d-%d", chromosome, loc.start, loc.end)))
     })


   observeEvent(input$geneModel, ignoreInit = TRUE, {
      model.name <- isolate(input$geneModel);
      if(model.name != " - "){
         printf("geneModel menu changed: %s", model.name)
         updateTabsetPanel(session, "trenaTabs", select="geneModelTab")
         }
      })

   output$geneModelTable <- renderDT({
      model.name <- input$geneModel;
      if(model.name != "no.model"){
         printf("    rendering DT, presumably because input$geneModel changes, model.name: %s", model.name);
         selection="single"
         tbl.model <- switch(model.name, "transcript.1.model"=tbl.model.1, "transcript.2.model"=tbl.model.2)
         printf("data.frame dimensions: %d, %d", nrow(tbl.model), ncol(tbl.model))
         if(nrow(tbl.model) > 10)
            tbl.model <- tbl.model[1:10,]
         rownames(tbl.model) <- NULL
         # later(function() {updateTabsetPanel(session, "trenaTabs", select="geneModelTab")}, 1);
         return(tbl.model)
         } # not "no.model"
     },
     selection = 'single',
     options=list(pageLength=25, dom="t")
     ) # renderDT

  } # server

#----------------------------------------------------------------------------------------------------
shinyOptions <- list()

if(Sys.info()[["nodename"]] == "trena.systemsbiology.net"){
   port <- 60031
   printf("running on trena, using port %d", port)
   shinyOptions <- list(host="0.0.0.0", port=port, launch.browser=FALSE)
   }
app <- shinyApp(ui=ui,server=server, options=shinyOptions)

