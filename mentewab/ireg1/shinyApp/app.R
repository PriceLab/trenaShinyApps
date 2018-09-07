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
totalColorCount <- 6
colors <- brewer.pal(totalColorCount, "Dark2")
colorNumber <- 0
state <- new.env(parent=emptyenv())
state$currentModel <- NULL
#----------------------------------------------------------------------------------------------------
load("data/tbls.dhs.RData")   # tbl.leaf.dhs tbl.buds.dhs
load("data/tbl.bindingSites.RData") #  "tbl.bindingSites"
#----------------------------------------------------------------------------------------------------
load("data/tbl.model.RData")  # tbl.model
#----------------------------------------------------------------------------------------------------
load("data/mtx.zinc.22810x42.RData")           # mtx
#----------------------------------------------------------------------------------------------------
geneSymbol <- "IREG1"
orf <- "AT2G38460"
tss.1 <- 16103164
targetGene <- orf
chromosome <- "chr2"
loc.start <- 16088281
loc.end <- 16117304
chrom.loc.string <- sprintf("%s:%d-%d", chromosome, loc.start - 1000, loc.end + 1000)
#----------------------------------------------------------------------------------------------------
ui <- fluidPage(

  tags$script(JS('setInterval(function(){console.log("keepAlive"); Shiny.setInputValue("foo", "var");}, 1 * 60 * 1000)')),
  titlePanel(sprintf("%s Trena Models ", geneSymbol)),

  sidebarLayout(
     sidebarPanel(
        width=3,
        selectInput("displayGenomicRegion", "Display Genomic Region:",
                    c("IREG1" = "gene",
                      "promoter" = "promoter",
                      "neighborhood" = "neighborhood"),
                      selected="neighborhood"),

        selectInput("addTrack", "Add Track:",
                   c("",
                     "bud DHS" = "bud.dhs",
                     "leaf DHS" = "leaf.dhs"),
                     selected="bud.dhs"),

        #sliderInput("fimo.motif.significance", "FIMO motif score",
        #            value = 2,
        #            step=0.1,
        #            min = 0,
        #            max = 12,
        #            round=-2),

        sliderInput("dhs.threshold", "DNase score",
                    value = 2,
                    step=0.1,
                    min = 0,
                    max = 10,
                    round=-2),

        actionButton("addBindingSites", "Add TF binding sites from regulatory model"),
        HTML("<br>"),
        actionButton("removeOptionalTracks", "Remove tracks")
        ),  # sidebarPanel

     mainPanel(
        tabsetPanel(type="tabs",
                    id="trenaTabs",
                    tabPanel(title="IGV",             value="igvTab",          igvShinyOutput('igvShiny')),
                    tabPanel(title="trena model",     value="geneModelTab",
                             mainPanel(h5("Row (transcription factor) selection will display gene expression scatter plot."),
                                       DTOutput("geneModelTable")
                                      )),
                    tabPanel(title="TF/IREG1 xy plot",  value="plotTab", plotOutput("xyPlot", height=800))
                    ) # tabsetPanel
      )  # mainPanel
    ) # sidebarLayout
  ) # fluidPage

#--------------------------------------------------------------------------------
server <- function(input, output, session) {

   dt.proxy <- dataTableProxy("geneModelTable")

   observeEvent(input$displayGenomicRegion, {
      regionName <- input$displayGenomicRegion
      regions <- list(gene="2:16,102,925-16,107,163",
                      promoter=sprintf("%s:%d-%d", chromosome, tss.1-2000, tss.1+499),
                      neighborhood="2:16,090,142-16,115,989"
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
      if(nchar(trackName) == 0) return()
      printf("addTrack: %s", trackName)
      updateTabsetPanel(session, "trenaTabs", select="igvTab");
      if(trackName == " ") return();
      tbl.dhs <- switch(trackName, "bud.dhs"=tbl.buds.dhs, "leaf.dhs"=tbl.leaf.dhs)
      myFunc <- function() loadBedGraphTrack(session, trackName, tbl.dhs, color="blue", trackHeight=50,
                                             autoscale=FALSE, min=0, max=10)
      later(myFunc, delay=0.5);
      later(function() {updateSelectInput(session, "addTrack", selected=character(0))}, 1)
      })

   observeEvent(input$addBindingSites, ignoreInit=TRUE, {
      modelName <- input$addBindingSites;
      printf("addBindingSites: %s", modelName);
      fimo.threshold <- isolate(input$fimo.motif.significance)
      dhs.threshold <-  isolate(input$dhs.threshold)
      printf("fimo: %f   dhs: %f", fimo.threshold, dhs.threshold)
      colors <- brewer.pal(8, "Dark2")
      totalColorCount <- length(colors)
      colorNumber <- 0
      for(tf in unique(tbl.bindingSites$geneSymbol)){
         tbl.tfbs <- subset(tbl.bindingSites, geneSymbol==tf)[, c("chrom", "start", "stop", "score")]
         colnames(tbl.tfbs) <- c("chr", "start", "end", "value")
         colorNumber <- (colorNumber %% totalColorCount) + 1
         color <- colors[colorNumber]
         loadBedGraphTrack(session, tf, tbl.tfbs, color=color, trackHeight=25, autoscale=FALSE, min=0, max=20)
         } # for tf
      #tbl.bs <- switch(modelName, "transcript.1.model" = tbl.bindingSites.1,
      #                            "transcript.2.model" = tbl.bindingSites.2,
      #                            "both.transcript.models" = rbind(tbl.bindingSites.1, tbl.bindingSites.2))
      #transcript.number <- switch(modelName, "transcript.1.model" = ".1",
      #                                       "transcript.2.model" = ".2",
      #                                       "both.transcript.models" = "")
      #tbl.bs <- tbl.bs[, c("chrom", "start", "stop", "score", "motif", "geneSymbol")]
      #colnames(tbl.bs) <- c("chr", "start", "end", "value", "sampleID", "geneSymbol")
      #printf("loadingBedGraphTrack, dim: %d %d", nrow(tbl.bs), ncol(tbl.bs))
      #myFunc <- function(){
      #   printf("entering myFunc, rows in tbl: %d", nrow(tbl.bs));
      #   tfs <- sort(unique(tbl.bs$geneSymbol))
      #   printf(" %d tfs", length(tfs))
      #   for(tf in tfs){
      #      colorNumber <<- (colorNumber + 1) %% totalColorCount
      #      tbl.bs.sub <- subset(tbl.bs, geneSymbol==tf)
      #      printf("new tf '%s', %d rows", tf, nrow(tbl.bs.sub))
      #      track.name <- sprintf("%s%s", tf, transcript.number)
      #      loadBedGraphTrack(session, track.name, tbl.bs.sub, color=colors[colorNumber], trackHeight=30)
      #      } # for tf
      #   }
      #  # myFunc <- function() loadBedGraphTrack(session, modelName, tbl.bs, color=colors[colorNumber], trackHeight=30)
      updateTabsetPanel(session, "trenaTabs", select="igvTab");
        # allow time for the igv tab to be displayed
      #later(myFunc, delay=1);
      })

   observeEvent(input$dhs.threshold, {
      dhs.threshold <- isolate(input$dhs.threshold)
      printf(" --- dhs threshold changed, redraw peaks: %f", dhs.threshold)
      x <- findRegionsAboveThreshold(tbl.buds.dhs$value, dhs.threshold, 10)
      tbl.peaks <- with(x, data.frame(chr=rep("chr2", length(starts)), start=starts, end=ends, stringsAsFactors=FALSE))
      printf("tbl.peaks rows: %d", nrow(tbl.peaks))
      tbl.peaks$start <- tbl.peaks$start + tbl.buds.dhs$start[1]
      tbl.peaks$end <- tbl.peaks$end + tbl.buds.dhs$start[1]
      updateTabsetPanel(session, "trenaTabs", select="igvTab")
      myFunc <- function() loadBedTrack(session, "dhs peaks", tbl.peaks, color="darkGreen", trackHeight=25)
      later(myFunc, delay=0.5);
      })

   output$igvShiny <- renderIgvShiny({
     printf("rendering igvShiny");
     igvShiny(list(genomeName="tair10",
                   initialLocus=sprintf("%s:%d-%d", chromosome, loc.start, loc.end)))
     })


   observeEvent(input$geneModel, ignoreInit = TRUE, {
      model.name <- isolate(input$geneModel);
      if(model.name != " "){
         printf("geneModel menu changed: %s", model.name)
         updateTabsetPanel(session, "trenaTabs", select="geneModelTab")
         }
      })

   observeEvent(input$removeOptionalTracks, ignoreInit = TRUE, {
      removeUserAddedTracks(session)
      })


   output$geneModelTable <- renderDT({
     printf("    rendering DT")
     selection="single"
         #tbl.model <- switch(model.name, "transcript.1.model"=tbl.model.1, "transcript.2.model"=tbl.model.2)
         #state$currentModel <- tbl.model
      printf("data.frame dimensions: %d, %d", nrow(tbl.model), ncol(tbl.model))
      if(nrow(tbl.model) > 10)
          tbl.model <- tbl.model[1:10,]
      rownames(tbl.model) <- NULL
         # later(function() {updateTabsetPanel(session, "trenaTabs", select="geneModelTab")}, 1);
      return(tbl.model)
      },
     selection = 'single',
     options=list(pageLength=25, dom="t")
     ) # renderDT

   observeEvent(input$geneModelTable_cell_clicked, {
      trueRowClick <- length(input$geneModelTable_cell_clicked) > 0
      printf("--- input$geneModelTable_cell_clicked: %s", trueRowClick)
      if(trueRowClick){
         selectedTableRow <- input$geneModelTable_row_last_clicked
         printf("selectedTableRow: %d", selectedTableRow)
         if(!is.null(selectedTableRow)){
            tf <- tbl.model$gene[selectedTableRow]
            printf("selected tf to plot: %s", tf)
            later(function(){selectRows(dt.proxy, NULL)}, 0.01);
            if(length(tf) > 0){
               updateTabsetPanel(session, "trenaTabs", select="plotTab");
               output$xyPlot = renderPlot(plotTfTargetGeneCorrelation(session, tf))
               } # if tf
            } # row not null
         } # if trueRowClick
      }) # geneModelTable_cell_clicked

} # server
#----------------------------------------------------------------------------------------------------
plotTfTargetGeneCorrelation <- function(session, tf)
{
   printf("want an xyplot of %s vs. FRD3 (%s)", tf, orf)
   plot(mtx[tf, ], mtx[orf, ],main="Gene Expression",
        xlab=tf, ylab=orf, col="red", pch=19) # , xlim=c(0,5), ylim=c(0,5))

} # plotTfTargetGeneCorrelation
#----------------------------------------------------------------------------------------------------
findRegionsAboveThreshold <- function(vec, threshold, minSpan)
{
   vec[vec < threshold] <- NA
   vec[vec >= threshold] <- 1
   vec.rle <- rle(vec)
   rle.segments.above.threshold <- which(vec.rle$lengths >= minSpan)
       # make all runs less than desired.minSpan are NA'd
   all.subMinimumSpans <- setdiff(1:length(vec.rle$values), rle.segments.above.threshold)
   vec.rle$values[all.subMinimumSpans] <- NA    # vec.rle$values[18] <- NA

   rle.region.count <- length(vec.rle$values)
   actual.index <- 1

   peak.starts <- vector("numeric", length(vec))
   peak.ends <- vector("numeric", length(vec))
   peak.count <- 0

   for(i in 1:rle.region.count){
      size <- vec.rle$length[i]
      value <- vec.rle$values[i]
      if(!is.na(value)){
         peak.count <- peak.count + 1
         region.start <- actual.index
         region.end   <- actual.index + size - 1
        #printf("peak found:  %d-%d", region.start, region.end)
         peak.starts[peak.count] <- region.start
         peak.ends[peak.count] <- region.end
         } # !is.na
    actual.index <- actual.index + size
    } # for i

   list(starts=peak.starts[seq_len(peak.count)], ends=peak.ends[seq_len(peak.count)])

} # findRegionsAboveThreshold
#------------------------------------------------------------------------------------------------------------------------
shinyOptions <- list()

if(Sys.info()[["nodename"]] == "trena.systemsbiology.net"){
   port <- 60031
   printf("running on trena, using port %d", port)
   shinyOptions <- list(host="0.0.0.0", port=port, launch.browser=FALSE)
   }
app <- shinyApp(ui=ui,server=server, options=shinyOptions)

