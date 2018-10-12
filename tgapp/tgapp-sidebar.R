library(shiny)
library(DT)
library(igvShiny)
library(htmlwidgets)
library(GenomicRanges)
library(later)
library(MotifDb)
library(RUnit)
library(RColorBrewer)
library(org.Hs.eg.db)
library(shiny)
library(shinydashboard)
library(igvShiny)
library(DT)

library(trenaSGM)
library(TrenaGene)
library(TrenaGenePlacentaData)
#------------------------------------------------------------------------------------------------------------------------
targetGene <- list(hugo="PSG1", ensembl="ENSG00000231924", combined="ENSG00000231924|PSG1")
genomeName <- "hg38"
dataDirectory <- system.file(package="TrenaGenePlacentaData", "extdata")

trenaGene <- TrenaGene(targetGene$hugo, genomeName, expressionDataDirectory=dataDirectory)
expression.data.set.names <- getExpressionMatrixNames(trenaGene)
tbl.enhancers <- getEnhancers(trenaGene)
tbl.dhs <- getEncodeDHS(trenaGene)
#------------------------------------------------------------------------------------------------------------------------
combineGeneNames <- function(ensemblIDs)
{
   tbl.geneNameMapper <- select(org.Hs.eg.db, keys=ensemblIDs, keytype="ENSEMBL", columns="SYMBOL")
   tbl.geneNameMapper$combined <- with(tbl.geneNameMapper, paste(ENSEMBL, SYMBOL, sep="|"))
   dim(tbl.geneNameMapper); #dim(mtx)
     # [1] 16831     3
     # [1] 16664   112
  dups <- which(duplicated(tbl.geneNameMapper$ENSEMBL)) # 167
  if(length(dups) > 0)
     tbl.geneNameMapper <- tbl.geneNameMapper[-dups,]
   return(tbl.geneNameMapper)

} # combineGeneNames
#------------------------------------------------------------------------------------------------------------------------
totalColorCount <- 6
colors <- brewer.pal(totalColorCount, "Dark2")
colorNumber <- 0
#----------------------------------------------------------------------------------------------------
ui <- fluidPage(

  tags$script(JS('setInterval(function(){console.log("keepAlive"); Shiny.setInputValue("foo", "var");}, 1 * 60 * 1000)')),
  titlePanel(sprintf("%s Trena Model & Disruptions", targetGene$hugo)),

  sidebarLayout(
     sidebarPanel(
        width=2,
        selectInput("displayGenomicRegion", "Display Genomic Region:",
                    c("",
                      "Tradtional Promoter" = "traditionalPromoter",
                      "Enhancers Region" = "enhancersRegion")),
        selectInput("addTrack", "Add Track:",
                   c("",
                     "Enhancers"="enhancers",
                     "DHS open chromatin - encode clustered"="dhs"
                     )),
        HTML("<br>"),
        actionButton("removeOptionalTracks", "Remove tracks")
        ),

     mainPanel(width=10,
        tabsetPanel(type="tabs",
                    id="trenaTabs",
                    tabPanel(title="IGV", value="igvTab",
                             fluidRow(
                                column(8, igvShinyOutput('igvShiny')),
                                column(4, sliderInput("obs", "Number of observations:", min = 1, max = 1000, value = 500))
                                )),
                    tabPanel(title="Create new trena model",  value="geneModelTab",
                            mainPanel(radioButtons("tfSelectionChoice", "Row (transcription factor) selection will display:",
                                                   c("XY plot" = "xyPlot",
                                                     "Binding sites" = "displayBindingSites"),
                                                   inline=TRUE),
                                      DTOutput("geneModelTable")
                                      ))
                   ) # tabsetPanel
      )  # mainPanel
    ) # sidebarLayout
  ) # fluidPage

#--------------------------------------------------------------------------------
server <- function(input, output, session) {


   dt.proxy <- dataTableProxy("geneModelTable")

   observeEvent(input$trackClick, {
      printf("browser snp click observed!")
      print(input$trackClick)
      motifName <- input$trackClick$id
      if(grepl("^tfbs-snp", motifName)){
         tokens <- strsplit(motifName, ":")[[1]]
         #pfm.name <- tokens[2]
         snp.name <- tokens[2]
         tbl.sub <- tbl.breaks[snp.name, c("seqMatch", "motifPos", "geneSymbol", "pctRef", "pctAlt")]
         #tbl.sub <- tbl.sub[, c("chrom", "start", "end", "strand", "sequence", "sequence.mut", "wtMutMatchDelta")]
         #colnames(tbl.sub)[7] <- "10^x match loss"
         #output$tfbsText <- renderText({sprintf("%s   %s", pfm.name, snp.name)})
         image.filename <- sprintf("png/%s.png", snp.name)
         printf("--- image file exists? %s, %s", image.filename, file.exists(image.filename))
         later(function(){
           output$snpMotifTable <- renderDT(tbl.sub, width="800px", options = list(scrollX = TRUE, dom='t'))
           output$logoImage <- renderImage({
                 list(src=image.filename,
                      contentType="image/png",
                      width=800,
                      height=800,
                      alt="png file not found")},
                 deleteFile=FALSE)
           }, 1)
         updateTabsetPanel(session, "trenaTabs", select="snpInfoTab");
         } # if tfbs-snp
      })

   observeEvent(input$removeOptionalTracks, ignoreInit = TRUE, {
      removeUserAddedTracks(session)
      })

   observeEvent(input$displayGenomicRegion, ignoreInit=TRUE,  {
      regionName <- isolate(input$displayGenomicRegion)
      printf("--- displayGenomicRegion change observed: '%s'", regionName)
      printf("    length: %d", nchar(regionName))
      if(nchar(regionName) > 0){
         updateTabsetPanel(session, "trenaTabs", select="igvTab");
         myFunc <- function(){
            showGenomicRegion(session, mapToChromLoc(regionName))
            } # myFunc
         later(myFunc, 0.2)
         later(function() {updateSelectInput(session, "displayGenomicRegion", selected=character(0))}, 1)
         } # if
      })

   observeEvent(input$addTrack, {
      newTrackName <- input$addTrack
      if(newTrackName == "") return()
      printf(" addTrack event: %s", newTrackName);
      displayTrack(session, newTrackName)
      later(function() {updateSelectInput(session, "addTrack", selected=character(0))}, 1)
      })

   output$igvShiny <- renderIgvShiny({
     igvShiny(list(genomeName=genomeName,
                   initialLocus=targetGene$hugo))
     })

  } # server

#----------------------------------------------------------------------------------------------------
displayTrack <- function(session, trackName)
{
   supportedTracks <- c("enhancers",
                        "dhs")

   printf("--- displayTrack('%s')", trackName)

   switch(trackName,
          enhancers = {displayEnhancersTrack(session)},
          dhs       = {displayEncodeDhsTrack(session)}
          )

    updateTabsetPanel(session, "trenaTabs", select="igvTab");


} # displayTrack
#------------------------------------------------------------------------------------------------------------------------
displayEnhancersTrack <- function(session)
{
   tbl.tmp <- tbl.enhancers[, c("chrom", "start", "end", "combinedScore")]
   loadBedGraphTrack(session, "GeneHancer", tbl.tmp, color="black", trackHeight=25, autoscale=FALSE, min=0, max=20)

} # displayEnhancersTrack
#------------------------------------------------------------------------------------------------------------------------
displayEncodeDhsTrack <- function(session)
{
   tbl.tmp <- tbl.dhs[, c("chrom", "chromStart", "chromEnd", "score")]
   colnames(tbl.tmp) <- c("chrom", "start", "end", "score")
   loadBedGraphTrack(session, "DHS", tbl.tmp, color="black", trackHeight=25, autoscale=TRUE)

} # displayEncodeDhsTrack
#------------------------------------------------------------------------------------------------------------------------
mapToChromLoc <- function(regionName)
{
   roi <- targetGene$hugo   # a safe fallback

   if(regionName == "traditionalPromoter"){
      tbl.transcripts <- getTranscriptsTable(trenaGene)
      chrom <- subset(tbl.transcripts, moleculetype=="gene")$chr
      strand <- subset(tbl.transcripts, moleculetype=="gene")$strand
      gene.start <- subset(tbl.transcripts, moleculetype=="gene")$start
      gene.end   <- subset(tbl.transcripts, moleculetype=="gene")$endpos
      upstream <- 5000
      downstream <- 5000
      if(strand == "+"){
         tss <- gene.start
         start <- tss - upstream
         end <- tss + downstream
         }
      if(strand == "-"){
         tss <- gene.end
         start <- tss - downstream
         end <- tss + upstream
         }
      roi <- sprintf("%s:%d-%d", chrom, start, end)
      }

   if(regionName == "enhancersRegion"){
      chrom <- tbl.enhancers$chrom[1]
      start <- min(tbl.enhancers$start) - 10000
      end   <- max(tbl.enhancers$end) + 10000
      roi <- sprintf("%s:%d-%d", chrom, start, end)
      }

   return(roi)

} # mapToChromLoc
#------------------------------------------------------------------------------------------------------------------------
buildFootprintModel <- function(upstream, downstream)
{
   tbl.gene <- subset(getTranscriptsTable(trenaGene), moleculetype=="gene")[1,]
   tss <- tbl.gene$start
   min.loc <- tss - upstream
   max.loc <- tss + downstream

   if(tbl.gene$strand == "-"){
      tss <- tbl.gene$endpos
      min.loc <- tss - downstream
      max.loc <- tss + upstream
      }

   chrom <-tbl.gene$chr
   tbl.regions <- data.frame(chrom=chrom, start=min.loc, end=max.loc, stringsAsFactors=FALSE)
   mtx <- loadExpressionData(trenaGene, "FilteredLengthScaledTPM8282018-vsn")

   build.spec <- list(title=sprintf("%s model %d", targetGene$hugo, 1),
                      type="footprint.database",
                      regions=tbl.regions,
                      geneSymbol=targetGene$hugo,
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20"),
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(identifierType="geneSymbol"),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.1,
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      orderModelByColumn="pearsonCoeff",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

     #------------------------------------------------------------
     # use the above build.spec: a small region, high correlation
     # required, MotifDb for motif/tf lookup
     #------------------------------------------------------------

   fpBuilder <- FootprintDatabaseModelBuilder(genomeName, targetGene$hugo, build.spec, quiet=TRUE)
   x <- build(fpBuilder)

} # buildFootprintModel
#------------------------------------------------------------------------------------------------------------------------
shinyOptions <- list()

if(Sys.info()[["nodename"]] == "trena.systemsbiology.net"){
   port <- 60031
   printf("running on trena, using port %d", port)
   shinyOptions <- list(host="0.0.0.0", port=port, launch.browser=FALSE)
   }
app <- shinyApp(ui=ui,server=server, options=shinyOptions)

