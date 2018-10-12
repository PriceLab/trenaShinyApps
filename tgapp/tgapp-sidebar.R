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
  #dim(tbl.geneNameMapper); dim(mtx)
   #rownames(mtx) <- tbl.geneNameMapper$combined
   return(tbl.geneNameMapper)

} # combineGeneNames
#------------------------------------------------------------------------------------------------------------------------
targetGene <- list(hugo="PSG1", ensembl="ENSG00000231924", combined="ENSG00000231924|PSG1")
genomeName <- "hg38"

dataDirectory <- system.file(package="TrenaGenePlacentaData", "extdata")

trenaGene <- TrenaGene(targetGene$hugo, genomeName, expressionDataDirectory=dataDirectory)
expression.data.set.names <- getExpressionMatrixNames(trenaGene)
tbl.enhancers <- getEnhancers(trenaGene)
tbl.dhs <- getEncodeDHS(trenaGene)


#----------------------------------------------------------------------------------------------------
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
                     "DHS open chromatin - encode clustered"="dhs",
                     "ChIP-seq"="chipSeq"
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

   #session$sendCustomMessage(type="showGenomicRegion", message=(list(region="FLT1")))
   dt.proxy <- dataTableProxy("geneModelTable")
   #later(function() {
   #   printf("about to select gene family");
   #   updateSelectInput(session, "displayGenomicRegion", selected="overview")
   #   }, 3)

    # define a reactive conductor. it returns a function, is
    # redefined with a change to dist or n
    # has reactive children  plot, summary, table, below, which mention it
    # thereby establishing the reactive chain.


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

   observeEvent(input$addAllBindingSitesButton, ignoreInit=TRUE,{
      addAllBindingSites(session)
     })

   observeEvent(input$removeOptionalTracks, ignoreInit = TRUE, {
      removeUserAddedTracks(session)
      })

   observeEvent(input$geneModel, {
      printf("geneModel event: %s", input$geneModel);
      updateTabsetPanel(session, "trenaTabs", select="geneModelTab");
      })

   observeEvent(input$showTargetGeneButton, {
     roi <- "chr6:41,152,337-41,213,272"
     # session$sendCustomMessage(type="showGenomicRegion", message=(list(region=roi)))
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

   observeEvent(input$findDisruptiveSNPsButton, {
      snpSignificanceThreshold <- isolate(session$input$IGAP.snp.significance)
      motifMatchThreshold <- isolate(session$input$fimo.motif.significance)
      bindingLossThreshold <- snpEffect.pval.threshold <- isolate(session$input$fimo.snp.effect)
      tbl.disruptions <- findDisruptiveSNPs(snpSignificanceThreshold, motifMatchThreshold, bindingLossThreshold)
      if(nrow(tbl.disruptions) > 0)
         displayDisruptiveSnps(session, tbl.disruptions)
      })

   observeEvent(input$findSnpsInModelButton, {
      trena.model <- isolate(input$model.trn)
      filters <- isolate(input$filters)
      snpShoulder <- isolate(input$snpShoulder)
      filter.result <- applyFilters(session, filters, trena.model, snpShoulder)
      })

   observeEvent(input$pwmMatchThreshold, {
      printf("  ** pwmMatchThreshold changed");
      currentPwmMatchThreshold <- input$pwmMatchThreshold
      })

   observeEvent(input$filters, {
      currentFilters <- input$filters
      print(paste(currentFilters, collapse=","))
      })

   observeEvent(input$model.trn, {
      currentModel <- input$model.trn
      print(paste(currentFilters, collapse=","))
      })

   observeEvent(input$roi, {
     currentValue = input$roi
     tss <- 88884466;   # get from MEF2C.data
     roi.string <- switch(currentValue,
                          tss_2kb = sprintf("chr5:%d-%d", tss-2000, tss+2000),
                          tss_5kb = sprintf("chr5:%d-%d", tss-5000, tss+5000),
                          studyRegion = "chr5:88391000-89322000")
     session$sendCustomMessage(type="roi", message=(list(roi=roi.string)))
     })

   observeEvent(input$showSNPsNearBindingSitesButton, {
      displaySNPsNearBindingSites(session)
      })

   observeEvent(input$geneModelTable_rows_selected, {
         selectedTableRow <- isolate(input$geneModelTable_rows_selected)
         tbl.model <- state$currentModel
         tf <- tbl.model$tf.hgnc[selectedTableRow]
         printf("%d %s", selectedTableRow, tf)
         tfSelectionAction <- isolate(input$tfSelectionChoice)
         printf("table row clicked: %s  - actionRequested: %s", tf, tfSelectionAction);
         if(tfSelectionAction == "xyPlot"){
            updateTabsetPanel(session, "trenaTabs", select="plotTab");
            #tokens <- strsplit(state$currentModelName, "\\.")[[1]]
            #expression.matrix.id <- tokens[length(tokens)]
            expression.matrix.id <- "tcx"
            output$xyPlot = renderPlot(plotTfTargetGeneCorrelation(session, tf, expression.matrix.id))
             }
         if(tfSelectionAction == "displayBindingSites"){
            displayBindingSites(session, tf)
            }
         }) # observe row selection event

   output$geneModelTable <- renderDT(
      {model.name <- input$geneModel;
       updateTabsetPanel(session, "trenaTabs", select="geneModelTab");
       printf("    rendering DT, presumably because input$geneModel changes, model.name: %s", model.name);
       selection="single"
       #tbl.model <- tbl.model[, c("tf.hgncall.models[[model.name]]$model
       printf("data.frame dimensions: %d, %d", nrow(tbl.model), ncol(tbl.model))
       #tbl.model <- switch(model.name,
       #   "Enhancers.TFClass" = tbl.model.enhancerAll.tfc,
       #   "Enhancers.MotifDb" = tbl.model.enhancerAll.mdb
       #   )
      #tbl.model <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE),]
      if(nrow(tbl.model) > 20)
         tbl.model <- tbl.model[1:20,]
      rownames(tbl.model) <- NULL
      state[["currentModel"]] <- tbl.model
      state[["currentModelName"]] <- model.name
      return(tbl.model)
      },
    selection = 'single',
    options=list(pageLength=25, dom="t"))

   output$igvShiny <- renderIgvShiny({
     igvShiny(list(genomeName=genomeName,
                   initialLocus=targetGene$hugo))
     })

  } # server

#----------------------------------------------------------------------------------------------------
displayTrack <- function(session, trackName)
{
   supportedTracks <- c("enhancers",
                        "dhs",
                        "chipSeq")

   printf("--- displayTrack('%s')", trackName)

   switch(trackName,
          enhancers = {printf("create enhancers track");
                       displayEnhancersTrack(session)},
          dhs       = {printf("create dhs track");
                       displayEncodeDhsTrack(session)},
                       # c("chrom", "chromStart", "chromEnd", "count", "score"))

          chipSeq   = {printf("create ChIP-seq track")}
          )

    updateTabsetPanel(session, "trenaTabs", select="igvTab");

#    if(trackName == "snps"){
#       tbl.tmp <- x$bedgraph[, c("chrom", "start", "end", "pScore")]
#       loadBedGraphTrack(session, "GWAS snp", tbl.tmp, color=x$color, trackHeight=25, autoscale=FALSE, min=0, max=10)
#       } # bedGraph for snps
#
#    if(trackName == "dhs"){
#       loadBedGraphTrack(session, "DHS", x$bedgraph, color=x$color, trackHeight=25, autoscale=TRUE)
#       } # bedGraph for snps
#
#
#   if(trackName %in% c("enhancer.snps", "breaking.snps", "chipSeq", "rs35349669", "spi1.footprints")){
#       loadBedTrack(session, trackName, x$bed, color=x$color, trackHeight=30)
#       } # bed format track

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
displayChipSeqByTF <- function(session)
{
   colors <- brewer.pal(8, "Dark2")
   totalColorCount <- length(colors)
   colorNumber <- 0

   for(gene in unique(tbl.csInModel$tf)){
      tbl.broad <- subset(tbl.csInModel, tf==gene)[, c("chr", "start", "end")]
      tbl.peak <- subset(tbl.csInModel, tf==gene)[, c("chr", "peakStart", "peakEnd")]
      colnames(tbl.peak) <- c("chr", "start", "end")
      colorNumber <- (colorNumber %% totalColorCount) + 1
      color <- colors[colorNumber]
      loadBedTrack(session, sprintf("%s ChIP", gene), tbl.broad, color=color, trackHeight=25)
      # loadBedTrack(session, sprintf("%s peak", gene), tbl.peak, color=color, trackHeight=25)
      } # for tf

} # displayChipSeqByTF
#------------------------------------------------------------------------------------------------------------------------
# tbl.pheno$Diagnosis distribution in temporal cortex data
#   AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_Covariates_Flowcell_1.csv
#
# cerebellum phenotypes not yet found  (11 jun 2018)
#
#                AD   84
#           Control   80
#  Pathologic Aging   30
#               PSP   84
#
plotTfTargetGeneCorrelation <- function(session, tf, expression.matrix.id)
{
   stopifnot(expression.matrix.id %in% c("cer", "tcx"))
   mtx <- switch(expression.matrix.id,
                 cer = mtx.cer,
                 tcx = mtx.tcx
                 )
   printf("want an xyplot of %s vs. INPP5D", tf)
   covariates.file <- "data/AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_Covariates_Flowcell_1.csv"
   tbl.pheno <- read.table(covariates.file, sep=",", as.is=TRUE, header=TRUE)

   if(expression.matrix.id == "tcx"){ # we have phenotype data
      samples.AD <- intersect(subset(tbl.pheno, Diagnosis=="AD")$ID, colnames(mtx))
      samples.ctl <- intersect(subset(tbl.pheno, Diagnosis=="Control")$ID, colnames(mtx))
      samples.pathAging <- intersect(subset(tbl.pheno, Diagnosis=="Pathologic Aging")$ID, colnames(mtx))
      samples.PSP  <- intersect(subset(tbl.pheno, Diagnosis=="PSP")$ID, colnames(mtx))
      samples.all <- colnames(mtx)

      color.choices <- list(AD="red",
                            pathAging="green",
                            PSP="blue",
                            control="black"
                            )

      colors <- rep("black", length(samples.all))
      colors[match(samples.AD, colnames(mtx))] <- color.choices$AD
      colors[match(samples.ctl, colnames(mtx))] <- color.choices$control
      colors[match(samples.pathAging, colnames(mtx))] <- color.choices$pathAging
      colors[match(samples.PSP, colnames(mtx))] <- color.choices$PSP

      plot(mtx[tf, samples.all], mtx[targetGene, samples.all],main="Gene Expression - Temporal Cortex",
           xlab=tf, ylab=targetGene, col=colors, pch=19, xlim=c(-3,3), ylim=c(-3,3))

      tbl.AD <- as.data.frame(t(mtx[c(tf, targetGene), samples.AD]))
      model.AD <- lm(sprintf("%s ~ %s", targetGene, tf), data=tbl.AD)
      abline(model.AD, col=color.choices$AD)

      tbl.control <- as.data.frame(t(mtx[c(tf, targetGene), samples.ctl]))
      model.control <- lm(sprintf("%s ~ %s", targetGene, tf), data=tbl.control)
      abline(model.control, col=color.choices$control)

      #browser()

      tbl.pathAging <- as.data.frame(t(mtx[c(tf, targetGene), samples.pathAging]))
      model.pathAging <- lm(sprintf("%s ~ %s", targetGene, tf), data=tbl.pathAging)
      abline(model.pathAging, col=color.choices$pathAging)

      tbl.PSP <- as.data.frame(t(mtx[c(tf, targetGene), samples.PSP]))
      model.PSP <- lm(sprintf("%s ~ %s", targetGene, tf), data=tbl.PSP)
      abline(model.PSP, col=color.choices$PSP)

      legend.names <- names(color.choices)
      for(i in seq_len(length(legend.names))){
         condition <- legend.names[i]
         rsq <- eval(parse(text=sprintf("summary(model.%s)$r.squared", condition)))
         legend.names [i] <- sprintf("%s (rSq=%5.3f)", condition, rsq)
         }
      legend (1.5, -1.8, legend.names, as.character(color.choices))
      }
   else{  # no phenotype data, use all samples
      tf.vals <- mtx[tf, ]
      targetGene.vals <- mtx[targetGene,]
      plot(tf.vals, targetGene.vals, main="Gene Expression - Cerebellum", xlab=tf, ylab=targetGene)
      }


} # plotTfTargetGeneCorrelation
#------------------------------------------------------------------------------------------------------------------------
displaySNPsNearBindingSites <- function(session)
{
   shoulder <- isolate(session$input$snpShoulder)
   snpSignificanceThreshold <- isolate(session$input$IGAP.snp.significance)
   motifMatchThreshold <- isolate(session$input$fimo.motif.significance)

   printf("--- displaySNPsNearBindingSites")
   printf("shoulder: %d", shoulder)
   printf("snpSig: %f", snpSignificanceThreshold)
   printf("motif match: %f", motifMatchThreshold)

   tbl.snpFiltered <- subset(tbl.snp, pScore >= snpSignificanceThreshold)
   if(nrow(tbl.snpFiltered) == 0){
      printf("no snps above pScore threshold")
      return()
      }

   tbl.tfbs <- subset(tbl.fimo, motif.pScore >= motifMatchThreshold)
   if(nrow(tbl.tfbs) == 0){
      printf("no binding sites above fimoe threshold")
      return()
      }

   tbl.tfbs <- tbl.tfbs[, c("chrom", "motifStart", "motifEnd", "motif.pScore", "motifName")]
   colnames(tbl.tfbs) <- c("chrom", "start", "end", "score", "motifName")
   tbl.snpsNearEnough <- data.frame()
   gr.snps <- GRanges(tbl.snp)
   gr.tfbs <- GRanges(data.frame(chrom=tbl.tfbs$chrom,
                                 start=tbl.tfbs$start - shoulder,
                                 end=tbl.tfbs$end+shoulder))
   tbl.ov <- as.data.frame(findOverlaps(gr.snps, gr.tfbs))
   if(nrow(tbl.ov) == 0){
      printf("no snps near binding sites")
      return()
      }

   colnames(tbl.ov) <- c("snp", "tfbs")
   tbl.snpsNearEnough <- tbl.snp[tbl.ov$snp, c("chrom", "start", "end", "rsid")]
   motifNames <- tbl.tfbs$motifName[tbl.ov$tfbs]
   tbl.snpsNearEnough$motifName <- motifNames
   rownames(tbl.snpsNearEnough) <- NULL
   tbl.snpsNearEnough <- unique(tbl.snpsNearEnough)
   tbl.snpsNearEnough$name <- with(tbl.snpsNearEnough, paste(rsid, motifName, sep=":"))
   tbl.snpsNearEnough <- tbl.snpsNearEnough[, c("chrom", "start", "end", "name")]

   loadBedTrack(session, "snps near binding site", tbl.snpsNearEnough, color="brown", trackHeight=25)

} # displaySNPsNearBindingSites
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

