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
#----------------------------------------------------------------------------------------------------

targetGene <- "INPP5D"
load("data/tbl.model.RData")                      # "tbl.model"
load("data/mtx.withDimers.cer.ros.tcx.RData")     # "mtx.cer" "mtx.ros" "mtx.tcx"
load("data/tbl.enhancers.RData")                  # "tbl.enhancers"

genome.name <- "hg38"

loc.min <- min(tbl.enhancers$start) - 2000
loc.max <-  max(tbl.enhancers$end) + 2000
chromosome <- tbl.enhancers$chr[1]

load("data/tbl.gwas.igap.snp.RData")              # "tbl.snp"
tbl.snp <- subset(tbl.snp, chrom==chromosome & start >= loc.min & end <= loc.max)
load("data/tbl.enhancer.gwas.igap.snp.RData")     # tbl.enhancer.snps
load("data/tbl.breaks.RData")                     # tbl.breaks
load("data/tbl.breaksBed.RData")                  # tbl.breaksBed
colnames(tbl.breaksBed)[4] <- "id"                # used by igvShiny.  TODO - rationalize this! (13 sep 2018)
load("data/tbl.fimo.5tfs.bindingSites.RData")     # tbl.fimo
load("data/tbl.chipSeqInModel.RData")             # tbl.chipSeqInModel
load("data/tbl.fp.RData")
load("data/tbl.dhs.RData")                        # tbl.dhs
addResourcePath("tmp", "tmp") # so the shiny webserver can see, and return track files for igv

model.names <- "fp + enhancers"
models <- list("fp + enhancers" = tbl.model)
#----------------------------------------------------------------------------------------------------
#load("enhancers.RData")
#load("snps.RData")               # tbl.snp, tbl.eqtl
#load("tbl.bindingSitesWithSnpsAndDeltaScores.RData")
#load("tbl.fimoHitsInEnhancers.RData")
#----------------------------------------------------------------------------------------------------
chromosome <- "chr2"
tss <- 233059795
currentShoulder <- 0
currentPwmMatchThreshold <- 90

state <- new.env(parent=emptyenv())
state[["currentModel"]] <- models[[1]]
state[["currentModelName"]] <- names(models)[1]

selectedTableRow <- 0

#state[["optionalTracks"]] <- c()
#state[["currentGenomicRegion"]] <- ""
#----------------------------------------------------------------------------------------------------
ui <- fluidPage(

  includeScript("message-handler.js"),
  tags$script(JS('setInterval(function(){console.log("keepAlive"); Shiny.setInputValue("foo", "var");}, 1 * 60 * 1000)')),
  tags$head(tags$link(rel = "stylesheet", type = "text/css",
                     href = "http://maxcdn.bootstrapcdn.com/font-awesome/4.2.0/css/font-awesome.min.css")
     ),

  titlePanel(sprintf("%s Trena Model & Disruptions", targetGene)),

  sidebarLayout(
     sidebarPanel(
        width=3,
        selectInput("displayGenomicRegion", "Display Genomic Region:",
                    c("All INPP5D enhancers" = "enhancers",
                      "INPP5D +/- 40kb" = "tssPlusMinus40kb",
                      "TSS-proximal interesting snps"="interestingSnpsNearTSS")),

        #selectInput("geneModel", "Choose model:", models),

        selectInput("addTrack", "Add Track:",
                   c("",
                     "Enhancers"="enhancers",
                     "DHS open chromatin - encode clustered"="dhs",
                     "ChIP-seq interleaved with footprints (4 TFs)"="chipSeqFootprints",
                     "ChIP-seq"="chipSeq",
                     "rs35349669"= "rs35349669",
                     "IGAP GWAS SNPs"="snps",
                     #"Footprints"="footprints"
                     #"SNPs in Enhancers"="enhancer.snps",
                     "Motif-breaking SNPs"="breaking.snps"
                     )),

        sliderInput("IGAP.snp.significance", "SNP score",
                    value = 4,
                    step=0.1,
                    min = 0,
                    max = 12),
        sliderInput("fimo.motif.significance", "FIMO motif score",
                    value = 2, #fivenum(-log10(tbl.bs$motif.pVal))[3],
                    step=0.1,
                    min = 0,
                    max = 12,
                    round=-2),
        actionButton("findDisruptiveSNPsButton", "SNPs which disrupt TF binding sites"),
        # sliderInput("fimo.snp.effect", "SNP binding affinity loss score",
        #              value = 0.5,
        #              step=0.1,
        #              min = 0,
        #              max = 12),
        HTML("<br>"),
        sliderInput("snpShoulder", "Proximity",
                    value = 10,
                    min = 0,
                    max = 100),
        actionButton("showSNPsNearBindingSitesButton", "SNPs near binding sites in enhancer regions"),
        HTML("<br>"),
        #actionButton("addAllBindingSitesButton", "Add TF binding sites from regulatory model"),
        HTML("<br>"),
        actionButton("removeOptionalTracks", "Remove tracks")
     ),

     mainPanel(
        tabsetPanel(type="tabs",
                    id="trenaTabs",
                   tabPanel(title="IGV",             value="igvTab",          igvShinyOutput('igvShiny')),
                   tabPanel(title="READ ME",         includeHTML("readme.html")),
                   tabPanel(title="trena model",     value="geneModelTab",
                            mainPanel(radioButtons("tfSelectionChoice", "Row (transcription factor) selection will display:",
                                                   c("XY plot" = "xyPlot",
                                                     "Binding sites" = "displayBindingSites"),
                                                   inline=TRUE),
                                      DTOutput("geneModelTable")
                                      )),
                   tabPanel(title="SNP", value="snpInfoTab",
                            mainPanel(
                               verbatimTextOutput("tfbsText"),
                               DTOutput("snpMotifTable"),
                               imageOutput("logoImage"),
                               width=12)
                            ),
                   tabPanel(title="TF/INPP5D xy plot",  value="plotTab", plotOutput("xyPlot", height=800))
                   ) # tabsetPanel
      )  # mainPanel
    ) # sidebarLayout
  ) # fluidPage

#--------------------------------------------------------------------------------
server <- function(input, output, session) {

   session$sendCustomMessage(type="showGenomicRegion", message=(list(region="FLT1")))
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
     session$sendCustomMessage(type="showGenomicRegion", message=(list(region=roi)))
     })

   observeEvent(input$displayGenomicRegion, {
      regionName <- input$displayGenomicRegion
      #printf("display event: %s", regionName)
      state$currentGenomicRegion <- regionName
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
     igvShiny(list(genomeName=genome.name,
                   initialLocus="chr2:233,027,632-233,264,596"))
     })

  } # server

#----------------------------------------------------------------------------------------------------
displayTrack <- function(session, trackName)
{
   supportedTracks <- c("enhancers",
                        "dhs",
                        "chipSeqFootprints",
                        "chipSeq",
                        "rs35349669",
                        "snps",
                        "footprints",
                        "enhancer.snps",
                        "breaking.snps")

   printf("--- displayTrack('%s')", trackName)

   if(!trackName %in% supportedTracks)
      return()

   if(trackName == "chipSeq"){
      displayChipSeqByTF(session)
      return()
      }

   if(trackName == "chipSeqFootprints"){
      displayInterleavedCHiPseqAndFootprints(session)
      return()
      }

   if(trackName == "footprints"){
      displayFootprintsByTF(session)
      return()
      }

   x <- switch(trackName,
      enhancers       = {list(bed=tbl.enhancers,        name="INPP5D enhancers",   color="black")},
      dhs             = {list(bedgraph=tbl.dhs,         name="DHS",                color="blue")},
      snps            = {list(bedgraph=tbl.snp,         name="GWAS snp",           color="darkRed")},
      enhancer.snps   = {list(bed=tbl.enhancer.snps,    name="Enhancer snp",       color="darkRed")},
      breaking.snps   = {list(bed=tbl.breaksBed,        name="Motif-breaking snp", color="maroon")},
      rs35349669      = {list(bed=data.frame(chr="chr2",
                                             start=233159830,
                                             end=233159830,
                                             name="rs35349669",
                                             stringsAsFactors=FALSE),
                              name="rs35349669",
                              color="purple")}
      )

    updateTabsetPanel(session, "trenaTabs", select="igvTab");

    if(trackName == "snps"){
       tbl.tmp <- x$bedgraph[, c("chrom", "start", "end", "pScore")]
       loadBedGraphTrack(session, "GWAS snp", tbl.tmp, color=x$color, trackHeight=25, autoscale=FALSE, min=0, max=10)
       } # bedGraph for snps

    if(trackName == "dhs"){
       loadBedGraphTrack(session, "DHS", x$bedgraph, color=x$color, trackHeight=25, autoscale=TRUE)
       } # bedGraph for snps

    if(trackName == "enhancers"){
       tbl.tmp <- tbl.enhancers[, c("chr", "start", "end", "score")]
       loadBedGraphTrack(session, "GeneHancer", tbl.tmp, color="black", trackHeight=25, autoscale=FALSE, min=0, max=20)
       } # bedGraph for snps

   if(trackName %in% c("enhancer.snps", "breaking.snps", "chipSeq", "rs35349669", "spi1.footprints")){
       loadBedTrack(session, trackName, x$bed, color=x$color, trackHeight=30)
       } # bed format track

} # displayTrack
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
displayInterleavedCHiPseqAndFootprints <- function(session)
{
   colors <- brewer.pal(8, "Dark2")
   totalColorCount <- length(colors)
   colorNumber <- 0

   for(gene in unique(tbl.csInModel$tf)){
      tbl.broad <- subset(tbl.csInModel, tf==gene)[, c("chr", "start", "end")]
      colorNumber <- (colorNumber %% totalColorCount) + 1
      color <- colors[colorNumber]
      loadBedTrack(session, sprintf("%s ChIP", gene), tbl.broad, color=color, trackHeight=25)
      tbl.bed <- subset(tbl.fp, tf==gene)[, c("chr", "start", "end")]
      loadBedTrack(session, sprintf("%s fp", gene), tbl.bed, color=color, trackHeight=25)
      } # for tf


} # displayInterleavedCHiPseqAndFootprints
#------------------------------------------------------------------------------------------------------------------------
displayFootprintsByTF <- function(session)
{
   colors <- brewer.pal(8, "Dark2")
   totalColorCount <- length(colors)
   colorNumber <- 0

   for(gene in unique(tbl.fp$tf)){
      tbl.bed <- subset(tbl.fp, tf==gene)[, c("chr", "start", "end")]
      colorNumber <- (colorNumber %% totalColorCount) + 1
      color <- colors[colorNumber]
      loadBedTrack(session, sprintf("%s fp", gene), tbl.bed, color=color, trackHeight=25)
      } # for tf

} # displayFootprintsByTF
#------------------------------------------------------------------------------------------------------------------------
findDisruptiveSNPs <- function(snpSignificanceThreshold, motifMatchThreshold, bindingLossThreshold)
{
   printf("--- findDisruptiveSNPs")
   printf("  snpSignificanceThreshold: %5.1f", snpSignificanceThreshold)
   printf("  motifMatchThreshold: %5.1f",      motifMatchThreshold)
   printf("  bindingLossThreshold: %5.1f",     bindingLossThreshold);
   browser()
   tbl.tmp <- subset(tbl.bs,  -log10(gwasPval) >= snpSignificanceThreshold |
                              -log10(eqtlPval) >= snpSignificanceThreshold)
   tbl.out <- subset(tbl.tmp,
                     -log10(motif.pVal )   >= motifMatchThreshold &
                     wtMutMatchDelta       >= bindingLossThreshold)
   printf("snps meet criteria: %d", nrow(tbl.out))
   tbl.out

} # findDisruptiveSNPs
#------------------------------------------------------------------------------------------------------------------------
test_findDisruptiveSNPs <- function()
{
   printf("--- test_findDisruptiveSNPs")
   checkTrue(exists("tbl.bs"))
   tbl.dsnps <- findDisruptiveSNPs(snpSignificanceThreshold=0, motifMatchThreshold=3, bindingLossThreshold=1)
   dim(tbl.dsnps)
   coi <- c("motif", "chrom", "start", "stop", "p.value", "rsid", "start.1", "wt", "mut.x", "snpPval.x", "wtMutMatchDelta")
   tbl.dsnps[, coi]

} # test_findDisruptiveSNPs
#------------------------------------------------------------------------------------------------------------------------
displayDisruptiveSnps <- function(session, tbl.disruptions)
{
   printf("--- displayDisruptiveSnps");

   coi <- c("chrom", "start", "end", "motifName", "tf", "motifScore", "rsid")
   tbl.condensed <- tbl.disruptions[, coi]
   tfs <- unique(tbl.condensed$tf)
   for(geneSymbol in tfs){
      tbl.tmp <- subset(tbl.condensed, tf==geneSymbol)
      tbl.tmp$motifName <- sprintf("tfbs-snp:%s:%s", tbl.tmp$motifName, tbl.tmp$rsid)
      temp.filename <- tempfile(tmpdir="./tmp", fileext=".bed")
      write.table(tbl.tmp, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
      session$sendCustomMessage(type="displayBedTrack",
                             message=(list(filename=temp.filename,
                                           trackName=sprintf("%s w/snp", geneSymbol),
                                           displayMode="EXPANDED",
                                           color="red",
                                           trackHeight=40)))
      } # for geneSymbol

} # displayDisruptiveSnps
#------------------------------------------------------------------------------------------------------------------------
calculateVisibleSNPs <- function(targetGene, roi, filters, snpShoulder)
{
   #return(sprintf("calculateVisibleSNPs(%s), %s, %s, %d", date(), targetGene, roi, snpShoulder))
   tbl.mayo.snps <- mef2c@misc.data$MAYO.eqtl.snps[, c("chrom", "start", "end", "score")]
   tbl.mayo.snps <- tbl.mayo.snps[order(tbl.mayo.snps$start, decreasing=FALSE),]

   tbl.igap.snps <- mef2c@misc.data[["IGAP.snpChip"]][, c("chrom", "start", "end", "score")]
   tbl.igap.snps <- tbl.igap.snps[order(tbl.igap.snps$start, decreasing=FALSE),]

   temp.filename <- sprintf("igvdata/tmp%d.bedGraph", as.integer(Sys.time()))
   write.table(tbl.igap.snps, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
   return(temp.filename)

} # calculateVisibleSNPs
#----------------------------------------------------------------------------------------------------
removeTrack <- function(session, trackName)
{
  session$sendCustomMessage(type="removeTrack", message=list(trackName=trackName))

} # removeTrack
#----------------------------------------------------------------------------------------------------
displayBindingSites <- function(session, target.tf)
{
   xyz <- "--- displayBindingSites, consult tbl.fimo..."
   #browser()
   if(!target.tf %in% tbl.fimo$tf){
      printf("tf %s not found in tbl.fimo", target.tf, state$currentModelName)
      return()
      }
   motif.pVal.threshold <- isolate(session$input$fimo.motif.significance)
   tbl.tfbs <- subset(tbl.fimo, tf==target.tf & -log10(motif.pVal) >= motif.pVal.threshold)
   tbl.tfbs <- tbl.tfbs[, c("chrom", "motifStart", "motifEnd", "motif.pScore")]
   temp.filename <- tempfile(tmpdir="tmp", fileext=".bed")
   write.table(tbl.tfbs, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
   updateTabsetPanel(session, "trenaTabs", select="igvTab");
   later(function(){
       session$sendCustomMessage(type="displayBedTrack",
                                   message=(list(filename=temp.filename,
                                                 trackName=target.tf,
                                                 color="darkGreen",
                                                 trackHeight=40)))
       }, 1)

} # displayBindingSites
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
addAllBindingSites <- function(session)
{
   printf("addBindingSites - disabled for now")

} # addAllTFbindingSites
#------------------------------------------------------------------------------------------------------------------------

shinyOptions <- list()

if(Sys.info()[["nodename"]] == "trena.systemsbiology.net"){
   port <- 60031
   printf("running on trena, using port %d", port)
   shinyOptions <- list(host="0.0.0.0", port=port, launch.browser=FALSE)
   }
app <- shinyApp(ui=ui,server=server, options=shinyOptions)

