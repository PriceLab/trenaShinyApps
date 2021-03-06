library(shiny)
library(DT)
library(igvShiny)
library(htmlwidgets)
library(GenomicRanges)
library(later)
library(MotifDb)
library(motifStack)
library(RUnit)
library(RColorBrewer)
#----------------------------------------------------------------------------------------------------
genome <- "hg38"
targetGene <- "CKAP2L"
chromosome <- "chr2"
tss <- 112764686
snp.rsid <- "rs7594852"
snp.bp <- 112764177
snp.loc <- sprintf("%s:%d-%d", chromosome, snp.bp, snp.bp)

colors <- brewer.pal(12, "Paired")
colors <- rep(colors, 100)
colorNumber <- 0
#----------------------------------------------------------------------------------------------------
load("../data/models.RData")
addResourcePath("tmp", "tmp") # so the shiny webserver can see, and return track files for igv

model.names <- names(models)
modelList <- list()
for(model.name in model.names){
   modelList <- c(modelList, eval(parse(text=sprintf('c("%s"="%s")', model.name, model.name))))
   }

load("../data/matrices.RData")
load("../data/tbl.pheno.RData")
# load("../data/tbl.fimo.55tfs.29406.bindingSites.RData")
load("../data/tbl.fimo.20tfs.1115.bindingSites.RData")
load("../data/tbl.dhs.RData")
load("../data/tbl.snp.RData")
#----------------------------------------------------------------------------------------------------
load("../data/enhancers.RData")
#load("snps.RData")               # tbl.snp, tbl.eqtl
#load("tbl.bindingSitesWithSnpsAndDeltaScores.RData")
#load("tbl.fimoHitsInEnhancers.RData")
#----------------------------------------------------------------------------------------------------

tss <- 112764686
snp <- "rs7594852"
snp.bp <- 112764177
snp.loc <- sprintf("%s:%d-%d", chromosome, snp.bp, snp.bp)

currentFilters <- list()
currentShoulder <- 0
currentPwmMatchThreshold <- 90
filter.result <- list(fp=data.frame(), snp=data.frame(), fpClean=data.frame())

state <- new.env(parent=emptyenv())
state[["currentModel"]] <- models[[1]]
state[["currentModelName"]] <- names(models)[1]
state[["optionalTracks"]] <- c()
state[["currentGenomicRegion"]] <- ""
#----------------------------------------------------------------------------------------------------
ui <- fluidPage(

  includeScript("message-handler.js"),
  tags$head(tags$link(rel = "stylesheet", type = "text/css",
                     href = "http://maxcdn.bootstrapcdn.com/font-awesome/4.2.0/css/font-awesome.min.css")
     ),

  titlePanel(sprintf("%s Trena Models", targetGene)),

  sidebarLayout(
     sidebarPanel(
        width=3,
        selectInput("displayGenomicRegion", "Display Genomic Region:",
                    c("CKAP2L region" = "ckap2l.1kb",
                      "Promoter +/- 5kb region"="ckap2l.promoter.10kb",
                      "Proximal enhancer region"="proximal.enhancers",
                      "Distal enhancer regions"="distal.enhancers",
                      "Region encompassing all enhancers"="all.enhancers")),

        selectInput("geneModel", "Choose model:", modelList),

        selectInput("addTrack", "Add Track:",
                   c(" - ",
                     "CKAP2L Enhancers"="enhancers",
                     "DHS (wgEncodeClustered)"="dhs",
                     "rs7594852"="rs7594852"
                     )),

       sliderInput("IGAP.snp.significance", "SNP score",
                    value = 2, #fivenum(-log10(tbl.bs$gwasPval))[3],
                    step=0.1,
                    min = 0,
                    max = 12),
       sliderInput("fimo.motif.significance", "FIMO motif score",
                    value = 2, #fivenum(-log10(tbl.bs$motif.pVal))[3],
                    step=0.1,
                    min = 0,
                    max = 12,
                    round=-2),
       sliderInput("fimo.snp.effect", "SNP binding affinity loss score",
                    value = 0.5,
                    step=0.1,
                    min = 0,
                    max = 12),
       # actionButton("findDisruptiveSNPsButton", "SNPs which disrupt TF binding sites"),
       HTML("<br><br><br>"),
       sliderInput("snpShoulder", "Proximity",
                   value = 5,
                   min = 0,
                   max = 100),
       actionButton("showSNPsNearBindingSitesButton", "Binding sites near SNPs")
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
                   tabPanel(title="SNP",             value="snpInfoTab",      mainPanel(
                                                                                 verbatimTextOutput("tfbsText"),
                                                                                 DTOutput("snpMotifTable"),
                                                                                 imageOutput("logoImage"))
                                                                                 ),
                   tabPanel(title="TF/CKAP2L plot",  value="plotTab", plotOutput("xyPlot", height=800))
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
         load("tbl.bindingSitesWithSnpsAndDeltaScores.RData")
         tokens <- strsplit(motifName, ":")[[1]]
         pfm.name <- tokens[2]
         snp.name <- tokens[3]
         tbl.sub <- subset(tbl.bs, motifName==pfm.name & rsid==snp.name)
         tbl.sub <- tbl.sub[, c("chrom", "start", "end", "strand", "sequence", "sequence.mut", "wtMutMatchDelta")]
         colnames(tbl.sub)[7] <- "10^x match loss"
         output$tfbsText <- renderText({sprintf("%s   %s", pfm.name, snp.name)})
         image.filename <- sprintf("png/%s.png", pfm.name)
         later(function(){
         output$snpMotifTable <- renderDT(tbl.sub, options = list(scrollX = TRUE))
         output$logoImage <- renderImage({
                 list(src=image.filename,
                      contentType="image/png",
                      width=800,
                      height=500,
                      alt="png file not found")},
                 deleteFile=FALSE)
         }, 1)
         updateTabsetPanel(session, "trenaTabs", select="snpInfoTab");
         } # if tfbs-snp
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
      printf("display event: %s", regionName)
      state$currentGenomicRegion <- regionName


      regions <- list(ckap2l.1kb="chr2:112,735,926-112,766,663",
                      ckap2l.promoter.10kb=sprintf("chr2:%d-%d", tss-5000, tss+4999),
                      proximal.enhancers="chr2:112,701,014-112,929,735",
                      distal.enhancers="chr2:112,010,865-112,096,266",
                      all.enhancers="chr2:111,958,695-113,011,015"
                      )
      new.region.name <- input$displayGenomicRegion
      if(!new.region.name %in% names(regions))
         return()
      roi <- regions[[new.region.name]]
      updateTabsetPanel(session, "trenaTabs", select="igvTab");
      myFunc <- function(){
         session$sendCustomMessage(type="showGenomicRegion", message=(list(region=roi)))
         } # myFunc
      later(myFunc, 1)
      later(function() {updateSelectInput(session, "displayGenomicRegion", selected=character(0))}, 1)
      })

   observeEvent(input$addTrack, {
      newTrackName <- input$addTrack
      printf(" addTrack event: %s", newTrackName);
      displayTrack(session, newTrackName)
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
      shoulder <- isolate(input$snpShoulder)
      snpSignificanceThreshold <- isolate(session$input$IGAP.snp.significance)
      motifMatchThreshold <- isolate(session$input$fimo.motif.significance)

      tbl.tfbs <- subset(tbl.fimo, -log10(motif.pVal) >= motifMatchThreshold)
      tbl.tfbs <- tbl.tfbs[, c("chrom", "motifStart", "motifEnd", "motif.pVal", "motifName", "tf")]
      tbl.tfbs$score <- -log10(tbl.tfbs$motif.pVal)
      tbl.tfbs <- tbl.tfbs[, c("chrom", "motifStart", "motifEnd", "score", "motifName", "tf")]
      colnames(tbl.tfbs) <- c("chrom", "start", "end", "score", "motifName", "tf")
      tbl.snpsNearEnough <- data.frame()
      gr.snps <- GRanges(tbl.snp)
      gr.tfbs <- GRanges(data.frame(chrom=tbl.tfbs$chrom,
                                    start=tbl.tfbs$start - shoulder,
                                    end=tbl.tfbs$end+shoulder))
      tbl.ov <- as.data.frame(findOverlaps(gr.snps, gr.tfbs))
      colnames(tbl.ov) <- c("snp", "tfbs")
      tbl.snpsNearEnough <- tbl.snp[tbl.ov$snp, c("chrom", "start", "end", "rsid")]
      tbl.tfbsNearEnough <- tbl.tfbs[tbl.ov$tfbs,]
      dups <- which(duplicated(tbl.tfbsNearEnough[, c("chrom", "start", "end", "tf")]))
      if(length(dups) > 0)
         tbl.tfbsNearEnough <- tbl.tfbsNearEnough[-dups,]
      tfs <- unique(tbl.tfbsNearEnough$tf)
      updateTabsetPanel(session, "trenaTabs", select="igvTab");
      for(this.tf in tfs){
         printf("tf near snp: %s", this.tf)
         temp.filename <- tempfile(tmpdir="tmp", fileext=".bed")
         tbl.tmp <- subset(tbl.tfbsNearEnough, tf==this.tf)[, c("chrom", "start", "end", "motifName", "score")]
         write.table(tbl.tmp, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
         session$sendCustomMessage(type="displayBedTrack",
                                   message=(list(filename=temp.filename,
                                                 trackName=sprintf("%s near snp", this.tf),
                                                 displayMode="EXPANDED",
                                                 color=colors[colorNumber],
                                                 trackHeight=40)))
         } # for this.tf
      })

   observeEvent(input$geneModelTable_cell_clicked, {
      trueRowClick <- length(input$geneModelTable_cell_clicked) > 0
      printf("--- input$geneModelTable_cell_clicked: %s", trueRowClick)
      if(trueRowClick){
         selectedTableRow <- input$geneModelTable_row_last_clicked
         if(!is.null(selectedTableRow)){
            tbl.model <- state$currentModel
            tf <- tbl.model$gene[selectedTableRow]
            later(function(){selectRows(dt.proxy, NULL)}, 0.01);
            if(length(tf) > 0){
               #browser()
               tfSelectionAction <- isolate(input$tfSelectionChoice)
               printf("table row clicked: %s  - actionRequested: %s", tf, tfSelectionAction);
               if(tfSelectionAction == "xyPlot"){
                  updateTabsetPanel(session, "trenaTabs", select="plotTab");
                  tokens <- strsplit(state$currentModelName, "\\.")[[1]]
                  expression.matrix.id <- tokens[length(tokens)]
                  output$xyPlot = renderPlot(plotTfTargetGeneCorrelation(session, tf, expression.matrix.id))
                  }
               if(tfSelectionAction == "displayBindingSites"){
                  displayBindingSites(session, tf)
                  }
               } # if tf
            } # row not null
         } # if trueRowClick
      }) # geneModelTable_cell_clicked

   #output$xyPlot = renderPlot(
   #    {print(" ---- renderPlot");
   #     selectedTableRow <- input$geneModelTable_row_last_clicked
   #     tbl.model <- state$currentModel
   #     browser()
   #     tf <- tbl.model$gene[selectedTableRow]
   #     printf("want an xyplot of %s vs. CKAP2L", tf)
   #     tbl.model <- state$currentModel
   #     tf <- tbl.model$gene[selectedTableRow]
   #     target <- "FLT1";
   #     tbl.pheno <- tbl.pheno[colnames(mtx),]
   #     males   <- which(tbl.pheno$sex_s == "male")
   #     females <- which(tbl.pheno$sex_s == "female")
   #     colors <- rep("blue", nrow(tbl.pheno))
   #     colors[females] <- "red"
   #     tbl.males <- as.data.frame(t(mtx[c(tf, target), males]))
   #     tbl.females <- as.data.frame(t(mtx[c(tf, target), females]))
   #     model.male <- lm(sprintf("%s ~ %s", target, tf), data=tbl.males)
   #     model.female <- lm(sprintf("%s ~ %s", target, tf), data=tbl.females)
   #     plot(as.numeric(mtx[tf,]), as.numeric(mtx[target,]), xlab=tf, ylab=target, col=colors)
   #     abline(model.male, col="blue")
   #     abline(model.female, col="red")
   #     printf("  male: ")
   #     print(summary(model.male))
   #     printf("  female: ")
   #     print(summary(model.female))
   #     })

   output$geneModelTable <- renderDT(
      {model.name <- input$geneModel;
       updateTabsetPanel(session, "trenaTabs", select="geneModelTab");
       printf("    rendering DT, presumably because input$geneModel changes, model.name: %s", model.name);
       selection="single"
       tbl.model <- models[[model.name]]$model
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
     igvShiny("hello shinyApp")
     })

  } # server

#----------------------------------------------------------------------------------------------------
displayTrack <- function(session, trackName)
{
   supportedTracks <- c("enhancers",
                        "dhs",
                        "rs7594852")

   printf("--- displayTrack('%s')", trackName)

   if(!trackName %in% supportedTracks)
      return()


   x <- switch(trackName,
      enhancers  = {list(bed=tbl.enhancers,  name="CKAP2L enhancers",   color="black")},
      dhs              = {list(bed=tbl.dhs,  name="DHS",                color="darkgray")},
      "rs7594852"      = {list(bed=tbl.snp,  name="rs7594852",          color="darkRed")},
      )

    temp.filename <- tempfile(tmpdir="tmp", fileext=".bed")
    write.table(x$bed, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
    updateTabsetPanel(session, "trenaTabs", select="igvTab");
    later(function(){
       session$sendCustomMessage(type="displayBedTrack",
                                 message=(list(filename=temp.filename,
                                               trackName=x$name,
                                               displayMode="EXPANDED",
                                               color=x$color,
                                               trackHeight=40)))
       }, 1)

} # displayTrack
#------------------------------------------------------------------------------------------------------------------------
findDisruptiveSNPs <- function(snpSignificanceThreshold, motifMatchThreshold, bindingLossThreshold)
{
   printf("--- findDisruptiveSNPs")
   printf("  snpSignificanceThreshold: %5.1f", snpSignificanceThreshold)
   printf("  motifMatchThreshold: %5.1f",      motifMatchThreshold)
   printf("  bindingLossThreshold: %5.1f",     bindingLossThreshold);
   #browser()
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
   # browser()
   xyz <- "--- displayBindingSites, consult tbl.fimo..."
   if(!target.tf %in% tbl.fimo$tf){
      printf("tf %s not found in tbl.fimo", target.tf, state$currentModelName)
      return()
      }
   motif.pVal.threshold <- isolate(session$input$fimo.motif.significance)
   tbl.tfbs <- subset(tbl.fimo, tf==target.tf & -log10(motif.pVal) >= motif.pVal.threshold)
   tbl.tfbs <- tbl.tfbs[, c("chrom", "motifStart", "motifEnd", "motif.pVal")]
   tbl.tfbs$score <- -log10(tbl.tfbs$motif.pVal)
   tbl.tfbs <- tbl.tfbs[, c("chrom", "motifStart", "motifEnd", "score")]
   temp.filename <- tempfile(tmpdir="tmp", fileext=".bed")
   write.table(tbl.tfbs, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
   updateTabsetPanel(session, "trenaTabs", select="igvTab");
   colorNumber <<- colorNumber + 1
   later(function(){
       session$sendCustomMessage(type="displayBedTrack",
                                   message=(list(filename=temp.filename,
                                                 trackName=target.tf,
                                                 color=colors[colorNumber],
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
   stopifnot(expression.matrix.id %in% c("cpm", "tmm"))
   mtx <- switch(expression.matrix.id,
                 cpm = mtx.cpm,
                 tmm = mtx.tmm
                 )
   printf("want an xyplot of %s vs. CKAP2l", tf)

   samples.male <- intersect(subset(tbl.pheno, Gender=="Male")$SequencingID, colnames(mtx))
   samples.female <- intersect(subset(tbl.pheno, Gender=="Female")$SequencingID, colnames(mtx))

   color.choices <- list(Male="blue",
                         Female="red"
                         )

   samples.all <- c(samples.male, samples.female)
   colors <- rep("black", length(samples.all))
   colors[match(samples.male, colnames(mtx))] <- color.choices$Male
   colors[match(samples.female, colnames(mtx))] <- color.choices$Female
   plot(mtx[tf, samples.all], mtx[targetGene, samples.all],main="Gene Expression",
        xlab=tf, ylab=targetGene, col=colors, pch=19) # , xlim=c(0,5), ylim=c(0,5))

   tbl.male <- as.data.frame(t(mtx[c(tf, targetGene), samples.male]))
   model.male <- lm(sprintf("%s ~ %s", targetGene, tf), data=tbl.male)
   print(model.male)
   abline(model.male, col=color.choices$Male)

   tbl.female <- as.data.frame(t(mtx[c(tf, targetGene), samples.female]))
   model.female <- lm(sprintf("%s ~ %s", targetGene, tf), data=tbl.female)
   print(model.female)
   abline(model.female, col=color.choices$Female)

#
#      tbl.control <- as.data.frame(t(mtx[c(tf, targetGene), samples.ctl]))
#      model.control <- lm(sprintf("%s ~ %s", targetGene, tf), data=tbl.control)
#      abline(model.control, col=color.choices$control)
#
#      #browser()
#
#      tbl.pathAging <- as.data.frame(t(mtx[c(tf, targetGene), samples.pathAging]))
#      model.pathAging <- lm(sprintf("%s ~ %s", targetGene, tf), data=tbl.pathAging)
#      abline(model.pathAging, col=color.choices$pathAging)
#
#      tbl.PSP <- as.data.frame(t(mtx[c(tf, targetGene), samples.PSP]))
#      model.PSP <- lm(sprintf("%s ~ %s", targetGene, tf), data=tbl.PSP)
#      abline(model.PSP, col=color.choices$PSP)
#
#      legend.names <- names(color.choices)
#      for(i in seq_len(length(legend.names))){
#         condition <- legend.names[i]
#         rsq <- eval(parse(text=sprintf("summary(model.%s)$r.squared", condition)))
#         legend.names [i] <- sprintf("%s (rSq=%5.3f)", condition, rsq)
#         }
#      legend (1.5, -1.8, legend.names, as.character(color.choices))
#      }
#   else{  # no phenotype data, use all samples
#      tf.vals <- mtx[tf, ]
#      targetGene.vals <- mtx[targetGene,]
#      plot(tf.vals, targetGene.vals, main="Gene Expression - Cerebellum", xlab=tf, ylab=targetGene)
#      }


} # plotTfTargetGeneCorrelation
#------------------------------------------------------------------------------------------------------------------------
#Sys.sleep(3)
#shinyApp(ui, server)
