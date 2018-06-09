library(shiny)
library(DT)
library(igvShiny)
library(htmlwidgets)
library(GenomicRanges)
library(later)
library(MotifDb)
library(motifStack)
library(RUnit)
#----------------------------------------------------------------------------------------------------
load("trem2-models.RData")
models <- all.models
load("dhs.RData") # tbl.dhs
addResourcePath("tmp", "tmp") # so the shiny webserver can see, and return track files for igv

model.names <- names(all.models)
modelList <- list()
for(model.name in model.names){
   modelList <- c(modelList, eval(parse(text=sprintf('c("%s"="%s")', model.name, model.name))))
   }
#----------------------------------------------------------------------------------------------------
load("enhancers.RData")
load("snps.RData")
load("tbl.bindingSitesWithSnpsAndDeltaScores.RData")
load("tbl.fimoHitsInEnhancers.RData")
#----------------------------------------------------------------------------------------------------
tss <- 41163186
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

  titlePanel("TREM2 Trena Models"),

  sidebarLayout(
     sidebarPanel(
        width=3,
        selectInput("displayGenomicRegion", "Display Genomic Region:",
                    c("TREM2 +/- 40kb" = "tssPlusMinus40kb",
                      "TREM2 2200 bp promoter"="promoter2200",
                      "TREM2 10kb promoter"="promoter10kb",
                      "TREM gene family"="overview")),

        selectInput("geneModel", "Choose model:", modelList),

        selectInput("addTrack", "Add Track:",
                   c(" - ",
                     "TREM2 Enhancers"="trem2.enhancers",
                     "TREML1 Enhancers"="treml1.enhancers",
                     "DHS (wgEncodeClustered)"="dhs",
                     "IGAP GWAS SNPs"="snps",
                     "Carrasquillo 2016 eqtl snps"="eqtl.snps"
                     )),

       sliderInput("IGAP.snp.significance", "-log10 (SNP pval)",
                    value = 2, #fivenum(-log10(tbl.bs$gwasPval))[3],
                    step=0.1,
                    min = 0,
                    max = 12),
       sliderInput("fimo.motif.significance", "-log10 (FIMO motif pval)",
                    value = 2, #fivenum(-log10(tbl.bs$motif.pVal))[3],
                    step=0.1,
                    min = 0,
                    max = 12,
                    round=-2),
       sliderInput("fimo.snp.effect", "-log10 (FIMO SNP effect) - delta qval",
                    value = 0.5,
                    step=0.1,
                    min = 0,
                    max = 12),
        actionButton("findDisruptiveSNPsButton", "Find SNPs which disrupt TF binding sites")
        ),

     mainPanel(
       tabsetPanel(type="tabs",
                   id="trenaTabs",
                   tabPanel(title="IGV",             value="igvTab",          igvShinyOutput('igvShiny')),
                   tabPanel(title="READ ME",         includeHTML("readme.html")),
                   tabPanel(title="trena model",     value="geneModelTab",    DTOutput("geneModelTable")),
                   tabPanel(title="SNP",             value="snpInfoTab",      mainPanel(
                                                                                 verbatimTextOutput("tfbsText"),
                                                                                 DTOutput("snpMotifTable"),
                                                                                 imageOutput("logoImage"))
                                                                                 ),
                   tabPanel(title="TF/TREM2 xy plot",  value="plotTab", plotOutput("xyPlot", height=800))
                   )
       ) # mainPanel
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
                      width=500,
                      height=500,
                      alt="alt text")},
                 deleteFile=FALSE)
         }, 0.5)
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
      regions <- list(overview="chr6:41,049,342-41,367,926",
                      promoter2200=sprintf("chr6:%d-%d", tss-200, tss+1999),
                      promoter10kb=sprintf("chr6:%d-%d", tss-5000, tss+4999),
                      tssPlusMinus40kb="chr6:41,122,760-41,200,161"
                      )
      new.region.name <- input$displayGenomicRegion
      if(!new.region.name %in% names(regions))
         return()
      roi <- regions[[new.region.name]]
      updateTabsetPanel(session, "trenaTabs", select="igvTab");
      myFunc <- function(){
         session$sendCustomMessage(type="showGenomicRegion", message=(list(region=roi)))
         } # myFunc
      later(myFunc, 5)
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

   observeEvent(input$geneModelTable_cell_clicked, {
      printf("--- input$geneModelTable_cell_clicked")
      selectedTableRow <- input$geneModelTable_row_last_clicked
      if(!is.null(selectedTableRow)){
         tbl.model <- state$currentModel
         tf <- tbl.model$gene[selectedTableRow]
         if(length(tf) > 0){
            printf("table row clicked: %s", tf );
            # updatePlot(session, tf, state$currentModelName);
            updateTabsetPanel(session, "trenaTabs", select="igvTab");
            later(function(){selectRows(dt.proxy, NULL)}, 0.5);
            later(function(){displayBindingSites(session, input, tf)}, 1);
            }
         } # row not null
      })

   output$xyPlot = renderPlot(
       {print(" ---- renderPlot");
        selectedTableRow <- input$geneModelTable_row_last_clicked
        tbl.model <- state$currentModel
        browser()
        tf <- tbl.model$gene[selectedTableRow]
        printf("want an xyplot of %s vs. TREM2", tf)
        tbl.model <- state$currentModel
        tf <- tbl.model$gene[selectedTableRow]
        target <- "FLT1";
        tbl.pheno <- tbl.pheno[colnames(mtx),]
        males   <- which(tbl.pheno$sex_s == "male")
        females <- which(tbl.pheno$sex_s == "female")
        colors <- rep("blue", nrow(tbl.pheno))
        colors[females] <- "red"
        tbl.males <- as.data.frame(t(mtx[c(tf, target), males]))
        tbl.females <- as.data.frame(t(mtx[c(tf, target), females]))
        model.male <- lm(sprintf("%s ~ %s", target, tf), data=tbl.males)
        model.female <- lm(sprintf("%s ~ %s", target, tf), data=tbl.females)
        plot(as.numeric(mtx[tf,]), as.numeric(mtx[target,]), xlab=tf, ylab=target, col=colors)
        abline(model.male, col="blue")
        abline(model.female, col="red")
        printf("  male: ")
        print(summary(model.male))
        printf("  female: ")
        print(summary(model.female))
        })

   output$geneModelTable <- renderDT(
      {model.name <- input$geneModel;
       printf("    rendering DT, presumably because input$geneModel changes, model.name: %s", model.name);
       selection="single"
       tbl.model <- all.models[[model.name]]$model
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
   supportedTracks <- c("trem2.enhancers",
                        "treml1.enhancers",
                        "dhs",
                        "snps",
                        "eqtl.snps")

   printf("--- displayTrack('%s')", trackName)

   if(!trackName %in% supportedTracks)
      return()


   x <- switch(trackName,
      trem2.enhancers  = {list(bed=tbl.enhancers.trem2,  name="TREM2 enhancers",   color="black")},
      treml1.enhancers = {list(bed=tbl.enhancers.treml1, name="TREML1 enahancers", color="black")},
      dhs              = {list(bed=tbl.dhs,              name="DHS",               color="darkgray")},
      snps             = {list(bed=tbl.snp,              name="GWAS snp",          color="darkRed")},
      eqtl.snps        = {list(bed=tbl.eqtl,             name="eqtl snp",          color="darkRed")}
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
createMotifLogo.pngs <- function(pfms)
{
   #load("tbl.bindingSitesWithSnpsAndDeltaScores.RData")
   #tfs <- unique(tbl.bs$tf)
   # pfms.1 <- query(MotifDb, "sapiens", c("hocomoco", "jaspar2018", "swissregulon"))
   #pfms.tf <- query(pfms.1, "", tfs)
   #length(pfms.tf)
   for(pfm.name in names(pfms.oi)){
      filename <- file.path("../shinyApp", "png", sprintf("%s.png", pfm.name))
      printf("filename: %s", filename)
      png(filename)
      plotMotifs(pfms.oi[pfm.name])
      dev.off()
      } # for i

} # createMotifLogo.pngs
#------------------------------------------------------------------------------------------------------------------------
removeTrack <- function(session, trackName)
{
  session$sendCustomMessage(type="removeTrack", message=list(trackName=trackName))

} # removeTrack
#----------------------------------------------------------------------------------------------------
displayBindingSites <- function(session, input, target.tf)
{
   #tbl.bindingSites <- all.models[[state$currentModelName]]$regulatoryRegions
   browser()
   #printf(" looking for binding sites in %s", state$currentModelName)
   if(!target.tf %in% tbl.fimoHitsInEnhancers$tf){
      printf("tf %s not found in tbl.fimoHitsInEnhancers", target.tf, state$currentModelName)
      return()
      }
      # two kinds of bindingSite tables: one from footprint database (has "fp_start"), one from
      # MotifMatcher (has "motifStart").  branch appropriately
   #if("fp_start" %in% colnames(tbl.bindingSites))
   #   tbl.tfbs <- subset(tbl.bindingSites, geneSymbol==tf)[, c("chrom", "fp_start", "fp_end", "shortMotif")]
   #if("motifStart" %in% colnames(tbl.bindingSites))
      #tbl.tfbs <- subset(tbl.bindingSites, geneSymbol==tf)[, c("chrom", "motifStart", "motifEnd", "shortMotif")]
     #threshold <- isolate(input$pwmMatchThreshold) / 100
     # printf("---  using slider threshold %f", threshold);
     #tbl.tfbs <- subset(tbl.bindingSites, motifRelativeScore >= threshold)[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")]
     # tbl.tfbs <- tbl.tfbs[order(tbl.tfbs$motifStart, decreasing=FALSE),]
   tbl.tfbs <- subset(tbl.fimoHitsInEnhancers, tf==target.tf & p.value < 1e-4)
   tbl.tfbs$chrom <- "chr6"
   tbl.tfbs <- tbl.tfbs[, c("chrom", "start", "end", "motif", "p.value")]
   printf("      found %d sites", nrow(tbl.tfbs))
   #printf("scores on binding sites: ")
   #print(fivenum(tbl.tfbs$motifRelativeScore))
   temp.filename <- tempfile(tmpdir="tmp", fileext=".bed")
   write.table(tbl.tfbs, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
   session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName=target.tf,
                                              color="darkRed",
                                              trackHeight=40)))



} # displayBindingSites
#------------------------------------------------------------------------------------------------------------------------
#Sys.sleep(3)
#shinyApp(ui, server)
