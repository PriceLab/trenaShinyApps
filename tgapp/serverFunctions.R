# serverFunctions.R
#------------------------------------------------------------------------------------------------------------------------
setupIgvAndTableToggling <- function(input)
{
   observeEvent(input$currentGenomicRegion, {
       newValue <- input$currentGenomicRegion
       state$chromLocRegion <- newValue
       # printf("currentGenomicRegion arrived, %s", newValue)
       })

   observeEvent(input$igvHideButton, {
     printf("igvHideButton: %s", input$igvHideButton)
     if(input$igvHideButton %% 2 == 1){
        printf("  --- hiding igv, widening dataTable")
        shinyjs::hide(id = "igvColumn")
        shinyjs::toggleClass("dataTableColumn", "col-sm-3")
        shinyjs::toggleClass("dataTableColumn", "col-sm-12")
      } else {
        printf("  --- showing igv, narrowing dataTable")
        shinyjs::toggleClass("dataTableColumn", "col-sm-12")
        shinyjs::toggleClass("dataTableColumn", "col-sm-3")
        shinyjs::show(id = "igvColumn")
        }
      })

   observeEvent(input$tableHideButton, {
      printf("tableHideButton: %s", input$tableHideButton)
      if(input$tableHideButton %% 2 == 1){
         shinyjs::hide(id = "modelSelectorColumn")
         shinyjs::hide(id = "dataTableColumn")
         shinyjs::toggleClass("igvColumn", "col-sm-9")
         shinyjs::toggleClass("igvColumn", "col-sm-12")
      } else {
         shinyjs::toggleClass("igvColumn", "col-sm-12")
         shinyjs::toggleClass("igvColumn", "col-sm-9")
         shinyjs::show(id = "dataTableColumn")
         shinyjs::show(id = "modelSelectorColumn")
         }
      })
} # setupIgvAndTableToggling
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
   browser()
   xyz <- "back from build"

} # buildFootprintModel
#------------------------------------------------------------------------------------------------------------------------
