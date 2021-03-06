# igap.R: support gene selection within the app
#------------------------------------------------------------------------------------------------------------------------
library(shiny)
library(shinydashboard)
library(shinyjs)
library(later)
library(igvShiny)
library(DT)
library(VariantAnnotation)
#------------------------------------------------------------------------------------------------------------------------
addResourcePath("tracks", "tracks")
#------------------------------------------------------------------------------------------------------------------------
library(trenaSGM)
library(TrenaProjectIGAP)
#------------------------------------------------------------------------------------------------------------------------
library(RColorBrewer)
colors <- brewer.pal(8, "Dark2")
totalColorCount <- length(colors)
colorNumber <- 0
#------------------------------------------------------------------------------------------------------------------------
state <- new.env(parent=emptyenv())
state$models <- list()
#------------------------------------------------------------------------------------------------------------------------
model.count <- 0   # for creating default model names
#------------------------------------------------------------------------------------------------------------------------
igap <- TrenaProjectIGAP();
setTargetGene(igap, getSupportedGenes(igap)[1])

expressionMatrixNames <- getExpressionMatrixNames(igap)
tbl.enhancers <- getEnhancers(igap)
tbl.dhs <- getEncodeDHS(igap)
footprintDatabases <- getFootprintDatabaseNames(igap)
tbl.transcripts <- getTranscriptsTable(igap)


#------------------------------------------------------------------------------------------------------------------------
# calculate the initial igv region, using the (typically? always? single transcripted designated as "gene"
#------------------------------------------------------------------------------------------------------------------------
# tbl.gene <- subset(tbl.transcripts, moleculetype=="gene")[1,]
# start <- tbl.gene$start
# end   <- tbl.gene$endpos
# chrom <- tbl.gene$chr
# padding <- round(0.25 * (end-start))
# initial.chromLoc <- sprintf("%s:%d-%d", chrom, start-padding, end+padding)
#------------------------------------------------------------------------------------------------------------------------
createSidebar <- function()
{
  dashboardSidebar(
    sidebarMenu(id="sidebarMenu",
       menuItem("IGV and Current Models", tabName = "igvAndTable"),
       menuItem("Build a new Model",      tabName = "buildModels")
      ),
    conditionalPanel(id="igvTableWidgets",
        condition = "input.sidebarMenu == 'igvAndTable'",
        actionButton(inputId = "igvHideButton", label = "Toggle IGV"),
        actionButton(inputId = "tableHideButton", label = "Toggle table"),
        selectInput("chooseGene", "Choose Gene:", c("", getSupportedGenes(igap))),

        selectInput("addTrack", "Add Track:",
                    c("",
                     "Enhancers"="enhancers",
                     "DHS open chromatin (encode, clustered)"="dhs"
                     )),
        selectInput("displayGenomicRegion", "Display Genomic Region:",
                    c("",
                      "Full Gene" = "fullGeneRegion",
                      "Full Enhancers Region" = "fullEnhancerRegion")),
        actionButton(inputId="removeUserAddedTracks", label="Remove Added Tracks")
        ) # conditionalPanel
     ) # dashboardSidebar

} # createSidebar
#------------------------------------------------------------------------------------------------------------------------
createBuildModelTab <- function()
{
   tab <- tabItem(tabName="buildModels",
             h4("Build trena regulatory model from DNase footprints in the genomic region currently displayed in IGV, with these constraints:"),
             br(),
             fluidRow(
                column(3, offset=1, checkboxGroupInput("footprintDatabases", "Footprint Databases",
                                            footprintDatabases, selected=footprintDatabases[1])),
                column(4, radioButtons("intersectWithRegions", "Intersect footprints with:",
                                             c("GeneHancer" = "genehancer",
                                               "Encode DHS" = "encodeDHS",
                                               "GeneHancer OR Encode DHS" = "geneHancerOrEncode",
                                               "GeneHancer AND Encode DHS" = "geneHancerPlusEncode",
                                               "Use all footprints in region" = "allDNAForFootprints"),
                                       selected="allDNAForFootprints")),
                column(3, radioButtons("motifMapping", "Map motifs to TFs via:",
                                       c("MotifDb", "TFClass", "MotifDb + TFClass"),
                                       selected="MotifDb"))
                ),
             fluidRow(
                column(6, selectInput("expressionSet", "With this gene expression dataset:",
                                      c("", expressionMatrixNames))),
                column(5, textInput("modelNameTextInput", "Model name", width=500))),
             fluidRow(column(width=1, offset=4, actionButton("buildModelButton", "Build model"))),
             br(),br(),
             wellPanel(style="overflow-y:scroll; height:200px", pre(id = "console"))
             ) # tabItem

   return(tab)

} # createBuildModelTab
#------------------------------------------------------------------------------------------------------------------------
createBody <- function()
{
   body <- #fluidRow(
      tabItems(
         tabItem(tabName="igvAndTable",
            fluidRow(
               column(width=9, id="igvColumn", igvShinyOutput('igvShiny', height="2000px")),
               column(width=3, id="dataTableColumn",
                      selectInput("modelSelector", NULL,  "- no models yet -"),
                      DTOutput("table"),
                      br(),
                      #wellPanel(
                        h4("TF selection will display:"),
                        selectInput("selectRowAction", NULL,  c("(no action)", "Footprints", "ChIP-seq hits",
                                                                sprintf("Plot expression, TF vs %s", getTargetGene(igap))))
                      ) # column
               ) # fluidRow
            ), # tabItem 1
         createBuildModelTab()
         ) # tabItems

   return(body)

} # createBody
#------------------------------------------------------------------------------------------------------------------------
ui <- dashboardPage(

  dashboardHeader(title=sprintf("trena %s", getTargetGene(igap))),
  createSidebar(),
  dashboardBody(
      tags$head(tags$style(HTML('
        .main-header .logo {
          font-family: "Georgia", Times, "Times New Roman", serif;
          font-weight: bold;
          font-size: 24px;
          }
        #igvShiny{
          background: white;
          border: 1px solid black;
          border-radius: 5px;
          margin: 5px;
          margin-right: 15px;
          }
        #table{
          border: 1px solid black;
          border-radius: 5px;
          }
       '))),
     createBody()),
  useShinyjs()
  ) # dashboardPage

#------------------------------------------------------------------------------------------------------------------------
server <- function(session, input, output){

   observeEvent(input$chooseGene, ignoreInit=TRUE, {
       newGene <- isolate(input$chooseGene)
       setTargetGene(igap, newGene)
       printf("igap targetGene: %s", getTargetGene(igap))
       shinyjs::html(selector=".logo", html=sprintf("trena %s", newGene), add=FALSE)
       showGenomicRegion(session, getGeneRegion(igap, flankingPercent=20))
       loadAndDisplayRelevantVariants(session, newGene)
       })


   output$igvShiny <- renderIgvShiny({
      options <- list(genomeName="hg38",
                      initialLocus=getGeneRegion(igap, flankingPercent=20),
                      displayMode="EXPANDED",
                      trackHeight=300)
      igvShiny(options) # , height=800)
      })

   output$table = DT::renderDataTable(data.frame(),
                                      width="800px",
                                      class='nowrap display',
                                      selection="single",
                                      extensions="FixedColumns",
                                      options=list(scrollX=TRUE,
                                                   scrollY="500px",
                                                   dom='t',
                                                   paging=FALSE,
                                                   autowWdth=FALSE,
                                                   fixedColumns=list(leftColumns=1)
                                                   ))

   observeEvent(input$table_rows_selected, {
      selectedTableRow <- isolate(input$table_rows_selected)
      dispatch.rowClickInModelTable(session, input, output, selectedTableRow)
      #current.model.name <- isolate(input$modelSelector)
      #if(current.model.name %in% ls(state$models)){
      #   tf.names <- state$models[[current.model.name]]$model$gene
      #   tf.name <- tf.names[selectedTableRow]
      #   printf("tf: %s", tf.name)
      #} else {
      #   printf("failed to find %s in current model %s",
      }) # observe row selection event

   observeEvent(input$currentGenomicRegion, {
      new.region <- isolate(input$currentGenomicRegion)
      #printf("new region: %s", new.region)
      state[["chromLocRegion"]] <- new.region
      })

   setupIgvAndTableToggling(session, input);
   setupAddTrack(session, input, output)
   setupDisplayRegion(session, input, output)
   setupBuildModel(session, input, output)

   #later(function(){
   #         state$chromLocRegion <- initial.chromLoc
   #         showGenomicRegion(session, initial.chromLoc)
   #         }, 2)

} # server
#------------------------------------------------------------------------------------------------------------------------
setupIgvAndTableToggling <- function(session, input)
{
   observeEvent(input$currentGenomicRegion, {
       newValue <- input$currentGenomicRegion
       state$chromLocRegion <- newValue
       # printf("currentGenomicRegion arrived, %s", newValue)
       })

   observeEvent(input$igvHideButton, {
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
      printf("--- calling redrawIgvWidget")
      redrawIgvWidget(session)
      })

   observeEvent(input$tableHideButton, {
      if(input$tableHideButton %% 2 == 1){
         shinyjs::hide(id = "modelSelectorColumn")
         shinyjs::hide(id = "dataTableColumn")
         shinyjs::toggleClass("igvColumn", "col-sm-9")
         shinyjs::toggleClass("igvColumn", "col-sm-12")
      } else {
         shinyjs::toggleClass("igvColumn", "col-sm-12")
         shinyjs::toggleClass("igvColumn", "col-sm-9")
         shinyjs::hide(id = "modelSelectorColumn")
         shinyjs::hide(id = "dataTableColumn")
         shinyjs::show(id = "dataTableColumn")
         shinyjs::show(id = "modelSelectorColumn")
         }
      printf("--- calling redrawIgvWidget")
      redrawIgvWidget(session)
      })
} # setupIgvAndTableToggling
#------------------------------------------------------------------------------------------------------------------------
mapToChromLoc <- function(regionName)
{
   roi <- getTargetGene(igap)   # a safe fallback

   if(regionName == "traditionalPromoter"){
      tbl.transcripts <- getTranscriptsTable(igap)
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
   tbl.gene <- subset(getTranscriptsTable(igap), moleculetype=="gene")[1,]
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
   mtx <- loadExpressionData(igap, "FilteredLengthScaledTPM8282018-vsn")

   build.spec <- list(title=sprintf("%s model %d", getTargetGene(igap), 1),
                      type="footprint.database",
                      regions=tbl.regions,
                      geneSymbol=getTargetGene(igap),
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

   fpBuilder <- FootprintDatabaseModelBuilder("hg38", getTargetGene(igap), build.spec, quiet=TRUE)
   x <- build(fpBuilder)
   browser()
   xyz <- "back from build"

} # buildFootprintModel
#------------------------------------------------------------------------------------------------------------------------
setupAddTrack <- function(session, input, output)
{
   observeEvent(input$addTrack, {
      newTrackName <- input$addTrack
      if(newTrackName == "") return()
      printf(" addTrack event: %s", newTrackName);
      displayTrack(session, newTrackName)
      later(function() {updateSelectInput(session, "addTrack", selected=character(0))}, 1)
      })

} # setupAddTrack
#------------------------------------------------------------------------------------------------------------------------
displayTrack <- function(session, trackName)
{
   supportedTracks <- c("enhancers",
                        "dhs")

   printf("--- displayTrack('%s')", trackName)

   switch(trackName,
          enhancers = {displayEnhancersTrack(session)},
          dhs       = {displayEncodeDhsTrack(session)}
          )

} # displayTrack
#------------------------------------------------------------------------------------------------------------------------
displayEnhancersTrack <- function(session)
{
   tbl.enhancers <- getEnhancers(igap)
   tbl.tmp <- tbl.enhancers[, c("chrom", "start", "end", "combinedScore")]
   loadBedGraphTrack(session, "GeneHancer", tbl.tmp, color="black", trackHeight=25, autoscale=FALSE, min=0, max=20)

} # displayEnhancersTrack
#------------------------------------------------------------------------------------------------------------------------
displayEncodeDhsTrack <- function(session)
{
   tbl.dhs <- getEncodeDHS(igap)
   tbl.tmp <- tbl.dhs[, c("chrom", "chromStart", "chromEnd", "score")]
   colnames(tbl.tmp) <- c("chrom", "start", "end", "score")
   loadBedGraphTrack(session, "DHS", tbl.tmp, color="black", trackHeight=25, autoscale=TRUE)

} # displayEncodeDhsTrack
#------------------------------------------------------------------------------------------------------------------------
setupDisplayRegion <- function(session, input, output)
{
   observeEvent(input$displayGenomicRegion, {
      requestedRegion <- input$displayGenomicRegion
      if(requestedRegion == "") return()
      printf(" displayRegion: %s", requestedRegion)
      margin <- 5000
      loc.string <-switch(requestedRegion,
                          fullEnhancerRegion = {getGeneEnhancersRegion(igap, 10)},
                          fullGeneRegion = {getGeneRegion(igap, 20)})
      showGenomicRegion(session, loc.string);
      later(function() {updateSelectInput(session, "displayGenomicRegion", selected=character(0))}, 1)
      })

   observeEvent(input$removeUserAddedTracks, {
      removeUserAddedTracks(session)
      })

} # setupDisplayRegion
#------------------------------------------------------------------------------------------------------------------------
setupBuildModel <- function(session, input, output)
{

   observeEvent(input$modelSelector, ignoreInit=TRUE, {
       modelName <- isolate(input$modelSelector)
       printf("--- new model selected: %s", modelName)
       new.table <- state$models[[modelName]]$model
       displayModel(session, input, output, new.table, modelName)
       })

   observeEvent(input$modelNameTextInput, {
       currentName <- isolate(input$modelNameTextInput)
       if(nchar(currentName) >= 1)
          shinyjs::enable("buildModelButton")
       else
          shinyjs::disable("buildModelButton")
       })

   observeEvent(input$buildModelButton, {
      getGenomicRegion(session)
      shinyjs::html(id="console", html="", add=FALSE)  # clear the console
      withCallingHandlers(
                          {buildModel(session, input, output);
                           model.count <- length(state$models)
                           new.model.name <- names(state$models)[model.count]
                           new.table <- state$models[[model.count]]$model
                           displayModel(session, input, output, new.table, new.model.name)
                           updateTabItems(session, "sidebarMenu", select="igvAndTable")
                           },
         message=function(m){
            #printf("got message for console");
            #print(m)
            shinyjs::html(id="console", html=m$message, add=TRUE)
            })
      #message(sprintf("back from buildModel, current models are %s", paste(names(state$models), collapse=", ")))
      #model.count <- length(state$models)
      #new.model.name <- names(state$models)[model.count]
      #new.table <- state$models[[model.count]]$model
      #output$buildInProgressTextArea <- renderText("build complete")
      #displayModel(session, input, output, new.table, new.model.name)
      #updateTabItems(session, "sidebarMenu", select="igvAndTable")
      })

} # setupBuildModel
#------------------------------------------------------------------------------------------------------------------------
buildModel <- function(session, input, output)
{
   model.name <- sprintf("trena.model.%s", input$modelNameTextInput)
   message(sprintf("about to build '%s'", model.name))
   footprint.database.names <- input$footprintDatabases
   tracks.to.intersect.with <- input$intersectWithRegions
   motifMapping <- isolate(input$motifMapping)
   if(tolower(motifMapping) == "motifdb + tfclass")
      motifMapping <- c("MotifDb", "TFClass")
   expressionMatrixName <- input$expressionSet
   full.roi <- state$chromLocRegion
   chrom.loc <- trena::parseChromLocString(full.roi)
   message(sprintf("  fpdb: %s", paste(footprintDatabases, collapse=", ")))
   message(sprintf("   roi: %s", full.roi))
   message(sprintf("   mtx: %s", expressionMatrixName))
   message(printf("  intersect with: %s", paste(tracks.to.intersect.with, collapse=",")))

   tbl.gene <- subset(tbl.transcripts, moleculetype=="gene")[1,]
   strand <- tbl.gene$strand
   tss <- tbl.gene$start
   if(strand == "-")
      tss <- tbl.gene$endpos

   run.trenaSGM(igap,
                model.name,
                chrom.loc$chrom, chrom.loc$start, chrom.loc$end,
                tss,
                expressionMatrixName,
                tracks.to.intersect.with,
                footprint.database.names,
                motifMapping)

} # buildModel
#------------------------------------------------------------------------------------------------------------------------
run.trenaSGM <- function(igap,
                         model.name,
                         chromosome, start.loc, end.loc,
                         tss,
                         expression.matrix.name,
                         tracks.to.intersect.with,
                         footprint.database.names,
                         motifMapping)
{
   message(sprintf("--- entering run.trenaSGM"))
      # no search for overlaps just yet: ignore "tracks.to.intersect.with"
   tbl.regions <- data.frame(chrom=chromosome, start=start.loc, end=end.loc, stringsAsFactors=FALSE)
   mtx <- getExpressionMatrix(igap, expression.matrix.name)
   geneSymbol <- getTargetGene(igap)
   message(sprintf("tbl.regions, width: %d", with(tbl.regions, end - start)))

   build.spec <- list(title=model.name,
                      type="footprint.database",
                      regions=tbl.regions,
                      geneSymbol=getTargetGene(igap),
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=footprint.database.names,
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(identifierType="geneSymbol"),
                      tfMapping=motifMapping,
                      tfPrefilterCorrelation=0.2,
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      orderModelByColumn="pearsonCoeff",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))
     save(build.spec, file=sprintf("%s.buildSpec.RData", model.name))

     #------------------------------------------------------------
     # use the above build.spec: a small region, high correlation
     # required, MotifDb for motif/tf lookup
     #------------------------------------------------------------

   fpBuilder <- FootprintDatabaseModelBuilder("hg38", getTargetGene(igap), build.spec, quiet=FALSE)
   x <- build(fpBuilder)
   message(sprintf("back from build, top 10 tfs:"))
   message(head(x$model, n=10))
   state$models[[model.name]] <- x

} # run.trenaSGM
#------------------------------------------------------------------------------------------------------------------------
# beautify the data.frame, display it in the UI DataTable, update the modelSelector pulldown
displayModel <- function(session, input, output, tbl.model, new.model.name)
{
   printf("--- entering displayModel: %s", new.model.name)
   print(dim(tbl.model))
   tf.names <- tbl.model$gene
   tbl.model <- tbl.model[, -1]
   if("lassoPValue" %in% colnames(tbl.model))
       tbl.model <- roundNumericColumns(tbl.model, 4, "lassoPValue")
   rownames(tbl.model) <- tf.names
   max.model.rows <- 20
   if(nrow(tbl.model) > max.model.rows)
      tbl.model <- tbl.model[1:max.model.rows,]

   output$table = DT::renderDataTable(tbl.model,
                                      width="800px",
                                      class='nowrap display',
                                      selection="single",
                                      extensions="FixedColumns",
                                      options=list(scrollX=TRUE,
                                                   scrollY="500px",
                                                   dom='t',
                                                   paging=FALSE,
                                                   autowWdth=FALSE,
                                                   fixedColumns=list(leftColumns=1)
                                                   ))


   updateSelectInput(session, "modelSelector",
                     choices=names(state$models),
                     selected=new.model.name)

} # displayModel
#------------------------------------------------------------------------------------------------------------------------
dispatch.rowClickInModelTable <- function(session, input, output, selectedTableRow)
{
      #current.model.name <- isolate(input$modelSelector)
      #if(current.model.name %in% ls(state$models)){
      #   printf("tf: %s", tf.name)

   current.model.name <- isolate(input$modelSelector)
   tf.names <- state$models[[current.model.name]]$model$gene
   if(length(tf.names) > 20) tf.names <- tf.names[1:20]
   tf.name <- tf.names[selectedTableRow]
   action.name     <- isolate(input$selectRowAction)
   expression.matrix.name <- isolate(input$expressionSet)
   #printf("%s of model %s, expression.set %s: %s", tf.name, current.model.name, expression.matrix.name, action.name)
   # browser()
   xyz <- "botton of dispatch.rowClick"

   if(action.name == "Footprints"){
      tbl.fp <- state$models[[current.model.name]]$regulatoryRegions
      tbl.fp.tf <- subset(tbl.fp, geneSymbol==tf.name)
      dups <- which(duplicated(tbl.fp.tf$loc))
      if(length(dups) > 0)
         tbl.fp.tf <- tbl.fp.tf[-dups,]
      tbl.tmp <- tbl.fp.tf[ c("chrom", "fp_start", "fp_end", "shortMotif")]
      colnames(tbl.tmp) <- c("chrom", "start", "end", "name")
      colorNumber <<- (colorNumber %% totalColorCount) + 1
      next.color <- colors[colorNumber]
      loadBedTrack(session, sprintf("FP-%s", tf.name), tbl.tmp, color=next.color, trackHeight=25)
      } # if footprints

   if(action.name == "ChIP-seq hits"){
      full.roi <- state$chromLocRegion
      chrom.loc <- trena::parseChromLocString(full.roi)
      if(!exists("tbl.chipSeq")){
         printf("=== retrieving chip-seq data from database")
         tbl.chipSeq <<- with(chrom.loc, getChipSeq(igap, chrom, start, end,  tf.names))
         printf("hits among all tfs in this model: %d", nrow(tbl.chipSeq))
         }
      tbl.hits <- subset(tbl.chipSeq, tf == tf.name)
      printf("ChIP-seq hits for %s: %d", tf.name, nrow(tbl.hits))
      if(nrow(tbl.hits) > 0){
         tbl.tmp <- tbl.hits[, c("chr", "start", "end", "name")]
         colorNumber <<- (colorNumber %% totalColorCount) + 1
         next.color <- colors[colorNumber]
         loadBedTrack(session, sprintf("Cs-%s", tf.name), tbl.tmp, color=next.color, trackHeight=25)
         }
      } # ChIP-seq hits


} # dispatch.rowClickInModelTable
#------------------------------------------------------------------------------------------------------------------------
display.chipseq.track <- function(session, input, output, tf)
{

} # display.chipseq.track
#------------------------------------------------------------------------------------------------------------------------
display.footprint.track <- function(session, input, output, tf)
{
   model.name <- isolate(input$modelSelector)

} # display.footprint.track
#------------------------------------------------------------------------------------------------------------------------
loadAndDisplayRelevantVariants <- function(session, newGene)
{
   variant.filenames <- getVariantDatasetNames(igap)
   short.names <- names(variant.filenames)
   file.paths <- unlist(variant.filenames, use.names=FALSE)
   vcf.file.indices <- grep(".vcf", short.names, ignore.case=TRUE)
   thisGene.indices <- grep(newGene, short.names, ignore.case=TRUE)
   final.indices <- intersect(vcf.file.indices, thisGene.indices)

   if(length(final.indices) == 0)
      return()

   current.gene.vcf.files <- file.paths[final.indices]
   for(i in seq_len(length(current.gene.vcf.files))){
      vcfFile <- current.gene.vcf.files[i]
      vcfData <- readVcf(vcfFile)
      loadVcfTrack(session, "vcf", vcfData)
      }

} # loadAndDisplayRelevantVariants
#------------------------------------------------------------------------------------------------------------------------
app <- shinyApp(ui, server)
