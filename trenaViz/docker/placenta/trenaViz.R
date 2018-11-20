# trenaViz.R:  build models for, and explore genomics of, a TrenaProject subclass
#------------------------------------------------------------------------------------------------------------------------
library(shiny)
library(shinydashboard)
library(shinyjs)
library(later)
library(igvShiny)
library(DT)
library(VariantAnnotation)
#------------------------------------------------------------------------------------------------------------------------
library(trenaSGM)
#------------------------------------------------------------------------------------------------------------------------
library(RColorBrewer)
colors <- brewer.pal(8, "Dark2")
totalColorCount <- length(colors)
colorNumber <- 0
#------------------------------------------------------------------------------------------------------------------------
MAX.TF.COUNT <- 50   # we display the top max.tf.cout TFs coming back from trena
#------------------------------------------------------------------------------------------------------------------------
state <- new.env(parent=emptyenv())
state$models <- list()
state$tbl.chipSeq <- NULL
#------------------------------------------------------------------------------------------------------------------------
model.count <- 0   # for creating default model names
#------------------------------------------------------------------------------------------------------------------------
# in interactive sessions, projectName should already be defined, eg, projectName <- "TrenaProjectIGAP"
# for non-interactive sessions, we require the projectName to be defined on the command line:
#  R -f trenaViz.R --args TrenaProjectIGAP

args <- commandArgs(trailingOnly=TRUE)

if(length(args) == 1)
   projectName <- args[1]

stopifnot(exists("projectName"))

library(projectName, character.only=TRUE)
initialization.command <- sprintf("trenaProject <- %s()", projectName)
eval(parse(text=initialization.command))
setTargetGene(trenaProject, getSupportedGenes(trenaProject)[1])
printf("targetGene: %s", getTargetGene(trenaProject))

state$chromLocRegion <- getGeneRegion(trenaProject, flankingPercent=20)
expressionMatrixNames <- getExpressionMatrixNames(trenaProject)
state$tbl.enhancers <- getEnhancers(trenaProject)
state$tbl.dhs <- getEncodeDHS(trenaProject)
footprintDatabases <- getFootprintDatabaseNames(trenaProject)
state$tbl.transcripts <- getTranscriptsTable(trenaProject)

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
       menuItem("Build a new Model",      tabName = "buildModels"),
       menuItem("Introductory video",     tabName = "video")
      ),
    conditionalPanel(id="igvTableWidgets",
        condition = "input.sidebarMenu == 'igvAndTable'",
        actionButton(inputId = "igvHideButton", label = "Toggle IGV"),
        actionButton(inputId = "tableHideButton", label = "Toggle table"),
        selectInput("chooseGene", "Choose Gene:", c("", getSupportedGenes(trenaProject))),
        trackList <-
        selectInput("addTrack", "Add Track:", additionalTracksOnOffer()),
        selectInput("displayGenomicRegion", "Display Genomic Region:",
                    c("",
                      "Full Gene" = "fullGeneRegion",
                      "Full Enhancers Region" = "fullEnhancerRegion")),
        actionButton(inputId="removeUserAddedTracks", label="Remove Added Tracks")
        ) # conditionalPanel
     ) # dashboardSidebar

} # createSidebar
#------------------------------------------------------------------------------------------------------------------------
additionalTracksOnOffer <- function()
{
   trackOfferings <- list("",
                          "Enhancers"="enhancers",
                          "DNase HS (encode, clustered)"="dhs")

   variantTrackNames <- getVariantDatasetNames(trenaProject)

   for(i in seq_len(length(variantTrackNames))){
      trackName <- variantTrackNames[i]
      if(grepl("gwas", trackName, ignore.case=TRUE)){
         new.list <- c(trackName)
         names(new.list) <- trackName
         trackOfferings <- c(trackOfferings, new.list)
         } # if 'gwas' in the name
      if(grepl(".variants", trackName, ignore.case=TRUE, fixed=TRUE)){
         new.list <- c(trackName)
         names(new.list) <- trackName
         trackOfferings <- c(trackOfferings, new.list)
         } # if '.variants' in the name
      } # for i

   return(trackOfferings)

} # additionalTracksOnOffer
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
                                               "GeneHancer AND Encode DHS" = "geneHancerAndEncode",
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
createVideoTab <- function()
{
   tab <- tabItem(tabName="video",
                  h4("Using trenaViz: a video tutorial"),
                  includeHTML("videoTutorial.html")
                  )
   return(tab)

} # createVideoTab
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
                                                                sprintf("Plot expression, TF vs target gene")))
                      ) # column
               ) # fluidRow
            ), # tabItem 1
         createBuildModelTab(),
         createVideoTab()
         ) # tabItems

   return(body)

} # createBody
#------------------------------------------------------------------------------------------------------------------------
ui <- dashboardPage(

  dashboardHeader(title=sprintf("trena %s", getTargetGene(trenaProject))),
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
          overflow: hidden;
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
       setTargetGene(trenaProject, newGene)
       showGenomicRegion(session, getGeneRegion(trenaProject, flankingPercent=20))
       printf("new targetGene: %s", getTargetGene(trenaProject))
       state$tbl.enhancers <- getEnhancers(trenaProject)
       state$tbl.dhs <- getEncodeDHS(trenaProject)
       state$tbl.transcripts <- getTranscriptsTable(trenaProject)
       shinyjs::html(selector=".logo", html=sprintf("trena %s", newGene), add=FALSE)
       loadAndDisplayRelevantVariants(session, newGene)
       })


   output$igvShiny <- renderIgvShiny({
      options <- list(genomeName="hg38",
                      initialLocus=getGeneRegion(trenaProject, flankingPercent=20),
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
      printf("--- calling redrawIgvWidget after igv toggle button")
      redrawIgvWidget(session)
      redrawModelDataTable()
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
      printf("--- calling redrawIgvWidget after table toggle button")
      redrawIgvWidget(session)
      redrawModelDataTable()
      })

} # setupIgvAndTableToggling
#------------------------------------------------------------------------------------------------------------------------
redrawModelDataTable <- function()
{
      # columns.adjust not actually needed, except perhaps if the column names have changed
      # jQuery.string <- "$('#table table.dataTable[id]').DataTable().columns.adjust().draw();"
      #

      # it is not at all clear where and how [id] is resolved by the js interpreter in the browser

  jQuery.string <- "$('#table table.dataTable[id]').DataTable().draw();"
  shinyjs::runjs(jQuery.string)

} # redrawModelDataTable
#------------------------------------------------------------------------------------------------------------------------
mapToChromLoc <- function(regionName)
{
   roi <- getTargetGene(trenaProject)   # a safe fallback

   if(regionName == "traditionalPromoter"){
      tbl.transcripts <- getTranscriptsTable(trenaProject)
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
      chrom <- state$tbl.enhancers$chrom[1]
      start <- min(state$tbl.enhancers$start) - 10000
      end   <- max(state$tbl.enhancers$end) + 10000
      roi <- sprintf("%s:%d-%d", chrom, start, end)
      }

   return(roi)

} # mapToChromLoc
#------------------------------------------------------------------------------------------------------------------------
buildFootprintModel <- function(upstream, downstream)
{
   tbl.gene <- subset(getTranscriptsTable(trenaProject), moleculetype=="gene")[1,]
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
   mtx <- loadExpressionData(trenaProject, "FilteredLengthScaledTPM8282018-vsn")

   build.spec <- list(title=sprintf("%s model %d", getTargetGene(trenaProject), 1),
                      type="footprint.database",
                      regions=tbl.regions,
                      geneSymbol=getTargetGene(trenaProject),
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

   fpBuilder <- FootprintDatabaseModelBuilder("hg38", getTargetGene(trenaProject), build.spec, quiet=TRUE)
   x <- build(fpBuilder)
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
   # supportedTracks <- c("enhancers",
   #                    "dhs",
   #                     "gwas",
   #                     "wgs.variants",
   #                     "dataframeVCF")

   printf("--- displayTrack('%s')", trackName)

   if(grepl(".variants", trackName, ignore.case=TRUE, fixed=TRUE))
      displayBedTrack(session, trackName)

   if(grepl("gwas", trackName, ignore.case=TRUE))
      displayBedTrack(session, trackName)

   if(trackName == "enhancers")
      displayEnhancersTrack(session)

   if(trackName == "dhs")
      displayEncodeDhsTrack(session)

   #switch(trackName,
   #       enhancers    = {displayEnhancersTrack(session)},
   #       dhs          = {displayEncodeDhsTrack(session)},
   #       gwas         = {displayGWASTrack(session, trackName)},
   #       wgs.variants = {displayGWASTrack(session, trackName)}
   #       )

} # displayTrack
#------------------------------------------------------------------------------------------------------------------------
displayEnhancersTrack <- function(session)
{
   tbl.enhancers <- getEnhancers(trenaProject)
   tbl.tmp <- tbl.enhancers[, c("chrom", "start", "end", "combinedScore")]
   loadBedGraphTrack(session, "GeneHancer", tbl.tmp, color="black", trackHeight=25, autoscale=FALSE, min=0, max=20)

} # displayEnhancersTrack
#------------------------------------------------------------------------------------------------------------------------
displayEncodeDhsTrack <- function(session)
{
   tbl.dhs <- getEncodeDHS(trenaProject)
   tbl.tmp <- tbl.dhs[, c("chrom", "chromStart", "chromEnd", "score")]
   colnames(tbl.tmp) <- c("chrom", "start", "end", "score")
   loadBedGraphTrack(session, "DHS", tbl.tmp, color="black", trackHeight=25, autoscale=TRUE)

} # displayEncodeDhsTrack
#------------------------------------------------------------------------------------------------------------------------
displayGWASTrack <- function(session, trackName)
{
   variantDatasetNames <- getVariantDatasetNames(trenaProject)
   printf("want to see trackName (%s) among variantDatasetNames", trackName)
   print(paste(variantDatasetNames, collapse=", "))

   if(trackName %in% variantDatasetNames){
      printf("  --- found %s to display", trackName)
      printf(" want to restrict variant table to this region: %s", state$chromLocRegion)
      loc <- parseChromLocString(state$chromLocRegion)
      tbl.variant <- getVariantDataset(trenaProject, trackName)
      if(trackName == "GWAS.snps"){
         tbl.variant <- tbl.variant[, c("chrom", "start", "end", "pScore")]
         colnames(tbl.variant) <- c("chrom", "start", "end", "value")
        }
      tbl.variant <- subset(tbl.variant, chrom==loc$chrom & start >= loc$start & end <= loc$end)
      printf("%d regions in %s", nrow(tbl.variant), trackName)
      showNotification(sprintf("%s: %d genomic features", trackName, nrow(tbl.variant)))
      if(nrow(tbl.variant) > 0)
         loadBedGraphTrack(session, trackName, tbl.variant, color="red", trackHeight=25, autoscale=TRUE) # , min=0, max=10)
      } # trackName found in variant data sets

} # displayGWASTrack
#------------------------------------------------------------------------------------------------------------------------
displayBedTrack <- function(session, trackName)
{
   variantDatasetNames <- getVariantDatasetNames(trenaProject)
   printf("want to see trackName (%s) among variantDatasetNames", trackName)
   print(paste(variantDatasetNames, collapse=", "))

   if(trackName %in% variantDatasetNames){
      printf("  --- found %s to display", trackName)
      printf(" want to restrict variant table to this region: %s", state$chromLocRegion)
      loc <- parseChromLocString(state$chromLocRegion)
      tbl.variant <- getVariantDataset(trenaProject, trackName)
      tbl.variant <- subset(tbl.variant, chrom==loc$chrom & start >= loc$start & end <= loc$end)
      printf("%d regions in %s", nrow(tbl.variant), trackName)
      showNotification(sprintf("%s: %d genomic features", trackName, nrow(tbl.variant)))
      if(nrow(tbl.variant) > 0)
         loadBedTrack(session, trackName, tbl.variant, color="red", trackHeight=25)
      } # trackName found in variant data sets

} # displayBedTrack
#------------------------------------------------------------------------------------------------------------------------
setupDisplayRegion <- function(session, input, output)
{
   observeEvent(input$displayGenomicRegion, {
      requestedRegion <- input$displayGenomicRegion
      if(requestedRegion == "") return()
      printf(" displayRegion: %s", requestedRegion)
      margin <- 5000
      loc.string <-switch(requestedRegion,
                          fullEnhancerRegion = {getGeneEnhancersRegion(trenaProject, 10)},
                          fullGeneRegion = {getGeneRegion(trenaProject, 20)})
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

   observe({
      x <- input$sidebarCollapsed;
      redrawIgvWidget(session)
      })

   observe({
       currentName <- input$modelNameTextInput
       expressionMatrixName <- input$expressionSet
       allInputsSpecified <- nchar(currentName) >= 1 & nchar(expressionMatrixName) > 2
       if(allInputsSpecified)
          shinyjs::enable("buildModelButton")
       else
          shinyjs::disable("buildModelButton")
       })

   observeEvent(input$buildModelButton, {
      getGenomicRegion(session)
      state$tbl.chipSeq <- NULL   # this will force a fresh database query after model is built
      shinyjs::html(id="console", html="", add=FALSE)  # clear the console
      tryCatch({
          withCallingHandlers({buildModel(session, input, output);
                               model.count <- length(state$models)
                               new.model.name <- names(state$models)[model.count]
                               new.table <- state$models[[model.count]]$model
                               displayModel(session, input, output, new.table, new.model.name)
                               updateTabItems(session, "sidebarMenu", select="igvAndTable")
                               },
             message=function(m){
             shinyjs::html(id="console", html=m$message, add=TRUE)
             })
          }, error=function(e){
               msg <- e$message
               print(msg)
               showModal(modalDialog(title="trena model building error", msg))
               }) # tryCatch
      }) # observe buildModelButton

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

   tbl.gene <- subset(state$tbl.transcripts, moleculetype=="gene")[1,]
   strand <- tbl.gene$strand
   tss <- tbl.gene$start
   if(strand == "-")
      tss <- tbl.gene$endpos

   run.trenaSGM(trenaProject,
                model.name,
                chrom.loc$chrom, chrom.loc$start, chrom.loc$end,
                tss,
                expressionMatrixName,
                tracks.to.intersect.with,
                footprint.database.names,
                motifMapping)

} # buildModel
#------------------------------------------------------------------------------------------------------------------------
run.trenaSGM <- function(trenaProject,
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
   #browser()
   tbl.regions <- buildRegionsTable(tracks.to.intersect.with, chromosome, start.loc, end.loc,
                                    state$tbl.enhancers, state$tbl.dhs)
   #tbl.regions <- data.frame(chrom=chromosome, start=start.loc, end=end.loc, stringsAsFactors=FALSE)
   mtx <- getExpressionMatrix(trenaProject, expression.matrix.name)
   geneSymbol <- getTargetGene(trenaProject)
   for(r in seq_len(nrow(tbl.regions))){
      message(sprintf("tbl.regions, width: %d", with(tbl.regions[r,], end - start)))
      }

   build.spec <- list(title=model.name,
                      type="footprint.database",
                      regions=tbl.regions,
                      geneSymbol=getTargetGene(trenaProject),
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
     # save(build.spec, file=sprintf("%s.buildSpec.RData", model.name))

     #------------------------------------------------------------
     # use the above build.spec: a small region, high correlation
     # required, MotifDb for motif/tf lookup
     #------------------------------------------------------------

   fpBuilder <- FootprintDatabaseModelBuilder("hg38", getTargetGene(trenaProject), build.spec, quiet=FALSE)
   x <- build(fpBuilder)
   message(sprintf("back from build, top 10 tfs:"))
   message(head(x$model, n=10))
   #save(x, file=sprintf("%s.results.RData", model.name))
   state$models[[model.name]] <- x

} # run.trenaSGM
#------------------------------------------------------------------------------------------------------------------------
buildRegionsTable <- function(tracks.to.intersect.with, chromosome, start.loc, end.loc, tbl.enhancers, tbl.dhs)
{
   if(tracks.to.intersect.with == "allDNAForFootprints")
      return(data.frame(chrom=chromosome, start=start.loc, end=end.loc, stringsAsFactors=FALSE))

   #browser()
   gr.region <- GRanges(seqnames=chromosome, IRanges(start.loc, end.loc))
   gr.enhancers <- GRanges(tbl.enhancers)
   gr.dhs <- GRanges(tbl.dhs)
     # strip off metadata so that the regions can be combined
   mcols(gr.region) <- NULL
   mcols(gr.enhancers) <- NULL
   mcols(gr.dhs) <- NULL
   gr.enhancers.or.dhs <- c(gr.enhancers, gr.dhs)
   gr.enhancers.and.dhs <- intersect(gr.enhancers, gr.dhs)

   #browser()
   tbl.out  <- switch(tracks.to.intersect.with,
                  genehancer = {
                     tbl.ov <- as.data.frame(findOverlaps(gr.enhancers, gr.region, type="any"))
                     indices <- unique(tbl.ov[,1])
                     as.data.frame(gr.enhancers[indices])
                     },
                  encodeDHS = {
                     tbl.ov <- as.data.frame(findOverlaps(gr.dhs, gr.region, type="any"))
                     indices <- unique(tbl.ov[,1])
                     as.data.frame(gr.dhs[indices])
                     },
                  geneHancerOrEncode = {
                     tbl.ov <- as.data.frame(findOverlaps(gr.enhancers.or.dhs, gr.region, type="any"))
                     indices <- unique(tbl.ov[,1])
                     as.data.frame(gr.enhancers.or.dhs[indices], row.names=NULL)
                     },
                  geneHancerAndEncode = {
                     tbl.ov <- as.data.frame(findOverlaps(gr.enhancers.and.dhs, gr.region, type="any"))
                     indices <- unique(tbl.ov[,1])
                     as.data.frame(gr.enhancers.and.dhs[indices], row.names=NULL)
                     })

   tbl.out <- tbl.out[, 1:3]
   colnames(tbl.out) <- c("chrom", "start", "end")
   tbl.out$chrom <- as.character(tbl.out$chrom)
   tbl.out$start <- as.numeric(tbl.out$start)
   tbl.out$end <- as.numeric(tbl.out$end)

   return(tbl.out)

} # buildRegionsTable
#------------------------------------------------------------------------------------------------------------------------
test_buildRegionsTable <- function()
{
   library(RUnit)
   variables.loaded <- load("buildRegionsTable.sampleData.RData")
   stopifnot(all(c("chromosome", "start.loc", "end.loc", "tbl.dhs", "tbl.enhancers") %in% variables.loaded))
   intersection.options <- c("genehancer", "encodeDHS", "geneHancerOrEncode",
                             "geneHancerAndEncode", "allDNAForFootprints")

   tbl.regions <- buildRegionsTable("allDNAForFootprints", chromosome, start.loc, end.loc, tbl.enhancers, tbl.dhs)
   checkEquals(colnames(tbl.regions), c("chrom", "start", "end"))
   checkEquals(unlist(lapply(tbl.regions, class), use.names=FALSE), c("character", "numeric", "numeric"))
   checkEquals(dim(tbl.regions), c(1, 3))

   tbl.regions <- buildRegionsTable("genehancer", chromosome, start.loc, end.loc, tbl.enhancers, tbl.dhs)
   checkEquals(colnames(tbl.regions), c("chrom", "start", "end"))
   checkEquals(unlist(lapply(tbl.regions, class), use.names=FALSE), c("character", "numeric", "numeric"))
   checkEquals(dim(tbl.regions), c(32, 3))

   tbl.regions <- buildRegionsTable("encodeDHS", chromosome, start.loc, end.loc, tbl.enhancers, tbl.dhs)
   checkEquals(colnames(tbl.regions), c("chrom", "start", "end"))
   checkEquals(unlist(lapply(tbl.regions, class), use.names=FALSE), c("character", "numeric", "numeric"))
   checkEquals(dim(tbl.regions), c(2541, 3))

   tbl.regions <- buildRegionsTable("geneHancerOrEncode", chromosome, start.loc, end.loc, tbl.enhancers, tbl.dhs)
   checkEquals(colnames(tbl.regions), c("chrom", "start", "end"))
   checkEquals(unlist(lapply(tbl.regions, class), use.names=FALSE), c("character", "numeric", "numeric"))
   checkEquals(dim(tbl.regions), c(2573, 3))

   tbl.regions <- buildRegionsTable("geneHancerAndEncode", chromosome, start.loc, end.loc, tbl.enhancers, tbl.dhs)
   checkEquals(colnames(tbl.regions), c("chrom", "start", "end"))
   checkEquals(unlist(lapply(tbl.regions, class), use.names=FALSE), c("character", "numeric", "numeric"))
   checkEquals(dim(tbl.regions), c(341, 3))

} # test_buildRegionsTable
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
   max.model.rows <- MAX.TF.COUNT
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
   if(length(tf.names) > MAX.TF.COUNT) tf.names <- tf.names[1:MAX.TF.COUNT]
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
      if(is.null(state$tbl.chipSeq)){
         showNotification("retrieving ChIP-seq data from database...")
         tbl.chipSeq <- with(chrom.loc, getChipSeq(trenaProject, chrom, start, end,  tf.names))
         #save(tbl.chipSeq, file="tbl.chipSeq.RData")
         state$tbl.chipSeq <- tbl.chipSeq
         hit.count <- nrow(state$tbl.chipSeq)
         printf("chipSeq hit.count: %d", hit.count)
         showNotification(sprintf("ChIP-seq hits across all TFs in this model: %d", hit.count))
         }
      tbl.hits <- subset(state$tbl.chipSeq, tf == tf.name)
      showNotification(sprintf("ChIP-seq hits for %s: %d", tf.name, nrow(tbl.hits)), type="message")
      if(nrow(tbl.hits) > 0){
         tbl.tmp <- tbl.hits[, c("chrom", "start", "endpos", "name")]
         colnames(tbl.tmp) <- c("chrom", "start", "end", "name")
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
   variant.filenames <- getVariantDatasetNames(trenaProject)
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
shinyOptions=list(launch.browser=FALSE)
if(Sys.info()[["nodename"]] == "riptide.local"){
   shinyOptions <- list(host="0.0.0.0", launch.browser=TRUE)
   }
app <- shinyApp(ui, server, options=shinyOptions)
#app <- shinyApp(ui, server)
runApp(app)

