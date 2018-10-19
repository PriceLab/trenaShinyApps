library(shiny)
library(shinydashboard)
library(shinyjs)
library(later)
library(igvShiny)
library(DT)
#------------------------------------------------------------------------------------------------------------------------
#source("serverFunctions.R")
#source("addTracks.R")
#source("displayRegion.R")
#source("buildModel.R")
#------------------------------------------------------------------------------------------------------------------------
library(trenaSGM)
library(TrenaGene)
library(TrenaGeneDataPlacenta)
library(TrenaGeneDataSkin)
#------------------------------------------------------------------------------------------------------------------------
state <- new.env(parent=emptyenv())
state$models <- list()
#------------------------------------------------------------------------------------------------------------------------
genomeName <- "hg38"
#targetGene <- list(hugo="PSG1", ensembl="ENSG00000231924", combined="ENSG00000231924|PSG1")
#dataDirectory <- system.file(package="TrenaGenePlacentaData", "extdata")
#psg1.data <- TrenaGeneDataPlacenta("PSG1")
#trenaGene <- TrenaGene(psg1.data)

targetGene <- list(hugo="PGF", ensembl="ENSG00000119630", combined="ENSG00000231924|PGF")
pgf.data <- TrenaGeneDataPlacenta("PGF")

targetGene <- list(hugo="COL1A1", ensembl="ENSG00000108821", combined="ENSG00000108821|COL1A1")
geneData <- TrenaGeneDataSkin(targetGene$hugo)
trenaGene <- TrenaGene(geneData)
expressionMatrixNames <- getExpressionMatrixNames(trenaGene)
tbl.enhancers <- getEnhancers(trenaGene)
tbl.dhs <- getEncodeDHS(trenaGene)
footprintDatabases <- getFootprintDatabaseNames(getGeneData(trenaGene))
tbl.transcripts <- getTranscriptsTable(trenaGene)

#------------------------------------------------------------------------------------------------------------------------
# calculate the initial igv region, using the (typically? always? single transcripted designated as "gene"
#------------------------------------------------------------------------------------------------------------------------
tbl.gene <- subset(tbl.transcripts, moleculetype=="gene")[1,]
start <- tbl.gene$start
end   <- tbl.gene$endpos
chrom <- tbl.gene$chr
padding <- round(0.25 * (end-start))
initial.chromLoc <- sprintf("%s:%d-%d", chrom, start-padding, end+padding)
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
        selectInput("addTrack", "Add Track:",
                    c("",
                     "Enhancers"="enhancers",
                     "DHS open chromatin (encode, clustered)"="dhs"
                     )),
        selectInput("displayGenomicRegion", "Display Genomic Region:",
                    c("",
                      "Full Gene" = "fullGeneRegion",
                      "Full Enhancers Region" = "fullEnhancerRegion"))
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
                column(2, offset=1, checkboxGroupInput("footprintDatabases", "Footprint Databases",
                                            footprintDatabases, selected=footprintDatabases[1])),
                column(2, radioButtons("intersectWithRegions", "Intersect footprints with:",
                                             c("GeneHancer" = "genehancer",
                                               "Encode DHS" = "encodeDHS",
                                               "GeneHancer OR Encode DHS" = "geneHancerOrEncode",
                                               "GeneHancer AND Encode DHS" = "geneHancerPlusEncode",
                                               "Use all footprints in region" = "allDNAForFootprints"),
                                            selected="allDNAForFootprints")),
                column(3, selectInput("expressionSet", "With this gene expression dataset:",
                                      expressionMatrixNames))),
             br(),
             fluidRow(column(width=3, offset=3, textInput("modelNameTextInput", "Model name", width=240))),
             fluidRow(column(width=1, offset=3, actionButton("buildModelButton", "Build model"))),
             fluidRow(column(width=1, offset=3, textOutput("buildInProgressTextArea")))
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
               column(width=9, id="igvColumn", igvShinyOutput('igvShiny', height="800px")),
               column(width=3, id="dataTableColumn",
                      selectInput("modelSelector", NULL,  "- no models yet -"),
                      DTOutput("table"),
                      br(),
                      wellPanel(
                        h4("TF selection displays:"),
                        selectInput("selectRowAction", NULL,  c("Footprints", "ChIP-seq hits",
                                                                sprintf("Plot expression, TF vs %s", targetGene$hugo))))
                      ) # column
               ) # fluidRow
            ), # tabItem 1
         createBuildModelTab()
         ) # tabItems

   return(body)

} # createBody
#------------------------------------------------------------------------------------------------------------------------
ui <- dashboardPage(

  dashboardHeader(title=sprintf("trena %s", targetGene$hugo)),
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

   output$igvShiny <- renderIgvShiny({
      options <- list(genomeName="hg38",
                      initialLocus=initial.chromLoc,
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
      printf("row selected: %d", selectedTableRow)
      #printf("car selected: %s", rownames(mtcars)[selectedTableRow])
      tfSelectionAction <- isolate(input$selectRowAction)
      printf("table row clicked, - actionRequested: %s", tfSelectionAction);
      current.model.name <- isolate(input$modelSelector)
      if(current.model.name %in% ls(state$models)){
         tf.names <- state$models[[current.model.name]]$model$gene
         tf.name <- tf.names[selectedTableRow]
         printf("tf: %s", tf.name)
         }
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

   later(function(){
            state$chromLocRegion <- initial.chromLoc
            showGenomicRegion(session, initial.chromLoc)
            }, 2)

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
setupDisplayRegion <- function(session, input, output)
{
   observeEvent(input$displayGenomicRegion, {
      requestedRegion <- input$displayGenomicRegion
      if(requestedRegion == "") return()
      printf(" displayRegion: %s", requestedRegion)
      margin <- 5000
      switch(requestedRegion,
             fullEnhancerRegion = {chrom <- tbl.enhancers$chrom[1];
                start <- min(tbl.enhancers$start) + margin;
                end <- max(tbl.enhancers$start) + -margin;
                loc.string <- sprintf("%s:%d-%d", chrom, start, end)
                },
             fullGeneRegion = {loc.string = initial.chromLoc})
      showGenomicRegion(session, loc.string);
      later(function() {updateSelectInput(session, "displayGenomicRegion", selected=character(0))}, 1)
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
      buildModel(session, input, output)
      printf("back from buildModel, current models are %s", paste(names(state$models), collapse=", "))
      model.count <- length(state$models)
      new.model.name <- names(state$models)[model.count]
      new.table <- state$models[[model.count]]$model
      output$buildInProgressTextArea <- renderText("build complete")
      displayModel(session, input, output, new.table, new.model.name)
      updateTabItems(session, "sidebarMenu", select="igvAndTable")
      })

} # setupBuildModel
#------------------------------------------------------------------------------------------------------------------------
buildModel <- function(session, input, output)
{
   model.name <- sprintf("trena.model.%s", input$modelNameTextInput)
   footprint.database.names <- input$footprintDatabases
   tracks.to.intersect.with <- input$intersectWithRegions
   expressionMatrixName <- input$expressionSet
   full.roi <- state$chromLocRegion
   chrom.loc <- trena::parseChromLocString(full.roi)
   printf("  fpdb: %s", paste(footprintDatabases, collapse=", "))
   printf("   roi: %s", full.roi)
   printf("   mtx: %s", expressionMatrixName)
   printf("  intersect with: %s", paste(tracks.to.intersect.with, collapse=","))

   tbl.gene <- subset(tbl.transcripts, moleculetype=="gene")[1,]
   strand <- tbl.gene$strand
   tss <- tbl.gene$start
   if(strand == "-")
      tss <- tbl.gene$endpos

   run.trenaSGM(trenaGene,
                model.name,
                chrom.loc$chrom, chrom.loc$start, chrom.loc$end,
                tss,
                expressionMatrixName,
                tracks.to.intersect.with,
                footprint.database.names)

} # buildModel
#------------------------------------------------------------------------------------------------------------------------
run.trenaSGM <- function(trenaGene,
                         model.name,
                         chromosome, start.loc, end.loc,
                         tss,
                         expression.matrix.name,
                         tracks.to.intersect.with,
                         footprint.database.names)
{
   printf("--- entering run.trenaSGM")
      # no search for overlaps just yet: ignore "tracks.to.intersect.with"
   tbl.regions <- data.frame(chrom=chromosome, start=start.loc, end=end.loc, stringsAsFactors=FALSE)
   mtx <- loadExpressionData(trenaGene, expression.matrix.name)
   geneSymbol <- getGeneSymbol(getGeneData(trenaGene))
   printf("tbl.regions, width: %d", with(tbl.regions, end - start))

   build.spec <- list(title=model.name,
                      type="footprint.database",
                      regions=tbl.regions,
                      geneSymbol=targetGene$hugo,
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=footprint.database.names,
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(identifierType="geneSymbol"),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.2,
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      orderModelByColumn="pearsonCoeff",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))
     save(build.spec, file=sprintf("%s.buildSpec.Rdata", model.name))

     #------------------------------------------------------------
     # use the above build.spec: a small region, high correlation
     # required, MotifDb for motif/tf lookup
     #------------------------------------------------------------

   fpBuilder <- FootprintDatabaseModelBuilder(genomeName, targetGene$hugo, build.spec, quiet=FALSE)
   x <- build(fpBuilder)
   printf("back from build, top 10 tfs:")
   print(head(x$model, n=10))
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

app <- shinyApp(ui, server)
