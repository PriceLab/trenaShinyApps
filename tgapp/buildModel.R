# buildModel.R
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
             fluidRow(column(width=1, offset=3, actionButton("buildModelButton", "Build model")))
             ) # tabItem

   return(tab)

} # createBuildModelTab
#------------------------------------------------------------------------------------------------------------------------
setupBuildModel <- function(session, input, output)
{
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
      tf.names <- new.table$gene
      new.table <- new.table[, -1]
      if("lassoPValue" %in% colnames(new.table))
          new.table <- roundNumericColumns(new.table, 2, "lassoPValue")
      rownames(new.table) <- tf.names
      max.model.rows <- 20
      if(nrow(new.table) > max.model.rows)
         new.table <- new.table[1:max.model.rows,]
      output$table = DT::renderDataTable(new.table,
                                         width="800px",
                                         class='nowrap display',
                                         extensions='FixedColumns',
                                         options=list(scrollX=TRUE, dom='t', pageLength=100,
                                                      fixedColumns=list(leftColumns=1)))
      browser()
      xyz <- 99
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
      # no search for overlaps just jet: ignore "tracks.to.intersect.with"
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
                      tfPrefilterCorrelation=0.1,
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      orderModelByColumn="pearsonCoeff",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))
     save(build.spec, file=sprintf("%s.buildSpec.Rdata", model.name))

     #------------------------------------------------------------
     # use the above build.spec: a small region, high correlation
     # required, MotifDb for motif/tf lookup
     #------------------------------------------------------------

   fpBuilder <- FootprintDatabaseModelBuilder(genomeName, targetGene$hugo, build.spec, quiet=TRUE)
   x <- build(fpBuilder)
   printf("back from build, top 10 tfs:")
   print(head(x$model, n=10))
   state$models[[model.name]] <- x


} # run.trenaSGM
#------------------------------------------------------------------------------------------------------------------------
test.run.trenaSGM <- function()
{
   printf("--- test.run.trenaSGM")
   run.trenaSGM(trenaGene,
                model.name,
                chromosome, start.loc, end.loc,
                tss,
                expression.matrix.name,
                tracks.to.intersect.with,
                footprint.database.names)


} # test.run.trenaSGM
#------------------------------------------------------------------------------------------------------------------------
