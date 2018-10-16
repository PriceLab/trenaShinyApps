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
             fluidRow(column(width=1, offset=3, actionButton("buildModelButton", "Build model")))
             ) # tabItem

   return(tab)

} # createBuildModelTab
#------------------------------------------------------------------------------------------------------------------------
setupBuildModel <- function(session, input, output)
{
   observeEvent(input$buildModelButton, {
      #browser()
      getGenomicRegion(session)
      buildModel(session, input, output)
      })

} # setupBuildModel
#------------------------------------------------------------------------------------------------------------------------
buildModel <- function(session, input, output)
{
   footprintDatabases <- input$footprintDatabases
   printf("  fpdb: %s", paste(footprintDatabases, collapse=", "))
   printf("   roi: %s", state$chromLocRegion)

} # buildModel
#------------------------------------------------------------------------------------------------------------------------

