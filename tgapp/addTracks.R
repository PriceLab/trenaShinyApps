# addTracks.R
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
