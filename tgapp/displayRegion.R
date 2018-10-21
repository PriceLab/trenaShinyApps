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
