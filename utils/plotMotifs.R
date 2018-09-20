library(motifStack)
#------------------------------------------------------------------------------------------------------------------------
plotMotifs <- function(motifs, group.size=5)
{
   motif.count <- length(motifs)

   stopifnot(motif.count > 0)

   if(motif.count == 1){
      pcm <- new("pcm", mat=motifs[[1]], name=names(motifs))
      plot(pcm)
      return(motif.count)
      }

   if(motif.count <= group.size){
      motifStack(lapply(names(motifs), function(mName) new("pfm", motifs[[mName]], name=mName)))
      return(motif.count)
      }

   group.count <- motif.count %/% group.size
   remainder <- motif.count %% group.size
   if(remainder > 0)
      group.count <- group.count + 1

   for(group in seq_len(group.count)){
      index.start <- 1 + ((group - 1)  * group.size)
      index.end   <- index.start + (group.size-1)
      if(group == group.count & remainder > 0)
         index.end <- index.start + (remainder - 1)
      printf("group %d, index.start:index.end %d:%d", group, index.start, index.end)
      motif.set <- motifs[index.start: index.end]
      if(length(motif.set) == 1){
         printf("--- final motif of length 1, duplicating in order to stack....")
         motif.set <- rep(motif.set, 2)
         }
      quartz()
      motifStack(lapply(names(motif.set), function(mName) new("pfm", motif.set[[mName]], name=mName)))
      }

} # plotMotifs
#------------------------------------------------------------------------------------------------------------------------

