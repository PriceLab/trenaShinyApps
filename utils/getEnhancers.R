getEnhancerRegions <- function(targetGene)
{
   base.dir <- "~/s/data/public/GeneHancer/v4.7"
   path <- file.path(base.dir, "enhancer_elite_ids.txt")
   tbl.elite <- read.table(path, sep="\t", as.is=TRUE, header=TRUE)            #  243281      6
   path <- file.path(base.dir, "enhancer_gene_scores.txt")
   tbl.geneScoresAll <- read.table(path, sep="\t", as.is=TRUE, header=TRUE)  #  934287      9
      # dim(tbl.elite); dim(tbl.geneScoresAll)
      # [1]  250733      7
      # [1] 1086552      9

   tbl.scores <- subset(tbl.geneScoresAll, symbol==targetGene)
   clusters <- tbl.scores$cluster_id
   tbl.enhancers <- subset(tbl.elite, cluster_id %in% clusters)[, c("chr", "enhancer_start", "enhancer_end")]
   colnames(tbl.enhancers) <- c("chrom", "start", "end")
   tbl.enhancers$chrom <- paste("chr", tbl.enhancers$chrom, sep="")
   rownames(tbl.enhancers) <- NULL

   tbl.enhancers

} # getEnhancerRegions
#------------------------------------------------------------------------------------------------------------------------
