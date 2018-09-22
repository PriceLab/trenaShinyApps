library(GenomicRanges)
library(igvR)
library(trena)
#------------------------------------------------------------------------------------------------------------------------
chromosome <- "chr1"

if(!exists("tbl.snps")){
   load("tbl.RData");
   tbl.snps <- tbl
   shoulder <- 10000
   region.min <- min(tbl.snps$bp.hg38) - shoulder
   region.max <- max(tbl.snps$bp.hg38) + shoulder
   }

if(!exists("igv")){
   igv <- igvR()
   Sys.sleep(10)
   setGenome(igv, "hg38")
   Sys.sleep(10)
   showGenomicRegion(igv, sprintf("%s:%d-%d", chromosome, region.min, region.max))
   Sys.sleep(10)
   }

tbl.wig <- data.frame(chr=rep("chr1", nrow(tbl.snps)),
                      start=tbl.snps$bp.hg38,
                      end=tbl.snps$bp.hg38,
                      value=-log10(tbl.snps$PVAL),
                      stringsAsFactors=FALSE)
track <- DataFrameQuantitativeTrack("snp", tbl.wig, color="red", autoscale=TRUE, trackHeight=25)
displayTrack(igv, track)


if(!exists("tbl.breaks")){
   print(load("results.all.snps.RData"))
   class(results)
   tbl.breaks <- as.data.frame(mcols(results))
   dim(tbl.breaks)
   table(tbl.breaks$effect)   # 372 strong, 118 weak
   #tbl.wig <- tbl.enhancers[, c("chr", "enhancer_start", "enhancer_end", "combined_score")]
   #colnames(tbl.wig) <- c("chr", "start", "end", "value")
   #track <- DataFrameQuantitativeTrack("enhancers", tbl.wig, color="black", autoscale=TRUE, trackHeight=25)
   #displayTrack(igv, track)
   }

if(!exists("tbl.models")){
   load("/Users/paul/github/trenaSGM/inst/batch/cory/AD/tbl.models.all.RData")
   }

if(!exists("tbl.enhancers")){
   load("tbl.enhancers.umo463f.RData")
   tbl.enhancers <- tbl.enhancers.umo463f
   genes.in.region <- unique(tbl.enhancers$geneSymbol)
   printf("gene.in.region: %d", length(genes.in.region))
       # only 36 120 genes.in.region have trena models.  30 are LINC or LOC, 6 MIR, 11 RP
       # 37 remain of possible interest as genes...
   genes.in.region.with.model <- intersect(genes.in.region, tbl.models$targetGene.hgnc)
   }

if(!exists("tbl.dhs")){
   load("tbl.dhs.RData")
   dim(tbl.dhs)
   tbl.wig <- tbl.dhs[, c("chrom", "chromStart", "chromEnd", "score")]
   colnames(tbl.wig) <- c("chr", "start", "end", "value")
   dim(tbl.wig)
   label <- "dhs"
   track <- DataFrameQuantitativeTrack(label, tbl.wig, color="darkGreen", autoscale=TRUE)
   displayTrack(igv, track)
   }


if(!exists("tbl.tss")){
   trena <- Trena("hg38")
   x <- lapply(genes.in.region.with.model, function(gene) getProximalPromoter(trena, gene, 2000, 1000))
   tbl.tss <- do.call(rbind, x)
   tbl.bed <- tbl.tss[, c(2:4, 1)]
   track <- DataFrameAnnotationTrack("genes", tbl.bed, color="purple", trackHeight=25)
   displayTrack(igv, track)
   }

if(!exists("tbl.tfs.by.gene")){
   x <- lapply(genes.in.region.with.model,
               function(gene) {
                  tbl.out <- subset(tbl.models, targetGene.hgnc==gene)[1:20,]
                  tbl.out <- tbl.out[order(abs(tbl.out$spearmanCoeff), decreasing=TRUE),]
                  tbl.out
                  })
   tbl.tfs.by.gene <- do.call(rbind, x)
   dim(tbl.tfs.by.gene)  # 720 15
   }

if(!exists("tbl.broken.tfs")){
   tbl.broken.tfs <- subset(tbl.breaks, geneSymbol %in% tbl.tfs.by.gene$tf.hgnc)
   tbl.broken.tfs$pctDelta <- tbl.broken.tfs$pctAlt - tbl.broken.tfs$pctRef
   }


# now, for each gene in region with a model:
#   1) get its tfs
#   2) get tf's breaks, if any
#   3) calculate distance to gene

if(!exists("tbl.broken.models")){
   uw.snp.locs <- tbl.snps$bp.hg38
   x <- list()
   for(gene in genes.in.region.with.model){
      #gene <- "ALDH9A1"
      tfs <- subset(tbl.tfs.by.gene, targetGene.hgnc==gene)$tf.hgnc
      tbl.gene <- subset(tbl.broken.tfs, geneSymbol %in% tfs & snpPos %in% uw.snp.locs)[, c("snpPos", "REF", "ALT",  "providerId", "geneSymbol", "pctRef", "pctAlt", "pctDelta", "seqMatch")]
      tbl.gene$posDelta = tbl.gene$snpPos - subset(tbl.tss, geneSymbol==gene)$start
      tbl.gene$targetGene <- rep(gene, nrow(tbl.gene))
      dups <- which(duplicated(tbl.gene[, c("snpPos", "REF", "ALT", "providerId")]))
      printf("length of dups: %d", length(dups))
      if(length(dups) > 0){
         printf("before dup elimination: %d", nrow(tbl.gene))
         tbl.gene <- tbl.gene[-dups,]
         printf("after dup elimination: %d", nrow(tbl.gene))
         }
      tbl.thisGenesModel <- subset(tbl.tfs.by.gene, targetGene.hgnc == gene)
      tbl.thisGenesModel <- tbl.thisGenesModel[order(tbl.thisGenesModel$rfScore, decreasing=TRUE),]
      tbl.thisGenesModel$rank <- seq_len(nrow(tbl.thisGenesModel))
      tbl.thisGenesModel <- subset(tbl.thisGenesModel, tf.hgnc %in% tbl.gene$geneSymbol)
      tf.order <- match(tbl.gene$geneSymbol, tbl.thisGenesModel$tf.hgnc)
      tbl.gene$rank <- tbl.thisGenesModel$rank[tf.order]
      x[[gene]] <- tbl.gene
      }
   tbl.brokenModels <- do.call(rbind, x)
   colnames(tbl.brokenModels)[grep("geneSymbol", colnames(tbl.brokenModels))] <- "TF"

   rownames(tbl.brokenModels) <- NULL
   preferred.column.order <- c("targetGene", "TF",  "pctDelta", "snpPos", "REF", "ALT", "providerId", "pctRef", "pctAlt", "seqMatch", "posDelta", "rank")
   tbl.brokenModels <- tbl.brokenModels[, preferred.column.order]
   }

