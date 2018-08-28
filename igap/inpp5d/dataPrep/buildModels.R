# buildModels.R
#------------------------------------------------------------------------------------------------------------------------
library(trenaSGM)
#------------------------------------------------------------------------------------------------------------------------
load("ratio.202.RData")  # ratio.202
load("../shinyApp/data/mtx.withDimers.cer.ros.tcx.RData")   # mtx.cer, mts.ros, mtx.tcx
#------------------------------------------------------------------------------------------------------------------------
genome <- "hg38"
targetGene <- "INPP5D"
chromosome <- "chr2"
tss <- 233059795
upstream <- 10000
downstream <- 10000
#------------------------------------------------------------------------------------------------------------------------
if(!exists("sgm"))
   sgm <- trenaSGM(genome, targetGene)
#------------------------------------------------------------------------------------------------------------------------
default.model <- function()
{
   start <- tss - downstream
   end   <- tss + upstream

   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)

   spec.0 <- list(title="fp.10kup.10kdown.04",
                  type="footprint.database",
                  geneSymbol=targetGene,
                  regions=tbl.regions,
                  tss=tss,
                  matrix=mtx.tcx,
                  db.host="khaleesi.systemsbiology.net",
                  databases=list("brain_hint_20"),
                  motifDiscovery="builtinFimo",
                  annotationDbFile=dbfile(org.Hs.eg.db),
                  tfPool=allKnownTFs(),
                  tfMapping="MotifDB",
                  tfPrefilterCorrelation=0.4,
                  orderModelByColumn="rfScore",
                  solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))


   shared.samples <- intersect(colnames(mtx.tcx), names(ratio.202))
   length(ratio.202)       # 276
   length(shared.samples)  # 252
   ncol(mtx.tcx)           # 264

   ratios <- ratio.202[shared.samples]
   mtx.tcx.trimmed <- mtx.tcx[, shared.samples]
   dim(mtx.tcx.trimmed)    # 16167 262

   inpp5d.0 <- mtx.tcx.trimmed["INPP5D",]
   inpp5d.1 <- inpp5d.0 * ratios
   inpp5d.2 <- inpp5d.0 * (1 - ratios)

      # check the arithmetic.  .1 and .2 should sum to .0
   stopifnot(max(fivenum(inpp5d.1 + inpp5d.2 - inpp5d.0)) < 1e-12)

   mtx.tcx.trimmed.transcript.1 <- mtx.tcx.trimmed
   mtx.tcx.trimmed.transcript.1["INPP5D",] <- inpp5d.1

   mtx.tcx.trimmed.transcript.2 <- mtx.tcx.trimmed
   mtx.tcx.trimmed.transcript.2["INPP5D",] <- inpp5d.2

      # hand check treatment of max ratio
   max.ratio <- fivenum(ratios)[5]
   checkEqualsNumeric(ratios[names(max.ratio)], 0.9066667, tol=1e-6)
   checkEqualsNumeric(mtx.tcx["INPP5D", names(max.ratio)], -0.3898564, tol=1e-6)
   checkEqualsNumeric(mtx.tcx.trimmed.transcript.1["INPP5D", names(max.ratio)], -0.3534698, tol=1e-6)
   checkEqualsNumeric(mtx.tcx.trimmed.transcript.2["INPP5D", names(max.ratio)], -0.0363866, tol=1e-6)

   spec.1 <- spec.0
   spec.2 <- spec.0

   spec.1$matrix <- mtx.tcx.trimmed.transcript.1
   spec.2$matrix <- mtx.tcx.trimmed.transcript.2

   strategies <- list(default=spec.0, transcript.1=spec.1, transcript.2=spec.2)
   x <- calculate(sgm, strategies)

   tbl.summary <- summarizeModels(sgm, orderBy="rfScore", maxTFpredictors=8)

} # default.model
#------------------------------------------------------------------------------------------------------------------------
