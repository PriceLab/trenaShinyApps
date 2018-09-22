library(FimoClient)
library(BSgenome.Hsapiens.UCSC.hg38)
#----------------------------------------------------------------------------------------------------
# we depend upon a properly loaded fimo server:
#    cd ~/github/fimoService/server/
#    make -f makefile.pshannon -n
#    python -i runFimoServer.py 5558 "/Users/paul/meme/bin/fimo" "../pfms/crem.human.meme"
#----------------------------------------------------------------------------------------------------
source("constants.R")
tbl.regions <- data.frame(chr=chromosome, start=snp.min.loc, end=snp.max.loc, stringsAsFactors=FALSE)
#----------------------------------------------------------------------------------------------------
FIMO_HOST <- "localhost"
FIMO_PORT <- 5558
fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=TRUE)

seq <- as.character(with(tbl.regions, getSeq(BSgenome.Hsapiens.UCSC.hg38, chrom, start, end)))
nchar(seq)  # 7M
tbl.hits <- requestMatch(fc, list(big=seq), 10e-3)
tbl.wig <- tbl.hits
tbl.wig$start <- tbl.wig$start + snp.min.loc
tbl.wig$stop <- tbl.wig$stop + snp.min.loc
tbl.wig$chr <- chromosome
tbl.wig <- tbl.wig[, c("chr", "start", "stop", "p.value")]
colnames(tbl.wig) <- c("chr", "start", "end", "value")
tbl.wig$value <- -log10(tbl.wig$value)
save(tbl.wig, file="tbl.bindingSites.crem.RData")
