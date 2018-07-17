# reduce the meaningless precsion in all of the numerical columns of trena models
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
roundNumericColumnsInDataframe <- function(tbl, digits, pvalColumnNames=NA)
{
  tbl.pvals <- data.frame()
  tbl.main <- tbl

  if(!(all(is.na(pvalColumnNames)))){
     pval.cols <- grep(pvalColumnNames, colnames(tbl))
     stopifnot(length(pval.cols) == length(pvalColumnNames))
     tbl.pvals <- tbl[, pval.cols, drop=FALSE]
     tbl.main <- tbl[, -pval.cols, drop=FALSE]
     }

  numeric_columns <- sapply(tbl.main, mode) == 'numeric'
  tbl.main[numeric_columns] <-  round(tbl.main[numeric_columns], digits)

  if(ncol(tbl.pvals) > 0){
     tbl.pvals <- apply(tbl.pvals, 2, function(col) as.numeric(formatC(col, format = "e", digits = 2)))
     }

  tbl.out <- cbind(tbl.main, tbl.pvals)[, colnames(tbl)]
  tbl.out

} # roundNumericColumnsInDataframe
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_roundNumericColumnsInDataframe()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_roundNumericColumnsInDataframe <- function()
{
   load("~/github/projects/priceLab/alison/flt1/consoleVersion/tbl.model.test.RData")
   tbl.fixed <- roundNumericColumnsInDataframe(tbl.model.test, 3, "lassoPValue")
   checkEquals(dim(tbl.fixed), dim(tbl.model.test))
   checkEquals(as.list(tbl.fixed[1,]),
               list(gene="MXI1",
                    betaLasso=0.395,
                    lassoPValue=1.64e-27,
                    pearsonCoeff=0.675,
                    rfScore=10.298,
                    betaRidge=0.106,
                    spearmanCoeff=0.593,
                    concordance=0.518,
                    pcaMax=2.766,
                    bindingSites=17))


} # test_roundNumericColumnsInDataframe
#------------------------------------------------------------------------------------------------------------------------
