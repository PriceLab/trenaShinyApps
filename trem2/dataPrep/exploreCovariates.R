covariates.file <- "AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_Covariates_Flowcell_1.csv"
tbl.co <- read.table(covariates.file, sep=",", as.is=TRUE, header=TRUE)
table(tbl.co$Diagnosis)
dim(tbl.co)
if(!exists("mtx.tcx"))
   load("mtx.withDimers.cer.ros.tcx.RData")
all(colnames(mtx.tcx) %in% tbl.co$ID)

patients.AD <- subset(tbl.co, Diagnosis=="AD")$ID
length(patients.AD)     # 84
patients.AD.tcx <- intersect(patients.AD, colnames(mtx.tcx))
length(patients.AD.tcx) # 80

patients.ctl <- subset(tbl.co, Diagnosis=="Control")$ID
length(patients.ctl)     # 80
patients.ctl.tcx <- intersect(patients.ctl, colnames(mtx.tcx))
length(patients.ctl.tcx) # 73

trem2.AD <- mtx.tcx["TREM2", patients.AD.tcx]
trem2.ctl <- mtx.tcx["TREM2", patients.ctl.tcx]
fivenum(trem2.AD)
fivenum(trem2.ctl)
t.test(trem2.AD, trem2.ctl)

tfec.AD <- mtx.tcx["TFEC", patients.AD.tcx]
tfec.ctl <- mtx.tcx["TFEC", patients.ctl.tcx]
fivenum(tfec.AD)
fivenum(tfec.ctl)
t.test(tfec.AD, tfec.ctl)  # 0.04
boxplot(tfec.AD, tfec.ctl)

tf <- "TFEC"

tbl.AD <- as.data.frame(t(mtx.tcx[c("TREM2", tf), patients.AD.tcx]))
tbl.ctl <- as.data.frame(t(mtx.tcx[c("TREM2", tf), patients.ctl.tcx]))

patients.both <- c(patients.AD.tcx, patients.ctl.tcx)
colors <- rep("blue", length(patients.both))
names(colors) <- patients.both
colors[patients.ctl.tcx] <- "red"
colors <- as.character(colors)

tf <- "TFEC"
tf <- "ELK3"
tf <- "SPI1"
target <- "TREM2"
plot(mtx.tcx[target, patients.both], mtx.tcx[tf, patients.both], col=colors)

model.AD <- lm(sprintf("%s ~ %s", target, tf), data=tbl.AD)
model.ctl <- lm(sprintf("%s ~ %s", target, tf), data=tbl.ctl)
abline(model.AD, col="red")
abline(model.ctl, col="blue")

colors <- rep("blue", nrow(tbl.pheno))
colors[females] <- "red"

tf <- "PLEK"
tf <- "LYL1"
tf <- "IKZF1"
tf <- "IRF5"

cor(mtx.tcx[tf,], mtx.tcx["TREM2",]) # 0.78
cor(mtx.tcx["TREM2", patients.AD.tcx], mtx.tcx[tf, patients.AD.tcx])   # 0.81
cor(mtx.tcx["TREM2", patients.ctl.tcx], mtx.tcx[tf, patients.ctl.tcx]) # 0.75

plot(mtx.tcx[tf,], mtx.tcx["TREM2",], col=colors)
tbl.AD <- as.data.frame(t(mtx.tcx[c("TREM2", tf), patients.AD.tcx]))
tbl.ctl <- as.data.frame(t(mtx.tcx[c("TREM2", tf), patients.ctl.tcx]))
model.AD <- lm(sprintf("%s ~ %s", target, tf), data=tbl.AD)
model.ctl <- lm(sprintf("%s ~ %s", target, tf), data=tbl.ctl)
abline(model.AD, col="red")
abline(model.ctl, col="blue")

