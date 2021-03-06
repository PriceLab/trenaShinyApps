library(shiny)
library(DT)
library(igvShiny)
library(htmlwidgets)
library(GenomicRanges)
library(later)
library(MotifDb)
library(motifStack)
#----------------------------------------------------------------------------------------------------
targetGene <- "RIN3"
load("rin3.models.RData")
models <- rin3.models
#load("small.model.RData")
addResourcePath("tmp", "tmp") # so the shiny webserver can see, and return track files for igv
#----------------------------------------------------------------------------------------------------
currentFilters <- list()
currentShoulder <- 0
currentPwmMatchThreshold <- 90
filter.result <- list(fp=data.frame(), snp=data.frame(), fpClean=data.frame())

state <- new.env(parent=emptyenv())
state[["currentModel"]] <- models[[1]]
state[["currentModelName"]] <- names(models)[1]
state[["optionalTracks"]] <- c()

#----------------------------------------------------------------------------------------------------
modelList <- c("proximal.promoter.fp"="proximal.promoter.fp")
modelList <- c()
model.names <- names(rin3.models)
for(model.name in model.names){
   modelList <- c(modelList, eval(parse(text=sprintf('c("%s"="%s")', model.name, model.name))))
   }
if(!exists("tbl.gwas.rin3")){
   load("tbl.gwas.rin3.RData")
   tbl.snps <- tbl.gwas.rin3[, c("CHR", "BP", "BP", "SNP", "P")]
   colnames(tbl.snps) <- c("chrom", "start", "end", "id", "pval")
   }
#----------------------------------------------------------------------------------------------------
ui <- fluidPage(

  includeScript("message-handler.js"),
  tags$head(tags$link(rel = "stylesheet", type = "text/css",
                     href = "http://maxcdn.bootstrapcdn.com/font-awesome/4.2.0/css/font-awesome.min.css")
     ),

  titlePanel("RIN3 Trena Models"),

  sidebarLayout(
     sidebarPanel(
        width=3,

       # sliderInput("pwmMatchThreshold",
       #             "PWM match threshold",
       #             value = 90,
       #             min = 75,
       #             max = 100),
       selectInput("displayGenomicRegion", "Display Genomic Region:",
                   c(" - " = "noRegionYet",
                     "RIN3 (entire gene)"="RIN3",
                     "RIN3 10kb promoter"="RIN3.10kb.promoter",
                     "Full RIN3 enhancers region"="all.RIN3.enhancers")),

        selectInput("geneModel", "Choose model:",
                    modelList
                    # c("fp.2000up.200down.cor02"="fp.2000up.200down.cor02",
                    #   "fp.10kup.10kdown.cor02"="fp.10kup.10kdown.cor02",
                    #   "fp.enhancers.cor02"="fp.enhancers.cor02",
                    #   "regions.motifMatching.5kup.5kdown.pwm80.cor02"="regions.motifMatching.5kup.5kdown.pwm80.cor02",
                    #   "regions.motifMatching.2kup.200down.pwm80.cor02"="regions.motifMatching.2kup.200down.pwm80.cor02",
                    #   "regions.motifMatching.2kup.200down.pwm90.cor02"="regions.motifMatching.2kup.200down.pwm90.cor02",
                    #   "regions.motifMatching.2kup.200down.pwm95.cor02"="regions.motifMatching.2kup.200down.pwm95.cor02",
                    #   "noDNA.allTFs.cor04"="noDNA.allTFs.cor04",
                    #   "fp.enhancers.cor02.5kup.5kdown"="fp.enhancers.cor02.5kup.5kdown",
                    #   "regions.motifMatching.enhancers.5kbup.5kbdown.pwm80.cor04"="regions.motifMatching.enhancers.5kbup.5kbdown.pwm80.cor04"
                    #   )
                    ),
       selectInput("addTrack", "Add Track:",
                   c(" - ",
                     "Enhancers"="enhancers",
                     "DHS (wgEncodeClustered)"="dhs",
                     "IGAP GWAS SNPs"="snps",
                     "open chromatin"="dhs")),
       sliderInput("IGAP.snp.significance", "-log10(GWAS SNP qval)",
                    value = 2,
                    min = 0,
                    max = 12),
       sliderInput("fimo.motif.significance", "-log10(FIMO motif pval)",
                    value = 5,
                    min = 0,
                    max = 12),
       sliderInput("fimo.snp.effect", "-log10(FIMO SNP effect) - delta qval",
                    value = 2,
                    min = 0,
                   max = 12),
        actionButton("findDisruptiveSNPsButton", "Find SNPs which disrupt TF binding sites")
        ), # sidebarPanel

     mainPanel(
       tabsetPanel(type="tabs",
                   id="trenaTabs",
                   tabPanel(title="IGV",             value="igvTab",       igvShinyOutput('igvShiny')),
                   tabPanel(title="READ ME",         includeHTML("readme.html")),
                   tabPanel(title="TRN",             value="snpTableTab",  DTOutput("geneModelTable")),
                   tabPanel(title="SNP",             value="snpInfoTab",   mainPanel(
                                                                              verbatimTextOutput("tfbsText"),
                                                                              DTOutput("snpMotifTable"),
                                                                              imageOutput("logoImage"))),
                   tabPanel(title="tf/target plot",  value="plotTab",      plotOutput("xyPlot", height=800))
                   )
       ) # mainPanel
    ) # sidebarLayout
  ) # fluidPage

#--------------------------------------------------------------------------------
server <- function(input, output, session) {

   session$sendCustomMessage(type="showGenomicRegion", message=(list(region="RIN3")))
   dt.proxy <- dataTableProxy("geneModelTable")

    # define a reactive conductor. it returns a function, is
    # redefined with a change to dist or n
    # has reactive children  plot, summary, table, below, which mention it
    # thereby establishing the reactive chain.

   observeEvent(input$trackClick, {
      printf("browser snp click observed!")
      print(input$trackClick)
      motifName <- input$trackClick$id
      if(grepl("^tfbs-snp", motifName)){
         load("tbl.bindingSitesWithSnpsAndDeltaScores.RData")
         tokens <- strsplit(motifName, ":")[[1]]
         pfm.name <- tokens[2]
         snp.name <- tokens[3]
         tbl.sub <- subset(tbl.bs, motif==pfm.name & rsid==snp.name)
         tbl.sub <- tbl.sub[, c("chrom", "start", "end", "strand", "sequence.wt", "sequence.mut", "wtMutMatchDelta")]
         colnames(tbl.sub)[7] <- "10^x match loss"
         output$tfbsText <- renderText({sprintf("%s   %s", pfm.name, snp.name)})
         image.filename <- sprintf("png/%s.png", pfm.name)
         output$snpMotifTable <- renderDT(tbl.sub, options = list(scrollX = TRUE))
         output$logoImage <- renderImage({
                 list(src=image.filename,
                      contentType="image/png",
                      width=500,
                      height=500,
                      alt="alt text")},
                 deleteFile=FALSE)
         #updateTabsetPanel(session, "trenaTabs", select="snpInfoTab");
         } # if tfbs-snp
      })

   observeEvent(input$geneModel, {
      printf("geneModel event: %s", input$geneModel);
      })

   observeEvent(input$findDisruptiveSNPsButton, {
      load("sequenceMatchesForTfsInModel.RData")   # calculated with fimo using only TFs in model
      motif.pval.threshold <- isolate(session$input$fimo.motif.significance)
      snpGWAS.pval.threshold <- isolate(session$input$IGAP.snp.significance)
      snpEffect.pval.threshold <- isolate(session$input$fimo.snp.effect)
      printf("motif: %f  snp: %f   disruption delta: %f", motif.pval.threshold, snpGWAS.pval.threshold, snpEffect.pval.threshold)
      #browser()
         # tbl.matches: of model's tfs' motifs to the full sequence around gene
      tbl.matches.filtered <- subset(tbl.matches, -log10(p.value) >= motif.pval.threshold)
      printf("  %d bs below p.value threshold", nrow(tbl.matches.filtered))
      tbl.snps.filtered <- tbl.snps #subset(tbl.snps, -log10(pval) >= snpGWAS.pval.threshold)
      gr.matches <- GRanges(tbl.matches.filtered)
      gr.snps <- GRanges(tbl.snps.filtered)
      tbl.ov <- as.data.frame(findOverlaps(gr.snps, gr.matches))
      colnames(tbl.ov) <- c("snp", "bindingSite")
      tbl.bs <- cbind(tbl.matches.filtered[tbl.ov$bindingSite,], tbl.snps.filtered[tbl.ov$snp,])
      state[["tbl.bs"]] <- tbl.bs
      tfs.with.snps <- unique(tbl.bs$tf)
      printf("--- tfs.with snps: %s", paste(tfs.with.snps, collapse=", "))
      for(tf.with.snp in tfs.with.snps){
         tbl.tf <- subset(tbl.bs, tf==tf.with.snp)[, c("chrom", "start", "end", "name", "q.value", "pval", "id")]
         tbl.tf$name <- sprintf("tfbs-snp:%s:%s", tbl.tf$name, tbl.tf$id)
         temp.filename <- tempfile(tmpdir="./tmp", fileext=".bed")
         #browser()
         write.table(tbl.tf, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
         session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName=tf.with.snp,
                                              displayMode="EXPANDED",
                                              color="blue",
                                              trackHeight=40)))
         } # for tf
      }) #


   observeEvent(input$displayGenomicRegion, {
      printf("display event: %s", input$displayGenomicRegion)
      regions <- list(RIN3="chr14:92,492,870-92,704,078",
                      RIN3.10kb.promoter="chr14:92,508,719-92,521,722",
                      all.RIN3.enhancers="chr14:91,990,979-93,109,884")
      new.region.name <- input$displayGenomicRegion
      roi <- regions[[new.region.name]]
      updateTabsetPanel(session, "trenaTabs", select="igvTab");
      myFunc <- function(){
         session$sendCustomMessage(type="showGenomicRegion", message=(list(region=roi)))
         } # myFunc
      later(myFunc, 1)
      })

   observeEvent(input$addTrack, {
      newTrackName <- input$addTrack
      displayTrack(session, newTrackName)
      })

      #regions <- list(RIN3="chr14:92,492,870-92,704,078",
      #                RIN3.10kb.promoter="chr14:92,508,719-92,521,722",
      #                all.RIN3.enhancers="chr14:91,990,979-93,109,884")
      #new.region.name <- input$displayGenomicRegion
      #roi <- regions[[new.region.name]]
      #session$sendCustomMessage(type="showGenomicRegion", message=(list(region=roi)))

   observeEvent(input$showTargetGeneButton, {
     roi <- "chr14:92,492,870-92,704,078"
     session$sendCustomMessage(type="showGenomicRegion", message=(list(region=roi)))
     })

   observeEvent(input$showEnhancerRegionButton, {
     roi <- "chr14:91,990,979-93,109,884"
     session$sendCustomMessage(type="showGenomicRegion", message=(list(region=roi)))
     })

   observeEvent(input$showSNPsButton, {
     trackName <- "snp"
     printf("-- loading snps")
     displayTrack(session, trackName)
     # optionalTracks <- state[["optionalTracks"]]
     # printf("entering button toggle, optionalTracks: ")
     # print(optionalTracks);
     # if(trackName %in% optionalTracks){
     #    removeTrack(session, trackName)
     #    optionalTracks <- optionalTracks[-grep(trackName, optionalTracks)]
     #    }
     # else{
     #    displayTrack(session, trackName)
     #    optionalTracks <- c(trackName, optionalTracks)
     #    }
     #  printf("--- updating optionalTracks: ")
     #  print(optionalTracks)
     #  state[["optionalTracks"]] <- optionalTracks
      })

   observeEvent(input$showEnhancersButton, {
     printf("--- show enhancers")
     displayTrack(session, "enhancers")
     #trackName <- "enhancers"
     #optionalTracks <- state[["optionalTracks"]]
     #printf("entering button toggle, optionalTracks: ")
     #print(optionalTracks);
     #if(trackName %in% optionalTracks){
     #   removeTrack(session, trackName)
     #   optionalTracks <- optionalTracks[-grep(trackName, optionalTracks)]
     #   }
     #else{
     #   displayTrack(session, trackName)
     #   optionalTracks <- c(trackName, optionalTracks)
     #   }
     # printf("--- updating optionalTracks: ")
     # print(optionalTracks)
     # state[["optionalTracks"]] <- optionalTracks
      })

   observeEvent(input$findSnpsInModelButton, {
      trena.model <- isolate(input$model.trn)
      filters <- isolate(input$filters)
      snpShoulder <- isolate(input$snpShoulder)
      filter.result <- applyFilters(session, filters, trena.model, snpShoulder)
      })

   observeEvent(input$pwmMatchThreshold, {
      printf("  ** pwmMatchThreshold changed");
      currentPwmMatchThreshold <- input$pwmMatchThreshold
      })

   observeEvent(input$tracks, {
     trackNames <- isolate(input$tracks)
     printf("currentTracks: %s", paste(trackNames, collapse=","))
     #newTracks <- setdiff(trackNames, currentTracks)
     newTracks <-  trackNames
     printf("newTracks: %s", paste(newTracks, collapse=","))
     if(length(newTracks) > 0){
        displayTrack(session, newTracks)
        }
     currentTracks <<- sort(unique(unlist(c(currentTracks, newTracks))))
     })

   observeEvent(input$filters, {
      currentFilters <- input$filters
      print(paste(currentFilters, collapse=","))
      })

   observeEvent(input$model.trn, {
      currentModel <- input$model.trn
      print(paste(currentFilters, collapse=","))
      })

   observeEvent(input$roi, {
     currentValue = input$roi
     tss <- 88884466;   # get from MEF2C.data
     roi.string <- switch(currentValue,
                          tss_2kb = sprintf("chr5:%d-%d", tss-2000, tss+2000),
                          tss_5kb = sprintf("chr5:%d-%d", tss-5000, tss+5000),
                          studyRegion = "chr5:88391000-89322000")
     session$sendCustomMessage(type="roi", message=(list(roi=roi.string)))
     })

   observeEvent(input$geneModelTable_cell_clicked, {
      printf("--- input$geneModelTable_cell_clicked")
      selectedTableRow <- input$geneModelTable_row_last_clicked
      if(!is.null(selectedTableRow)){
         tbl.model <- state$currentModel
         tf <- tbl.model$gene[selectedTableRow]
         if(length(tf) > 0){
            printf("table row clicked: %s", tf );
            # updatePlot(session, tf, state$currentModelName);
            updateTabsetPanel(session, "trenaTabs", select="igvTab");
            printf(" hoping igvTab will be raised...")
            # service(2);
            later(function(){selectRows(dt.proxy, NULL)}, 1);
            displayBindingSites(session, input, tf) # , state$currentModelName);
            }
         } # row not null
      })

   output$xyPlot = renderPlot(
       {print(" ---- renderPlot");
        selectedTableRow <- input$geneModelTable_row_last_clicked
        tbl.model <- state$currentModel
        tf <- tbl.model$gene[selectedTableRow]
        printf("want an xyplot of %s vs. FLT1", tf)
        tbl.model <- state$currentModel
        tf <- tbl.model$gene[selectedTableRow]
        target <- "FLT1";
        tbl.pheno <- tbl.pheno[colnames(mtx),]
        males   <- which(tbl.pheno$sex_s == "male")
        females <- which(tbl.pheno$sex_s == "female")
        colors <- rep("blue", nrow(tbl.pheno))
        colors[females] <- "red"
        tbl.males <- as.data.frame(t(mtx[c(tf, target), males]))
        tbl.females <- as.data.frame(t(mtx[c(tf, target), females]))
        model.male <- lm(sprintf("%s ~ %s", target, tf), data=tbl.males)
        model.female <- lm(sprintf("%s ~ %s", target, tf), data=tbl.females)
        plot(as.numeric(mtx[tf,]), as.numeric(mtx[target,]), xlab=tf, ylab=target, col=colors)
        abline(model.male, col="blue")
        abline(model.female, col="red")
        printf("  male: ")
        print(summary(model.male))
        printf("  female: ")
        print(summary(model.female))
        })

  output$summary <- renderPrint({
    #print(input$targetGene)
    #print(input$roi)
    #print(input$filters)
    #print(input$pwmMatchThreshold)
    #print(input$model.trn)
    })

   output$geneModelTable <- renderDT(
      {model.name <- input$geneModel;
       printf("    rendering DT, presumably because input$geneModel changes, model.name: '%s'", model.name);
       printf("    names of models: %s", paste(names(models), collapse=", "))
       selection="single"
       print(lapply(models[[model.name]], dim))
       tbl.model <- models[[model.name]]$model
       printf("data.frame dimensions: %d, %d", nrow(tbl.model), ncol(tbl.model))
       #tbl.model <- switch(model.name,
       #   "Enhancers.TFClass" = tbl.model.enhancerAll.tfc,
       #   "Enhancers.MotifDb" = tbl.model.enhancerAll.mdb
       #   )
      #tbl.model <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE),]
      if(nrow(tbl.model) > 20)
         tbl.model <- tbl.model[1:20,]
      rownames(tbl.model) <- NULL
      state[["currentModel"]] <- tbl.model
      state[["currentModelName"]] <- model.name
      return(tbl.model)
      },
    selection = 'single',
    options=list(pageLength=25, dom="t"))

   output$igvShiny <- renderIgvShiny({
     igvShiny("hello shinyApp")
     })

  } # server

#----------------------------------------------------------------------------------------------------
calculateVisibleSNPs <- function(targetGene, roi, filters, snpShoulder)
{
   #return(sprintf("calculateVisibleSNPs(%s), %s, %s, %d", date(), targetGene, roi, snpShoulder))
   tbl.mayo.snps <- mef2c@misc.data$MAYO.eqtl.snps[, c("chrom", "start", "end", "score")]
   tbl.mayo.snps <- tbl.mayo.snps[order(tbl.mayo.snps$start, decreasing=FALSE),]

   tbl.igap.snps <- mef2c@misc.data[["IGAP.snpChip"]][, c("chrom", "start", "end", "score")]
   tbl.igap.snps <- tbl.igap.snps[order(tbl.igap.snps$start, decreasing=FALSE),]

   temp.filename <- sprintf("igvdata/tmp%d.bedGraph", as.integer(Sys.time()))
   write.table(tbl.igap.snps, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
   return(temp.filename)

} # calculateVisibleSNPs
#----------------------------------------------------------------------------------------------------
removeTrack <- function(session, trackName)
{
  session$sendCustomMessage(type="removeTrack", message=list(trackName=trackName))

} # removeTrack
#----------------------------------------------------------------------------------------------------
displayTrack <- function(session, trackName)
{
   if("snps" %in% trackName){
      temp.filename <- tempfile(tmpdir="./tmp", fileext=".bed")
      printf("-- writing %d gwas snps to %s", nrow(tbl.snps), temp.filename)
      write.table(tbl.snps, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
      file.exists(temp.filename)
      session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName="snp",
                                              displayMode="EXPANDED",
                                              color="red",
                                              trackHeight=40)))
      } # snps

   if(trackName == "snpsPlusBindingSites"){
      # printf("find and display snps in binding sites")
      # load("sequenceMatchesForTfsInModel.RData")
      # motif.qval.threshold <- isolate(session$input$fimo.motif.significance)
      # snpGWAS.pval.threshold <- isolate(session$input$IGAP.snp.significance)
      # snpEffect.pval.threshold <- isolate(session$input$fimo.snp.effect)
      # printf("motif: %f  snp: %f   disruption delta: %f", motif.qval.threshold, snpGWAS.pval.threshold, snpEffect.pval.threshold)
      # browser()
      # tbl.matches <- subset(tbl.matches, q.value < 0.05)
      # printf("  %d bs below q.value threshold", nrow(tbl.matches))
      # gr.matches <- GRanges(tbl.matches)
      # gr.snps <- GRanges(tbl.snps)
      # tbl.ov <- as.data.frame(findOverlaps(gr.snps, gr.matches))
      # colnames(tbl.ov) <- c("snp", "bindingSite")
      # tbl.bs <- cbind(tbl.matches[tbl.ov$bindingSite,], tbl.snps[tbl.ov$snp,])
      # tfs.with.snps <- unique(tbl.bs$tf)
      # printf("--- tfs.with snps: %s", paste(tfs.with.snps, collapse=", "))
      # for(tf.with.snp in tfs.with.snps){
      #    tbl.tf <- subset(tbl.bs, tf==tf.with.snp)[, c("chrom", "start", "end", "name", "q.value", "pval")]
      #    tbl.tf$name <- sprintf("tfbs-snp:%s", tbl.tf$name)
      #    temp.filename <- tempfile(tmpdir="./tmp", fileext=".bed")
      #    write.table(tbl.tf, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
      #    session$sendCustomMessage(type="displayBedTrack",
      #                           message=(list(filename=temp.filename,
      #                                         trackName=tf.with.snp,
      #                                         displayMode="EXPANDED",
      #                                         color="blue",
      #                                         trackHeight=40)))
      #    } # for tf
       } # snpsInBindingSites

   if("enhancers" %in% trackName){
      load("tbl.enhancers.RData")
      temp.filename <- tempfile(tmpdir="./tmp", fileext=".bed")
      write.table(tbl.enhancers, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
      printf(" enhancers, custom message, displayBedTrack: %s", temp.filename);
      session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName="enhancers",
                                              color="black",
                                              displayMode="SQUISHED",
                                              trackHeight=40)))
      } # enhancers

   if(trackName == "dhs"){
      load("tbl.dhs.RData")
      temp.filename <- tempfile(tmpdir="./tmp", fileext=".bed")
      write.table(tbl.dhs, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
      printf(" dhs, custom message, displayBedTrack: %s", temp.filename);
      session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName="DHS",
                                              color="gray",
                                              displayMode="SQUISHED",
                                              trackHeight=30)))
      } # dhs

   if("tss" %in% trackName) {
      tss <- mef2c@misc.data$TSS
      tbl.tss <- data.frame(chrom="chr5", start=tss, end=tss, stringsAsFactors=FALSE)
      temp.filename <- sprintf("igvdata/tmp%d.bed", as.integer(Sys.time()))
      write.table(tbl.tss, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
      session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName="TSS",
                                              color="blue",
                                              trackHeight=40)))
     } # tss

    if("fp" %in% trackName){
      roi <- list(chrom="chr5", start=88391000, end=89322000)
      tbl.fp <- getFootprints(mef2c, roi)[, c("chrom", "start", "end", "shortMotifName", "score")]
      temp.filename <- sprintf("igvdata/tmp%d.bed", as.integer(Sys.time()))
      write.table(tbl.fp, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
      session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName="Footprints",
                                              color="gray",
                                              trackHeight=40)))
      } # fp

    if("tfbs.model.1" %in% trackName){
      tbl.motifs <- mef2c@misc.data[["allDNA-jaspar2018-human-mouse-motifs"]]
      tfs <- getModels(mef2c)[["mef2c.cory.wgs.cer.tfClass"]]$gene[1:5]
      for(tf in tfs){
        tbl.tfbs <- get.tf.bindingSites(tf)
        temp.filename <- sprintf("igvdata/tmp%d.bed", as.integer(Sys.time()))
        write.table(tbl.tfbs, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
        session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName=tf,
                                              color="darkRed",
                                              trackHeight=40)))
        } # for tf
      } # tfbs.model.1

    if("tfbs.model.2" %in% trackName){
      tbl.motifs <- mef2c@misc.data[["allDNA-jaspar2018-human-mouse-motifs"]]
      tfs <- getModels(mef2c)[["mef2c.cory.wgs.tcx.tfClass"]]$gene[1:5]
      for(tf in tfs){
        tbl.tfbs <- get.tf.bindingSites(tf)
        temp.filename <- sprintf("igvdata/tmp%d.bed", as.integer(Sys.time()))
        write.table(tbl.tfbs, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
        session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName=tf,
                                              color="darkRed",
                                              trackHeight=40)))
        } # for tf
      } # tfbs.model.2

    if("tfbs.model.3" %in% trackName){
      tbl.motifs <- mef2c@misc.data[["allDNA-jaspar2018-human-mouse-motifs"]]
      tfs <- getModels(mef2c)[["mef2c.cory.wgs.ros.tfClass"]]$gene[1:5]
      for(tf in tfs){
        tbl.tfbs <- get.tf.bindingSites(tf)
        temp.filename <- sprintf("igvdata/tmp%d.bed", as.integer(Sys.time()))
        write.table(tbl.tfbs, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
        session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName=tf,
                                              color="darkRed",
                                              trackHeight=40)))
        } # for tf
      } # tfbs.model.3

} # displayTrack
#------------------------------------------------------------------------------------------------------------------------
displayBindingSites <- function(session, input, tf)
{
   tbl.bindingSites <- models[[state$currentModelName]]$regulatoryRegions
   if(!tf %in% tbl.bindingSites$geneSymbol){
      printf("tf %s not in regulatory regions for model '%s'", tf, state$currentModelName)
      return()
      }
      # two kinds of bindingSite tables: one from footprint database (has "fp_start"), one from
      # MotifMatcher (has "motifStart").  branch appropriately
   if("fp_start" %in% colnames(tbl.bindingSites)){
      tbl.tfbs <- subset(tbl.bindingSites, geneSymbol==tf)[, c("chrom", "fp_start", "fp_end", "shortMotif")]
      window.start <- min(tbl.tfbs$fp_start)
      window.end   <- max(tbl.tfbs$fp_end)
      window.chrom <- tbl.tfbs$chrom[1]
      }
   if("motifStart" %in% colnames(tbl.bindingSites)){
      tbl.tfbs <- subset(tbl.bindingSites, geneSymbol==tf)[, c("chrom", "motifStart", "motifEnd", "shortMotif")]
      window.start <- min(tbl.tfbs$motifStart)
      window.end   <- max(tbl.tfbs$motifEnd)
      window.chrom <- tbl.tfbs$chrom[1]
      }
     #threshold <- isolate(input$pwmMatchThreshold) / 100
     # printf("---  using slider threshold %f", threshold);
     #tbl.tfbs <- subset(tbl.bindingSites, motifRelativeScore >= threshold)[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")]
     # tbl.tfbs <- tbl.tfbs[order(tbl.tfbs$motifStart, decreasing=FALSE),]
   printf("      found %d sites over %d bases", nrow(tbl.tfbs), window.end - window.start)

   #printf("scores on binding sites: ")
   #print(fivenum(tbl.tfbs$motifRelativeScore))
   temp.filename <- tempfile(tmpdir="tmp", fileext=".bed")
   write.table(tbl.tfbs, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
   updateTabsetPanel(session, "trenaTabs", select="igvTab");

   myFunc <- function(){
      genomicRegion <- sprintf("%s:%d-%d", window.chrom, window.start-500, window.end+500)
      session$sendCustomMessage(type="showGenomicRegion", message=(list(region=genomicRegion)))
      session$sendCustomMessage(type="displayBedTrack",
                                  message=(list(filename=temp.filename,
                                                trackName=tf,
                                                color="darkRed",
                                                trackHeight=40)))
     } # myFunc

   later(myFunc, 1)

} # displayBindingSites
#------------------------------------------------------------------------------------------------------------------------
Sys.sleep(3)
shinyApp(ui, server)
