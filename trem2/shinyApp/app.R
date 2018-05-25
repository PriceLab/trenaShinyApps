library(shiny)
library(DT)
library(igvShiny)
library(htmlwidgets)
library(GenomicRanges)
#----------------------------------------------------------------------------------------------------


  # "tbl.model.enhancerAll.mdb"  "tbl.model.enhancerAll.tfc"
  # "tbl.motifs.enhancerAll.mdb" "tbl.motifs.enhancerAll.tfc"
  # "tbl.enhancer"



load("trem2-models.RData")
models <- all.models
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
ui <- fluidPage(

  includeScript("message-handler.js"),
  tags$head(tags$link(rel = "stylesheet", type = "text/css",
                     href = "http://maxcdn.bootstrapcdn.com/font-awesome/4.2.0/css/font-awesome.min.css")
     ),

  titlePanel("TREM2 Trena Models"),

  sidebarLayout(
     sidebarPanel(
        width=3,
        selectInput("geneModel", "Choose model:",
                    c("fp.2000up.200down.cor02"="fp.2000up.200down.cor02",
                      "fp.10kup.10kdown.cor02"="fp.10kup.10kdown.cor02",
                      "fp.enhancers.cor02"="fp.enhancers.cor02",
                      "regions.motifMatching.5kup.5kdown.pwm80.cor02"="regions.motifMatching.5kup.5kdown.pwm80.cor02",
                      "regions.motifMatching.2kup.200down.pwm80.cor02"="regions.motifMatching.2kup.200down.pwm80.cor02",
                      "regions.motifMatching.2kup.200down.pwm90.cor02"="regions.motifMatching.2kup.200down.pwm90.cor02",
                      "regions.motifMatching.2kup.200down.pwm95.cor02"="regions.motifMatching.2kup.200down.pwm95.cor02",
                      "noDNA.allTFs.cor04"="noDNA.allTFs.cor04",
                      "fp.enhancers.cor02.5kup.5kdown"="fp.enhancers.cor02.5kup.5kdown",
                      "Summary Model (of above 9)" = "summary"
                      )),

        sliderInput("pwmMatchThreshold",
                    "PWM match threshold",
                    value = 90,
                    min = 75,
                    max = 100),
        actionButton("showTargetGeneButton", "Show TREM2"),
        actionButton("showSNPsButton", "SNPs"),
        actionButton("showEnhancersButton", "Enhancers")
        ),

     mainPanel(
       tabsetPanel(type="tabs",
                   id="trenaTabs",
                   tabPanel(title="IGV",     value="igvTab",       igvShinyOutput('igvShiny')),
                   tabPanel(title="Summary", value="summaryTab",   verbatimTextOutput("summary")),
                   tabPanel(title="TRN",     value="snpTableTab",  DTOutput("geneModelTable")),
                   tabPanel(title="tf/target plot",  value="plotTab",      plotOutput("xyPlot", height=800))
                   )
       ) # mainPanel
    ) # sidebarLayout
  ) # fluidPage

#--------------------------------------------------------------------------------
server <- function(input, output, session) {

   session$sendCustomMessage(type="showGenomicRegion", message=(list(region="FLT1")))

    # define a reactive conductor. it returns a function, is
    # redefined with a change to dist or n
    # has reactive children  plot, summary, table, below, which mention it
    # thereby establishing the reactive chain.

   observeEvent(input$geneModel, {
      printf("geneModel event: %s", input$geneModel);
      })

   observeEvent(input$showTargetGeneButton, {
     roi <- "chr6:41,152,337-41,213,272"
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
         #browser()
         if(length(tf) > 0){
            printf("table row clicked: %s", tf );
            # updatePlot(session, tf, state$currentModelName);
            updateTabsetPanel(session, "trenaTabs", select="igvTab");
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
    print(input$targetGene)
    print(input$roi)
    print(input$filters)
    print(input$pwmMatchThreshold)
    print(input$model.trn)
    })

   output$geneModelTable <- renderDT(
      {model.name <- input$geneModel;
       printf("    rendering DT, presumably because input$geneModel changes, model.name: %s", model.name);
       selection="single"
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
      #state[["currentModel"]] <- tbl.model
      #state[["currentModelName"]] <- model.name
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
   if("snp" %in% trackName){
       tbl.snps <- read.table("snps.tsv", stringsAsFactors=FALSE, header=TRUE)
      temp.filename <- tempfile(tmpdir="./tmp", fileext=".bed")
      write.table(tbl.snps, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
      file.exists(temp.filename)
      session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName="snp",
                                              displayMode="EXPANDED",
                                              color="red",
                                              trackHeight=40)))
      } # snps

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

    if("dhs" %in% trackName){
      tbl.dhs <- mef2c@misc.data$tbl.dhs
      temp.filename <- sprintf("igvdata/tmp%d.bed", as.integer(Sys.time()))
      write.table(tbl.dhs, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
      session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName="DHS",
                                              color="darkGreen",
                                              trackHeight=40)))
      } # dhs

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
   printf("--- displayBindingSites: %s, %s", tf, tf %in% names(tbls.bindingSites))
   tbl.bindingSites <- tbls.bindingSites[[tf]]
   threshold <- isolate(input$pwmMatchThreshold) / 100
   printf("---  using slider threshold %f", threshold);
   tbl.tfbs <- subset(tbl.bindingSites, motifRelativeScore >= threshold)[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")]
   tbl.tfbs <- tbl.tfbs[order(tbl.tfbs$motifStart, decreasing=FALSE),]
   printf("      found %d sites", nrow(tbl.tfbs))
   #printf("scores on binding sites: ")
   #print(fivenum(tbl.tfbs$motifRelativeScore))
   temp.filename <- temp.filename <- tempfile(tmpdir="tmp", fileext=".bedGraph")
   write.table(tbl.tfbs, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
   session$sendCustomMessage(type="displayBedGraphTrack",
                             message=(list(filename=temp.filename,
                                           trackName=sprintf("%s-%2d%%", tf, threshold * 100),
                                           color="green",
                                           minValue=0,
                                           maxValue=1,
                                           trackHeight=30)))

} # displayBindingSites
#------------------------------------------------------------------------------------------------------------------------
Sys.sleep(3)
shinyApp(ui, server)
