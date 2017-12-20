library(shiny)
library(ggplot2)
library(gplots)
library(data.table)
library(plyr)
library(pathview)
library(stringr)
library(biomformat)
library(KEGGREST)
library(shinythemes)
library(matrixStats)
library(png)  # For writePNG function
library(devtools)
library(d3heatmap)
library(gplots)
library(RColorBrewer)
library(rhandsontable)
library(limma)
source('mg.key.R')
load("pathview.Rdata")
options(shiny.maxRequestSize = 10 * 1024 ^ 3)
#mdt$Group <- paste(gsub('is', 'IGBS', gsub('sws', 'SS', mdt$Origin)), mdt$Source)
#mdt$Group <- gsub('inflow inflow', 'SW', gsub(' inoculum', '', mdt$Group))
ko.path.name <- ko.path.name[-grep('ko01100', ko.path.name$ko), ]
kegg <- kegg[-grep('ko01100', kegg$ko), ]
loadRdata <- function(fname) {
  load(file = fname)
  list <-
    list(
      "funtaxall" = funtaxall,
      "mdt" = mdt,
      "kegg" = kegg,
      "ko.path.name" = ko.path.name,
      "d.kres" = d.kres
    )
  return(list)
}

source("global.R")
mergeMetagenomes <- function(funtax, newName, prevNames){
  indC<-which(names(funtax)%in%c(prevNames))
  funtax$y <- rowSums(funtax[,..indC])
  names(funtax)[names(funtax) == 'y'] <- newName
  return(funtax)
}
Intfuntax <- function(result2, t1, tn, f1, fn, t2=NULL, f2=NULL){
  if(t1!="toplevel"){
    idx<-unique(unlist(sapply(tn,grep,x=result2[,get(t1)])))
    result2 <- result2[idx]
  }
  if (f1 != "toplevel") {
    result2 <- result2[grep(fn, result2[, get(f1)])]
  }
  if (!is.null(t2) & !is.null(f2)) {
    result2 <- ddply(result2, c(t2, f2), numcolwise(sum))
  } else{
    if (!is.null(t2)) {
      result2 <- ddply(result2, t2, numcolwise(sum))
    }
    if (!is.null(f2)) {
      result2 <- ddply(result2, f2, numcolwise(sum))
    }
  }
  if (dim(result2)[1] == 0) {
    result2 <- NULL
  } else {
    indR <-
      which(!is.na(result2[, 1]) & !is.na(result2[, 2]) &
              result2[, 2] != "")
    if (length(indR) > 0) {
      result2 <- result2[indR, ]
    } else {
      result2 <- NULL
    }
  }
  return(result2)
}
make2d <- function(funtax) {
  obj <-
    matrix(nrow = length(unique(funtax[, 2])), ncol = length(unique(funtax[, 1])))
  colnames(obj) <- unique(funtax[, 1])
  rownames(obj) <- unique(funtax[, 2])
  for (x in 1:nrow(funtax)) {
    obj[as.character(funtax[x, 2]), as.character(funtax[x, 1])] <-
      funtax[x, 3]
  }
  return(obj)
}
getSpecieFromAbundMD5 <-
  function(taxall,
           tx = tx,
           sp = SpName,
           aggregate = FALSE) {
    drops <-
      c("domain",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species" ,
        "usp")
    drops <- drops[drops != tx]
    dd.res <- taxall[,!(names(taxall) %in% drops), with = FALSE]
    idx <- unique(unlist(sapply(sp, grep, x = dd.res[, get(tx)])))
    d.res <- dd.res[idx]
    if (aggregate &
        dim(d.res)[1] > 1)
      d.res <- aggregate(. ~ ufun, as.data.frame(d.res)[, -c(1, 2)], FUN = sum)
    return(d.res)
  }
getSpecieFromAbundMD5_2 <- function(d.bm,
                                    sp = SpName,
                                    aggregate = FALSE) {
  idx <- unique(unlist(sapply(sp, grep, x = d.bm$usp)))
  d.res <- d.bm[idx]
  if (aggregate &
      dim(d.res)[1] > 1)
    d.res <- aggregate(. ~ ufun, as.data.frame(d.res)[, -c(1, 3)], FUN = sum)
  return(d.res)
}
getSpecieFromAbundMD5_2<-function(d.bm,sp=SpName,aggregate=FALSE){
  idx<-unique(unlist(sapply(sp,grep,x=d.bm$usp)))
  d.res<-d.bm[idx]
  if(aggregate&dim(d.res)[1]>1) d.res<-aggregate(.~ufun,as.data.frame(d.res)[,-c(1,3)],FUN = sum)
  return(d.res)
}
plotHeatmap<-function(obj,n,norm=TRUE,log=TRUE,fun=sd,...){
  mat = returnAppropriateObj(obj, norm, log)
  otusToKeep = which(rowSums(mat) > 0)
  if (length(otusToKeep) == 0) {
    return(matrix(NA, ncol = 1, nrow = 1))
  } else if (length(otusToKeep) == 1) {
    otuStats <- fun(mat)
  } else{
    otuStats = apply(mat[otusToKeep,], 1, fun)
  }
  otuIndices = otusToKeep[order(otuStats, decreasing = TRUE)[1:min(c(length(otusToKeep), n, dim(mat)[1]))]]
  mat2 = mat[otuIndices,]
  return(mat2)
}
returnAppropriateObj <- function(obj, norm, log) {
  if (class(obj) != 'matrix' |
      any(dim(obj) == 0)) {
    return(matrix(NA, ncol = 1, nrow = 1))
    #stop('Obj should be a matrix')
  }
  res <- avearrays(obj)
  
  if (log) {
    res <- log2(res + 1)
  }
  if (norm) {
    res <- scale(res)
  }
  return(res)
}
get_ko_data <- function(funtax, taxon, metagenomes) {
  d <- getSpecieFromAbundMD5_2(funtax, sp = taxon, aggregate = FALSE)
  indC <- c(which(names(d) == 'md5'), match(metagenomes, names(d)))
  if (dim(d)[1] == 0) {
    adk5 <- NULL
  } else{
    d5.1 <- d[, list(m5 = unlist(str_split(md5, ','))), by = .(usp, ufun, md5)]
    d5 <- merge(d5.1, d[, ..indC], by = 'md5')
    dk5 <-
      unique(merge(
        d5,
        d.kres,
        all = FALSE,
        by.x = 'm5',
        by.y = 'md5'
      )[, -c('md5', '.id')])
    adk5 <-
      aggregate(. ~ ko, as.data.frame(dk5[, -c('m5', 'usp', 'ufun', 'annotation')]), FUN =
                  sum)
  }
  return(adk5)
}

pathImage <- function(funtax, sp.li, mgm, pathwi, kostat, nms) {
  withProgress(
    message = paste(
      "Drawing KEGG pathway",
      pathwi,
      "for",
      sp.li,
      ".",
      "Please wait!"
    ),
    value = 10,
    {
      adk5 <- get_ko_data(funtax, sp.li, mgm)
      rownames(adk5) <- adk5$ko
      adk5 <- adk5[, -1]
      indRow <- match(rownames(adk5), rownames(kostat))
      indCol <- match(colnames(adk5), colnames(kostat))
      if (any(is.na(indCol)) | any(is.na(indRow))) {
        showModal(
          modalDialog(
            title = titleForNonMatchKostat,
            textForNonMatchKostat,
            easyClose = TRUE,
            footer = NULL
          )
        )
      } else {
        kostat <- kostat[indRow, indCol]
        diffr <- kostat - adk5
        if (any(diffr < 0)) {
          showModal(
            modalDialog(
              title = titleForNonMatchKostat,
              textForNonMatchKostat,
              easyClose = TRUE,
              footer = NULL
            )
          )
        } else {
          adk5 <- adk5 / kostat * 100
          obj <- as.matrix(adk5)
          colnames(obj) <- nms
          obj <- avearrays(obj)
          idx <- match(kegg$K[kegg$ko == paste0('ko', pathwi)], rownames(obj))
          idx <- idx[!is.na(idx)]
          #save(obj,pathwi,sp.li,file=paste0('dump.',pathwi,'.',sp.li,'.Rdata'))
          if (!is.null(pathwi)) {
            pathview(
              gene.data = obj,
              pathway.id = pathwi,
              species = "ko",
              out.suffix = paste0(sp.li, ".ko"),
              kegg.native = T,
              limit = list(gene = range(as.vector(obj[idx, ])), cpd =
                             1)
            )
            mg.key(
              fname = paste0(
                getwd(),
                "/",
                "ko",
                pathwi,
                ".",
                sp.li,
                ".ko.multi.png"
              ),
              names = nms,
              node.size = node.size,
              cex = 0.25,
              lwd = 0.25,
              tstep = 0.05,
              rstep = 0.05
            )
          }
        }
      }
    }
  )
}
filter_stats <- function(funtax, taxon, metagenomes, sd_cutoff) {
  adk5 <- get_ko_data(funtax, taxon, metagenomes)
  if (is.null(adk5)) {
    dk7 <- NULL
  } else{
    indM <- which(names(adk5) %in% c(metagenomes))
    adk5 <- as.data.table(adk5)
    dk6 <-
      data.frame(ID = adk5[, "ko"],
                 Means = rowMeans(adk5[, ..indM]),
                 SD = rowSds(as.matrix(adk5[, ..indM])))
    dk6idx <- order(dk6$SD, decreasing = TRUE)[which((dk6$Means != 0))]
    cutoff <- as.integer(max(2, sd_cutoff * length(dk6idx) / 100))
    dk7 <- adk5[dk6idx[1:cutoff], ]
  }
  return(dk7)
}
getpathfromKO <- function(KO) {
  temp <- kegg[K == KO]
  temp <- gsub('ko', '', paste0(unlist(temp[, "ko"]), collapse = ","))
}
pathwayHeatmap <- function(funtax, sp.lis, mgms, ko_sd) {
  adk5 <- filter_stats(funtax, sp.lis, mgms, ko_sd)
  if (is.null(adk5)) {
    a8 <- NULL
  } else{
    lastcol <- ncol(adk5) + 1
    for (y in 1:nrow(adk5)) {
      adk5[y, "pathwayID"] <- getpathfromKO(adk5[y, "ko"])
    }
    indC <- which(names(adk5) %in% c('ko', mgms))
    adk5 <- as.data.table(adk5)
    adk6 <- adk5[, list(pat = unlist(str_split(pathwayID, ','))), by = .(ko)]
    a6 <- merge(adk6, adk5[, ..indC], by = 'ko')
    a7 <- aggregate(. ~ pat, as.data.frame(a6[, -c('ko')]), FUN = sum)
    pnameIdx <- match(paste0('ko', a7$pat), ko.path.name$ko)
    a7 <- a7[!is.na(pnameIdx), ]
    pnameIdx <- pnameIdx[!is.na(pnameIdx)]
    rownames(a7) <- ko.path.name$name[pnameIdx]
    a8 <- as.matrix(a7[, -1])
    rownames(a8) <- ko.path.name$name[pnameIdx]
  }
  return(a8)
}
getpathsfromKOs <- function(KOs) {
  setkey(kegg, K)
  temp <- kegg[KOs]
  unlist(str_split(gsub('ko', '', paste0(
    unlist(temp[, "ko"]), collapse = ","
  )), ','))
}
getPathwayList <- function(funtax, sp.li, mgm, ko_sd) {
  dk7 <- filter_stats(funtax, sp.li, mgm, ko_sd)
  kos <- unique(dk7[, "ko"])
  listko <- getpathsfromKOs(unique(dk7[, "ko"]))
  ko.path.name$ko <- gsub('ko', '', ko.path.name$ko)
  pathandnames <-
    as.matrix(ko.path.name[which(ko.path.name$ko %in% listko), ])
  return(pathandnames)
}

chooseDends <- function(res) {
  if (dim(res)[1] == 1 & dim(res)[2] == 1) {
    dend <- 'none'
  } else if (dim(res)[1] == 1) {
    dend <- 'column'
  } else if (dim(res)[2] == 1) {
    dend <- 'row'
  } else{
    dend <- 'both'
  }
  return(dend)
}
kostatgen <- function(funtaxall, d.kres) {
  funtaxall <-
    funtaxall[, -c(
      "usp",
      "species",
      "genus",
      "family",
      "order",
      "class",
      "phylum",
      "domain",
      "ufun",
      "FUN2",
      "FUN3",
      "FUN4"
    )]
  d5.1 <- funtaxall[, list(m5 = unlist(str_split(md5, ','))), by = .(md5)]
  d5 <- merge(d5.1, funtaxall, by = 'md5')
  dk5 <-
    unique(merge(
      d5,
      d.kres,
      all = FALSE,
      by.x = 'm5',
      by.y = 'md5'
    )[, -c('md5', '.id')])
  adk5 <-
    aggregate(. ~ ko, as.data.frame(dk5[, -c('m5', 'annotation')]), FUN = sum)
  rownames(adk5) <- adk5$ko
  adk5 <- adk5[, -1]
  indF <-  which(adk5 == 0, arr.ind = TRUE)
  adk5[indF] <- 10 ^ (-5)
  return(adk5)
}

ui <- fluidPage(
  titlePanel(maintitle),
  sidebarPanel(
    "Note: values in cells are log2 of read abundance.\n",
    conditionalPanel(condition = "input.conditionedPanels==2 || input.conditionedPanels==3 || input.conditionedPanels==4 || input.conditionedPanels==5",
                     uiOutput("Mgall")),
    conditionalPanel(condition = "input.conditionedPanels==1",
                     uiOutput("mg1")),
    conditionalPanel(
      condition = "input.conditionedPanels==1 || input.conditionedPanels==2 || input.conditionedPanels==3 || input.conditionedPanels==4 || input.conditionedPanels==5",
      selectInput(
        inputId = "tl1",
        label = taxone,
        choices = tax1n,
        selected = tax1selected,
        selectize = FALSE
      )
    ),
    conditionalPanel(
      condition = "input.conditionedPanels==1 || input.conditionedPanels==2 || input.conditionedPanels==3 || input.conditionedPanels==4",
      uiOutput("taxNames")
    ),
    conditionalPanel(
      condition = "input.conditionedPanels==5",
      uiOutput("taxNamesKEGG")
    ),
    conditionalPanel(
      condition = "input.conditionedPanels==1 || input.conditionedPanels==3",
      selectInput(
        inputId = "tl2",
        label = taxtwo,
        choices = tax2n,
        selected = tax2selected
      )
    ),
    conditionalPanel(
      condition = "input.conditionedPanels==1 || input.conditionedPanels==2 || input.conditionedPanels==3",
      selectInput(
        inputId = "fl1",
        label = funcone,
        choices = func1n,
        selected = func1selected,
        selectize = FALSE
      ),
      uiOutput("funNames")
    ),
    conditionalPanel(
      condition = "input.conditionedPanels==1 || input.conditionedPanels==2",
      selectInput(
        inputId = "fl2",
        label = functwo,
        choices = func2n,
        selected = func2selected
      )
    ),
    conditionalPanel(
      condition = "input.conditionedPanels==1",
      sliderInput(
        "pix1",
        "Heatmap height (px)",
        value = 20,
        min = 10,
        max = 100
      )
    ),
    conditionalPanel(
      condition = "input.conditionedPanels==2",
      sliderInput(
        "pix2",
        "Heatmap height (px)",
        value = 20,
        min = 10,
        max = 100
      )
    ),
    conditionalPanel(
      condition = "input.conditionedPanels==3",
      sliderInput(
        "pix3",
        "Heatmap height (px)",
        value = 20,
        min = 10,
        max = 100
      )
    ),
    conditionalPanel(
      condition = "input.conditionedPanels==4",
      actionButton("goButton", "GO"),
      sliderInput(
        "pix4",
        "Heatmap height (px)",
        value = 20,
        min = 10,
        max = 100
      )
    ),
    conditionalPanel(condition = "input.conditionedPanels==5",
                     actionButton("path", "GO"),
                     uiOutput("PathwayID")),
    conditionalPanel(condition = "input.conditionedPanels==1 ||input.conditionedPanels==2 || input.conditionedPanels==3 || input.conditionedPanels==4",
                     textInput("filename", "Enter file name")),
    conditionalPanel(
      condition = "input.conditionedPanels==1 ||input.conditionedPanels==2 || input.conditionedPanels==3 || input.conditionedPanels==4",
      radioButtons(
        inputId = "var3",
        label = "Select the file format",
        choices = list("png", "pdf")
      )
    ),
    conditionalPanel(condition = "input.conditionedPanels==1",
                     uiOutput("downLink1")),
    conditionalPanel(condition = "input.conditionedPanels==2",
                     uiOutput("downLink2")),
    conditionalPanel(condition = "input.conditionedPanels==3",
                     uiOutput("downLink3")),
    conditionalPanel(condition = "input.conditionedPanels==4",
                     uiOutput("downLink4")),
    conditionalPanel(
      condition = "input.conditionedPanels==5 || input.conditionedPanels==4",
      sliderInput(
        "ko_sd",
        "SD cutoff for KO terms",
        value = 20,
        min = 0,
        max = 100,
        step = 0.1
      )
    ),
    conditionalPanel(condition = "input.conditionedPanels==5",
                     actionButton("down5", "Download KEGG map")),
    conditionalPanel(
      condition = "input.conditionedPanels==6",
      fileInput('InFile', 'Upload previously saved Rdata file.'),
      actionButton("loadRdata", "Upload Rdata"),
      textInput(
        "Rdataname",
        "Enter file name for Rdata being saved (with '.Rdata' in the end"
      ),
      actionButton("saveRdata", "Save current Rdata")
    ),
    conditionalPanel(
      condition = "input.conditionedPanels==7",
      h3("Your Metadata"),
      selectInput(
        "colName",
        ColNameSelectorTitle,
        choices = colnames(mdt),
        selected = colName,
        selectize = FALSE
      )
    ),
    width = 3
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel(titletab1, uiOutput("dynamic1") , value = 1),
      tabPanel(titletab2, uiOutput("dynamic2"), value = 2),
      tabPanel(titletab3, uiOutput("dynamic3"), value = 3),
      tabPanel(titletab4, uiOutput("dynamic4"), value = 4),
      tabPanel(
        titletab5,
        imageOutput("Pathway", width = "100%", height = "400px"),
        value = 5
      ),
      id = "conditionedPanels",
      tabPanel(
        titletab6,
        selectInput(
          inputId = "set_taxlevel1" ,
          label = taxone,
          choices = tax1n,
          selected = tax1selected
        ),
        selectInput(
          inputId = "set_taxlevel2",
          label = taxtwo,
          choices = tax2n,
          selected = tax2selected
        ),
        selectInput(
          inputId = "set_funlevel1",
          label = funcone,
          choices = func1n,
          selected = func1selected
        ),
        selectInput(
          inputId = "set_funlevel2",
          label = functwo,
          choices = func2n,
          selected = func2selected
        ),
        selectInput(
          inputId = "colorPalette",
          label = labelColPal,
          choices = rownames(brewer.pal.info),
          selected = currentPalette
        ),
        actionButton("showAllCols", "Show All Color Palettes"),
        plotOutput("paletteOutput"),
        actionButton("save_changes", "Save Changes"),
        value = 6
      ),
      tabPanel(
        titletab7,
        rHandsontableOutput("hot"),
        div(
          class = 'row',
          div(h3("Edit Your Metadata"), class = "col-sm-5", uiOutput("ui_newcolname")),
          div(
            class = "col-sm-4",
            h3("Select the type of a new column"),
            radioButtons(
              "newcolumntype",
              "Type of a new column",
              c("integer", "double", "character")
            )
          ),
          div(class = "col-sm-3")
        ),
        actionButton("addcolumn", "Create a new column"),
        actionButton("addcolumn2", "Save a new column"),
        actionButton("saveBtn", "Save"),
        value = 7
      ) #dataTableOutput("table1")
    ),
    width = 9
  )
)

server <- function(input, output, session) {
  #Reactive values
  mgall <- reactive({
    input$mgall
  })
  mg1   <- reactive({
    input$mg1
  })
  tl1   <- reactive({
    input$tl1
  })
  tl2   <- reactive({
    input$tl2
  })
  fl1   <- reactive({
    input$fl1
  })
  fl2   <- reactive({
    input$fl2
  })
  taxnames <-
    reactiveValues(tn = as.character(funtaxall$genus[(nrow(funtaxall) / 4)]),
                   tnKEGG = as.character(funtaxall$genus[(nrow(funtaxall) / 4)]),
                   fn = as.character(funtaxall$FUN4[(nrow(funtaxall) / 3)]))
  tn    <- reactive({
    if (tl1() == 'toplevel') {
      taxnames$tn <- 'toplevel'
    } else{
      taxnames$tn <- input$tn
    }
    return(taxnames$tn)
  })
  
  tnKEGG<- reactive({
    if (tl1() == 'toplevel') {
      taxnames$tnKEGG <- 'toplevel'
    } else{
      taxnames$tnKEGG <- input$tnKEGG
    }
    return(taxnames$tnKEGG)
  })
  
  fn    <- reactive({
    input$fn
  })
  ko_sd <- reactive({
    input$ko_sd
  })
  numrow1 <- reactiveValues(plot1 = 0)
  numrow2 <- reactiveValues(plot2 = 0)
  numrow3 <- reactiveValues(plot3 = 0)
  numrow4 <- reactiveValues(plot4 = 0)
  downHeat1 <- reactiveValues(is = TRUE)
  downHeat2 <- reactiveValues(is = TRUE)
  downHeat3 <- reactiveValues(is = TRUE)
  downHeat4 <- reactiveValues(is = TRUE)
  set_taxlevel1 <- reactive({
    input$set_taxlevel1
  })
  set_taxlevel2 <- reactive({
    input$set_taxlevel2
  })
  set_funlevel1 <- reactive({
    input$set_funlevel1
  })
  set_funlevel2 <- reactive({
    input$set_funlevel2
  })
  colorPalette <- reactive({
    input$colorPalette
  })
  colPal <- reactive({
    brewer.pal(9, input$colorPalette)
  })
  colName <- reactive({
    input$colName
  })
  ######################################
  #FUNCTIONS
  plotInput1 <- function(pdf) {
    tl1 <- tl1()
    tl2 <- tl2()
    tn  <- tn()
    fl1 <- fl1()
    fl2 <- fl2()
    fn  <- fn()
    mg1 <- mg1()
    colPal <- colPal()
    keepcols <- which(names(funtaxall) %in% c(tl1, tl2, fl1, fl2, mg1))
    funtax <- funtaxall[, ..keepcols]
    funtax <- Intfuntax(funtax, tl1, tn, fl1, fn, t2 = tl2, f2 = fl2)
    if (is.null(funtax)) {
      showModal(
        modalDialog(
          title = titleIntfuntax,
          textIntfuntaxInPlot1,
          easyClose = TRUE,
          footer = NULL
        )
      )
    } else{
      obj <- make2d(funtax)
      obj[is.na(obj)] <- 0
      rowmean <- data.frame(Means = rowMeans(obj))
      colmean <- data.frame(Means = colMeans(obj))
      obj <- obj[which(rowmean$Means != 0), which(colmean$Means != 0)]
      if (dim(obj)[1] > 1) {
        res <- plotHeatmap(obj, 50, trace = "none", norm = FALSE)
      } else{
        res <- returnAppropriateObj(obj, norm = FALSE, log = TRUE)
      }
      res[is.na(res)] <- 0
      main <- makePlot1Title(tl1, tl2, tn, fl1, fl2, fn, mg1)
      x <-
        heatmap.2(
          res,
          dendrogram = chooseDends(res),
          col = colPal,
          main = 'F vs. T heatmap',
          sepcolor = "black",
          sepwidth = c(0.05, 0.05),
          key = TRUE,
          keysize = 0.75,
          key.par = list(cex = 0.7),
          symkey = FALSE,
          density.info = "none",
          cexRow = 1,
          cexCol = 1,
          cex.main = 1,
          margins = c(20, 30),
          trace = "none",
          srtCol = 50
        )
      if (pdf) {
        op = par(mar = c(0, 0, 0, 0), family = "mono")
        width <- 1100
        height <- 1800
        plot(
          c(0, width),
          c(0, height),
          type = "n",
          xlab = "",
          ylab = "",
          xaxs = "i",
          yaxs = "i"
        )
        main <- gsub('\t', '    ', main)
        x <- 100 + strwidth(main, units = 'user', cex = 1) / 2
        y <- height - 100 - strheight(main, units = 'user', cex = 1) / 2
        text(x, y, labels = main, adj = 0)
        par(op)
        
      }
    }
  }
  plotInput2 <- function(pdf) {
    tl1 <- tl1()
    tn  <- tn()
    fl1 <- fl1()
    fl2 <- fl2()
    fn  <- fn()
    mg2 <- mgall()
    colPal <- colPal()
    keepcols <- which(names(funtaxall) %in% c(tl1, fl1, fl2, mg2))
    funtax <- funtaxall[, ..keepcols]
    funtax <- Intfuntax(funtax, tl1, tn, fl1, fn, f2 = fl2)
    if (is.null(funtax)) {
      showModal(
        modalDialog(
          title = titleIntfuntax,
          textIntfuntaxInPlot2,
          easyClose = TRUE,
          footer = NULL
        )
      )
    } else{
      obj <- funtax[, -c(1)]
      rownames(obj) <- funtax[, fl2]
      colnames(obj) <- as.character(mdt[c(colnames(obj)), colName()])
      dk6 <- data.frame(Means = rowMeans(obj))
      obj <- obj[which(dk6$Means != 0), ]
      obj <- as.matrix(obj)
      main <- makePlot2Title(tl1, tn, fl1, fl2, fn, mg2)
      if (dim(obj)[1] > 1) {
        res <- plotHeatmap(obj, 50, trace = "none", norm = FALSE)
      } else{
        res <- returnAppropriateObj(obj, norm = FALSE, log = TRUE)
      }
      res[is.na(res)] <- 0
      x <-
        heatmap.2(
          res,
          dendrogram = chooseDends(res),
          col = colPal,
          sepcolor = "black",
          main = plot2Title,
          sepwidth = c(0.05, 0.05),
          key = TRUE,
          keysize = 0.75,
          key.par = list(cex = 0.7),
          symkey = FALSE,
          density.info = "none",
          cexRow = 1,
          cexCol = 1,
          cex.main = 0.5,
          margins = c(20, 30),
          trace = "none",
          srtCol = 50
        )
      if (pdf) {
        op = par(mar = c(0, 0, 0, 0), family = "mono")
        width <- 1100
        height <- 1800
        plot(
          c(0, width),
          c(0, height),
          type = "n",
          xlab = "",
          ylab = "",
          xaxs = "i",
          yaxs = "i"
        )
        main <- gsub('\t', '    ', main)
        x <- 100 + strwidth(main, units = 'user', cex = 1) / 2
        y <- height - 100 - strheight(main, units = 'user', cex = 1) / 2
        text(x, y, labels = main, adj = 0)
        par(op)
        
      }
    }
  }
  plotInput3 <- function(pdf) {
    tl1 <- tl1()
    tl2 <- tl2()
    tn  <- tn()
    fl1 <- fl1()
    fn  <- fn()
    mg3 <- mgall()
    colPal <- colPal()
    keepcols <- which(names(funtaxall) %in% c(tl1, tl2, fl1, mg3))
    funtax <- funtaxall[, ..keepcols]
    funtax <- Intfuntax(funtax, tl1, tn, fl1, fn, t2 = tl2)
    if (is.null(funtax)) {
      showModal(
        modalDialog(
          title = titleIntfuntax,
          textIntfuntaxInPlot3,
          easyClose = TRUE,
          footer = NULL
        )
      )
    } else{
      obj <- funtax[, -c(1)]
      rownames(obj) <- funtax[, tl2]
      colnames(obj) <- as.character(mdt[c(colnames(obj)), colName()])
      dk6 <- data.frame(Means = rowMeans(obj))
      obj <- obj[which(dk6$Means != 0), ]
      obj <- as.matrix(obj)
      if (dim(obj)[1] > 1) {
        res <- plotHeatmap(obj, 50, trace = "none", norm = FALSE)
      } else{
        res <- returnAppropriateObj(obj, norm = FALSE, log = TRUE)
      }
      res[is.na(res)] <- 0
      main <- makePlot3Title(tl1, tl2, tn, fl1, fn, mg3)
      x <-
        heatmap.2(
          res,
          dendrogram = chooseDends(res),
          col = colPal,
          main = plot3Title,
          sepcolor = "black",
          sepwidth = c(0.05, 0.05),
          key = TRUE,
          keysize = 0.75,
          key.par = list(cex = 0.7),
          symkey = FALSE,
          density.info = "none",
          cexRow = 1,
          cexCol = 1,
          margins = c(20, 30),
          trace = "none",
          srtCol = 50
        )
      if (pdf) {
        op = par(mar = c(0, 0, 0, 0), family = "mono")
        width <- 1100
        height <- 1800
        plot(
          c(0, width),
          c(0, height),
          type = "n",
          xlab = "",
          ylab = "",
          xaxs = "i",
          yaxs = "i"
        )
        main <- gsub('\t', '    ', main)
        x <- 100 + strwidth(main, units = 'user', cex = 1) / 2
        y <- height - 100 - strheight(main, units = 'user', cex = 1) / 2
        text(x, y, labels = main, adj = 0)
        par(op)
        
      }
    }
  }
  plotInput4 <- function(pdf) {
    tl1 <- tl1()
    tn  <- tn()
    mgall <- mgall()
    ko_sd <- ko_sd()
    colPal <- colPal()
    keepcols <- which(names(funtaxall) %in% c(tl1, "ufun", "md5", mgall))
    funtax <- funtaxall[, ..keepcols]
    names(funtax)[names(funtax) == tl1] <- 'usp'
    obj <- pathwayHeatmap(funtax, tn, mgall, ko_sd)
    if (is.null(obj)) {
      showModal(
        modalDialog(
          title = titleIntfuntax,
          textDimErrorPlot4,
          easyClose = TRUE,
          footer = NULL
        )
      )
    } else{
      colnames(obj) <- as.character(mdt[c(colnames(obj)), colName()])
      mat3 <-
        plotHeatmap(obj,
                    100,
                    norm = FALSE,
                    log = TRUE,
                    trace = "none")
      mat3[is.na(mat3)] <- 0
      main <- makePlot4Title(tl1, tn, ko_sd, mgall)
      x <-
        heatmap.2(
          mat3,
          dendrogram = chooseDends(mat3),
          col = colPal,
          main = plot4Title,
          sepcolor = "black",
          sepwidth = c(0.05, 0.05),
          key = TRUE,
          keysize = 0.75,
          key.par = list(cex = 0.7),
          symkey = FALSE,
          density.info = "none",
          cexRow = 1,
          cexCol = 1,
          margins = c(20, 30),
          trace = "none",
          srtCol = 50
        )
      if (pdf) {
        op = par(mar = c(0, 0, 0, 0), family = "mono")
        width <- 1100
        height <- 1800
        plot(
          c(0, width),
          c(0, height),
          type = "n",
          xlab = "",
          ylab = "",
          xaxs = "i",
          yaxs = "i"
        )
        main <- gsub('\t', '    ', main)
        x <- 100 + strwidth(main, units = 'user', cex = 1) / 2
        y <- height - 100 - strheight(main, units = 'user', cex = 1) / 2
        text(x, y, labels = main, adj = 0)
        par(op)
        
      }
    }
  }
  plotInput5 <- function() {
    sp.li <- tnKEGG()
    tl1 <- tl1()
    pathwi <- pathw()
    mgall <- mgall()
    keepcols <-
      which(names(funtaxall) %in% c(tl1, "ufun", "md5", mgall))
    funtax <- funtaxall[, ..keepcols]
    funtax <- Intfuntax(funtax, tl1, sp.li, 'toplevel', NULL)
    names <-
      as.character(mdt[match(mgall, rownames(mdt)), colName()])
    x <-
      pathImage(funtax, sp.li, mgall, pathwi, kostat, nms = names)
  }  
  makePlot1Title <- function(tl1, tl2, tn, fl1, fl2, fn, mg1) {
    if (fl1 == 'toplevel')
      fn <- 'root'
    if (tl1 == 'toplevel')
      tn <- 'root'
    main <-
      paste(
        plot1Title,
        '\n',
        'Metagenome dimension:\n',
        '\tSelection:\t',
        mg1,
        '\n',
        'Functional dimension:\n',
        '\tSelection:\n\t\tlevel:\t"',
        fl1,
        '"\n\t\tvalues:\n',
        paste0('\t\t\t"', fn, '"\n'),
        '\tAggregation:\n\t\tlevel:\t"',
        fl2,
        '"\n',
        'Taxonomy dimension:\n',
        '\tSelection:\n\t\tlevel:\t"',
        tl1,
        '"\n\t\tvalues:\n',
        paste('\t\t\t"', tn, '"\n', collapse = ', '),
        '\tAggregation:\n\t\tlevel:\t"',
        tl2,
        '"\n'
      )
    return(main)
  }
  makePlot2Title <- function(tl1, tn, fl1, fl2, fn, mg2) {
    if (fl1 == 'toplevel')
      fn <- 'root'
    if (tl1 == 'toplevel')
      tn <- 'root'
    main <-
      paste(
        plot2Title,
        '\n',
        'Metagenome dimension:\n',
        '\tSelection:\t',
        paste0(mg2, collapse = ', '),
        '\n',
        '\tAggregation:\t',
        colName(),
        '\n',
        'Functional dimension:\n',
        '\tSelection:\n\t\tlevel:\t"',
        fl1,
        '"\n\t\tvalues:\n',
        paste0('\t\t\t"', fn, '"\n'),
        '\tAggregation:\n\t\tlevel:\t"',
        fl2,
        '"\n',
        'Taxonomy dimension:\n',
        '\tSelection:\n\t\tlevel:\t"',
        tl1,
        '"\n\t\tvalues:\n',
        paste0('\t\t\t"', tn, '"', collapse = ',\n')
      )
    return(main)
  }
  makePlot3Title <- function(tl1, tl2, tn, fl1, fn, mg3) {
    if (fl1 == 'toplevel')
      fn <- 'root'
    if (tl1 == 'toplevel')
      tn <- 'root'
    main <-
      paste(
        plot3Title,
        '\n',
        'Metagenome dimension:\n',
        '\tSelection:\t',
        paste0(mg3, collapse = ', '),
        '\n',
        '\tAggregation:\t',
        colName(),
        '\n',
        'Functional dimension:\n',
        '\tSelection:\n\t\tlevel:\t"',
        fl1,
        '"\n\t\tvalues:\n',
        paste0('\t\t\t"', fn, '"\n'),
        'Taxonomy dimension:\n',
        '\tSelection:\n\t\tlevel:\t"',
        tl1,
        '"\n\t\tvalues:\n',
        paste0('\t\t\t"', tn, '"', collapse = ',\n'),
        '\n\tAggregation:\n\t\tlevel:\t"',
        tl2,
        '"\n'
      )
    return(main)
  }
  makePlot4Title <- function(tl1, tn, ko_sd, mgall) {
    if (tl1 == 'toplevel')
      tn <- 'root'
    main <-
      paste(
        plot4Title,
        '\n',
        'Metagenome dimension:\n',
        '\tSelection:\t',
        mgall,
        '\n',
        '\tAggregation:\t',
        colName(),
        '\n',
        'Taxonomy dimension:\n',
        '\tSelection:\n\t\tlevel:\t"',
        tl1,
        '"\n\t\tvalues:\n',
        paste0('\t\t\t"', tn, '"', collapse = ',\n')
      )
    return(main)
  }
#########################################  
#  
  output$Mgall <- renderUI({
    selectInput(
      inputId = "mgall",
      label = metagenomeone,
      choices = setNames(c(colnames(funtaxall)[-c(1:13)]), mdt[, colName()]),
      selected = c(colnames(funtaxall)[metagenome1selected]),
      selectize = TRUE,
      multiple = TRUE
    )
  })
  output$mg1 <- renderUI({
    selectInput(
      inputId = "mg1",
      label = metagenometwo,
      choices = setNames(c(colnames(funtaxall)[-c(1:13)]), mdt[, colName()]),
      selected = c(colnames(funtaxall)[metagenome2selected]),
      selectize = FALSE
    )
  })
  output$taxNames <- renderUI({
    x <- input$tl1
    if (x != "toplevel") {
      selectInput(
        inputId = "tn",
        label = taxthree,
        multiple = TRUE,
        choices = as.vector(unique(funtaxall[, get(x)])),
        selected = taxnames$tn
      )
    }
  })
  output$taxNamesKEGG <- renderUI({
    x <- input$tl1
    if (x != "toplevel") {
      selectInput(
        inputId = "tnKEGG",
        label = taxthree,
        multiple = FALSE,
        choices = as.vector(unique(funtaxall[, get(x)])),
        selected = taxnames$tn
      )
    }
  })
  output$funNames <- renderUI({
    y <- input$fl1
    if (y != "toplevel") {
      selectInput(
        inputId = "fn",
        label = functhree,
        choices = as.vector(unique(funtaxall[, get(y)])),
        selected = taxnames$fn
      )
    }
  })
  output$downLink1 <- renderUI({
    if (downHeat1$is == TRUE) {
      downloadButton(outputId = "down1", label = "Download the heatmap")
    }
  })
  output$down1 <- downloadHandler(
    filename =  function() {
      paste(input$filename, input$var3, sep = ".")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      if (input$var3 == "png") {
        png(file,
            width = 3000,
            height = 2300,
            pointsize = 35) # open the png device
        plotInput1(pdf = FALSE)
      } else {
        pdf(file, width = 15, height = 15) # open the pdf device
        plotInput1(pdf = TRUE)
      }
      dev.off()
    }
  )
  output$plot1 <- renderD3heatmap({
    tl1 <- tl1()
    tl2 <- tl2()
    tn  <- tn()
    fl1 <- fl1()
    fl2 <- fl2()
    fn  <- fn()
    mg1 <- mg1()
    colPal <- colPal()
    if (is.null(fn) | is.null(tn)) {
      downHeat1$is <- FALSE
      return()
    }
    keepcols <-
      which(names(funtaxall) %in% c(tl1, tl2, fl1, fl2, mg1))
    funtax <- funtaxall[, ..keepcols]
    funtax <-
      Intfuntax(funtax, tl1, tn, fl1, fn, t2 = tl2, f2 = fl2)
    if (is.null(funtax)) {
      showModal(
        modalDialog(
          title = titleIntfuntax,
          textIntfuntaxInPlot1,
          easyClose = TRUE,
          footer = NULL
        )
      )
    } else{
      obj <- make2d(funtax)
      obj[is.na(obj)] <- 0
      rowmean <- data.frame(Means = rowMeans(obj))
      colmean <- data.frame(Means = colMeans(obj))
      idxM <- which(rowmean$Means != 0)
      obj <- obj[idxM, which(colmean$Means != 0)]
      if (length(idxM) > 1) {
        res <- plotHeatmap(obj, 50, trace = "none", norm = FALSE)
      } else{
        res <- returnAppropriateObj(obj, norm = FALSE, log = TRUE)
      }
      res[is.na(res)] <- 0
      numrow1$plot1 <- dim(res)[1]
      if (dim(res)[1] > 1 & dim(res)[2] > 1) {
        downHeat1$is <- TRUE
        d3heatmap(
          res,
          dendrogram = chooseDends(res),
          xaxis_height = 220,
          yaxis_width = 270,
          yaxis_font_size = "10px",
          xaxis_font_size = "10px",
          scalecolors = colPal
        )
      } else{
        downHeat1$is <- FALSE
        showModal(
          modalDialog(
            title = titleForDimErrorPopup,
            textForDimErrorPopup,
            easyClose = TRUE,
            footer = NULL
          )
        )
      }
    }
  })
  output$dynamic1 <- renderUI({
    d3heatmapOutput("plot1", height = paste0(numrow1$plot1 * input$pix1 + 220, "px"))
  })
  
  output$downLink2 <- renderUI({
    if (downHeat2$is == TRUE) {
      downloadButton(outputId = "down2", label = "Download the heatmap")
    }
  })
  output$down2 <- downloadHandler(
    filename =  function() {
      paste(input$filename, input$var3, sep = ".")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      if (input$var3 == "png") {
        png(file,
            width = 3000,
            height = 2300,
            pointsize = 35) # open the png device
        plotInput2(FALSE)
      } else{
        pdf(file, width = 15, height = 15) # open the pdf device
        plotInput2(TRUE)
      }
      dev.off()
    }
  )
  output$plot2 <- renderD3heatmap({
    tl1 <- tl1()
    tn  <- tn()
    fl1 <- fl1()
    fl2 <- fl2()
    fn  <- fn()
    mg2 <- mgall()
    colPal <- colPal()
    if (is.null(fn) | is.null(tn)) {
      downHeat2$is <- FALSE
      return()
    }
    keepcols <- which(names(funtaxall) %in% c(tl1, fl1, fl2, mg2))
    funtax <- funtaxall[, ..keepcols]
    funtax <- Intfuntax(funtax, tl1, tn, fl1, fn, f2 = fl2)
    if (is.null(funtax)) {
      showModal(
        modalDialog(
          title = titleIntfuntax,
          textIntfuntaxInPlot2,
          easyClose = TRUE,
          footer = NULL
        )
      )
    } else{
      obj <- funtax[, -c(1)]
      rownames(obj) <- funtax[, fl2]
      colnames(obj) <-
        as.character(mdt[c(colnames(obj)), colName()])
      dk6 <- data.frame(Means = rowMeans(obj))
      obj <- obj[which(dk6$Means != 0), ]
      obj <- as.matrix(obj)
      if (dim(obj)[1] > 1) {
        res <- plotHeatmap(obj, 50, trace = "none", norm = FALSE)
      } else {
        res <- returnAppropriateObj(obj, norm = FALSE, log = TRUE)
      }
      res[is.na(res)] <- 0
      numrow2$plot2 <- dim(res)[1]
      main <- makePlot2Title(tl1, tn, fl1, fl2, fn, mg2)
      if (dim(res)[1] > 1 & dim(res)[2] > 1) {
        downHeat2$is <- TRUE
        d3heatmap(
          res,
          dendrogram = chooseDends(res),
          xaxis_height = 220,
          yaxis_width = 270,
          yaxis_font_size = "10px",
          xaxis_font_size = "10px",
          scalecolors = colPal
        )
      } else {
        downHeat2$is <- FALSE
        showModal(
          modalDialog(
            title = titleForDimErrorPopup,
            textForDimErrorPopup,
            easyClose = TRUE,
            footer = NULL
          )
        )
      }
    }
  })
  output$dynamic2 <- renderUI({
    d3heatmapOutput("plot2", height = paste0(numrow2$plot2 * input$pix2 + 220, "px"))
  })
  output$downLink3 <- renderUI({
    if (downHeat3$is == TRUE) {
      downloadButton(outputId = "down3", label = "Download the heatmap")
    }
  })
  output$down3 <- downloadHandler(
    filename =  function() {
      paste(input$filename, input$var3, sep = ".")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      if (input$var3 == "png") {
        png(file,
            width = 3000,
            height = 2300,
            pointsize = 35) # open the png device
        plotInput3(FALSE)
      } else{
        pdf(file, width = 15, height = 15) # open the pdf device
        plotInput3(TRUE)
      }
      dev.off()
    }
  )
  output$plot3 <- renderD3heatmap({
    tl1 <- tl1()
    tl2 <- tl2()
    tn  <- tn()
    fl1 <- fl1()
    fn  <- fn()
    mg3 <- mgall()
    colPal <- colPal()
    if (is.null(fn) | is.null(tn)) {
      downHeat3$is <- FALSE
      return()
    }
    keepcols <- which(names(funtaxall) %in% c(tl1, tl2, fl1, mg3))
    funtax <- funtaxall[, ..keepcols]
    funtax <- Intfuntax(funtax, tl1, tn, fl1, fn, t2 = tl2)
    if (is.null(funtax)) {
      showModal(
        modalDialog(
          title = titleIntfuntax,
          textIntfuntaxInPlot3,
          easyClose = TRUE,
          footer = NULL
        )
      )
    } else {
      obj <- funtax[,-c(1)]
      rownames(obj) <- funtax[, tl2]
      colnames(obj) <-
        as.character(mdt[c(colnames(obj)), colName()])
      dk6 <- data.frame(Means = rowMeans(obj))
      obj <- obj[which(dk6$Means != 0),]
      obj <- as.matrix(obj)
      main <- makePlot3Title(tl1, tl2, tn, fl1, fn, mg3)
      if (dim(obj)[1] > 1) {
        res <- plotHeatmap(obj, 50, trace = "none", norm = FALSE)
      } else {
        res <- returnAppropriateObj(obj, norm = FALSE, log = TRUE)
      }
      res[is.na(res)] <- 0
      numrow3$plot3 <- dim(res)[1]
      if (dim(res)[1] > 1 & dim(res)[2] > 1) {
        downHeat3$is <- TRUE
        d3heatmap(
          res,
          dendrogram = chooseDends(res),
          xaxis_height = 220,
          yaxis_width = 270,
          yaxis_font_size = "10px",
          xaxis_font_size = "10px",
          scalecolors = colPal
        )
      } else {
        downHeat3$is <- FALSE
        showModal(
          modalDialog(
            title = titleForDimErrorPopup,
            textForDimErrorPopup,
            easyClose = TRUE,
            footer = NULL
          )
        )
      }
    }
  })
  output$dynamic3 <- renderUI({
    d3heatmapOutput("plot3", height = paste0(numrow3$plot3 * input$pix3 + 220, "px"))
  })
  output$downLink4 <- renderUI({
    if (downHeat4$is == TRUE) {
      downloadButton(outputId = "down4", label = "Download the heatmap")
    }
  })
  output$down4 <- downloadHandler(
    filename =  function() {
      paste(input$filename, input$var3, sep = ".")
    },
    #content is a function with argument file. content writes the plot to the device
    content = function(file) {
      if (input$var3 == "png") {
        png(file,
            width = 3000,
            height = 3000,
            pointsize = 35) # open the png device
        plotInput4(FALSE)
      } else {
        w <- h <- 15
        if (numrow4$plot4 > 50) {
          h <- h * 2
        }
        pdf(file, width = w, height = h) # open the pdf device
        plotInput4(TRUE)
      }
      dev.off()
    }
  )
  observeEvent(input$path, {
    output$PathwayID <- renderUI({
      tn <- tnKEGG()
      if (!is.null(tn)) {
        tl1 <- tl1()
        mgall <- mgall()
        ko_sd <- ko_sd()
        keepcols <-
          which(names(funtaxall) %in% c(tl1, "ufun", "md5", mgall))
        funtax <- funtaxall[, ..keepcols]
        if (tl1 == 'toplevel') {
          keepcols <- c('toplevel', keepcols)
          funtax$toplevel <- factor('toplevel')
          sp.li <- 'toplevel'
        }
        funtax <- Intfuntax(funtax, tl1, tn, 'toplevel', NULL)
        if (is.null(funtax)) {
          showModal(
            modalDialog(
              title = titleIntfuntaxPath,
              textIntfuntaxPath,
              easyClose = TRUE,
              footer = NULL
            )
          )
        } else {
          names(funtax)[names(funtax) == tl1] <- 'usp'
          pathandnames <-
            as.matrix(getPathwayList(
              funtax,
              sp.li =  tn,
              mgm =  mgall,
              ko_sd = ko_sd
            ))
          selectInput(
            inputId = "PathwayID",
            label = "Input Pathway ID",
            choices = setNames(as.vector(pathandnames[, "ko"]), pathandnames[, "name"])
          )
        }
      }
    })
  })
  
  pathw <- reactive({
    input$PathwayID
  })
  
  observeEvent(input$goButton, {
    output$plot4 <- renderD3heatmap({
      tl1 <- tl1()
      tn  <- tn()
      mgall <- mgall()
      ko_sd <- ko_sd()
      colPal <- colPal()
      if (is.null(tn)) {
        downHeat4$is <- FALSE
        return()
      }
      keepcols <-
        which(names(funtaxall) %in% c(tl1, "ufun", "md5", mgall))
      funtax <- funtaxall[, ..keepcols]
      names(funtax)[names(funtax) == tl1] <- 'usp'
      obj <- pathwayHeatmap(funtax, tn, mgall, ko_sd)
      if (is.null(obj)) {
        showModal(
          modalDialog(
            title = titleIntfuntax,
            textDimErrorPlot4,
            easyClose = TRUE,
            footer = NULL
          )
        )
      } else{
        colnames(obj) <- as.character(mdt[c(colnames(obj)), colName()])
        mat3 <-
          plotHeatmap(obj,
                      100,
                      norm = FALSE,
                      log = TRUE,
                      trace = "none")
        mat3[is.na(mat3)] <- 0
        numrow4$plot4 <- dim(mat3)[1]
        main <- makePlot4Title(tl1, tn, ko_sd, mgall)
        if (dim(mat3)[1] > 1 & dim(mat3)[2] > 1) {
          downHeat4$is <- TRUE
          d3heatmap(
            mat3,
            dendrogram = chooseDends(mat3),
            xaxis_height = 220,
            yaxis_width = 270,
            yaxis_font_size = "10px",
            xaxis_font_size = "10px",
            scalecolors = colPal
          )
        } else{
          downHeat4$is <- FALSE
          showModal(
            modalDialog(
              title = titleForDimErrorPopup,
              textForDimErrorPopup,
              easyClose = TRUE,
              footer = NULL
            )
          )
        }
      }
    })
  })
  output$dynamic4 <- renderUI({
    d3heatmapOutput("plot4", height = paste0(numrow4$plot4 * input$pix4 + 220, "px"))
  })
  
  observeEvent(input$down5, {
    showModal(modalDialog(
      title = textForDownloadingMetadata,
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  kostat <- kostatgen(funtaxall, d.kres)
  output$Pathway <- renderImage({
    sp.li <- tnKEGG()
    tl1 <- tl1()
    pathwi <- pathw()
    mgall <- mgall()
    if (length(sp.li) > 1) {
      showModal(
        modalDialog(
          title = titleSevTaxInImage,
          textSevTaxInImage,
          easyClose = TRUE,
          footer = NULL
        )
      )
    }
    sp.li <- sp.li[[1]]
    if (!is.null(sp.li)) {
      keepcols <-
        which(names(funtaxall) %in% c(tl1, "ufun", "md5", mgall))
      funtax <- funtaxall[, ..keepcols]
      if (tl1 == 'toplevel') {
        keepcols <- c('toplevel', keepcols)
        funtax$toplevel <- factor('toplevel')
        sp.li <- 'toplevel'
        #save(funtax,file='funtax.tmp.Rdata')
      }
      funtax <- Intfuntax(funtax, tl1, sp.li, 'toplevel', NULL)
      if (is.null(funtax)) {
        showModal(
          modalDialog(
            title = titleIntfuntaxPath,
            textIntfuntaxPath,
            easyClose = TRUE,
            footer = NULL
          )
        )
      } else {
        names(funtax)[names(funtax) == tl1] <- 'usp'
        names <-
          as.character(mdt[match(mgall, rownames(mdt)), colName()])
        pathImage(funtax, sp.li, mgall, pathwi, kostat, names)
      }
    }
    list(
      src = paste0(getwd(), "/", "ko", pathwi, ".", sp.li, ".ko.multi.leg.png"),
      contentType = 'png',
      alt = "Press GO to select Pathway!"
    )
  }, deleteFile = FALSE)
  
  
  
  ######Metadata&Settings
  
  observeEvent(input$loadRdata, {
    inFile <- input$InFile
    listA <- loadRdata(inFile$datapath)
    funtaxall <<- listA$funtaxall
    mdt <<- listA$mdt
    values[["hot"]] <<- mdt
    d.kres <<- listA$d.kres
    kegg <<- listA$kegg
    ko.path.name <<- listA$ko.path.name
    showModal(
      modalDialog(
        title = "Congratulations!",
        "Other Rdata was successfully uploaded",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  
  observeEvent(input$saveRdata, {
    if (!is.null(values[["hot"]])) {
      #write.table(values[["hot"]], "mtcarsss")
      mdt <<- values[["hot"]]
      save(funtaxall,
           mdt,
           d.kres,
           ko.path.name,
           kegg,
           file = input$Rdataname)
      showModal(
        modalDialog(
          title = "Congratulations!",
          "Rdata was successfully saved",
          easyClose = TRUE,
          footer = NULL
        )
      )
    }
  })
  
  
  output$paletteOutput <- renderPlot({
    display.brewer.pal(8, input$colorPalette)
  }, height = 200, width = 500)
  observeEvent(input$showAllCols, {
    showModal(
      modalDialog(
        renderPlot({
          display.brewer.all()
        }, height = 700, width = 500),
        title = "All color palettes",
        easyClose = TRUE,
        footer = NULL,
        size = "l"
      )
    )
  })
  
  values = reactiveValues()
  
  data = reactive({
    if (is.null(input$hot)) {
      hot = mdt
    } else {
      hot = hot_to_r(input$hot)
    }
    
    # this would be used as a function input
    values[["hot"]] = hot
    hot
  })
  
  observeEvent(input$saveBtn, {
    showModal(
      modalDialog(
        title = titleForSavingMetadata,
        textForSavingMetadata,
        easyClose = TRUE,
        footer = NULL
      )
    )
    if (!is.null(values[["hot"]])) {
      #write.table(values[["hot"]], "mtcarsss")
      mdt <<- values[["hot"]]
      save(mdt, funtaxall, d.kres, kegg, ko.path.name, file = "pathview.Rdata")
    }
  })
  
  output$hot <- renderRHandsontable({
    if (pressed$a == FALSE) {
      DF <- values[["DF"]]
      if (!is.null(DF)) {
        rhandsontable(DF, useTypes = FALSE, stretchH = "all")
      }
    } else {
      DF = data()
      if (!is.null(DF))
        rhandsontable(DF, useTypes = TRUE, stretchH = "all")
    }
  })
  
  observe({
    input$saveBtn
    if (!is.null(input$hot)) {
      DF = hot_to_r(input$hot)
    } else {
      if (is.null(values[["DF"]]))
        DF <- mdt
      else
        DF <- values[["DF"]]
    }
    values[["DF"]] <- DF
  })
  
  output$ui_newcolname <- renderUI({
    textInput("newcolumnname",
              "Name a new column",
              sprintf("newcol%s", 1 + ncol(values[["DF"]])))
  })
  
  pressed <- reactiveValues(a = TRUE)
  
  observeEvent(input$addcolumn, {
    pressed$a <- FALSE
    DF <- isolate(values[["DF"]])
    values[["previous"]] <- DF
    newcolumn <-
      eval(parse(text = sprintf(
        '%s(nrow(DF))', isolate(input$newcolumntype)
      )))
    values[["DF"]] <-
      setNames(cbind(DF, newcolumn, stringsAsFactors = FALSE),
               c(names(DF), isolate(input$newcolumnname)))
  })
  
  observeEvent(input$addcolumn2, {
    pressed$a <- TRUE
  })
  
  observeEvent(input$save_changes, {
    tax1selected <- set_taxlevel1()
    tax2selected <- set_taxlevel2()
    func1selected <- set_funlevel1()
    func2selected <- set_funlevel2()
    currentPalette <- colorPalette()
    save(tax1selected,
         tax2selected,
         func1selected,
         func2selected,
         currentPalette,
         file = "Settings.Rdata")
  })
}
shinyApp(ui = ui, server = server)
