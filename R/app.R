library(shiny)
library(pander)
library(knitr)
library(ggplot2)
library(mmnet)
library(RCurl)
library(metagenomeSeq)
library(gplots)
library(xtable)
library(data.table)
library(RJSONIO)
library(plyr)
library(pathview)
library(stringr)
library(biomformat)
library(KEGGREST)
library(shinythemes)
library(matrixStats)
library(png)  # For writePNG function
library(devtools)
install_github("Alanocallaghan/d3heatmap") #It has color key/color bar
library(d3heatmap)
library(gplots)
library(RColorBrewer)
library(rhandsontable)
library(limma)
load("pathview.Rdata")
options(shiny.maxRequestSize=10*1024^3)
mdt$Group<-paste(gsub('is','IGBS',gsub('sws','SS',mdt$Origin)),mdt$Source)
mdt$Group<-gsub('inflow inflow','SW',gsub(' inoculum','',mdt$Group))

ko.path.name<-ko.path.name[-grep('ko01100',ko.path.name$ko),]
kegg <- kegg[-grep('ko01100', kegg$ko),]
loadRdata <- function(fname){
  load(file = fname)
  return(funtaxall)
}

funtaxall <- loadRdata("pathview.Rdata")
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
  if(f1!="toplevel"){
    result2 <- result2[grep(fn, result2[,get(f1)])]
  }
  if(!is.null(t2) & !is.null(f2)){
    result2<- ddply(result2, c(t2,f2), numcolwise(sum))
  }else{
    if(!is.null(t2)){
      result2<- ddply(result2, t2, numcolwise(sum))
    }
    if(!is.null(f2)){
      result2<- ddply(result2, f2, numcolwise(sum))
    }}
  result2 <- result2[which(!is.na(result2[,1])& !is.na(result2[,2])& result2[,2]!= ""),]
  return(result2)
}
make2d <- function(funtax){
  obj <- matrix(nrow = length(unique(funtax[,2])), ncol = length(unique(funtax[,1])))
  colnames(obj) <- unique(funtax[,1])
  rownames(obj) <- unique(funtax[,2])
  for (x in 1:nrow(funtax)){
    obj[as.character(funtax[x,2]), as.character(funtax[x,1])]<- funtax[x,3]
  }
  return(obj)
}
getSpecieFromAbundMD5<-function(taxall, tx=tx, sp = SpName, aggregate=FALSE){
  drops <- c("domain","phylum", "class", "order", "family", "genus", "species" ,"usp" )
  drops <- drops[drops!= tx]
  dd.res <- taxall[ , !(names(taxall) %in% drops), with = FALSE]
  idx<-unique(unlist(sapply(sp,grep,x=dd.res[,get(tx)])))
  d.res<-dd.res[idx]
  if(aggregate&dim(d.res)[1]>1) d.res<-aggregate(.~ufun,as.data.frame(d.res)[,-c(1,2)],FUN = sum)
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
  if(length(otusToKeep)==0){
    return(matrix(NA,ncol = 1,nrow = 1))
  }else if(length(otusToKeep)==1){
    otuStats<-fun(mat)
  }else{  
    otuStats = apply(mat[otusToKeep, ], 1, fun)
  }
  otuIndices = otusToKeep[order(otuStats, decreasing = TRUE)[1:min(c(length(otusToKeep),n,dim(mat)[1]))]]
  mat2 = mat[otuIndices, ]
  return(mat2)
}
returnAppropriateObj <- function(obj, norm, log) {
  if (class(obj) != 'matrix' |
      any(dim(obj) == 0)){
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
  d<-getSpecieFromAbundMD5_2(funtax,sp = taxon,aggregate = FALSE)
  indC<-c(which(names(d)=='md5'),match(metagenomes,names(d)))
  d5.1<-d[,list(m5=unlist(str_split(md5,','))),by=.(usp,ufun,md5)]
  d5<-merge(d5.1,d[,..indC],by='md5')
  dk5<-unique(merge(d5,d.kres,all=FALSE,by.x='m5',by.y='md5')[,-c('md5','.id')])
  adk5<-aggregate(.~ko,as.data.frame(dk5[,-c('m5', 'usp', 'ufun', 'annotation')]),FUN=sum)
}

pathImage<-function(funtax, sp.li, mgm, pathwi, kostat,nms) {
  withProgress(message = paste("Drawing KEGG pathway", pathwi, "for", sp.li, ".", "Please wait!"), value = 10, {
  adk5<-get_ko_data(funtax, sp.li, mgm)
  rownames(adk5)<-adk5$ko
  adk5<-adk5[,-1]
  indRow <- match(rownames(adk5),rownames(kostat))
  indCol <- match(colnames(adk5),colnames(kostat))
  if(any(is.na(indCol))|any(is.na(indRow))){
    showModal(modalDialog(
      title = titleForNonMatchKostat, textForNonMatchKostat, easyClose = TRUE, footer = NULL
    ))
  } else {
    kostat <- kostat[indRow, indCol]
    ind0 <- apply(kostat,2,function(.x){which(.x==0)})
    if(any(sapply(1:length(ind0), function(i){any(adk5[ind0[[i]],names(ind0)[i]]!=0)}))){
      showModal(modalDialog(
        title = titleForNonMatchZerosKostat, textForNonMatchZerosKostat, easyClose = TRUE, footer = NULL
      ))
    } else {
      sapply(1:length(ind0), function(i){kostat[ind0[[i]],names(ind0)[i]] <<- 10^(-5)})
      cat('pathImg',pathwi,class(adk5),dim(adk5),apply(adk5,2,max),nms,'\n')
      adk5 <- adk5/kostat*100
      cat(class(adk5),dim(adk5),colnames(adk5),'\n',mgm,'\n')
      obj<-as.matrix(adk5)
      cat('pathImg!',class(obj),dim(obj),apply(obj,2,max),'\n')
      colnames(obj)<-nms
      obj<-avearrays(obj)
      cat('pathImg!!',class(obj),dim(obj),apply(obj,2,max),'\n')
      #obj<-log10(avearrays(obj)+1)
      idx<-match(kegg$K[kegg$ko==paste0('ko',pathwi)],rownames(obj))
      idx<-idx[!is.na(idx)]
      cat('pathImg!!!',length(idx),length(which(is.na(idx))),'\n')
      cat('---\t',head(kegg$K[kegg$ko==paste0('ko',pathwi)]),'\n')
      cat('---\t',head(rownames(obj)),'\n')
      cat('---\t',head(rownames(adk5)),'\n')
      cat('pathImg!V',length(idx),dim(obj),apply(obj[idx,],2,max),'\n')
      save(obj,pathwi,sp.li,file=paste0('dump.',pathwi,'.',sp.li,'.Rdata'))
      pathview(gene.data = obj, pathway.id = pathwi,
               species = "ko", out.suffix = paste0(sp.li,".ko"), kegg.native = T,
               limit = list(gene=range(as.vector(obj[idx,])),cpd=1))
    }
  }
  })
}
filter_stats <- function(funtax, taxon, metagenomes, sd_cutoff) {
  adk5 <- get_ko_data(funtax, taxon, metagenomes)
  indM <- which(names(adk5)%in%c(metagenomes))
  adk5 <- as.data.table(adk5)
  dk6 <- data.frame(ID = adk5[,"ko"], Means=rowMeans(adk5[,..indM]), SD=rowSds(as.matrix(adk5[,..indM])))
  dk6idx<-order(dk6$SD,decreasing = TRUE)[which((dk6$Means!=0))]
  cutoff<-as.integer(max(2,sd_cutoff*length(dk6idx)/100))
  dk7 <- adk5[dk6idx[1:cutoff],]
}
getpathfromKO <- function(KO){
  temp <- kegg[K == KO]
  temp <- gsub('ko','',paste0(unlist(temp[,"ko"]), collapse = ","))
}
pathwayHeatmap<-function(funtax,sp.lis, mgms, ko_sd) {
  adk5 <- filter_stats(funtax, sp.lis, mgms, ko_sd)
  lastcol<- ncol(adk5)+1
  for (y in 1:nrow(adk5)){adk5[y,"pathwayID"] <- getpathfromKO(adk5[y,"ko"])}
  indC<-which(names(adk5)%in%c('ko', mgms))
  adk5 <- as.data.table(adk5)
  adk6<-adk5[,list(pat=unlist(str_split(pathwayID,','))),by=.(ko)]
  a6<-merge(adk6,adk5[,..indC],by='ko')
  a7<-aggregate(.~pat,as.data.frame(a6[,-c('ko')]),FUN=sum)
  pnameIdx<-match(paste0('ko',a7$pat),ko.path.name$ko)
  a7<-a7[!is.na(pnameIdx),]
  pnameIdx<-pnameIdx[!is.na(pnameIdx)]
  rownames(a7)<-ko.path.name$name[pnameIdx]
  a8<- as.matrix(a7[,-1])
  rownames(a8) <- ko.path.name$name[pnameIdx]
  return(a8)
}
getpathsfromKOs <- function(KOs){
  setkey(kegg, K)
  temp <- kegg[KOs]
  unlist(str_split(gsub('ko','',paste0(unlist(temp[,"ko"]), collapse = ",")), ','))
}
getPathwayList <- function(funtax, sp.li, mgm, ko_sd) {
  dk7 <- filter_stats(funtax, sp.li, mgm, ko_sd)
  kos<- unique(dk7[,"ko"])
  listko <- getpathsfromKOs(unique(dk7[,"ko"]))
  ko.path.name$ko <- gsub('ko','',ko.path.name$ko)
  pathandnames <- as.matrix(ko.path.name[which(ko.path.name$ko %in% listko),])
  return(pathandnames)
}

chooseDends <- function(res){
  if(dim(res)[1]==1&dim(res)[2]==1){
    dend <- 'none'
  } else if(dim(res)[1]==1){
    dend <- 'column'
  } else if(dim(res)[2]==1){
    dend <- 'row'
  } else{
    dend <- 'both'
  }
  return(dend)
}
kostat <- function(funtaxall, d.kres){
  funtaxall <- funtaxall[,-c("usp","species", "genus", "family", "order", "class", "phylum", "domain","ufun", "FUN2", "FUN3", "FUN4")]
  d5.1<-funtaxall[,list(m5=unlist(str_split(md5,','))),by=.(md5)]
  d5<-merge(d5.1,funtaxall,by='md5')
  dk5<-unique(merge(d5,d.kres,all=FALSE,by.x='m5',by.y='md5')[,-c('md5','.id')])
  adk5<-aggregate(.~ko,as.data.frame(dk5[,-c('m5', 'annotation')]),FUN=sum)
  rownames(adk5)<- adk5$ko
  adk5 <- adk5[,-1]
}
ui <- fluidPage(
  titlePanel(maintitle),
  sidebarPanel(
    conditionalPanel(condition = "input.conditionedPanels==2 || input.conditionedPanels==3 || input.conditionedPanels==4 || input.conditionedPanels==5",
                     uiOutput("Mgall")),
    conditionalPanel(condition = "input.conditionedPanels==1",
                     uiOutput("mg1")),
    conditionalPanel(condition = "input.conditionedPanels==1 || input.conditionedPanels==2 || input.conditionedPanels==3 || input.conditionedPanels==4 || input.conditionedPanels==5",
                     selectInput(inputId = "tl1", label = taxone, choices = tax1n, selected = tax1selected, selectize = FALSE),
                     uiOutput("taxNames")),
    conditionalPanel(condition = "input.conditionedPanels==1 || input.conditionedPanels==3",
                     selectInput(inputId = "tl2", label = taxtwo, choices = tax2n, selected = tax2selected)),
    conditionalPanel(condition = "input.conditionedPanels==1 || input.conditionedPanels==2 || input.conditionedPanels==3",
                     selectInput(inputId = "fl1", label = funcone, choices = func1n, selected = func1selected, selectize = FALSE),
                     uiOutput("funNames")),
    conditionalPanel(condition = "input.conditionedPanels==1 || input.conditionedPanels==2",
                     selectInput(inputId = "fl2", label = functwo, choices = func2n, selected = func2selected)),
    conditionalPanel(condition = "input.conditionedPanels==1",
                     sliderInput("pix1", "height", value = 20, min = 10, max = 100)),
    conditionalPanel(condition = "input.conditionedPanels==2",
                     sliderInput("pix2", "height", value = 20, min = 10, max = 100)),
    conditionalPanel(condition = "input.conditionedPanels==3",
                     sliderInput("pix3", "height", value = 20, min = 10, max = 100)),
    conditionalPanel(condition = "input.conditionedPanels==4",
                     actionButton("goButton", "GO"),
                     sliderInput("pix4", "height", value = 20, min = 10, max = 100)
                     ),
    conditionalPanel(condition = "input.conditionedPanels==5",
                     actionButton("path", "GO"),
                     uiOutput("PathwayID")
                     ),
    conditionalPanel(condition = "input.conditionedPanels==1 ||input.conditionedPanels==2 || input.conditionedPanels==3 || input.conditionedPanels==4",
                     textInput("filename","Enter file name")
                     ),
    conditionalPanel(condition = "input.conditionedPanels==1 ||input.conditionedPanels==2 || input.conditionedPanels==3 || input.conditionedPanels==4",
                     radioButtons(inputId = "var3", label = "Select the file type", choices = list("png", "pdf"))
                     ),
    conditionalPanel(condition = "input.conditionedPanels==1",
                     uiOutput("downLink1")
                     ),
    conditionalPanel(condition = "input.conditionedPanels==2",
                     uiOutput("downLink2")
                     ),
    conditionalPanel(condition = "input.conditionedPanels==3",
                     uiOutput("downLink3")
                     ),
    conditionalPanel(condition = "input.conditionedPanels==4",
                     uiOutput("downLink4")
                     ),
    conditionalPanel(condition = "input.conditionedPanels==5 || input.conditionedPanels==4",
                     sliderInput("ko_sd", "SD cutoff for KO terms", value = 20, min = 0, max = 100,step = 0.1)
                     ),
    conditionalPanel(condition = "input.conditionedPanels==5",
                     actionButton("down5", "Download KEGG map")
                     ),
    conditionalPanel(condition = "input.conditionedPanels==6",
                     fileInput('InFile', 'Upload previously saved Rdata file.'),
                     actionButton("loadRdata", "Upload Rdata"),
                     textInput("Rdataname","Enter file name for Rdata being saved (with '.Rdata' in the end"),
                     actionButton("saveRdata", "Save current Rdata")
                     ),
    conditionalPanel(condition = "input.conditionedPanels==7",
                     h3("Your Metadata"),
                     selectInput("colName", ColNameSelectorTitle, choices = colnames(mdt), selected = colName, selectize = FALSE)
                     ),
    width = 3),
  
  mainPanel(
    tabsetPanel(
      tabPanel("F&T", uiOutput("dynamic1") , value = 1), 
      tabPanel("F&M", uiOutput("dynamic2"), value = 2),
      tabPanel("T&M", uiOutput("dynamic3"), value = 3),
      tabPanel("Pathway Abundance Heatmap", uiOutput("dynamic4"), value = 4),
      tabPanel("KEGG Pathway Map", imageOutput("Pathway",width = "100%", height = "400px"), value = 5),
      id = "conditionedPanels", 
      tabPanel("Settings", selectInput(inputId = "set_taxlevel1" , label = set_taxone, choices = tax1n, selected = tax1selected),
               selectInput(inputId = "set_taxlevel2", label = set_taxtwo, choices = tax2n, selected = tax2selected),
               selectInput(inputId = "set_funlevel1", label = set_funcone, choices = func1n, selected = func1selected),
               selectInput(inputId = "set_funlevel2", label = set_functwo, choices = func2n, selected = func2selected),
               selectInput(inputId = "colorPalette", label = "Choose color palette for heatmaps", choices = rownames(brewer.pal.info), selected = currentPalette),
               actionButton("showAllCols", "Show All Color Palettes"),
               plotOutput("paletteOutput"),
               actionButton("save_changes", "Save Changes"), value = 6),
      tabPanel("Metadata",rHandsontableOutput("hot"),
               div(class='row', div(h3("Edit Your Metadata"), class="col-sm-5", uiOutput("ui_newcolname")), div(class="col-sm-4", h3("Select the type of a new column"), 
               radioButtons("newcolumntype", "Type of a new column", c("integer", "double", "character"))),
               div(class="col-sm-3")
               ), 
               actionButton("addcolumn", "Press here to create a coumn"),
               actionButton("addcolumn2", "Save a new column"),
               actionButton("saveBtn", "Save the metadata"), value = 7) #dataTableOutput("table1")
    ), width = 9)
)

server <- function(input, output, session) {
  
  #Settings
  # themes, colorPalette,
  set_taxlevel1 <-reactive({input$set_taxlevel1})
  set_taxlevel2 <-reactive({input$set_taxlevel2})
  set_funlevel1 <-reactive({input$set_funlevel1})
  set_funlevel2 <-reactive({input$set_funlevel2})
  colorPalette <-reactive({input$colorPalette})
  colPal <- reactive({brewer.pal(9, input$colorPalette)})
  colName <- reactive({input$colName})
  
  output$Mgall <- renderUI({
    selectInput(inputId = "mgall", label = metagenomeone, choices = setNames(c(colnames(funtaxall)[-c(1:13)]), mdt[,colName()]), selected = c(colnames(funtaxall)[metagenome1selected]), selectize = TRUE, multiple = TRUE)
  })
  output$mg1 <- renderUI({
    selectInput(inputId = "mg1", label = metagenometwo, choices = setNames(c(colnames(funtaxall)[-c(1:13)]), mdt[,colName()]), selected = c(colnames(funtaxall)[metagenome2selected]), selectize = FALSE)
  })  
  observeEvent(input$loadRdata, {
    inFile <- input$InFile
    funtaxall <<- loadRdata(inFile$datapath)
  })
  
  observeEvent(input$saveRdata, {
    save(funtaxall, mdt, file = input$Rdataname)
  })
  output$paletteOutput <- renderPlot({
    display.brewer.pal(8,input$colorPalette)
  }, height = 200, width = 500)
  observeEvent(input$showAllCols,{
    showModal(modalDialog(
      renderPlot({display.brewer.all()}, height = 700, width = 500),title = "All color palettes", easyClose = TRUE, footer = NULL, size = "l"))
  })
    mgall <-reactive({input$mgall})
    mg1   <-reactive({input$mg1})
    tl1   <-reactive({input$tl1})
    tl2   <-reactive({input$tl2})
    fl1   <-reactive({input$fl1})
    fl2   <-reactive({input$fl2})
    taxnames<-reactiveValues(tn=as.character(funtaxall$genus[(nrow(funtaxall)/2)]))
    tn    <-reactive({
#      isolate({
      if(tl1()=='toplevel'){
        cat('toplevel\n')
        taxnames$tn<-'toplevel'
      }else{
        taxnames$tn<-input$tn
      }
#      })
      return(taxnames$tn)
    })
    fn    <-reactive({input$fn})
    ko_sd <-reactive({input$ko_sd})
    numrow1 <- reactiveValues(plot1 =0)
    numrow2 <- reactiveValues(plot2 =0)
    numrow3 <- reactiveValues(plot3 =0)
    numrow4 <- reactiveValues(plot4 =0)
    downHeat1 <- reactiveValues(is =TRUE)
    downHeat2 <- reactiveValues(is =TRUE)
    downHeat3 <- reactiveValues(is =TRUE)
    downHeat4 <- reactiveValues(is =TRUE)
    
    plotInput1 <- function(){
      tl1 <- tl1()
      tl2 <- tl2()
      tn  <- tn()
      fl1 <- fl1()
      fl2 <- fl2()
      fn  <- fn()
      mg1 <- mg1()
      colPal <- colPal() 
      keepcols<-which(names(funtaxall)%in%c(tl1, tl2, fl1, fl2, mg1))
      funtax <- funtaxall[,..keepcols]
      funtax <- Intfuntax(funtax,tl1,tn,fl1,fn,t2 = tl2,f2 = fl2)
      obj <- make2d(funtax)
      obj[is.na(obj)] <- 0
      rowmean <- data.frame(Means=rowMeans(obj))
      colmean <- data.frame(Means=colMeans(obj))
      obj <- obj[which(rowmean$Means!=0),which(colmean$Means!=0)]
      if(dim(obj)[1]>1){
        res<-plotHeatmap(obj,50,trace = "none",norm=FALSE)
      }else{
        res<-returnAppropriateObj(obj,norm = FALSE,log = TRUE)
      }
      res[is.na(res)] <- 0 
      x <- heatmap.2(res,dendrogram = chooseDends(res), col = colPal, sepcolor="black", sepwidth=c(0.05,0.05), key=TRUE, keysize=0.75, key.par = list(cex=0.7), symkey=FALSE, density.info="none",cexRow=1,cexCol=1,margins=c(20,30),trace="none",srtCol=50)
    }
   output$downLink1 <- renderUI({
     if(downHeat1$is==TRUE){
       downloadButton(outputId = "down1", label = "Download the heatmap")
     }
   })
    output$down1 <- downloadHandler(
      filename =  function() {
        paste(input$filename, input$var3, sep=".")
      },
      # content is a function with argument file. content writes the plot to the device
      content = function(file) {
        if(input$var3 == "png")
          png(file, width = 3000, height = 2300, pointsize = 35) # open the png device
        else
          pdf(file, width = 15, height = 15) # open the pdf device
          plotInput1()
          dev.off()
      })
      
    output$taxNames <- renderUI({x <- input$tl1
    if(x!="toplevel"){
    selectInput(inputId = "tn", label = taxthree, multiple=TRUE, choices = as.vector(unique(funtaxall[,get(x)])), selected = taxnames$tn)
    }})
    output$funNames <- renderUI({y <- input$fl1
    if(y!="toplevel"){
      selectInput(inputId = "fn", label = functhree, choices = as.vector(unique(funtaxall[,get(y)])), selected = as.character(funtaxall$FUN4[(nrow(funtaxall)/2)]))
    }})
    
    output$plot1 <- renderD3heatmap({
      tl1 <- tl1()
      tl2 <- tl2()
      tn  <- tn()
      fl1 <- fl1()
      fl2 <- fl2()
      fn  <- fn()
      mg1 <- mg1()
      colPal <- colPal()
      if(is.null(fn)|is.null(tn)){
        downHeat1$is <- FALSE
        return()
      }
        keepcols <- which(names(funtaxall) %in% c(tl1, tl2, fl1, fl2, mg1))
        funtax <- funtaxall[, ..keepcols]
        funtax <- Intfuntax(funtax, tl1, tn, fl1, fn, t2 = tl2, f2 = fl2)
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
    })
    output$dynamic1 <- renderUI({
    d3heatmapOutput("plot1", height = paste0(numrow1$plot1*input$pix1+220, "px"))
  })

  plotInput2 <- function(){ 
    tl1 <- tl1()
    tn  <- tn()
    fl1 <- fl1()
    fl2 <- fl2()
    fn  <- fn()
    mg2 <- mgall()
    colPal <- colPal()
    keepcols<-which(names(funtaxall)%in%c(tl1, fl1, fl2, mg2))
    funtax <- funtaxall[,..keepcols]
    funtax <- Intfuntax(funtax,tl1,tn,fl1,fn,f2 = fl2)
    obj <- funtax[,-c(1)]
    rownames(obj)<-funtax[,fl2]
    colnames(obj)<-as.character(mdt[c(colnames(obj)), colName()])
    dk6 <- data.frame(Means=rowMeans(obj))
    obj <- obj[which(dk6$Means!=0),]
    obj <- as.matrix(obj)
    if(dim(obj)[1]>1){
      res<-plotHeatmap(obj,50,trace = "none",norm=FALSE)
    }else{
      res<-returnAppropriateObj(obj,norm = FALSE,log = TRUE)
    }
    res[is.na(res)] <- 0
    x <- heatmap.2(res,dendrogram = chooseDends(res), col = colPal, sepcolor="black", sepwidth=c(0.05,0.05), key=TRUE, keysize=0.75, key.par = list(cex=0.7), symkey=FALSE, density.info="none",cexRow=1,cexCol=1,margins=c(20,30),trace="none",srtCol=50)
  }
  output$downLink2 <- renderUI({
    if(downHeat2$is==TRUE){
      downloadButton(outputId = "down2", label = "Download the heatmap")
    }
  })
  output$down2 <- downloadHandler(
    filename =  function() {
      paste(input$filename, input$var3, sep=".")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      if(input$var3 == "png")
        png(file, width = 3000, height = 2300, pointsize = 35) # open the png device
      else
        pdf(file, width = 15, height = 15) # open the pdf device
      plotInput2()
      dev.off()
    })
  
    output$plot2 <- renderD3heatmap({
      tl1 <- tl1()
      tn  <- tn()
      fl1 <- fl1()
      fl2 <- fl2()
      fn  <- fn()
      mg2 <- mgall()
      colPal <- colPal()
      if(is.null(fn)|is.null(tn)){
        downHeat2$is <- FALSE
        return()
        }
      keepcols<-which(names(funtaxall)%in%c(tl1, fl1, fl2, mg2))
      funtax <- funtaxall[,..keepcols]
      funtax <- Intfuntax(funtax,tl1,tn,fl1,fn,f2 = fl2)
      obj <- funtax[,-c(1)]
      rownames(obj)<-funtax[,fl2]
      colnames(obj)<-as.character(mdt[c(colnames(obj)), colName()])
      dk6 <- data.frame(Means=rowMeans(obj))
      obj <- obj[which(dk6$Means!=0),]
      obj <- as.matrix(obj)
      if(dim(obj)[1]>1){
        res<-plotHeatmap(obj,50,trace = "none",norm=FALSE)
      }else{
        res<-returnAppropriateObj(obj,norm = FALSE,log = TRUE)
      }
      res[is.na(res)] <- 0
      numrow2$plot2 <- dim(res)[1]
      if(dim(res)[1]>1 & dim(res)[2]>1){
        downHeat2$is <- TRUE
        d3heatmap(res,dendrogram = chooseDends(res), xaxis_height = 220, yaxis_width = 270, yaxis_font_size = "10px", xaxis_font_size = "10px", scalecolors = colPal)
      }else{
        downHeat2$is <- FALSE
        showModal(modalDialog(
          title = titleForDimErrorPopup, textForDimErrorPopup, easyClose = TRUE, footer = NULL
        ))
      }
    })
  output$dynamic2 <- renderUI({
    d3heatmapOutput("plot2", height = paste0(numrow2$plot2*input$pix2+220, "px"))
  })
  
  plotInput3 <- function(){ 
    tl1 <- tl1()
    tl2 <- tl2()
    tn  <- tn()
    fl1 <- fl1()
    fn  <- fn()
    mg3 <- mgall()
    colPal <- colPal()
    keepcols<-which(names(funtaxall)%in%c(tl1, tl2, fl1, mg3))
    funtax <- funtaxall[,..keepcols]
    funtax <- Intfuntax(funtax,tl1,tn,fl1,fn,t2 = tl2)
    obj <- funtax[,-c(1)]
    rownames(obj)<-funtax[,tl2]
    colnames(obj)<-as.character(mdt[c(colnames(obj)), colName()])
    dk6 <- data.frame(Means=rowMeans(obj))
    obj <- obj[which(dk6$Means!=0),]
    obj <- as.matrix(obj)
    if(dim(obj)[1]>1){
      res<-plotHeatmap(obj,50,trace = "none",norm=FALSE)
    }else{
      res<-returnAppropriateObj(obj,norm = FALSE,log = TRUE)
    }
    res[is.na(res)] <- 0
    x <- heatmap.2(res,dendrogram = chooseDends(res), col = colPal, sepcolor="black", sepwidth=c(0.05,0.05), key=TRUE, keysize=0.75, key.par = list(cex=0.7), symkey=FALSE, density.info="none",cexRow=1,cexCol=1,margins=c(20,30),trace="none",srtCol=50)
  }
  output$downLink3 <- renderUI({
    if(downHeat3$is==TRUE){
      downloadButton(outputId = "down3", label = "Download the heatmap")
    }
  })
  output$down3 <- downloadHandler(
    filename =  function() {
      paste(input$filename, input$var3, sep=".")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      if(input$var3 == "png")
        png(file, width = 3000, height = 2300, pointsize = 35) # open the png device
      else
        pdf(file, width = 15, height = 15) # open the pdf device
      plotInput3()
      dev.off()
    })
  
  
    output$plot3 <- renderD3heatmap({
      tl1 <- tl1()
      tl2 <- tl2()
      tn  <- tn()
      fl1 <- fl1()
      fn  <- fn()
      mg3 <- mgall()
      colPal <- colPal()
      if(is.null(fn)|is.null(tn)){
        downHeat3$is <- FALSE
        return()
      }
      keepcols<-which(names(funtaxall)%in%c(tl1, tl2, fl1, mg3))
      funtax <- funtaxall[,..keepcols]
      funtax <- Intfuntax(funtax,tl1,tn,fl1,fn,t2 = tl2)
      obj <- funtax[,-c(1)]
      rownames(obj)<-funtax[,tl2]
      colnames(obj)<-as.character(mdt[c(colnames(obj)), colName()])
      dk6 <- data.frame(Means=rowMeans(obj))
      obj <- obj[which(dk6$Means!=0),]
      obj <- as.matrix(obj)
      if(dim(obj)[1]>1){
        res<-plotHeatmap(obj,50,trace = "none",norm=FALSE)
      }else{
        res<-returnAppropriateObj(obj,norm = FALSE,log = TRUE)
      }
      res[is.na(res)] <- 0
      numrow3$plot3 <- dim(res)[1]
      if(dim(res)[1]>1 & dim(res)[2]>1){
        downHeat3$is <- TRUE
        d3heatmap(res,dendrogram = chooseDends(res), xaxis_height = 220, yaxis_width = 270, yaxis_font_size = "10px", xaxis_font_size = "10px", scalecolors = colPal) 
      }else{
        downHeat3$is <- FALSE
        showModal(modalDialog(
          title = titleForDimErrorPopup, textForDimErrorPopup, easyClose = TRUE, footer = NULL
        ))
      }
    })
  output$dynamic3 <- renderUI({
    d3heatmapOutput("plot3", height = paste0(numrow3$plot3*input$pix3+220, "px"))
  })
  
  
  plotInput4 <- function(){ 
    tl1 <- tl1()
    tn  <- tn() 
    mgall <-mgall()
    ko_sd <- ko_sd()
    colPal <- colPal()
    keepcols<-which(names(funtaxall)%in%c(tl1,"ufun","md5", mgall))
    funtax <- funtaxall[,..keepcols]
    names(funtax)[names(funtax) == tl1] <- 'usp'
    obj<-pathwayHeatmap(funtax, tn, mgall, ko_sd)
    colnames(obj)<-as.character(mdt[c(colnames(obj)), colName()])
    mat3 <- plotHeatmap(obj,100,norm = FALSE, log = TRUE,trace = "none")
    mat3[is.na(mat3)] <- 0
    x <- heatmap.2(mat3,dendrogram = chooseDends(mat3), col = colPal, sepcolor="black", sepwidth=c(0.05,0.05), key=TRUE, keysize=0.75, key.par = list(cex=0.7), symkey=FALSE, density.info="none",cexRow=1,cexCol=1,margins=c(20,30),trace="none",srtCol=50)
  }
  output$downLink4 <- renderUI({
    if(downHeat4$is==TRUE){
      downloadButton(outputId = "down4", label = "Download the heatmap")
    }
  })
  output$down4 <- downloadHandler(
    filename =  function() {
      paste(input$filename, input$var3, sep=".")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      if(input$var3 == "png")
        png(file, width = 3000, height = 3000, pointsize = 35) # open the png device
      else
        w<-h<-15
        if(numrow4$plot4>50){
          h<-h*2
        }
        pdf(file, width = w, height = h) # open the pdf device
      plotInput4()
      dev.off()
    })
  
  observeEvent(input$path, {
    output$PathwayID <- renderUI({
      tn <- tn()
      if(!is.null(tn)){
      tl1 <- tl1()
      mgall <- mgall()
      cat('mgall',mgall,'tn',class(tn),'tl1',tl1,'\n')
      ko_sd <- ko_sd()
      keepcols<-which(names(funtaxall)%in%c(tl1,"ufun","md5", mgall))
      funtax <- funtaxall[,..keepcols]
      if(tl1=='toplevel'){
        keepcols<-c('toplevel',keepcols)
        funtax$toplevel<-factor('toplevel')
        sp.li<-'toplevel'
        #save(funtax,file='funtax.tmp.Rdata')
      }
      funtax <- Intfuntax(funtax,tl1,tn,'toplevel',NULL)
      names(funtax)[names(funtax) == tl1] <- 'usp'
      pathandnames <- as.matrix(getPathwayList(funtax, sp.li =  tn, mgm =  mgall, ko_sd = ko_sd))
      selectInput(inputId = "PathwayID", label = "Input Pathway ID", choices = setNames(as.vector(pathandnames[, "ko"]), pathandnames[,"name"]))
    }})})

  pathw <- reactive({input$PathwayID})

  observeEvent(input$goButton, {
  output$plot4 <- renderD3heatmap({
    tl1 <- tl1()
    tn  <- tn() 
    mgall <-mgall()
    ko_sd <- ko_sd()
    colPal <- colPal()
    if(is.null(tn)){
      downHeat4$is <- FALSE
      return()
    }
    keepcols<-which(names(funtaxall)%in%c(tl1,"ufun","md5", mgall))
    funtax <- funtaxall[,..keepcols]
    names(funtax)[names(funtax) == tl1] <- 'usp'
    obj<-pathwayHeatmap(funtax, tn, mgall, ko_sd)
    colnames(obj)<-as.character(mdt[c(colnames(obj)), colName()])
    mat3 <- plotHeatmap(obj,100,norm = FALSE, log = TRUE,trace = "none")
    mat3[is.na(mat3)] <- 0
    numrow4$plot4 <- dim(mat3)[1]
    if(dim(mat3)[1]>1 & dim(mat3)[2]>1){
      downHeat4$is <- TRUE
      d3heatmap(mat3,dendrogram = chooseDends(mat3), xaxis_height = 220, yaxis_width = 270, yaxis_font_size = "10px", xaxis_font_size = "10px", scalecolors = colPal)
    }else{
      downHeat4$is <- FALSE
      showModal(modalDialog(
        title = titleForDimErrorPopup, textForDimErrorPopup, easyClose = TRUE, footer = NULL
      ))
    }
  })})
  output$dynamic4 <- renderUI({
    d3heatmapOutput("plot4", height = paste0(numrow4$plot4*input$pix4+220, "px"))
  })
    

  plotInput5 <- function(){
    sp.li<- tn()
    tl1 <- tl1()
    pathwi<- pathw()
    mgall <-mgall()
    keepcols<-which(names(funtaxall)%in%c(tl1,"ufun","md5", mgall))
    funtax <- funtaxall[,..keepcols]
    funtax <- Intfuntax(funtax,tl1,sp.li,'toplevel',NULL)
    x <- pathImage(funtax, sp.li, mgall, pathwi,nms = names)
    }
  
  observeEvent(input$down5,{
    showModal(modalDialog(
      title = textForDownloadingMetadata, 
      easyClose = TRUE,
      footer = NULL
    ))
  })

  kostat <- kostat(funtaxall, d.kres)
  output$Pathway <- renderImage({
    sp.li<- tn()
    tl1 <- tl1()
    pathwi<- pathw()
    mgall <-mgall()
    if(!is.null(sp.li)){
      keepcols<-which(names(funtaxall)%in%c(tl1,"ufun","md5", mgall))
    funtax <- funtaxall[,..keepcols]
    if(tl1=='toplevel'){
      keepcols<-c('toplevel',keepcols)
      funtax$toplevel<-factor('toplevel')
      sp.li<-'toplevel'
      #save(funtax,file='funtax.tmp.Rdata')
    }
    funtax <- Intfuntax(funtax,tl1,sp.li,'toplevel',NULL)
    names(funtax)[names(funtax) == tl1] <- 'usp'
    names(funtax)[names(funtax) == tl1] <- 'usp'
    cat(mgall,'-',dim(funtax),'-',names(funtax),'\n-',as.character(mdt[match(mgall,rownames(mdt)), colName()]),'\n')
    names<-as.character(mdt[match(mgall,rownames(mdt)), colName()])
    pathImage(funtax, sp.li, mgall, pathwi, kostat,names)
    }
    list(src = paste0(getwd(),"/","ko", pathwi, ".", sp.li, ".ko.multi.png"),
         contentType = 'png',
         alt = "Press GO to select Pathway!")
  }, deleteFile = FALSE)
  
  
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
  
  observeEvent(input$saveBtn,{
      showModal(modalDialog(
      title = titleForSavingMetadata, textForSavingMetadata, 
      easyClose = TRUE,
      footer = NULL
    ))
    if (!is.null(values[["hot"]])) {
      #write.table(values[["hot"]], "mtcarsss")
      mdt <<- values[["hot"]]
      save(mdt, funtaxall, d.kres, kegg, ko.path.name, file = "pathview.Rdata")
    }
  })
  
  output$hot <- renderRHandsontable({
    if (pressed$a==FALSE){
      DF <- values[["DF"]] 
      if (!is.null(DF)){
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
    textInput("newcolumnname", "Name a new column", sprintf("newcol%s", 1+ncol(values[["DF"]])))
  })

  pressed <- reactiveValues(a=TRUE)
  
  observeEvent(input$addcolumn, {
    pressed$a <- FALSE
    DF <- isolate(values[["DF"]])
    values[["previous"]] <- DF
    newcolumn <- eval(parse(text=sprintf('%s(nrow(DF))', isolate(input$newcolumntype))))
    values[["DF"]] <- setNames(cbind(DF, newcolumn, stringsAsFactors=FALSE), c(names(DF), isolate(input$newcolumnname)))
  })
  
  observeEvent(input$addcolumn2, {
    pressed$a <- TRUE
  })
  
  observeEvent(input$save_changes, {
    tax1selected <-set_taxlevel1()
    tax2selected <-set_taxlevel2()
    func1selected <-set_funlevel1()
    func2selected <-set_funlevel2()
    currentPalette <- colorPalette()
    save(tax1selected, tax2selected, func1selected, func2selected, currentPalette, file = "Settings.Rdata")
  })
}
shinyApp(ui = ui, server = server)