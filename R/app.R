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
options(shiny.maxRequestSize=10*1024^3)
load("mdt.Rdata")
load("keggmappings.Rdata")
ko.path.name<-ko.path.name[-grep('ko01100',ko.path.name$ko),]
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
    result2 <- result2[grep(tn, result2[,get(t1)])]
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
  d.res<-dd.res[grep(sp,dd.res[,get(tx)])]
  if(aggregate&dim(d.res)[1]>1) d.res<-aggregate(.~ufun,as.data.frame(d.res)[,-c(1,2)],FUN = sum)
  return(d.res)
}
getSpecieFromAbundMD5_2<-function(d.bm,sp=SpName,aggregate=FALSE){
  d.res<-d.bm[grep(sp,d.bm$usp),]
  if(aggregate&dim(d.res)[1]>1) d.res<-aggregate(.~ufun,as.data.frame(d.res)[,-c(1,3)],FUN = sum)
  return(d.res)
}
plotHeatmap<-function(obj,n,norm=TRUE,log=TRUE,fun=sd,...){
  mat = returnAppropriateObj(obj, norm, log)
  otusToKeep = which(rowSums(mat) > 0)
  if(length(otusToKeep)==1){
    otuStats<-fun(mat)
  }else{  
    otuStats = apply(mat[otusToKeep, ], 1, fun)
  }
  otuIndices = otusToKeep[order(otuStats, decreasing = TRUE)[1:min(c(length(otusToKeep),n,dim(mat)[1]))]]
  mat2 = mat[otuIndices, ]
}
returnAppropriateObj<-function(obj, norm, log){
  if(class(obj)!='matrix') stop('Obj should be a matrix')
  res<-obj
  if(log){
    res<-log2(res+1)
  }
  if(norm){
    res<-scale(res)
  }
  return(res)
}
get_ko_data <- function(funtax, taxon, metagenomes) {
  d<-getSpecieFromAbundMD5_2(funtax,sp = taxon,aggregate = FALSE)
  indC<-which(names(d)%in%c('md5',metagenomes))
  d5.1<-d[,list(m5=unlist(str_split(md5,','))),by=.(usp,ufun,md5)]
  d5<-merge(d5.1,d[,..indC],by='md5')
  dk5<-unique(merge(d5,d.kres,all=FALSE,by.x='m5',by.y='md5')[,-c('md5','.id')])
  adk5<-aggregate(.~ko,as.data.frame(dk5[,-c('m5', 'usp', 'ufun', 'annotation')]),FUN=sum)
}
pathImage<-function(funtax, sp.li, mgm, pathwi) {
  adk5<-get_ko_data(funtax, sp.li, mgm)
  rownames(adk5)<-adk5$ko
  adk5<-adk5[,-1]
  pathview(gene.data = adk5, pathway.id = pathwi,
           species = "ko", out.suffix = paste0(sp.li,".ko"), kegg.native = T,
           limit = list(gene=range(as.vector(as.matrix(adk5))),cpd=1))
}
filter_stats <- function(funtax, taxon, metagenomes, sd_cutoff) {
  adk5 <- get_ko_data(funtax, taxon, metagenomes)
  indM <- which(names(adk5)%in%c(metagenomes))
  adk5 <- as.data.table(adk5)
  dk6 <- data.frame(ID = adk5[,"ko"], Means=rowMeans(adk5[,..indM]), SD=rowSds(as.matrix(adk5[,..indM])))
  dk7 <- adk5[which((dk6$Means!=0) & (dk6$SD>sd_cutoff)),]
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
  cat(unlist(head(a7,1)),names(a7),'\n')
  cat(pnameIdx,'\n')
  pnameIdx<-pnameIdx[!is.na(pnameIdx)]
  cat('heatmap',dim(adk5),dim(adk6),dim(a7),length(pnameIdx),'\n')
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
  getpathsfromKOs(unique(dk7[,"ko"]))
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

ui <- fluidPage(
  titlePanel(maintitle),
  sidebarPanel(
    conditionalPanel(condition = "input.conditionedPanels==2 || input.conditionedPanels==3 || input.conditionedPanels==4 || input.conditionedPanels==5",
                     selectInput(inputId = "mgall", label = metagenomeone, choices = setNames(c(colnames(funtaxall)[-c(1:13)]), mdt[,colName]), selected = c(colnames(funtaxall)[metagenome1selected]), selectize = TRUE, multiple = TRUE)),
    conditionalPanel(condition = "input.conditionedPanels==1",
                     selectInput(inputId = "mg1", label = metagenometwo, choices = setNames(c(colnames(funtaxall)[-c(1:13)]), mdt[,colName]), selected = c(colnames(funtaxall)[metagenome2selected]), selectize = FALSE)),
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
    conditionalPanel(condition = "input.conditionedPanels==1 ||input.conditionedPanels==2 || input.conditionedPanels==3 || input.conditionedPanels==4 ||input.conditionedPanels==5 ",
                     radioButtons(inputId = "var3", label = "Select the file type", choices = list("png", "pdf"))
                     ),
    conditionalPanel(condition = "input.conditionedPanels==1",
                     downloadButton(outputId = "down1", label = "Download the heatmap")
                     ),
    conditionalPanel(condition = "input.conditionedPanels==2",
                     downloadButton(outputId = "down2", label = "Download the heatmap")
    ),
    conditionalPanel(condition = "input.conditionedPanels==3",
                     downloadButton(outputId = "down3", label = "Download the heatmap")
    ),
    conditionalPanel(condition = "input.conditionedPanels==4",
                     downloadButton(outputId = "down4", label = "Download the heatmap")
    ),
    conditionalPanel(condition = "input.conditionedPanels==5 || input.conditionedPanels==4",
                     sliderInput("ko_sd", "SD cutoff for KO terms", value = 2, min = 0, max = 20)
                     ),
    conditionalPanel(condition = "input.conditionedPanels==5",
                     downloadButton(outputId = "down5", label = "Download KEGG map")
    ),
    conditionalPanel(condition = "input.conditionedPanels==6",
                     fileInput('InFile', 'Upload previously saved Rdata file.'),
                     actionButton("loadRdata", "Upload Rdata"),
                     textInput("Rdataname","Enter file name for Rdata being saved (with '.Rdata' in the end"),
                     actionButton("saveRdata", "Save current Rdata")
    ),
    conditionalPanel(condition = "input.conditionedPanels==7",
                     h3("Table options"),
                     radioButtons("useType", "Use Data Types", c("FALSE", "TRUE")),
                     h3("Download the metadata"), 
                     div(class='row', 
                         div(class="col-sm-6", 
                             actionButton("save", "Download")),
                         div(class="col-sm-6",
                             radioButtons("fileType", "File type", c("txt", "pdf")))
                     )
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
               selectInput(inputId = "colorPalette", label = "Choose color palette for heatmaps", choices = rownames(brewer.pal.info[which(brewer.pal.info$category=="seq"),]), selected = currentPalette),
               plotOutput("paletteOutput"), 
               actionButton("save_changes", "Save Changes"), value = 6),
      tabPanel("Metadata",rHandsontableOutput("hot"), h3("Add a new column"),
               div(class='row', 
                   div(class="col-sm-5", 
                       uiOutput("ui_newcolname"),
                       actionButton("addcolumn", "Add")),
                   div(class="col-sm-4", 
                       radioButtons("newcolumntype", "Type", c("integer", "double", "character"))),
                   div(class="col-sm-3")
               ), value = 7) #dataTableOutput("table1")
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
  
  observeEvent(input$loadRdata, {
    inFile <- input$InFile
    funtaxall <<- loadRdata(inFile$datapath)
  })
  
  observeEvent(input$saveRdata, {
    save(funtaxall, mdt, file = input$Rdataname)
  })
  output$paletteOutput <- renderPlot({
    display.brewer.all(type = "seq")
  })
  
    mgall <-reactive({input$mgall})
    mg1   <-reactive({input$mg1})
    tl1   <-reactive({input$tl1})
    tl2   <-reactive({input$tl2})
    fl1   <-reactive({input$fl1})
    fl2   <-reactive({input$fl2})
    tn    <-reactive({
      if(tl1()=='toplevel'){
        cat('toplevel\n')
        return('toplevel')
      }else{
        return(input$tn)
      }
    })
    fn    <-reactive({input$fn})
    ko_sd <-reactive({input$ko_sd})
    numrow1 <- reactiveValues(plot1 =0)
    numrow2 <- reactiveValues(plot2 =0)
    numrow3 <- reactiveValues(plot3 =0)
    numrow4 <- reactiveValues(plot4 =0)
    
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
    selectInput(inputId = "tn", label = taxthree, choices = as.vector(unique(funtaxall[,get(x)])), selected = as.character(funtaxall$genus[(nrow(funtaxall)/2)])) 
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
      numrow1$plot1 <- dim(res)[1]
      d3heatmap(res,dendrogram = chooseDends(res), xaxis_height = 220, yaxis_width = 280, yaxis_font_size = "10px", xaxis_font_size = "10px", scalecolors = colPal)
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
    colnames(obj)<-as.character(mdt[c(colnames(obj)), colName])
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
      keepcols<-which(names(funtaxall)%in%c(tl1, fl1, fl2, mg2))
      funtax <- funtaxall[,..keepcols]
      funtax <- Intfuntax(funtax,tl1,tn,fl1,fn,f2 = fl2)
      obj <- funtax[,-c(1)]
      rownames(obj)<-funtax[,fl2]
      colnames(obj)<-as.character(mdt[c(colnames(obj)), colName])
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
      d3heatmap(res,dendrogram = chooseDends(res), xaxis_height = 220, yaxis_width = 270, yaxis_font_size = "10px", xaxis_font_size = "10px", scalecolors = colPal)
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
    colnames(obj)<-as.character(mdt[c(colnames(obj)), colName])
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
      keepcols<-which(names(funtaxall)%in%c(tl1, tl2, fl1, mg3))
      funtax <- funtaxall[,..keepcols]
      funtax <- Intfuntax(funtax,tl1,tn,fl1,fn,t2 = tl2)
      obj <- funtax[,-c(1)]
      rownames(obj)<-funtax[,tl2]
      colnames(obj)<-as.character(mdt[c(colnames(obj)), colName])
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
      
      d3heatmap(res,dendrogram = chooseDends(res), xaxis_height = 220, yaxis_width = 270, yaxis_font_size = "10px", xaxis_font_size = "10px", scalecolors = colPal) 
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
    colnames(obj)<-as.character(mdt[c(colnames(obj)), colName])
    mat3 <- plotHeatmap(obj,100,norm = FALSE, log = FALSE,trace = "none")
    mat3[is.na(mat3)] <- 0
    x <- heatmap.2(mat3,dendrogram = chooseDends(mat3), col = colPal, sepcolor="black", sepwidth=c(0.05,0.05), key=TRUE, keysize=0.75, key.par = list(cex=0.7), symkey=FALSE, density.info="none",cexRow=1,cexCol=1,margins=c(20,30),trace="none",srtCol=50)
  }
  
  output$down4 <- downloadHandler(
    filename =  function() {
      paste(input$filename, input$var3, sep=".")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      if(input$var3 == "png")
        png(file, width = 3000, height = 3000, pointsize = 35) # open the png device
      else
        pdf(file, width = 15, height = 15) # open the pdf device
      plotInput4()
      dev.off()
    })
  
  observeEvent(input$path, {
    output$PathwayID <- renderUI({
      tn <- tn()
      tl1 <- tl1()
      mgall <- mgall()
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
      selectInput(inputId = "PathwayID", label = "Input Pathway ID", as.vector(getPathwayList(funtax, sp.li =  tn, mgm =  mgall, ko_sd = ko_sd)))
    })})

  pathw <- reactive({input$PathwayID})

  observeEvent(input$goButton, {
  output$plot4 <- renderD3heatmap({
    tl1 <- tl1()
    tn  <- tn() 
    mgall <-mgall()
    ko_sd <- ko_sd()
    colPal <- colPal()
    keepcols<-which(names(funtaxall)%in%c(tl1,"ufun","md5", mgall))
    funtax <- funtaxall[,..keepcols]
    names(funtax)[names(funtax) == tl1] <- 'usp'
    obj<-pathwayHeatmap(funtax, tn, mgall, ko_sd)
    colnames(obj)<-as.character(mdt[c(colnames(obj)), colName])
    mat3 <- plotHeatmap(obj,100,norm = FALSE, log = FALSE,trace = "none")
    mat3[is.na(mat3)] <- 0
    numrow4$plot4 <- dim(mat3)[1]
    d3heatmap(mat3,dendrogram = chooseDends(mat3), xaxis_height = 220, yaxis_width = 270, yaxis_font_size = "10px", xaxis_font_size = "10px", scalecolors = colPal)
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
    names(funtax)[names(funtax) == tl1] <- 'usp'
    x <- pathImage(funtax, sp.li, mgall, pathwi)
    }
  
  output$down5 <- downloadHandler(
    filename =  function() {
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      plotInput5()
    })
  
  output$Pathway <- renderImage({
    sp.li<- tn()
    tl1 <- tl1()
    pathwi<- pathw()
    mgall <-mgall()
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
    pathImage(funtax, sp.li, mgall, pathwi)
    cat(paste0(getwd(),"/","ko", pathwi, ".", sp.li, ".ko.multi.png"))
    list(src = paste0(getwd(),"/","ko", pathwi, ".", sp.li, ".ko.multi.png"),
         contentType = 'png',
         alt = "Press GO to select Pathway!")
  }, deleteFile = FALSE)
  
  
  
  values <- reactiveValues()
  
  data1 <- reactive({
    if (!is.null(input$hot)) {
      DF <- hot_to_r(input$hot)
    } else {
      if (is.null(values[["DF"]]))
        DF <- mdt
      else
        DF <- values[["DF"]]
    }
    values[["DF"]] <- DF
    DF
  })
  
  output$hot <- renderRHandsontable({
    DF = data1()
    if (!is.null(DF))
      rhandsontable(DF, useTypes = as.logical(input$useType), stretchH = "all")
  })
  output$ui_newcolname <- renderUI({
    textInput("newcolumnname", "Name", sprintf("newcol%s", 1+ncol(values[["DF"]])))
  })
  observeEvent(input$addcolumn, {
    DF <- isolate(values[["DF"]])
    values[["previous"]] <- DF
    newcolumn <- eval(parse(text=sprintf('%s(nrow(DF))', isolate(input$newcolumntype))))
    values[["DF"]] <- setNames(cbind(DF, newcolumn, stringsAsFactors=FALSE), c(names(DF), isolate(input$newcolumnname)))
  })
  #output$table1 <- renderDataTable(as.matrix(mdt))
  
  observeEvent(input$save, {
    fileType <- isolate(input$fileType)
    finalDF <- isolate(values[["DF"]])
    if(fileType == "txt"){
      dput(finalDF, file=file.path(outdir, sprintf("%s.txt", outfilename)))
    }
    else{
      saveRDS(finalDF, file=file.path(outdir, sprintf("%s.rds", outfilename)))
    }
  }
  )
  
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