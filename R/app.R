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
library(png)  # For writePNG function
library(devtools)
#install_github("Alanocallaghan/d3heatmap") #It has color key/color bar
library(d3heatmap)
load("pathview.Rdata")
load("keggmappings.Rdata")
source("global.R")

Intfuntax <- function(funtax, t1, tn, f1, fn, t2=NULL, f2=NULL){
  result2 <- funtax[grep(tn, funtax[,get(t1)])]
  result2 <- result2[grep(fn, result2[,get(f1)])]
  if(!is.null(t2)&!is.null(f2)){
    result2<- ddply(result2, c(t2,f2), numcolwise(sum))
  }else{
    if(!is.null(t2)){
      result2<- ddply(result2, t2, numcolwise(sum))
    }
    if(!is.null(f2)){
      result2<- ddply(result2, f2, numcolwise(sum))
    }}
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
getSpecieFromAbund<-function(d.bm,sp = sp, tx=tx, fun=fun, fN=fN, aggregate=FALSE){
  es<-d.bm[grep(sp,d.bm[,get(tx)])]
  d.res <- es[grep(fN,d.bm[,get(fun)])]
  if(aggregate&dim(d.res)[1]>1) d.res<-aggregate(.~ufun,as.data.frame(d.res)[,-1],FUN = sum)
  return(d.res)
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
  otuIndices = otusToKeep[order(otuStats, decreasing = TRUE)[1:min(c(n,dim(mat)[1]))]]
  mat2 = mat[otuIndices, ]
}
plotSP<-function(d3,sp,tx,tx2,fN, fun){
  d<-getSpecieFromAbund(d3,sp = sp, tx = tx, fN=fN, fun=fun, aggregate = FALSE)
  d.sp<- ddply(d, tx2, numcolwise(sum))
  obj<-as.matrix(d.sp[,-1])
  rownames(obj)<-d.sp[,tx2]
  if(dim(obj)[1]>1){
    res<-plotHeatmap(obj,30,trace = "none", col = heatmapCols,norm=FALSE)
  }else{
    res<-returnAppropriateObj(obj,norm = FALSE,log = TRUE)
  }
  return(res)
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
pathImage<-function(funtax, sp.li, mgm, pathwi) {
  cat("pathImage\n")
  cat(mgm, "\n")
  cat(sp.li, "\n")
  cat(pathwi, "\n")
  d<-getSpecieFromAbundMD5_2(funtax,sp = sp.li,aggregate = FALSE)
  cat(names(d), "\n")
  indC<-which(names(d)%in%c('md5',mgm))
  d5.1<-d[,list(m5=unlist(str_split(md5,','))),by=.(usp,ufun,md5)]
  d5<-merge(d5.1,d[,..indC],by='md5')
  cat(names(d5),"\n")
  dk5<-unique(merge(d5,d.kres,all=FALSE,by.x='m5',by.y='md5')[,-c('md5','.id')])
  indC<-which(names(dk5)%in%c('ko',mgm))
  adk5<-aggregate(.~ko,as.data.frame(dk5[,..indC]),FUN=sum)
  rownames(adk5)<-adk5$ko
  adk5<-adk5[,-1]
  pathview(gene.data = adk5, pathway.id = pathwi,
           species = "ko", out.suffix = paste0(sp.li,".ko"), kegg.native = T,
           limit = list(gene=range(as.vector(as.matrix(adk5))),cpd=1))
}

getpathfromKO <- function(KO){
  temp <- kegg[K == KO]
  temp <- gsub('ko','',paste0(unlist(temp[,"ko"]), collapse = ","))
}
pathwayHeatmap<-function(funtax,sp.lis, mgms) {
  cat("pathwayHeatmap\n")
  cat(mgms, "\n")
  cat(sp.lis, "\n")
  d<-getSpecieFromAbundMD5_2(funtax,sp = sp.lis,aggregate = FALSE)
  cat(names(d), "\n")
  indC<-which(names(d)%in%c('md5',mgms))
  d5.1<-d[,list(m5=unlist(str_split(md5,','))),by=.(usp,ufun,md5)]
  d5<-merge(d5.1,d[,..indC],by='md5')
  cat(names(d5),"\n")
  dk5<-unique(merge(d5,d.kres,all=FALSE,by.x='m5',by.y='md5')[,-c('md5','.id')])
  indC<-which(names(dk5)%in%c('ko',mgms))
  adk5<-aggregate(.~ko,as.data.frame(dk5[,-c('m5', 'usp', 'ufun', 'annotation')]),FUN=sum)
  lastcol<- ncol(adk5)+1
  for (y in 1:nrow(adk5)){adk5[y,lastcol] <- getpathfromKO(adk5[y,"ko"])}
  colnames(adk5)[lastcol] <- "pathwayID"
  #for (y in 1:nrow(adk5)){adk5[y] <- mutate(adk5, pathwayID = paste(as.character(gsub('^path:ko','',matrix(keggLink("pathway", adk5[y,"ko"]), ncol=2, byrow=TRUE)[,2])), collapse = ","))}
  indC<-which(names(adk5)%in%c('ko', mgms))
  adk5 <- as.data.table(adk5)
  adk6<-adk5[,list(pat=unlist(str_split(pathwayID,','))),by=.(ko)]
  a6<-merge(adk6,adk5[,..indC],by='ko')
  a7<-aggregate(.~pat,as.data.frame(a6[,-c('ko')]),FUN=sum)
  rownames(a7)<-a7$pat
  a8<- as.matrix(a7[,-1])
  rownames(a8) <- a7$pat
  return(a8)
}

getpathsfromKOs <- function(KOs){
  setkey(kegg, K)
  temp <- kegg[KOs]
  unlist(str_split(gsub('ko','',paste0(unlist(temp[,"ko"]), collapse = ",")), ','))
  
}
getPathwayList <- function(funtax, sp.li, mgm) {
  cat(mgm, "\n")
  cat(sp.li, "\n")
  d<-getSpecieFromAbundMD5_2(funtax,sp = sp.li,aggregate = FALSE)
  cat(names(d), "\n")
  indC<-which(names(d)%in%c('md5',mgm))
  d5.1<-d[,list(m5=unlist(str_split(md5,','))),by=.(usp,ufun,md5)]
  d5<-merge(d5.1,d[,..indC],by='md5')
  cat(names(d5),"\n")
  dk5<-unique(merge(d5,d.kres,all=FALSE,by.x='m5',by.y='md5')[,-c('md5','.id')])
  kos<- unique(dk5[,"ko"])
  #getpathsfromKOs(unique(dk5[,"ko"]))
  #Since, if num of KOs are more than 500, (function kegglink)it shows error 403. That is why we are doing following:
  if(nrow(kos)>300){
    unikos <- NULL
    kossep <- NULL
    i <- ceiling(nrow(kos)/300)
    for (x in 1:i){
      kossep[[x]]<-kos[((x-1)*300+1):(x*300)]
      unikos[[x]]<-unique(unlist(str_split(paste(as.character(gsub('^path:ko','',matrix(keggLink("pathway", kossep[[x]]$ko), ncol=2, byrow=TRUE)[,2]))),',')))
    }
    unikos<-as.character(unique(unlist(unikos)))
  }else{
    eloop <- NULL
    eloop<-paste(as.character(gsub('^path:ko','',matrix(keggLink("pathway", kos$ko), ncol=2, byrow=TRUE)[,2])))
    unikos <-unique(unlist(str_split(eloop,',')))
  }
}

ui <- fluidPage(
  titlePanel(maintitle),
  shinythemes::themeSelector(),
  navbarPage(theme = "simplex", ""),
  sidebarPanel(
    conditionalPanel(condition = "input.conditionedPanels==2 || input.conditionedPanels==3 || input.conditionedPanels==4 || input.conditionedPanels==5",
                     selectInput(inputId = "mgall", label = metagenomeone, choices = metagenome1n, selected = metagenome1selected, selectize = TRUE, multiple = TRUE)
    ),#setNames(rownames(mdt), mdt[,"MGN"])
    conditionalPanel(condition = "input.conditionedPanels==1",
                     selectInput(inputId = "mg1", label = metagenometwo, choices = metagenome2n, selected = metagenome2selected, selectize = FALSE)),
    conditionalPanel(condition = "input.conditionedPanels==1 || input.conditionedPanels==2 || input.conditionedPanels==3 || input.conditionedPanels==4",
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
                     sliderInput("pix1", "height", value = 400, min = 100, max = 1000)),
    conditionalPanel(condition = "input.conditionedPanels==2",
                     sliderInput("pix2", "height", value = 400, min = 100, max = 1000)),
    conditionalPanel(condition = "input.conditionedPanels==3",
                     sliderInput("pix3", "height", value = 400, min = 100, max = 1000)),
    conditionalPanel(condition = "input.conditionedPanels==4",
                     actionButton("goButton", "GO"),
                     sliderInput("pix4", "height", value = 400, min = 100, max = 1000)
    ),
    conditionalPanel(condition = "input.conditionedPanels==5",
                     actionButton("path", "GO"),
                     uiOutput("PathwayID")),
    # downloadButton('downloadData', 'Download'),
    width = 3),
  
  mainPanel(
    tabsetPanel(
      tabPanel("F&T", uiOutput("dynamic1"), value = 1), 
      tabPanel("F&M", uiOutput("dynamic2"), value = 2),
      tabPanel("T&M", uiOutput("dynamic3"), value = 3),
      tabPanel("Pathway Abundance Heatmap", uiOutput("dynamic4"), value = 4),
      tabPanel("KEGG Pathway Map", imageOutput("Pathway",width = "100%", height = "400px"), value = 5),
      id = "conditionedPanels", 
      tabPanel("Metadata", dataTableOutput("table1"))
    ), width = 9)
)

server <- function(input, output) {
  
    mgall <-reactive({input$mgall})
    mg1   <-reactive({input$mg1})
    tl1   <-reactive({input$tl1})
    tl2   <-reactive({input$tl2})
    fl1   <-reactive({input$fl1})
    fl2   <-reactive({input$fl2})
    tn    <-reactive({input$tn})
    fn    <-reactive({input$fn})
    
    # output$downloadData <- downloadHandler(
    #   filename = "pv.out", content = pv.out , contentType = 'image/png'
    # )
    
    output$taxNames <- renderUI({x <- input$tl1
    selectInput(inputId = "tn", label = taxthree, choices = as.vector(unique(funtaxall[,get(x)])), selected = as.character(funtaxall$genus[(nrow(funtaxall)/2)])) 
    })
    
    output$funNames <- renderUI({y <- input$fl1
    selectInput(inputId = "fn", label = functhree, choices = as.vector(unique(funtaxall[,get(y)])), selected = as.character(funtaxall$FUN4[(nrow(funtaxall)/2)]))
    })
    
    output$plot1 <- renderD3heatmap({
      tl1 <- tl1()
      tl2 <- tl2()
      tn  <- tn()
      fl1 <- fl1()
      fl2 <- fl2()
      fn  <- fn()
      mg1 <- mg1()
      keepcols<-which(names(funtaxall)%in%c(tl1, tl2, fl1, fl2, mg1))
      funtax <- funtaxall[,..keepcols]
      funtax <- Intfuntax(funtax,tl1,tn,fl1,fn,t2 = tl2,f2 = fl2)
      obj <- make2d(funtax)
      obj[is.na(obj)] <- 0
      if(dim(obj)[1]>1){
        res<-plotHeatmap(obj,30,trace = "none", col = heatmapCols,norm=FALSE)
      }else{
        res<-returnAppropriateObj(obj,norm = FALSE,log = TRUE)
      }
      d3heatmap(obj)
    })
  output$dynamic1 <- renderUI({
    d3heatmapOutput("plot1", height = paste0(input$pix1, "px"))
  })
  
    output$plot2 <- renderD3heatmap({
      tl1 <- tl1()
      tn  <- tn()
      fl1 <- fl1()
      fl2 <- fl2()
      fn  <- fn()
      mg2 <- mgall()
      keepcols<-which(names(funtaxall)%in%c(tl1, fl1, fl2, mg2))
      funtax <- funtaxall[,..keepcols]
      funtax <- Intfuntax(funtax,tl1,tn,fl1,fn,f2 = fl2)
      obj <- as.matrix(funtax[,-c(1)])
      rownames(obj)<-funtax[,fl2]
      colnames(obj)<-as.character(mdt[c(gsub('mgm','',mg2)), 3])
      if(dim(obj)[1]>1){
        res<-plotHeatmap(obj,30,trace = "none", col = heatmapCols,norm=FALSE)
      }else{
        res<-returnAppropriateObj(obj,norm = FALSE,log = TRUE)
      }
      d3heatmap(obj)
    })
  output$dynamic2 <- renderUI({
    d3heatmapOutput("plot2", height = paste0(input$pix2, "px"))
  })
  
    output$plot3 <- renderD3heatmap({
      tl1 <- tl1()
      tl2 <- tl2()
      tn  <- tn()
      fl1 <- fl1()
      fn  <- fn()
      mg3 <- mgall()
      keepcols<-which(names(funtaxall)%in%c(tl1, tl2, fl1, mg3))
      funtax <- funtaxall[,..keepcols]
      funtax <- Intfuntax(funtax,tl1,tn,fl1,fn,t2 = tl2)
      obj <- as.matrix(funtax[,-1])
      rownames(obj)<-funtax[,tl2]
      colnames(obj)<-as.character(mdt[c(gsub('mgm','',mg3)), 3])
      if(dim(obj)[1]>1){
        res<-plotHeatmap(obj,30,trace = "none", col = heatmapCols,norm=FALSE)
      }else{
        res<-returnAppropriateObj(obj,norm = FALSE,log = TRUE)
      }
      d3heatmap(res)
    })
  output$dynamic3 <- renderUI({
    d3heatmapOutput("plot3", height = paste0(input$pix3, "px"))
  })
  
  observeEvent(input$path, {
    output$PathwayID <- renderUI({
      tn <- tn()
      tl1 <- tl1()
      mgall <- mgall()
      keepcols<-which(names(funtaxall)%in%c(tl1,"ufun","md5", mgall))
      funtax <- funtaxall[,..keepcols]
      colnames(funtax) <- c("usp","md5","ufun", mgall)
      selectInput(inputId = "PathwayID", label = "Input Pathway ID", as.vector(getPathwayList(funtax, sp.li =  tn, mgm =  mgall)))
    })})

  #sp.lis<- reactive({input$SpecieNames})
  #sp.l<- reactive({input$SpecieN})
  pathw <- reactive({input$PathwayID})

  observeEvent(input$goButton, {
  output$plot4 <- renderD3heatmap({
    tl1 <- tl1()
    tn  <- tn() 
    mgall <-mgall()
    keepcols<-which(names(funtaxall)%in%c(tl1,"ufun","md5", mgall))
    funtax <- funtaxall[,..keepcols]
    colnames(funtax) <- c("usp","md5","ufun", mgall)
    obj<-pathwayHeatmap(funtax, tn, mgall)
    colnames(obj)<-as.character(mdt[c(gsub('mgm','',mgall)), 3])
    mat3 <- plotHeatmap(obj,100,norm = FALSE, log = FALSE,trace = "none", col = heatmapCols)
    d3heatmap(mat3)
  })})
  output$dynamic4 <- renderUI({
    d3heatmapOutput("plot4", height = paste0(input$pix4, "px"))
  })
    
  output$Pathway <- renderImage({
    sp.li<- tn()
    tl1 <- tl1()
    pathwi<- pathw()
    mgall <-mgall()
    keepcols<-which(names(funtaxall)%in%c(tl1,"ufun","md5", mgall))
    funtax <- funtaxall[,..keepcols]
    colnames(funtax) <- c("usp","md5","ufun", mgall)
    pathImage(funtax, sp.li, mgall, pathwi)
    cat(paste0(getwd(),"/","ko", pathwi, ".", sp.li, ".ko.multi.png"))
    list(src = paste0(getwd(),"/","ko", pathwi, ".", sp.li, ".ko.multi.png"),
         contentType = 'png',
         # width = "100%", 
         # height = "400px",
         alt = "Press GO to select Pathway!")
  }, deleteFile = FALSE)
  
  output$table1 <- renderDataTable(as.matrix(mdt))
}
shinyApp(ui = ui, server = server)