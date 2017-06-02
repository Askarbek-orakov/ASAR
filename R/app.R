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
library(d3heatmap)
library(KEGGREST)
library(png)  # For writePNG function
load("pathview.Rdata") 
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
  mat2 = mat[otuIndices, ]}

plotSP<-function(d3,sp,tx,tx2,fN, fun, fun2){
  d<-getSpecieFromAbund(d3,sp = sp, tx = tx, fN=fN, fun=fun, aggregate = FALSE)
  d.sp<- ddply(d, tx2, numcolwise(sum))
  obj<-as.matrix(d.sp[,-1])
  rownames(obj)<-d.sp[,tx2]
  colnames(obj)<-mdt$MGN
  if(dim(obj)[1]>1){
    res<-plotHeatmap(obj,100,trace = "none", col = heatmapCols,norm=FALSE)
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
# CompareMGonKO <- function(ko,mg1, mg2,tx = tx, sp = SpName){
#   taxon<-funtaxall[grep(sp,funtaxall[,get(tx)])]
#   drops <- c("usp","ufun","md5") 
#   d5<-taxon[,m5:=unlist(str_split(md5,',')),by=.(usp,ufun,md5)]
#   # find function to leave all other than by.. , probably "..."
#   dk5<-unique(merge(d5,d.kres,all=FALSE,by.x='m5',by.y='md5')[,.(m5,usp,ufun,mgm4714659.3,mgm4714661.3,mgm4714663.3,mgm4714665.3,mgm4714667.3,mgm4714669.3,mgm4714671.3,mgm4714673.3,mgm4714675.3,mgm4714677.3,mgm4714679.3,ko)])
#   mg1total <-dk5[dk5$ko == get(ko),]
#   #filter by ko and sum abundances for each metagenome and present as heatmap.
#   
# }

pathImage<-function(sp.li, mgm, pathwi) {
  cat("pathImage\n")
  cat(mgm, "\n")
  cat(sp.li, "\n")
  cat(pathwi, "\n")
  d<-getSpecieFromAbundMD5_2(d.bm,sp = sp.li,aggregate = FALSE)
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
  #plotPNG(func = pathview(gene.data = adk5, pathway.id = pathwi,species = "ko", out.suffix = paste0(sp.li,".ko"), kegg.native = T, limit = list(gene=range(as.vector(as.matrix(adk5))),cpd=1)), filename = tempfile(fileext = ".png"), width = 400, height = 400, res = 72)
  return()
}

pathwayHeatmap<-function(sp.lis, mgms) {
  cat("pathwayHeatmap\n")
  cat(mgms, "\n")
  cat(sp.lis, "\n")
  d<-getSpecieFromAbundMD5_2(d.bm,sp = sp.lis,aggregate = FALSE)
  cat(names(d), "\n")
  indC<-which(names(d)%in%c('md5',mgms))
  d5.1<-d[,list(m5=unlist(str_split(md5,','))),by=.(usp,ufun,md5)]
  d5<-merge(d5.1,d[,..indC],by='md5')
  cat(names(d5),"\n")
  dk5<-unique(merge(d5,d.kres,all=FALSE,by.x='m5',by.y='md5')[,-c('md5','.id')])
  indC<-which(names(dk5)%in%c('ko',mgms))
  adk5<-aggregate(.~ko,as.data.frame(dk5[,-c('m5', 'usp', 'ufun', 'annotation')]),FUN=sum)
  lastcol<- ncol(adk5)+1
  for (y in 1:nrow(adk5)){adk5[y,lastcol] <- paste(as.character(gsub('^path:ko','',matrix(keggLink("pathway", adk5[y,"ko"]), ncol=2, byrow=TRUE)[,2])), collapse = ",")}
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

getPathwayList <- function(sp.li, mgm) {
  cat(mgm, "\n")
  cat(sp.li, "\n")
  d<-getSpecieFromAbundMD5_2(d.bm,sp = sp.li,aggregate = FALSE)
  cat(names(d), "\n")
  indC<-which(names(d)%in%c('md5',mgm))
  d5.1<-d[,list(m5=unlist(str_split(md5,','))),by=.(usp,ufun,md5)]
  d5<-merge(d5.1,d[,..indC],by='md5')
  cat(names(d5),"\n")
  dk5<-unique(merge(d5,d.kres,all=FALSE,by.x='m5',by.y='md5')[,-c('md5','.id')])
  kos <- unique(dk5[,"ko"])
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

# funtree <- read.delim("subsys.txt", header = FALSE, quote = "")
# funtree <- funtree[,-5]
# colnames(funtree) <- c("FUN4", "FUN3", "FUN2", "FUN1")
# funtaxall <- merge(taxall, funtree, by.x = 'ufun', by.y = 'FUN1')[,.(usp,species,genus,family,order,class,phylum,domain,md5,ufun,FUN2,FUN3,FUN4,mgm4714659.3,mgm4714661.3,mgm4714663.3,mgm4714665.3,mgm4714667.3,mgm4714669.3,mgm4714671.3,mgm4714673.3,mgm4714675.3,mgm4714677.3,mgm4714679.3)]

ui <- fluidPage(
  titlePanel("METAGENOMIC ANALYSIS by ASAR"),
  sidebarPanel(
  selectInput(inputId = "taxlevel", label = "Choose Taxonomic Level",c("strain" = "usp", "species" = "species", "genus" = "genus", "family" = "family", "order" = "order", "class" = "class", "phylum" = "phylum", "domain" = "domain"), selected = "usp", selectize = FALSE),
  actionButton("do", "GO"),
  uiOutput("taxNames"),
  p("For Taxonomic Content analysis I want taxon chosen above to be separated by :"),
  selectInput(inputId = "taxlevel2", label = "Choose Another Taxonomic Level",c("strain" = "usp", "species" = "species", "genus" = "genus", "family" = "family", "order" = "order", "class" = "class", "phylum" = "phylum"), selected = "usp"),
  width = 3),
  sidebarPanel(
  selectInput(inputId = "funlevel", label = "Choose Functional Level",c("level 1" = "ufun", "level 2" = "FUN2", "level 3" = "FUN3", "level 4" = "FUN4"), selected = "FUN4", selectize = FALSE),
  actionButton("fundo", "GO"),
  uiOutput("funNames"),
  p("Aggregate selected function by next functional level:"),
  selectInput(inputId = "funlevel2", label = "Choose functional Level",c("level 1" = "ufun", "level 2" = "FUN2", "level 3" = "FUN3", "level 4" = "FUN4"), selected = "ufun")
  , width = 3),
  sidebarPanel(
    selectInput(inputId = "SpecieNames", "Choose Specie", as.vector(unique(d.bm[,"usp"])), selected = NULL),
    selectInput(inputId = "Metagenome", label = "Select Multiple Metagenome Samples", choices = c(colnames(d.bm[,-c(1:3)])), selected = NULL, selectize = TRUE, multiple = TRUE),
    actionButton("goButton", "Go!")
    ,  width = 3),
  
  sidebarPanel(
    selectInput(inputId = "Metagenomes", label = "Select Maximum of 5 Metagenome Samples(increase in sample number slows down the process)", choices = c(colnames(d.bm[,-c(1:3)])), selected = NULL, selectize = TRUE, multiple = TRUE),
    selectInput(inputId = "SpecieN", "Choose Specie", as.vector(unique(d.bm[,"usp"])), selected = NULL),
    actionButton("path", "GO"),
    uiOutput("PathwayID")
    , width = 3),

  mainPanel(
    tabsetPanel(
      tabPanel("Functional Heatmap", d3heatmapOutput("plot1", width = "100%", height = "1500px")), 
      tabPanel("Functional Table", tableOutput("table1")), 
      tabPanel("Taxonomic Content Heatmap", d3heatmapOutput("plot2", width = "100%", height = "1500px")),
      tabPanel("Pathway Abundance Heatmap", d3heatmapOutput("plot3",width = "100%", height = "400px")),
      tabPanel("KEGG Pathway Map", imageOutput("Pathway",width = "100%", height = "400px"))
      ), width = 9)
)
server <- function(input, output) {
  observeEvent(input$do, { 
    output$taxNames <- renderUI({x <- input$taxlevel
    selectInput(inputId = "SpecieName", label = "Input Specie Name", as.vector(unique(taxall[,get(x)])))
    })})
  observeEvent(input$fundo, { 
    output$funNames <- renderUI({y <- input$funlevel
    selectInput(inputId = "FunctionName", label = "Input Function Name", as.vector(unique(funtaxall[,get(y)])))
    })})
  observeEvent(input$path, {
    output$PathwayID <- renderUI({y <- input$SpecieN
    x <- input$Metagenomes
    selectInput(inputId = "PathwayID", label = "Input Pathway ID", as.vector(getPathwayList(sp.li =  y, mgm =  x)))
    })})
  
  txa <- reactive({input$taxlevel})
  txa2 <- reactive({input$taxlevel2})
  spName <- reactive({input$SpecieName})
  
  fun <- reactive({input$funlevel})
  fun2 <- reactive({input$funlevel2})
  funName <- reactive({input$FunctionName})
  
  sp.lis<- reactive({input$SpecieNames})
  mgms <-reactive({input$Metagenome})
  
  sp.l<- reactive({input$SpecieN})
  mg <-reactive({input$Metagenomes})
  pathw <- reactive({input$PathwayID})
  
  output$plot1 <- renderD3heatmap({SpName <- spName()
  tx <- txa()
  d<-getSpecieFromAbundMD5(taxall,tx = tx, sp = SpName,aggregate = TRUE)
  obj<-as.matrix(d[,-1])
  rownames(obj)<-d$ufun
  colnames(obj)<-mdt$MGN
  mat2 <- plotHeatmap(obj,100,norm = FALSE,trace = "none", col = heatmapCols)
  d3heatmap(mat2)
  }) 
  
  output$table1 <- renderTable({SpName <- spName()
  tx <- txa()
  d<-getSpecieFromAbundMD5(taxall,tx=tx, sp = SpName,aggregate = TRUE)
  obj<-as.matrix(d[,-1])
  rownames(obj)<-d$ufun
  colnames(obj)<-mdt$MGN
  mat2<-plotHeatmap(obj,100,trace = "none", col = heatmapCols)
  })
  
  output$plot2 <- renderD3heatmap({SpName <- spName()
  fun <- fun()
  FunName <- funName()
  tx <- txa()
  tx2 <- txa2()
  drops <- c("usp", "species", "genus", "family", "order", "class", "phylum", "domain", "md5", "ufun", "FUN2", "FUN3", "FUN4")
  drops <- drops[drops!= tx]
  drops <- drops[drops!= tx2]
  drops <- drops[drops!= fun]
  plot <- plotSP(funtaxall[ , !(names(funtaxall) %in% drops), with = FALSE], sp = SpName, tx = tx, tx2 = tx2, fun = fun, fN = FunName)
  d3heatmap(plot)})
  
  output$plot3 <- renderD3heatmap({
  input$goButton
  sp.lis <- sp.lis()
  mgms <-mgms()
  obj<-pathwayHeatmap(sp.lis, mgms)
  mat3 <- plotHeatmap(obj,100,norm = TRUE,trace = "none", col = heatmapCols)
  d3heatmap(mat3)
  }) 
  
  output$Pathway <- renderImage({
    sp.li<- sp.l()
    pathwi<- pathw()
    mgm <-mg()
    pathImage(sp.li, mgm, pathwi)
    cat(paste0(getwd(),"/","ko", pathwi, ".", sp.li, ".ko.multi.png"))
    list(src = paste0(getwd(),"/","ko", pathwi, ".", sp.li, ".ko.multi.png"),
         contentType = 'png',
         # width = "100%", 
         # height = "400px",
         alt = "Please wait... We are generating KEGG MAP and saving it in your working directory!")
  }, deleteFile = FALSE)
  
}
shinyApp(ui = ui, server = server)