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

Intfuntax <- function(funtax, t1, tn, f1, fn, t2=NULL, f2=NULL){
  result <- funtax[grep(tn, funtax[,get(t1)])]
  result2 <- result[grep(fn, result[,get(f1)])]
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
  mat2 = mat[otuIndices, ]}

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
    conditionalPanel(condition = "input.conditionedPanels==1",
                     selectInput(inputId = "tl1_1", label = "Choose taxlevel 1",c("strain" = "usp", "species" = "species", "genus" = "genus", "family" = "family", "order" = "order", "class" = "class", "phylum" = "phylum", "domain" = "domain"), selected = "usp", selectize = FALSE),
                     actionButton("do1", "GO"),
                     uiOutput("taxNames1"),
                     selectInput(inputId = "tl1_2", label = "Aggregation taxlevel 2",c("strain" = "usp", "species" = "species", "genus" = "genus", "family" = "family", "order" = "order", "class" = "class", "phylum" = "phylum"), selected = "usp"),
                     
                     selectInput(inputId = "fl1_1", label = "Choose funLevel 1",c("level 1" = "ufun", "level 2" = "FUN2", "level 3" = "FUN3", "level 4" = "FUN4"), selected = "FUN4", selectize = FALSE),
                     actionButton("fundo1", "GO"),
                     uiOutput("funNames1"),
                     selectInput(inputId = "fl1_2", label = "Choose funLevel 2",c("level 1" = "ufun", "level 2" = "FUN2", "level 3" = "FUN3", "level 4" = "FUN4"), selected = "ufun"),
                     
                     selectInput(inputId = "mg1", label = "Choose one metagenome sample", choices = c(colnames(d.bm[,-c(1:3)])), selected = NULL, selectize = FALSE),
                     
                     actionButton("bh1", "Build heatmap!")
    ),
    conditionalPanel(condition = "input.conditionedPanels==2",
                     selectInput(inputId = "mg2", label = "Choose one/several metagenome sample", choices = c(colnames(d.bm[,-c(1:3)])), selected = NULL, selectize = TRUE, multiple = TRUE),
                     
                     selectInput(inputId = "fl2_1", label = "Choose funLevel 1",c("level 1" = "ufun", "level 2" = "FUN2", "level 3" = "FUN3", "level 4" = "FUN4"), selected = "FUN4", selectize = FALSE),
                     actionButton("fundo2", "GO"),
                     uiOutput("funNames2"),
                     selectInput(inputId = "fl2_2", label = "Choose funLevel 2",c("level 1" = "ufun", "level 2" = "FUN2", "level 3" = "FUN3", "level 4" = "FUN4"), selected = "ufun"),
                     
                     selectInput(inputId = "tl2_1", label = "Choose taxlevel",c("strain" = "usp", "species" = "species", "genus" = "genus", "family" = "family", "order" = "order", "class" = "class", "phylum" = "phylum", "domain" = "domain"), selected = "usp", selectize = FALSE),
                     actionButton("do2", "GO"),
                     uiOutput("taxNames2"),
                     
                     actionButton("bh2", "Build heatmap!")
    ),
    conditionalPanel(condition = "input.conditionedPanels==3",
                     selectInput(inputId = "mg3", label = "Choose one/several metagenome sample", choices = c(colnames(d.bm[,-c(1:3)])), selected = NULL, selectize = TRUE, multiple = TRUE),
                     
                     selectInput(inputId = "tl3_1", label = "Choose taxlevel 1",c("strain" = "usp", "species" = "species", "genus" = "genus", "family" = "family", "order" = "order", "class" = "class", "phylum" = "phylum", "domain" = "domain"), selected = "usp", selectize = FALSE),
                     actionButton("do3", "GO"),
                     uiOutput("taxNames3"),
                     selectInput(inputId = "tl3_2", label = "Aggregate taxlevel 2",c("strain" = "usp", "species" = "species", "genus" = "genus", "family" = "family", "order" = "order", "class" = "class", "phylum" = "phylum"), selected = "usp"),
                     
                     selectInput(inputId = "fl3_1", label = "Choose funlevel",c("level 1" = "ufun", "level 2" = "FUN2", "level 3" = "FUN3", "level 4" = "FUN4"), selected = "FUN4", selectize = FALSE),
                     actionButton("fundo3", "GO"),
                     uiOutput("funNames3"),
                     
                     actionButton("bh3", "Build heatmap!")
    ),
    conditionalPanel(condition = "input.conditionedPanels==4",
                     selectInput(inputId = "SpecieNames", "Choose Specie", as.vector(unique(d.bm[,"usp"]))),
                     selectInput(inputId = "Metagenome", label = "Select Multiple Metagenome Samples", choices = c(colnames(d.bm[,-c(1:3)])), selected = NULL, selectize = TRUE, multiple = TRUE),
                     actionButton("goButton", "GO")
    ),
    
    conditionalPanel(condition = "input.conditionedPanels==5",
                     selectInput(inputId = "Metagenomes", label = "Select Maximum of 5 Metagenome Samples(increase in sample number slows down the process)", choices = c(colnames(d.bm[,-c(1:3)])), selected = NULL, selectize = TRUE, multiple = TRUE),
                     selectInput(inputId = "SpecieN", "Choose Specie", as.vector(unique(d.bm[,"usp"])), selected = NULL),
                     actionButton("path", "GO"),
                     uiOutput("PathwayID")
    ),
    width = 3),
  
  mainPanel(
    tabsetPanel(
      tabPanel("F&T", d3heatmapOutput("plot1"), value = 1), 
      tabPanel("F&M", d3heatmapOutput("plot2", width = "100%", height = "1500px"), value = 2),
      tabPanel("T&M", d3heatmapOutput("plot3", width = "100%", height = "1500px"), value = 3),
      tabPanel("Pathway Abundance Heatmap", d3heatmapOutput("plot4",width = "100%", height = "1500px"), value = 4),
      tabPanel("KEGG Pathway Map", imageOutput("Pathway",width = "100%", height = "400px"), value = 5),
      id = "conditionedPanels"
    ), width = 9)
)

server <- function(input, output) {
  observeEvent(input$do1, { 
    output$taxNames1 <- renderUI({x <- input$tl1_1
    selectInput(inputId = "tn1", label = "Select taxon", as.vector(unique(funtaxall[,get(x)])))
    })})
  observeEvent(input$fundo1, { 
    output$funNames1 <- renderUI({y <- input$fl1_1
    selectInput(inputId = "fn1", label = "Select function", as.vector(unique(funtaxall[,get(y)])))
    })})
  
  observeEvent(input$do2, { 
    output$taxNames2 <- renderUI({x <- input$tl2_1
    selectInput(inputId = "tn2", label = "Select taxon", as.vector(unique(funtaxall[,get(x)])), selected = NULL)
    })})
  observeEvent(input$fundo2, { 
    output$funNames2 <- renderUI({y <- input$fl2_1
    selectInput(inputId = "fn2", label = "Select function", as.vector(unique(funtaxall[,get(y)])), selected = NULL)
    })})
  
  observeEvent(input$do3, { 
    output$taxNames3 <- renderUI({x <- input$tl3_1
    selectInput(inputId = "tn3", label = "Select taxon", as.vector(unique(funtaxall[,get(x)])), selected = NULL)
    })})
  observeEvent(input$fundo3, { 
    output$funNames3 <- renderUI({y <- input$fl3_1
    selectInput(inputId = "fn3", label = "Select function", as.vector(unique(funtaxall[,get(y)])), selected = NULL)
    })})
  
  observeEvent(input$bh1, {
    output$plot1 <- renderD3heatmap({
      tl1_1 <- input$tl1_1
      tl1_2 <- input$tl1_2
      tn1   <- input$tn1
      fl1_1 <- input$fl1_1
      fl1_2 <- input$fl1_2
      fn1   <- input$fn1
      mg1   <- input$mg1
      keepcols<-which(names(funtaxall)%in%c(tl1_1, tl1_2, fl1_1, fl1_2, mg1))
      funtax <- funtaxall[,..keepcols]
      funtax <- Intfuntax(funtax,tl1_1,tn1,fl1_1,fn1,t2 = tl1_2,f2 = fl1_2)
      obj <- make2d(funtax)
      obj[is.na(obj)] <- 0
      if(dim(obj)[1]>1){
        res<-plotHeatmap(obj,30,trace = "none", col = heatmapCols,norm=FALSE)
      }else{
        res<-returnAppropriateObj(obj,norm = FALSE,log = TRUE)
      }
      d3heatmap(obj,Rowv = FALSE,Colv=FALSE)
    })})
  
  observeEvent(input$bh2, {
    output$plot2 <- renderD3heatmap({
      tl2_1 <- input$tl2_1
      tn2   <- input$tn2
      fl2_1 <- input$fl2_1
      fl2_2 <- input$fl2_2
      fn2   <- input$fn2
      mg2   <- input$mg2
      keepcols<-which(names(funtaxall)%in%c(tl2_1, fl2_1, fl2_2, mg2))
      funtax <- funtaxall[,..keepcols]
      funtax <- Intfuntax(funtax,tl2_1,tn2,fl2_1,fn2,f2 = fl2_2)
      obj <- as.matrix(funtax[,-c(1)])
      rownames(obj)<-funtax[,fl2_2]
      colnames(obj)<-as.character(mg2)
      if(dim(obj)[1]>1){
        res<-plotHeatmap(obj,30,trace = "none", col = heatmapCols,norm=FALSE)
      }else{
        res<-returnAppropriateObj(obj,norm = FALSE,log = TRUE)
      }
      d3heatmap(obj,Rowv = FALSE,Colv=FALSE)
    })})
  observeEvent(input$bh3, {
    output$plot3 <- renderD3heatmap({
      tl3_1 <- input$tl3_1
      tl3_2 <- input$tl3_2
      tn3   <- input$tn3
      fl3_1 <- input$fl3_1
      fn3   <- input$fn3
      mg3   <- input$mg3
      keepcols<-which(names(funtaxall)%in%c(tl3_1, tl3_2, fl3_1, mg3))
      funtax <- funtaxall[,..keepcols]
      funtax <- Intfuntax(funtax,tl3_1,tn3,fl3_1,fn3,t2 = tl3_2)
      obj <- as.matrix(funtax[,-1])
      rownames(obj)<-funtax[,tl3_2]
      colnames(obj)<-as.character(mg3)
      if(dim(obj)[1]>1){
        res<-plotHeatmap(obj,30,trace = "none", col = heatmapCols,norm=FALSE)
      }else{
        res<-returnAppropriateObj(obj,norm = FALSE,log = TRUE)
      }
      d3heatmap(res,Rowv = FALSE,Colv=FALSE)
    })
  })
  observeEvent(input$path, {
    output$PathwayID <- renderUI({y <- input$SpecieN
    x <- input$Metagenomes
    selectInput(inputId = "PathwayID", label = "Input Pathway ID", as.vector(getPathwayList(sp.li =  y, mgm =  x)))
    })})
  
  # txa <- reactive({input$taxlevel})
  # txa2 <- reactive({input$taxlevel2})
  # spName <- reactive({input$SpecieName})
  # 
  # fun <- reactive({input$funlevel})
  # fun2 <- reactive({input$funlevel2})
  # funName <- reactive({input$FunctionName})
  # 
  sp.lis<- reactive({input$SpecieNames})
  mgms <-reactive({input$Metagenome})
  
  sp.l<- reactive({input$SpecieN})
  mg <-reactive({input$Metagenomes})
  pathw <- reactive({input$PathwayID})
  
  # output$plot1 <- renderD3heatmap({SpName <- spName()
  # tx <- txa()
  # d<-getSpecieFromAbundMD5(taxall,tx = tx, sp = SpName,aggregate = TRUE)
  # obj<-as.matrix(d[,-1])
  # rownames(obj)<-d$ufun
  # colnames(obj)<-mdt$MGN
  # mat2 <- plotHeatmap(obj,100,norm = FALSE,trace = "none", col = heatmapCols)
  # d3heatmap(mat2)
  # }) 
  
  # output$table1 <- renderTable({SpName <- spName()
  # tx <- txa()
  # d<-getSpecieFromAbundMD5(taxall,tx=tx, sp = SpName,aggregate = TRUE)
  # obj<-as.matrix(d[,-1])
  # rownames(obj)<-d$ufun
  # colnames(obj)<-mdt$MGN
  # mat2<-plotHeatmap(obj,100,trace = "none", col = heatmapCols)
  # })
  
  # output$plot2 <- renderD3heatmap({SpName <- spName()
  # fun <- fun()
  # fun2 <- fun2()
  # FunName <- funName()
  # tx <- txa()
  # tx2 <- txa2()
  # drops <- c("usp", "species", "genus", "family", "order", "class", "phylum", "domain", "md5", "ufun", "FUN2", "FUN3", "FUN4")
  # drops <- drops[drops!= tx]
  # drops <- drops[drops!= tx2]
  # drops <- drops[drops!= fun]
  # plot <- plotSP(funtaxall[ , !(names(funtaxall) %in% drops), with = FALSE], sp = SpName, tx = tx, tx2 = tx2, fun = fun, fN = FunName)
  # d3heatmap(plot,Rowv = FALSE,Colv=FALSE)
  # })
  # 
  output$plot4 <- renderD3heatmap({
    input$goButton
    sp.lis <- sp.lis()
    mgms <-mgms()
    obj<-pathwayHeatmap(sp.lis, mgms)
    mat3 <- plotHeatmap(obj,100,norm = FALSE, log = FALSE,trace = "none", col = heatmapCols)
    d3heatmap(mat3,Rowv = FALSE,Colv=FALSE)
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