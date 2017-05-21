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
CompareMGonKO <- function(ko,mg1, mg2,tx = tx, sp = SpName){
  taxon<-funtaxall[grep(sp,funtaxall[,get(tx)])]
  drops <- c("usp","ufun","md5") 
  d5<-taxon[,m5:=unlist(str_split(md5,',')),by=.(usp,ufun,md5)]
  # find function to leave all other than by.. , probably "..."
  dk5<-unique(merge(d5,d.kres,all=FALSE,by.x='m5',by.y='md5')[,.(m5,usp,ufun,mgm4714659.3,mgm4714661.3,mgm4714663.3,mgm4714665.3,mgm4714667.3,mgm4714669.3,mgm4714671.3,mgm4714673.3,mgm4714675.3,mgm4714677.3,mgm4714679.3,ko)])
  mg1total <-dk5[dk5$ko == get(ko),]
  #filter by ko and sum abundances for each metagenome and present as heatmap.
  
}

koTaxaMetagenome<-function(sp.li, mgm, kon) {
  d<-getSpecieFromAbundMD5(d.bm,sp = sp.li,aggregate = FALSE)
  d5<-d[,list(m5=unlist(str_split(md5,',')),mgm),by=.(usp,ufun,md5)]
  dk5<-unique(merge(d5,d.kres,all=FALSE,by.x='m5',by.y='md5')[,.(m5,usp,ufun,mgm,ko)])
  adk5<-aggregate(.~ko,as.data.frame(dk5[,-c(1:3)]),FUN=sum)
  rownames(adk5)<-adk5$ko
  pv.out <- pathview(gene.data = adk5[,c(1,8,3,10)+1], pathway.id = gsub('^K','',kon),
                     species = "ko", out.suffix = paste0(sp.l,".ko.plankton"), kegg.native = T,
                     limit = list(gene=range(as.vector(as.matrix(adk5[,c(1,8,3,10)+1]))),cpd=1))
}

funtree <- read.delim("subsys.txt", header = FALSE, quote = "")
funtree <- funtree[,-5]
colnames(funtree) <- c("FUN4", "FUN3", "FUN2", "FUN1")
funtaxall <- merge(taxall, funtree, by.x = 'ufun', by.y = 'FUN1')[,.(usp,species,genus,family,order,class,phylum,domain,md5,ufun,FUN2,FUN3,FUN4,mgm4714659.3,mgm4714661.3,mgm4714663.3,mgm4714665.3,mgm4714667.3,mgm4714669.3,mgm4714671.3,mgm4714673.3,mgm4714675.3,mgm4714677.3,mgm4714679.3)]

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
    selectInput(inputId = "SpecieN", "Choose Specie", as.vector(unique(taxall[,"genus"])), selected = NULL),
    selectInput(inputId = "Metagenomes", label = "Select Multiple Metagenome Samples", choices = c(colnames(d.bm[,-c(1:3)])), selected = NULL, selectize = TRUE, multiple = TRUE),
    selectInput(inputId = "KONames", "Choose KEGG pathway", as.vector(unique(adk5[,"ko"])), selected = "KO0001", selectize = FALSE),
    actionButton("KOnames", "GO"),
    uiOutput("KONames")
    ),

  mainPanel(
    tabsetPanel(
      tabPanel("Functional Heatmap", d3heatmapOutput("plot1", width = "100%", height = "1500px")), 
      tabPanel("Functional Table", tableOutput("table1")), 
      tabPanel("Taxonomic Content Heatmap", d3heatmapOutput("plot2", width = "100%", height = "1500px")),
      tabPanel("KEGG Pathway Map", imageOutput("image1",width = "100%", height = "400px"))
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
  
  txa <- reactive({input$taxlevel})
  txa2 <- reactive({input$taxlevel2})
  spName <- reactive({input$SpecieName})
  
  fun <- reactive({input$funlevel})
  fun2 <- reactive({input$funlevel2})
  funName <- reactive({input$FunctionName})
  
  sp.l<- reactive({input$SpecieN})
  kn<- reactive({input$KONames})
  mg <-reactive({input$Metagenomes})
  
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
  
  output$image1 <- renderImage({
    outfile <- tempfile(fileext = ".png")
    sp.li<- sp.l()
    kon<- kn()
    mgm <-mg()
    ima <- koTaxaMetagenome(sp.l, mg, kn)
  })
  
}
shinyApp(ui = ui, server = server)