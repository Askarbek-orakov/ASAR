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
load("d.bm.Rdata")
tax.df.from.biome <- function(){
  dat <- read_biom("mgm.biome")
  tax <- dat$rows
  d.tax <- lapply(X = tax, FUN = function(ex){list(ex$id, ex$metadata$taxonomy$strain, ex$metadata$taxonomy$species, ex$metadata$taxonomy$genus, ex$metadata$taxonomy$family, ex$metadata$taxonomy$order, ex$metadata$taxonomy$class, ex$metadata$taxonomy$phylum, ex$metadata$taxonomy$domain)})
  dd.tax <- lapply(d.tax, function(x) {
    x[sapply(x, is.null)] <- NA
    return(x)
  })
  tax.df <- data.frame(matrix(unlist(dd.tax), nrow=length(unlist(dd.tax))/9, byrow=T))
  colnames(tax.df) <- list("id", "strain", "species", "genus", "family", "order", "class", "phylum", "domain")
  rownames(tax.df) <- tax.df$id
  tax.df <- tax.df[,-1]
  return(tax.df)
}
getSpecieFromAbund<-function(d.bm,sp = SpName,aggregate=FALSE){
  d.res<-d.bm[grep(sp,d.bm$usp),]
  if(aggregate&dim(d.res)[1]>1) d.res<-aggregate(.~ufun,as.data.frame(d.res)[,-1],FUN = sum)
  return(d.res)
}
getSpecieFromAbundMD5<-function(d.bm, tx=tx, sp = SpName, aggregate=FALSE){
  drops <- c("usp", "species", "genus", "family", "order", "class", "phylum", "domain")
  drops <- drops[drops!= tx]
  dd.res <- taxall[ , !(names(taxall) %in% drops), with = FALSE]
  d.res<-dd.res[grep(sp,dd.res[,get(tx)])]
  #if(get(tx)=="usp") {d.res<-taxall[grep(sp,taxall[,get(tx)]),-(2:8)]
  #} else if(get(tx)=="species") {d.res<-taxall[grep(sp,taxall[,get(tx)]),-c(1,3,4,5,6,7,8)]
  #} else if(get(tx)=="genus") {d.res<-taxall[grep(sp,taxall[,get(tx)]),-c(1,2,4,5,6,7,8)]
  #} else if(get(tx)=="family") {d.res<-taxall[grep(sp,taxall[,get(tx)]),-c(1,2,3,5,6,7,8)]
  #} else if(get(tx)=="order") {d.res<-taxall[grep(sp,taxall[,get(tx)]),-c(1,2,3,4,6,7,8)]
  #} else if(get(tx)=="class") {d.res<-taxall[grep(sp,taxall[,get(tx)]),-c(1,2,3,4,5,7,8)]
  #} else if(get(tx)=="phylum") {d.res<-taxall[grep(sp,taxall[,get(tx)]),-c(1,2,3,4,5,6,8)]
  #} else if(get(tx)=="domain") {d.res<-taxall[grep(sp,taxall[,get(tx)]),-c(1,2,3,4,5,6,7)]}
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

plotSP<-function(d3,sp){
  d<-getSpecieFromAbund(d3,sp = sp,aggregate = FALSE)
  d.sp<-aggregate(.~usp,as.data.frame(d)[,-2],FUN = sum)
  obj<-as.matrix(d.sp[,-1])
  rownames(obj)<-d.sp$usp
  colnames(obj)<-mdt$MGN
  if(dim(obj)[1]>1){
    res<-plotHeatmap(obj,50,trace = "none", col = heatmapCols,main=paste(sp,'\n log-transformed'),norm=FALSE)
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
heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)
tax.df <- tax.df.from.biome()
taxall<- merge(d.bm,tax.df,all=FALSE,by.x='usp',by.y='strain')[,.(usp,species,genus,family,order,class,phylum,domain,md5,ufun,mgm4714659.3,mgm4714661.3,mgm4714663.3,mgm4714665.3,mgm4714667.3,mgm4714669.3,mgm4714671.3,mgm4714673.3,mgm4714675.3,mgm4714677.3,mgm4714679.3)]

ui <- fluidPage(
  selectInput(inputId = "taxlevel", label = "choose taxonomic level",c("strain" = "usp", "species" = "species", "genus" = "genus", "family" = "family", "order" = "order", "class" = "class", "phylum" = "phylum", "domain" = "domain"), selected = "usp"),
  actionButton("do", "GO"),
  uiOutput("taxNames"),
  #textInput(inputId = "SpecieName", label = "Input Specie Name", value = "Geobacter"),
  mainPanel(
    d3heatmapOutput("plot1", width = "150%"),
    tableOutput("table1"),
    d3heatmapOutput("plot2", width = "150%")
  )
)
#as.list(unique(taxall$x))
server <- function(input, output) {
  observeEvent(input$do, { 
    output$taxNames <- renderUI({x <- input$taxlevel
    selectInput(inputId = "SpecieName", label = "Input Specie Name", as.vector(unique(taxall[,get(x)])))
    })})
  txa <- reactive({input$taxlevel})
 
  spName <- reactive({input$SpecieName})
  output$plot1 <- renderD3heatmap({SpName <- spName()
  tx <- txa()
  d<-getSpecieFromAbundMD5(taxall,tx = tx, sp = SpName,aggregate = FALSE)
  obj<-as.matrix(d[,-(1:3)])
  rownames(obj)<-d$ufun
  colnames(obj)<-mdt$MGN
  mat2 <- plotHeatmap(obj,100,norm = FALSE,trace = "none", col = heatmapCols)
  d3heatmap(mat2)
  }) 
  output$table1 <- renderTable({SpName <- spName()
  tx <- txa()
  d<-getSpecieFromAbundMD5(taxall,tx=tx, sp = SpName,aggregate = FALSE)
  obj<-as.matrix(d[,-(1:3)])
  rownames(obj)<-d$ufun
  colnames(obj)<-mdt$MGN
  mat2<-plotHeatmap(obj,100,trace = "none", col = heatmapCols)
  })
  output$plot2 <- renderD3heatmap({SpName <- spName()
  plot <- plotSP(d.bm[,-3],sp = SpName)
  d3heatmap(plot)})
}
shinyApp(ui = ui, server = server)