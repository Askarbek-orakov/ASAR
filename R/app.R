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
ui <- fluidPage(
  #selectInput(inputId = "SpecieName", label = "Input Specie Name", as.list(unique(d.bm$usp))),
 textInput(inputId = "SpecieName", label = "Input Specie Name", value = "Geobacter"),
 mainPanel(
   d3heatmapOutput("plot1", width = "150%"),
   tableOutput("table1"),
   d3heatmapOutput("plot2", width = "150%")
 )
)
getSpecieFromAbund<-function(d.bm,sp = SpName,aggregate=FALSE){
  d.res<-d.bm[grep(sp,d.bm$usp),]
  if(aggregate&dim(d.res)[1]>1) d.res<-aggregate(.~ufun,as.data.frame(d.res)[,-1],FUN = sum)
  return(d.res)
}
getSpecieFromAbundMD5<-function(d.bm, sp = SpName, aggregate=FALSE){
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


server <- function(input, output) {
  
  spName <- reactive({input$SpecieName})
   output$plot1 <- renderD3heatmap({SpName <- spName()
    d<-getSpecieFromAbundMD5(d.bm,sp = SpName,aggregate = FALSE)
    obj<-as.matrix(d[,-(1:3)])
    rownames(obj)<-d$ufun
    colnames(obj)<-mdt$MGN
    mat2 <- plotHeatmap(obj,100,norm = FALSE,trace = "none", col = heatmapCols)
    d3heatmap(mat2)
    }) 
  output$table1 <- renderTable({SpName <- spName()
                                d<-getSpecieFromAbundMD5(d.bm,sp = SpName,aggregate = FALSE)
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