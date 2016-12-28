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
ui <- fluidPage(
  
 textInput(inputId = "SpecieName", label = "Input Specie Name", value = "Geobacter"),
 mainPanel(
   tableOutput("table1"),
   plotOutput("plot1"),
   plotOutput("plot2")
 )
  
)

getSpecieFromAbundMD5<-function(d.bm,sp='Geobacter',aggregate=FALSE){
  d.res<-d.bm[grep(sp,d.bm$usp),]
  if(aggregate&dim(d.res)[1]>1) d.res<-aggregate(.~ufun,as.data.frame(d.res)[,-c(1,3)],FUN = sum)
  return(d.res)
}
plotHeatmap<-function(obj,n,norm=TRUE,log=TRUE,fun=sd,...){
  #uses returnAppropriateObj function to get normalized and log(2)-ed obj matrix as mat
  mat = returnAppropriateObj(obj, norm, log)
  #filters to otusToKeep rows which have total sum of more than 0 
  otusToKeep = which(rowSums(mat) > 0)
  #applies sd(standard deviation) function to otusToKeep rows.
  if(length(otusToKeep)==1){
    otuStats<-fun(mat)
  }else{  
    otuStats = apply(mat[otusToKeep, ], 1, fun)
  }
  #otuStats in otusToKeep are ordered in increasing order. ??[1:min(c(n,dim(mat)[1]
  otuIndices = otusToKeep[order(otuStats, decreasing = TRUE)[1:min(c(n,dim(mat)[1]))]]
  #copies otuIndices rows from mat to mat2
  mat2 = mat[otuIndices, ]
  #creates a heatmap with clustered dendrogram
  heatmap.2(mat2, hclustfun = function(.x)hclust(.x,method = 'ward.D2'),srtCol=45,
            key.title=NA,
            key.xlab=NA,
            key.par=list(mgp=c(1.5, 0.5, 0),
                         mar=c(1, 2.5, 1, 0.5)),
            lmat=rbind( c(0, 3, 4), c(2,1,0 ) ),
            lwid=c(1.5, 4, 1.5 ),...)
  
  #srtRow=45,...)
  invisible(mat2)
}
plotSP<-function(d3,sp){
  #filters rows with given species from d3 to d
  d<-getSpecieFromAbund(d3,sp = sp,aggregate = FALSE)
  # aggregates usp column by d values except column 2 by summing
  d.sp<-aggregate(.~usp,as.data.frame(d)[,-2],FUN = sum)
  #takes d.sp without first column, i.e. only matrix body
  obj<-as.matrix(d.sp[,-1])
  #names to columns and rows of obj givem from specie names from d.sp and from metagenome IDs from mdt
  rownames(obj)<-d.sp$usp
  colnames(obj)<-mdt$MGN
  #checks if obj has more than one dimensions and plots Heatmap
  if(dim(obj)[1]>1){
    res<-plotHeatmap(obj,50,trace = "none", col = heatmapCols,main=paste(sp,'\n log-transformed'),norm=FALSE)
  }else{
    res<-returnAppropriateObj(obj,norm = FALSE,log = TRUE)
  }
  invisible(res)
}
showTable<-function(obj,main,landscape=TRUE){
  mat<-as.matrix(obj)
  #order mat rows in decreasing order by rowSums
  mat<-mat[order(rowSums(mat),decreasing = TRUE),]
  addtorow          <- list()
  addtorow$pos      <- list()
  addtorow$pos[[1]] <- c(0)
  addtorow$command  <- c(paste("\\hline \n",
                               "\\endhead \n",
                               "\\hline \n",
                               "\\multicolumn{3}{l}{\\footnotesize Continued on next page} \n",
                               "\\endfoot \n",
                               "\\endlastfoot \n",sep=""))
  #gives prints
  if(landscape){
    cat(sprintf("\\newpage\n \\begin{landscape} \n\\begin{center}\n\\captionof{table}{Summary %s }\n\\scriptsize",main))
  }else{
    cat(sprintf("\\newpage \n\\begin{center}\n\\captionof{table}{Summary %s }\n\\scriptsize",main))
  }
  #column names as vector are passed to ind
  ind<-1:length(colnames(mat))
  #alig is string of 'p{5cm}' and 'r' legth(ind) times
  alig<-c('p{5cm}',rep('r',length(ind)))
  
  print(xtable(mat[,order(colnames(mat))[ind]], 
               align = paste(alig,collapse = ''),digits = 2)
        ,size="small",include.colnames = TRUE,
        tabular.environment="longtable",
        floating=FALSE,include.rownames=TRUE,
        add.to.row = addtorow,hline.after=c(-1))
  if(landscape){
    cat("\\end{center}\n \\end{landscape}")
  }else{
    cat("\\end{center}\n ")
  }
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

server <- function(input, output) {
  SpecieName <- reactive(input$SpecieName)
  d<-getSpecieFromAbundMD5(d.bm,sp = SpecieName,aggregate = FALSE)
  #d5<-d[,list(m5=unlist(str_split(md5,',')),mgm4714659.3,mgm4714661.3,mgm4714663.3,mgm4714665.3,mgm4714667.3,mgm4714669.3,mgm4714671.3,mgm4714673.3,mgm4714675.3,mgm4714677.3,mgm4714679.3),by=.(usp,ufun,md5)]
  #dk5<-unique(merge(d5,d.kres,all=FALSE,by.x='m5',by.y='md5')[,.(m5,usp,ufun,mgm4714659.3,mgm4714661.3,mgm4714663.3,mgm4714665.3,mgm4714667.3,mgm4714669.3,mgm4714671.3,mgm4714673.3,mgm4714675.3,mgm4714677.3,mgm4714679.3,ko)])
  #adk5<-aggregate(.~ko,as.data.frame(dk5[,-c(1:3)]),FUN=sum)
  #rownames(adk5)<-adk5$ko
  
  obj<-as.matrix(d[,-(1:3)])
  rownames(obj)<-d$ufun
  colnames(obj)<-mdt$MGN
  res<-plotHeatmap(obj,100,trace = "none", col = heatmapCols,main = c(SpecieName, ' functions, \nnorm and log-transformed'))
  output$table1 <- renderTable(showTable(res,SpecieName))
  output$plot1 <- renderPlot(plotHeatmap(obj,100,norm = FALSE,trace = "none", col = heatmapCols,main=c(SpecieName, ' functions \nlog-transformed')))
  output$plot2 <- renderPlot(plotSP(d.bm[,-3],sp = SpecieName))
  
}
shinyApp(ui = ui, server = server)