#'Download annotations for all samples and check download completeness
#' 
#'Gets IDs from file$MG.RAST.ID to put it into myGetMgrastAnnotation function to get annotations and then checks each for download completess
#'@param file mdt by default, which should have MGRAST IDs in row MG.RAST.ID
#'@return ind file with indexes of unloaded samples
#'@export
getAllAnnots <- function() {
  ids<-mdt$MG.RAST.ID
  annot<-lapply(ids,myGetMgrastAnnotation,webkey=key)
  ind<-which(!sapply(annot,function(.x)grepl("Download\\s+complete", .x[nrow(.x), 1])))
  ind
}

#'Plots graphs for single specific species
#'
#'@details 1) gets annoations of reads for single species in "d" file
#'@details 2) splits rows with several "md5"s  to separate rows with unique "m5" values in "d5" file
#'@details 3) confers unique rows resulted from merging of "d5" and "d.kres" files at "m5" and "md5" columns respectively. "d.kres" file adds up kegg orthology IDs by each unique "m5" value
#'@details 4) "adk5" file has enrichments on each KEGG Orthology ID for each metagenome sample
#'@details 5) "gdk5" has the same information as "adk5", but for 5 samples only.
#'@details 6) outputs pathway maps for each pathway for each metagenome sample
#'@details 7) generates general heatmap and table of expressions of all kegg functions for all samples for species entered
#'@param 
#'@return mapped KEGG pathway maps
#'@return Heatmap and table of all functions and samples in given specie
#'@export
extractFigures <- function(SpecieName) {
  #gets annotations(usp, ufun, md5, metagenomeID) of reads for single species in "d" file, does not aggregate to generic functions
  d<-getSpecieFromAbundMD5(d.bm,sp = SpecieName,aggregate = FALSE) 
  #  colnam <- as.list(colnames(d.bm[,-c(1:3)]))
  d5<-d[,list(m5=unlist(str_split(md5,',')),mgm4714659.3,mgm4714661.3,mgm4714663.3,mgm4714665.3,mgm4714667.3,mgm4714669.3,mgm4714671.3,mgm4714673.3,mgm4714675.3,mgm4714677.3,mgm4714679.3),by=.(usp,ufun,md5)]
  dk5<-unique(merge(d5,d.kres,all=FALSE,by.x='m5',by.y='md5')[,.(m5,usp,ufun,mgm4714659.3,mgm4714661.3,mgm4714663.3,mgm4714665.3,mgm4714667.3,mgm4714669.3,mgm4714671.3,mgm4714673.3,mgm4714675.3,mgm4714677.3,mgm4714679.3,ko)])
  adk5<-aggregate(.~ko,as.data.frame(dk5[,-c(1:3)]),FUN=sum)
  rownames(adk5)<-adk5$ko
  gdk5<-adk5[,2:6]
  names(gdk5)<-c('MFC1_p','MFC3_p','MFC1_a','MFC3_a','inflowSW')
  gdk5[,5]<-adk5[,1+which(mdt$MFCID==names(gdk5)[5])]
  for(i in 1:4) gdk5[,i]<-rowSums(adk5[,1+which(mdt$MFCID==names(gdk5)[i])])
  
  
  for(pid in pwlist){
    pv.out <- pathview(gene.data = adk5[,c(1,8,3,10)+1], pathway.id = pid,
                       species = "ko", out.suffix = paste0(SpecieName,".ko.plankton"), kegg.native = T,
                       limit = list(gene=range(as.vector(as.matrix(adk5[,c(1,8,3,10)+1]))),cpd=1))
    pv.out <- pathview(gene.data = log2(adk5[,c(1,8,3,10)+1]+1), pathway.id = pid,
                       species = "ko", out.suffix = paste0(SpecieName,".log2.ko.plankton"), kegg.native = T,
                       limit = list(gene=range(as.vector(as.matrix(log2(adk5[,c(1,8,3,10)+1]+1)))),cpd=1))
    pv.out <- pathview(gene.data = adk5[,c(4,11,9,6)+1], pathway.id = pid,
                       species = "ko", out.suffix = paste0(SpecieName,".ko.anode"), kegg.native = T,
                       limit = list(gene=range(as.vector(as.matrix(adk5[,c(4,11,9,6)+1]))),cpd=1))
    pv.out <- pathview(gene.data = log2(adk5[,c(4,11,9,6)+1]+1), pathway.id = pid,
                       species = "ko", out.suffix = paste0(SpecieName,".log2.ko.anode"), kegg.native = T,
                       limit = list(gene=range(as.vector(as.matrix(log2(adk5[,c(4,11,9,6)+1]+1)))),cpd=1))
    pv.out <- pathview(gene.data = gdk5, pathway.id = pid,
                       species = "ko", out.suffix = paste0(SpecieName,".ko"), kegg.native = T,
                       limit = list(gene=range(as.vector(as.matrix(gdk5)))))
    pv.out <- pathview(gene.data = log2(gdk5+1), pathway.id = pid,
                       species = "ko", out.suffix = paste0(SpecieName,".log2.ko"), kegg.native = T,
                       limit = list(gene=range(as.vector(as.matrix(log2(gdk5+1))))))
  }
  obj<-as.matrix(d[,-(1:3)])
  rownames(obj)<-d$ufun
  colnames(obj)<-mdt$MGN
  res<-plotHeatmap(obj,100,trace = "none", col = heatmapCols,main = c(SpecieName, ' functions, \nnorm and log-transformed'))
  showTable(res,SpecieName)
  plotHeatmap(obj,100,norm = FALSE,trace = "none", col = heatmapCols,main=c(SpecieName, ' functions \nlog-transformed'))
  plotSP(d.bm[,-3],sp = SpecieName)
}
#'Looks for one bacterial species from the table in file "d.bm" and if present copies whole row of the species to the table in file "d.res".
#'
#'@param file name of file that contains table of bacterial species, function, md5 and sum of bacteria present in ten metagenomes for bacterial species.
#'@param sp bacterial species.
#'@param aggregate=FALSE 
#'@return file "d.res" that contains species name, function, mdm5 and sum of bacteria present in ten metagenomes for one bacterial species
#'@export
getSpecieFromAbundMD5<-function(d.bm,sp='Geobacter',aggregate=FALSE){
  d.res<-d.bm[grep(sp,d.bm$usp),]
  if(aggregate&dim(d.res)[1]>1) d.res<-aggregate(.~ufun,as.data.frame(d.res)[,-c(1,3)],FUN = sum)
  return(d.res)
}

#'return appropriate object
#'
#'Checks if "obj" is matrix. Then normalizes and/or takes log2 when their corresponding values are TRUE
#'@param obj should be data of matrix class
#'@param norm logical value for normalization
#'@param log logical value for taking logarithm of 2
#'@return processed "res" file
#'@export

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
#' Plot a Heatmap
#' 
#' @export 
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
#'plot Heatmap for species
#'
#'@export
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
#'@example analyzeMatches("mgm4714659")
#'problem: should we leave paste0("ghead -n -1 ./",.x) without .path inclusion?
analyzeMatches <- function(metagenomeID, .path = ".") {
  #makes a list of files for that metagenomeID
  mlist<-dir(path = .path,pattern = c(metagenomeID, ".*"))
  #extract all table except last row
  mannot<-lapply(mlist,function(.x){fread(paste0("ghead -n -1 ./",.x),sep="\t",header = TRUE)})
  #create table with number of rows and columns equal to number of tables imported
  matches<-matrix(0,nrow = length(mannot),ncol = length(mannot))
  #determine number of rows for each table
  for(i in 1:length(mannot)) matches[i,i]<-dim(mannot[[i]])[1]
  # calculate numbers of md5sums between tables in both directions
  for(i in 2:length(mannot)-1){
    for(j in (i+1):length(mannot)){
      matches[i,j]<-length(which(!is.na(match(mannot[[i]]$`hit m5nr id (md5sum)`,mannot[[j]]$`hit m5nr id (md5sum)`))))
      matches[j,i]<-length(which(!is.na(match(mannot[[j]]$`hit m5nr id (md5sum)`,mannot[[i]]$`hit m5nr id (md5sum)`))))
    }
  }
  #give names to rows and columns according to table names.
  rownames(matches)<-mlist
  colnames(matches)<-mlist
  # export tables of matches correlated to itself with 'cor'
  xtable(matches)
  xtable(cor(matches))
}
#devtools::install_github("biomformat", "joey711") NEEDED
tax.df.from.biome <- function(){
  dat <- read_biom("mgm.biome")
  tax <- dat$rows
  d.tax <- lapply(X = tax, FUN = function(ex){list(ex$id, ex$metadata$ncbi_tax_id, ex$metadata$taxonomy$strain, ex$metadata$taxonomy$species, ex$metadata$taxonomy$genus, ex$metadata$taxonomy$family, ex$metadata$taxonomy$order, ex$metadata$taxonomy$class, ex$metadata$taxonomy$phylum, ex$metadata$taxonomy$domain)})
  dd.tax <- lapply(d.tax, function(x) {
         x[sapply(x, is.null)] <- NA
         return(x)
     })
  tax.df <- data.frame(matrix(unlist(dd.tax), nrow=length(unlist(dd.tax))/10, byrow=T))
  colnames(tax.df) <- list("id", "ncbi_tax_id", "strain", "species", "genus", "family", "order", "class", "phylum", "domain")
  return(tax.df)
}
