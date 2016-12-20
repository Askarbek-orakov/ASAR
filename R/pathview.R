#'Download annotations for all samples and check download completeness
#' 
#'Gets IDs from file$MG.RAST.ID to put it into myGetMgrastAnnotation function to get annotations and then checks each for download completess
#'@param file mdt by default, which should have MGRAST IDs in row MG.RAST.ID
#'@return ind file with indexes of unloaded samples
#'@export
getAllAnnots <- function(file = mdt) {
  ids<-file$MG.RAST.ID
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
#'@return formatted annotation file for metagenome sample annotation from input file
#'@export
extractFigures <- function(SpecieName) {
  d<-getSpecieFromAbundMD5(d.bm,sp = SpecieName,aggregate = FALSE)
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
