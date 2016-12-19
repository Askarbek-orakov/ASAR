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

#'
#'
#'Takes in Filename of metagenome sample annotation data to create formatted "anno" file. 
#'@param file filename containing anotation of your sample
#'@return formatted annotation file for metagenome sample annotation from input file
#'@export
mergeAnnots<-function(f,s){
  keycols<-names(f)[1:12]
  setkeyv(f,keycols)
  setkeyv(s,keycols)
  d.merge<-merge(f,s,all=TRUE,suffixes = c('.fun','.sp'))
  names(d.merge)<-gsub('semicolon separated list of annotations.','',names(d.merge))
  return(d.merge)  
}

#'Analysis of maches
#'
#'Creates a list called "mlist" that contains functional and taxanomical analysis of one metagenome (mgm4714659.*). 
#'@param file filename containing anotation of your sample
#'@return formatted annotation file for metagenome sample annotation from input file
#'@export