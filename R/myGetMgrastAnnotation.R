#'Get MG-RAST annotation
#'
#'Takes in ID of metagenome sample at MG-RAST and webkey to access that and generates a formatted "anno" file containing that annotation.
#'@param MetagenomeID Id of your metagenome for which you get annotation from MG-RAST.
#'@param evalue = 5 by default.
#'@param identity = 60 (%) by default.
#'@param length = 15 (bp) by default.
#'@param resource has two parameters: source = "KO" by default and type = "ontology" by default.
#'@param webkey authentication key to access annotation file from MG-RAST.
#'@return Formatted annotation file for sample having ID of "MetagenomeID".
#'@export
myGetMgrastAnnotation<-function(MetagenomeID, evalue = 5, identity = 60, length = 15, 
                                resource = c(source = "KO", type = "ontology"), webkey){
  MetagenomeID <- paste0("mgm", MetagenomeID)
  server.resource <- "http://api.metagenomics.anl.gov/1/annotation/similarity/"
  server.resource <- paste0(server.resource, MetagenomeID)
  message(paste("Loading the annotations form MG-RAST of", 
                MetagenomeID), domain = NA)
  message("The time spent in this step is proportional to the total amount of remote data...")
  param <- list(source = resource["source"], type = resource["type"], 
                evalue = evalue, identity = identity, length = length, 
                auth = webkey)
  anno <- tryCatch(getForm(server.resource, .params = param, 
                           .opts = list(noprogress = FALSE, 
                                        progressfunction = function(down,up) {
                                          cat(paste("\r", 
                                                    "loading", 
                                                    paste0(round(down[2]/1024^2,2), "MB"),
                                                    "..."))
                                          flush.console()
                                        })), error = function(e) {
                                          msg <- conditionMessage(e)
                                          structure(msg, class = "try-error")
                                        })
  if (inherits(anno, "try-error")) {
    warning(anno)
    return(FALSE)
  }
  invalid.source <- which(grepl("Invalid\\s+ontology\\s+source", 
                                anno))
  if (length(invalid.source)) 
    stop("invalid ontology source")
  if (length(which(grepl("insufficient\\s+permissions", 
                         anno)))) 
    stop("invalid webkey")
  anno <- read.delim(textConnection(anno), header = FALSE, 
                     sep = "\t", stringsAsFactor = F)
  return(anno)
  cat("\n", MetagenomeID, "annotation data loading completed")
}
#'Get Annotation from File
#'
#'Takes in File name of metagenome sample annotation data to create formatted "anno" file. 
#'@param file file name containing anotation of your sample.
#'@return Formatted annotation file for metagenome sample annotation from input file.
#'@export
getAnnotationFromFile<-function(file){
  l<-readLines(file)
  if(!grepl("Download\\s+complete", l[length(l)])){
    stop('Download is incomplete')
  }
  anno<-read.delim(file, header = FALSE, fill=TRUE,  sep = "\t", stringsAsFactor = F)
  anno[nrow(anno)+1,1]<-'Download complete'
  return(anno)
}

#'Merge functional and taxanomical analysis data of a read.
#'
#'For each read it merges files of functional and taxonomical analysis from MG-RAST, so that it results one table containing md5, function and taxonomy.
#'@param f a list of SEED function from files for different metagenomes that end as "*.3.fseed" as MG-RAST output.
#'@param s a list of SEED Orthologs from files for different metagenomes that end as "*.3.seed" as MG-RAST output.
#'@return table in file "d.merge" consisting of md5, functional and taxonomical data per one read.
#'@export
mergeAnnots<-function(f,s){
  keycols<-names(f)[1:12]
  setkeyv(f,keycols)
  setkeyv(s,keycols)
  d.merge<-merge(f,s,all=TRUE,suffixes = c('.fun','.sp'))
  names(d.merge)<-gsub('semicolon separated list of annotations.','',names(d.merge))
  return(d.merge)  
}
#'
#'
#'
#'@param .
#'@param .
#'@return 
#'@export
expandNamesDT<-function(d.merge){
  d<-d.merge[,.(`query sequence id`,`hit m5nr id (md5sum)`,fun,sp)]
  names(d)<-c('id','md5sum','fun','sp')
  #dt[ , list( pep = unlist( strsplit( pep , ";" ) ) ) , by = pro ]
  unique(d[ , list(ufun = gsub('(\\]|\\[)','',unlist( strsplit( fun, "\\]; *\\[" ) )) ,fun,sp,ab=.N) , by = md5sum ])->d.ufun
  unique(d.ufun[ , list(fun,usp=gsub('(\\]|\\[)','',unlist( strsplit( sp, "\\]; *\\[" ) )) ,sp,ab) , by = .(md5sum,ufun) ])->d.uspfun
  return(d.uspfun)
}

#'Calculates number of reads. 
#'
#'
#'@param 
#'@param 
#'@return 
#'@export
getAbundanceFromDT<-function(d.ab){
  d.ab<-d.ab[,sum:=sum(ab),by=.(usp,ufun)]
  d.ab<-unique(d.ab[,.(usp,ufun,sum)])
  return(d.ab)
}
#'Calculates number of sequances. 
#'
#'
#'@param 
#'@param 
#'@return table containing all reads, all annotations, function, species, 
#'@export
getAbundanceMD5FromDT<-function(d.ab){
  d.ab<-d.ab[,.(sum=sum(ab),md5=paste(md5sum,collapse = ',')),by=.(usp,ufun)]
  d.ab<-unique(d.ab[,.(usp,ufun,sum,md5)])
  return(d.ab)
}

getSpecieFromAbund<-function(d.bm,sp='Geobacter',aggregate=FALSE){
  d.res<-d.bm[grep(sp,d.bm$usp),]
  if(aggregate&dim(d.res)[1]>1) d.res<-aggregate(.~ufun,as.data.frame(d.res)[,-1],FUN = sum)
  return(d.res)
}
getFunctionFromAbund<-function(d.bm,fun='protein',aggregate=FALSE){
  d.res<-d.bm[grep(fun,d.bm$ufun),]
  if(aggregate&dim(d.res)[1]>1) d.res<-aggregate(.~usp,as.data.frame(d.res)[,-2],FUN = sum)
  return(d.res)
}

#'Looks for one bacterial species from the table in file "d.bm" and if present copies whole row of the species to the table in file "d.res".
#'
#' 
#'@param file name of file that contains table of bacterial species, function, md5 and sum of bacteria present in ten metagenomes for bacterial species.
#'@param sp bacterial species.
#'@param aggregate 
#'@return file "d.res" that contains species name, function, mdm5 and sum of bacteria present in ten metagenomes for one bacterial species
#'@export
getSpecieFromAbundMD5<-function(d.bm,sp='Geobacter',aggregate=FALSE){
  d.res<-d.bm[grep(sp,d.bm$usp),]
  if(aggregate&dim(d.res)[1]>1) d.res<-aggregate(.~ufun,as.data.frame(d.res)[,-c(1,3)],FUN = sum)
  return(d.res)
}

#'@return tabel containing bacterial species, corresponding function and 

getSpecieFromAbundMD5<-function(d.bm,sp='Geobacter',aggregate=FALSE){
  d.res<-d.bm[grep(sp,d.bm$usp),]
  if(aggregate&dim(d.res)[1]>1) d.res<-aggregate(.~ufun,as.data.frame(d.res)[,-c(1,3)],FUN = sum)
  return(d.res)
}
getFunctionFromAbundMD5<-function(d.bm,fun='protein',aggregate=FALSE){
  d.res<-d.bm[grep(fun,d.bm$ufun),]
  if(aggregate&dim(d.res)[1]>1) d.res<-aggregate(.~usp,as.data.frame(d.res)[,-c(2,3)],FUN = sum)
  return(d.res)
}









#make d.res contains all data about metagenome (ids function taxonomy abundance)
# d.kres contains 
{r make.d.res}
d.res<-ldply(.data = res,
             .fun = function(.x){
               ab<-.x$ab;
               ab$mgid=.x$name;
               return(ab)
             })
names(d.res)
d.kres<-unique(ldply(.data = kres,
                     .fun = function(.x){
                       ab<-.x$ab;
                       return(ab)
                     }))
names(d.kres)