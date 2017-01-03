#'Get MG-RAST annotation
#' 
#'Takes in ID of metagenome sample at MG-RAST and webkey to access that and generates a formatted "anno" file containing that annotation.
#'@usage myGetMgrastAnnotation(MetagenomeID, evalue = 5, identity = 60, length = 15, resource = c(source = "KO", type = "ontology"), webkey)
#'@param MetagenomeID Id of your metagenome for which you get annotation from MG-RAST.
#'@param evalue = 5 by default.
#'@param identity = 60 (%) by default.
#'@param length = 15 (bp) by default.
#'@param resource has two parameters: source = "KO" by default and type = "ontology" by default.
#'@param webkey authentication key to access annotation file from MG-RAST.
#'@return Formatted annotation file for sample having ID of "MetagenomeID".
#'@example myGetMgrastAnnotation('mgm4714675.3.ko', evalue = 5, identity = 60, length = 15, resource = c(source = "KO", type = "ontology"), webkey)
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
#'@useage getAnnotationFromFile(file)
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

#'Merge functional and taxonomical analysis data of a read.
#'
#'For each read it merges files of functional and taxonomical analysis from MG-RAST, so that it results one table containing ids, md5, function and taxonomy.
#'@useage mergeAnnots(f,s)
#'@param f a list of SEED function from files for different metagenomes that end as "*.3.fseed" as MG-RAST output.
#'@param s a list of SEED Orthologs from files for different metagenomes that end as "*.3.seed" as MG-RAST output.
#'@return table in file "d.merge" consisting of ids, md5, functional and taxonomical data per one read.
#'@export
mergeAnnots<-function(f,s){
  keycols<-names(f)[1:12]
  setkeyv(f,keycols)
  setkeyv(s,keycols)
  d.merge<-merge(f,s,all=TRUE,suffixes = c('.fun','.sp'))
  names(d.merge)<-gsub('semicolon separated list of annotations.','',names(d.merge))
  return(d.merge)  
}

#' Separate functional analysis data into a new table in a file removing duplicates.
#' 
#'Checks md5sum and removes duplicates and copies functional analysis data into the file "d.uspfun".
#'@useage expandNamesDT(x)
#'@param x file that is output by function "mergeAnnots" and from which functional analysis data is going to be extracted. 
#'@details Copies coloumns, namely `query sequence id`,`hit m5nr id (md5sum)`,fun and sp with all rows from file "d.merge" to the file "d".
#'@details Changes names of coloumns `query sequence id`,`hit m5nr id (md5sum)`,fun and sp to 'id','md5sum','fun'and 'sp' respectively.
#'@details Checks coloumn 'md5sum' and removes duplicates by transfering rows to "d.ufun" and removes quare brackets. 
#'@details Checks coloumn 'md5sum' and 'ufun' and removes duplicates by transfering rows to "d.uspfun" and removes quare brackets. 
#'@return table in the file "d.uspfun" which consists of ids, md5sum, function and species name.  
#'@seealso @seealso \code{\link{mergeAnnots}}
#'@export
#'@example expandNamesDT(d.merge)
expandNamesDT<-function(x){
  d<-x[,.(`query sequence id`,`hit m5nr id (md5sum)`,fun,sp)]
  names(d)<-c('id','md5sum','fun','sp')
  #dt[ , list( pep = unlist( strsplit( pep , ";" ) ) ) , by = pro ]
  unique(d[ , list(ufun = gsub('(\\]|\\[)','',unlist( strsplit( fun, "\\]; *\\[" ) )) ,fun,sp,ab=.N) , by = md5sum ])->d.ufun
  unique(d.ufun[ , list(fun,usp=gsub('(\\]|\\[)','',unlist( strsplit( sp, "\\]; *\\[" ) )) ,sp,ab) , by = .(md5sum,ufun) ])->d.uspfun
  return(d.uspfun)
}
#'Calculates number of reads for bacterial species and function.
#' 
#'First calculates sum of reads ('ab') for every group in bacterial species (usp) and function (ufun).
#'Columns 'species' (usp), 'functions' (ufun) and 'sum' (sum) of reads are returned as a data.table and duplicated rows by all columns are removed.
#'@useage getAbundanceFromDT(d.ab)
#'@param d.ab file that should be input for the calculation.
#'@return file that was input, but adding sum of reads in the table.
#'@export
getAbundanceFromDT<-function(d.ab){
  d.ab<-d.ab[,sum:=sum(ab),by=.(usp,ufun)]
  d.ab<-unique(d.ab[,.(usp,ufun,sum)])
  return(d.ab)
}
#'Calculates number of sequences.  
#'
#'First calculates sum of sequences ('ab') for every group in bacterial species (usp) and function (ufun) and renames coloumn as 'sum'. Then names 'md5sum as a 'md5' and separates elements of it by comma.
#'Columns 'species' (usp), 'functions' (ufun), 'sum' (sum) and md5 of sequences are returned as a data.table and duplicated rows by all columns are removed.
#'@useage getAbundanceMD5FromDT(d.ab)
#'@param d.ab file that should be input for the calculation.
#'@return  file that was input, but adding sum of sequences and sum of md5 in the table.
#'@export
getAbundanceMD5FromDT<-function(d.ab){
  d.ab<-d.ab[,.(sum=sum(ab),md5=paste(md5sum,collapse = ',')),by=.(usp,ufun)]
  d.ab<-unique(d.ab[,.(usp,ufun,sum,md5)])
  return(d.ab)
}

#'Collecting table for certain biological species
#'
#'Looks for certain biological species from the table in input file and if present copies whole row to the table in file "d.res".
#'@useage getSpecieFromAbund(file,sp,aggregate=FALSE)
#'@param file name of file that contains table of bacterial species, function, md5 and sum of bacteria present in metagenomes for bacterial species.
#'@param sp biological species present in metagenome. 
#'@param aggregate=FALSE 
#'@return file "d.res" that contains species name, function, mdm5 and sum of bacteria present in ten metagenomes for one certain biological function.
#'@example getSpecieFromAbund(d.bm,sp='Geobacter',aggregate=FALSE)
#'@export
getSpecieFromAbund<-function(file, sp,aggregate=FALSE){
  d.res<-file[grep(sp,file$usp),]
  if(aggregate&dim(d.res)[1]>1) d.res<-aggregate(.~ufun,as.data.frame(d.res)[,-1],FUN = sum)
  return(d.res)
}
#'Collecting table for certain biological function 
#'
#'Looks for certain biological function from the table in input file and if present copies whole row to the table in file "d.res".
#'@useage getFunctionFromAbund(file,fun,aggregate=FALSE)
#'@param file name of file that contains table of bacterial species, function, md5 and sum of bacteria present in metagenomes for bacterial species.
#'@param fun biological function either by KEGG or SEED. 
#'@param aggregate=FALSE 
#'@return file "d.res" that contains species name, function, mdm5 and sum of bacteria present in ten metagenomes for one certain biological function.
#'@example getFunctionFromAbund(d.bm,fun='protein',aggregate=FALSE)
#'@export
getFunctionFromAbund<-function(file,fun,aggregate=FALSE){
  d.res<-file[grep(fun,file$ufun),]
  if(aggregate&dim(d.res)[1]>1) d.res<-aggregate(.~usp,as.data.frame(d.res)[,-2],FUN = sum)
  return(d.res)
}

#'Collecting table for certain biological species. 
#'
#'Looks for one bacterial species from the table in input file and if present copies whole row of the species to the table in file "d.res".
#'@useage getSpecieFromAbundMD5<-function(file, sp, aggregate=FALSE)
#'@param file name of file that contains table of bacterial species, function, md5 and sum of bacteria present in ten metagenomes for bacterial species.
#'@param sp bacterial species.
#'@param aggregate=FALSE 
#'@return file "d.res" that contains species name, function, mdm5 and sum of bacteria present in ten metagenomes for one bacterial species
#'@export
#'@example getSpecieFromAbundMD5(d.bm,sp='Geobacter',aggregate=FALSE)
getSpecieFromAbundMD5<-function(file, sp, aggregate=FALSE){
  d.res<-file[grep(sp,file$usp),]
  if(aggregate&dim(d.res)[1]>1) d.res<-aggregate(.~ufun,as.data.frame(d.res)[,-c(1,3)],FUN = sum)
  return(d.res)
}

#'Collecting table for certain biological function
#'
#'Looks for certain biological function from the table in input file  and if present copies whole row to the table in file "d.res".
#'@useage getFunctionFromAbundMD5(file,fun,aggregate=FALSE)
#'@param file name of file that contains table of bacterial species, function, md5 and sum of bacteria present in ten metagenomes for bacterial species.
#'@param fun function.
#'@param aggregate=FALSE 
#'@return file "d.res" that contains species name, function, mdm5 and sum of bacteria present in ten metagenomes for one certain function.
#'@example getFunctionFromAbundMD5(d.bm,fun='protein',aggregate=FALSE)
#'@export
getFunctionFromAbundMD5<-function(file,fun,aggregate=FALSE){
  d.res<-file[grep(fun,file$ufun),]
  if(aggregate&dim(d.res)[1]>1) d.res<-aggregate(.~usp,as.data.frame(d.res)[,-c(2,3)],FUN = sum)
  return(d.res)
}