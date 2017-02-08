## THESE ARE FUNCTIONS TO PREPARE "Pathview.RData" file
## if you have required files for each sample, running this code should create pathview.Rdata file containing files needed to visualize paths and make tables and heatmaps.

#'Separate functional analysis data into a new table in a file removing duplicates.
#' 
#'Checks md5sum and removes duplicates and copies functional analysis data into the file "d.uspfun".
#'@usage expandNamesDT(x)
#'@param x file that is output by function "mergeAnnots" and from which functional analysis data is going to be extracted. 
#'@details Copies coloumns, namely `query sequence id`,`hit m5nr id (md5sum)`,fun and sp with all rows from file "d.merge" to the file "d".
#'@details Changes names of coloumns `query sequence id`,`hit m5nr id (md5sum)`,fun and sp to 'id','md5sum','fun'and 'sp' respectively.
#'@details Checks coloumn 'md5sum' and removes duplicates by transfering rows to "d.ufun" and removes quare brackets. 
#'@details Checks coloumn 'md5sum' and 'ufun' and removes duplicates by transfering rows to "d.uspfun" and removes quare brackets. 
#'@return table in the file "d.uspfun" which consists of ids, md5sum, function and species name.  
#'@seealso \code{\link{mergeAnnots}}
#'@export
#'@example expandNamesDT(d.merge)
expandNamesDT<-function(d.merge){
  d<-d.merge[,.(`query sequence id`,`hit m5nr id (md5sum)`,fun,sp)]
  names(d)<-c('id','md5sum','fun','sp')
  #dt[ , list( pep = unlist( strsplit( pep , ";" ) ) ) , by = pro ]
  unique(d[ , list(ufun = gsub('(\\]|\\[)','',unlist( strsplit( fun, "\\]; *\\[" ) )) ,fun,sp,ab=.N) , by = md5sum ])->d.ufun
  unique(d.ufun[ , list(fun,usp=gsub('(\\]|\\[)','',unlist( strsplit( sp, "\\]; *\\[" ) )) ,sp,ab) , by = .(md5sum,ufun) ])->d.uspfun
  return(d.uspfun)
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

#'Load metadata of metagenome samples
#'
#'Takes in a file containing metadata and assigns source and origin values depending on MetagenomeID and transfers it to the file "mdt" as output.
#'@usage load.metadata(file)
#'@param file file containing metadata of selected samples ususally exported from MG-RAST
#'@return formatted tab-delimited metadata table called "mdt"
#'@export
#'@example load.metadata("jobs.tsv")
#'command "mdt <- load.metadata("jobs.tsv")" should be run
load.metadata <- function(file) {
  mdt<-read.delim(file)
  mdt$Name<-sub('^[0-9]+_([^_]+)_[^_]+_([^_]+).*$','\\1_\\2',mdt$Metagenome.Name)
  mdt$EID<-sub('^([^_]+)_.*$','\\1',mdt$Name)
  mdt$MFCID<-NA
  mdt$MFCID[mdt$EID=='S5'|mdt$EID=='S6']<-'MFC1_p'
  mdt$MFCID[mdt$EID=='S7'|mdt$EID=='S8']<-'MFC3_p'
  mdt$MFCID[mdt$EID=='S9'|mdt$EID=='S10']<-'MFC1_a'
  mdt$MFCID[mdt$EID=='S11'|mdt$EID=='S12']<-'MFC3_a'
  mdt$MFCID[mdt$EID=='S9'|mdt$EID=='S10']<-'MFC1_a'
  mdt$MFCID[mdt$EID=='S11'|mdt$EID=='S12']<-'MFC3_a'
  mdt$MFCID[mdt$EID=='S1']<-'inflowSW'
  mdt$MFCID[mdt$EID=='S2']<-'rIS'
  mdt$MFCID[mdt$EID=='S3']<-'rSW'
  mdt$MFCID[mdt$EID=='S4']<-'sMiz'
  mdt$Source<-NA
  samplesToKeep<-grep('_a',mdt$MFCID)
  mdt$Source[samplesToKeep]<-'anode'
  samplesToKeep<-grep('_p',mdt$MFCID)
  mdt$Source[samplesToKeep]<-'plankton'
  samplesToKeep<-grep('^r',mdt$MFCID)
  mdt$Source[samplesToKeep]<-'inoculum'
  samplesToKeep<-grep('^inflow',mdt$MFCID)
  mdt$Source[samplesToKeep]<-'inflow'
  mdt$Origin<-NA
  samplesToKeep<-grep('S(3|5|9|6|10)',mdt$Metagenome.Name)
  mdt$Origin[samplesToKeep]<-'sws'
  samplesToKeep<-grep('S(2|7|11|8|12)',mdt$Metagenome.Name)
  mdt$Origin[samplesToKeep]<-'is'
  mdt$Origin[mdt$MFCID=="inflowSW"]<-'inflow'
  rownames(mdt) <- as.character(mdt[, 1])
  return(mdt)
}
#'Loading files output by MG-RAST with functional analysis by SEED.
#' 
#' After entering location of a file (path) and name of the file (pattern), this function reads in information from the file as a list.
#' @usage load.fdata.from.file(path)
#' @param path location of a file that should be input.
#' @return list called "fannot" 
#' @export
#' @example load.fdata.from.file(path = ".")
#' command ">fannot <- load.fdata.from.file()" should be run
load.fdata.from.file <- function(path = ".") {
  flist<-dir(path = path, pattern = "*.3.fseed$")
  cat(paste(flist,collapse = "\n"))
  fannot<-lapply(flist,function(.x){fread(paste0('ghead -n -1 ./', .x),sep='\t',header = TRUE)})
  return(fannot)
}

#'Loading file KEGG Orthology output by MG-RAST.
#' 
#' After entering location of a file (path) and name of the file (pattern), this function reads in information from the file as a list.
#' @usage load.kodata.from.file (path)
#' @param path location of a file that should be input.
#' @return list called "ko" 
#' @export
#' @example load.kodata.from.file()
#'@return list "ko"
#' command ">ko <- load.kodata.from.file()" should be run
load.kodata.from.file <- function(path = '.') {
  klist<-dir(path = path, pattern = '^m.*.ko$')
  cat(paste(klist,collapse = '\n'))
  ko<-lapply(klist,function(.x){fread(paste0('ghead -n -1 ./', .x),sep='\t',header = TRUE)})
}
#'Loading files output by MG-RAST with functional and taxonomical analysis by SEED.
#' 
#' After entering location of a file (path) and name of the file (pattern), this function reads in information from the file as a list.
#' @usage load.sdata.from.file (path)
#' @param path location of a file that should be input.
#' @return list called "ko" 
#' @export
#' @example load.sdata.from.file()
#'@return list "sannot"
#' command ">sannot <- load.sdata.from.file()" should be run
load.sdata.from.file <- function(path = '.') {
  slist<-dir(path = path,pattern = '*.3.seed$')
  flist<-dir(path = path, pattern = "*.3.fseed$")
  cat(paste(slist,collapse = '\n'))
  if(length(slist)!=length(flist)) stop('Length of functional and specie annotation should match\n')
  sannot<-lapply(slist,function(.x){fread(paste0('ghead -n -1 ./', .x),sep='\t',header = TRUE)})
}

#' 
#' 
#' 
#' @details 
#' @usage command "kres.res <- our.merge()" should be run
#' @param path destination of .fseed files. "." by default.
#' @return list of kres and res dataframes
#' @example our.merge()

our.merge <- function(path = ".") {
  flist<-dir(path = path, pattern = "*.3.fseed$")
  nms<-gsub('.fseed$','',flist)
  res<-list()
  kres<-list()
  #fannot <- load.fdata.from.file()
  #sannot <- load.sdata.from.file()
  #ko <- load.kodata.from.file()
  for(i in 1:length(fannot)){
    f<-fannot[[i]]
    #f$`query sequence id`<-gsub('\\|KO$','',f$`query sequence id`)
    s<-sannot[[i]]
    #s$`query sequence id`<-gsub('\\|SEED$','',s$`query sequence id`)
    d.k<-unique(ko[[i]][,list(`hit m5nr id (md5sum)`,`semicolon separated list of annotations`)])
    #creates 3rd column with accession number itself only
    d.k1<-d.k[,list(`semicolon separated list of annotations`,ko=unlist(gsub('accession=\\[K([0-9]+)\\].*','K\\1',unlist(str_split(`semicolon separated list of annotations`,';'))))),by=.(`hit m5nr id (md5sum)`)]
    names(d.k1)<-c('md5','annotation','ko')
    #kres is ready
    kres[[nms[i]]]<-list(ab=d.k1,name=nms[i])
    
    keycols<-names(f)[1:12]
    setkeyv(f,keycols)
    setkeyv(s,keycols)
    d.merge<-merge(f,s,all=TRUE,suffixes = c('.fun','.sp'))
    names(d.merge)<-gsub('semicolon separated list of annotations.','',names(d.merge))
    d.uspfun<-expandNamesDT(d.merge)
    d.ab<-getAbundanceMD5FromDT(d.uspfun)
    #res is ready
    res[[nms[i]]]<-list(ab=d.ab,name=nms[i])
    
    #d.m.t<-table(d.merge$fun,d.merge$sp)
    cat(paste(i,nms[i],'\n'))
  }
  return(list(kres=kres, res=res))
}

#' 
#' 
#' 
#' @details 
#' @usage 
#' @param 
#' @return 
#' @example 
#'@return 
#command "d.res <- make.d.res(kres.res)" should be run
make.d.res <- function(.list) {
  d.res<-ldply(.data = .list$res,
               .fun = function(.x){
                 ab<-.x$ab;
                 ab$mgid=.x$name;
                 return(ab)
               })
}
#' 
#' 
#' 
#' @details 
#' @usage 
#' @param 
#' @return 
#' @example 
#'@return 
#command "d.kres <- make.d.kres(kres.res)" should be run
make.d.kres <- function(.list){
  d.kres<-unique(ldply(.data = .list$kres,
                       .fun = function(.x){
                         ab<-.x$ab;
                         return(ab)
                       }))
}
#' 
#' 
#' 
#' @details 
#' @usage 
#' @param 
#' @return 
#' @example 
#'@return 
our.aggregate <- function() {
  dcast(setDT(d.res),usp+ufun+md5 ~ mgid,value.var = 'sum',fill = 0)->d.bm
  
  d.sp<-aggregate(.~usp,as.data.frame(d.bm)[,-c(2,3)],FUN = sum)
  d.fun<-aggregate(.~ufun,as.data.frame(d.bm)[,-c(1,3)],FUN = sum)
  return(d.bm)
}

#'Extreacting taxonomical levels from the file "mgm.biom"
#'
#'This function creates new table with colounms of eight taxonomical levels. 
#' @details 
#' @usage tax.df.from.biome()
#' @param 
#' @return table "tax.df" containing eight taxonomical levels.
#' @export
#'@example tax.df.from.biome()
#devtools::install_github("biomformat", "joey711") NEEDED
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
## END OF FUNCTIONS FOR PATHVIEW.RDATA



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



#'Listing formatted annotation as a string.
#'
#'If one fragment of DNA matches more than one function or taxon, MG-RAST returns list separated by semicolon. 
#'extractOTU replaces square brackets by space and reads string as a string formatting it for future work. 
#'@usage extractOTU(.x) 
#'@param .x input list of annotations.
#'@return List of formatted annotations as a string.
#'@export
#'@seealso \code{\link{Read.ko}} and \code{\link{read.seed.sp}} that apply extractOTU to splitted data frames.
extractOTU<-function(.x){
  .res<-data.frame(otu=gsub('(\\[|\\])','',
                            unlist(
                              strsplit(
                                as.character(
                                  .x$semicolon.separated.list.of.annotations),
                                ';'))))
}

#' Takes in JSON format as R friendly-format and saves in "Biom.RData" file
#' 
#' This function is to make JSON format avalable in R and save in .RData format for further use.
#' Biom additonally gives for each strain its taxonomic lineage.
#' @usage read.biome(file)
#' @param file file to be read in.
#' @return Biom.RData file containing R friendly data from JSON file.
#' @export
#' @example read.biome('mgm.biome')
read.biome <- function(file) {
  json_res<-fromJSON(file)
  biom(json_res)->biom_res
  biom_df<-biom_data(biom_res)
  save(biom_df,biom_res,file = 'Biom.RData')
}
#' Read KO annotation and  apply function "extractOTU" to splitted data frame. 
#' 
#' To the file for KEGG Orthology this function applies function "extractOTU" to get properly formatted file for the future work.
#' Functions "read.seed.sp" and "Read.ko" basically are same, while former reads file for KEGG Orthology and latter reads file for taxonmy by SEED.
#' @details To each splitted data frame applies function "extractOTU" and returns results in a data frame.
#' @details Takes data-frame arguments and combines by rows.
#' @usage Read.ko(ids)
#' @param ids list of IDs to read
#' @return file called "SEED.RData" 
#' @seealso \code{\link{extractOTU}} and \code{\link{read.seed.sp}}
#' @example Read.ko('mgm4714675.3.ko') 
#' @export
Read.ko <- function(ids) {
  for(sid in ids){
    d1<-fread(c(sid, '.3.ko'))
    dd<-ddply(.data = d1[1:11],.(query.sequence.id,hit.m5nr.id..md5sum.,alignment.length.,e.value),.fun = extractOTU)
    seed.df<-rbind(seed.df,cbind(data.frame(mgrast.id=rep(sid,dim(dd)[1])),dd))
    save(seed.df,file = 'SEED.RData')
  }
}

#' Read SEED species and apply function "extractOTU" to splitted data frame.
#'  
#' To the file for taxonmy by SEED this function applies function "extractOTU" to get properly formatted file for the future work.
#' Functions "read.seed.sp" and "Read.ko" basically are same, while former reads file for KEGG Orthology and latter reads file for taxonmy by SEED.
#' @details Reads taxonmy by SEED for metagenomes.
#' @details To each splitted data frame (query.sequence.id,hit.m5nr.id..md5sum.,alignment.length.,e.value) applies function "extractOTU" and returns results in a data frame.
#' @details Takes data-frame arguments and combines by rows.
#' @usage read.seed.sp(ids, path= '.')
#' @param ids name of the file
#' @param path adress where the file is located with that name 
#' @return file called "SEED.RData"
#' @seealso \code{\link{extractOTU}} and \code{\link{Read.ko}}
#' @example read.seed.sp('mgm4714675.3.seed', 'tmp')
#' @export
read.seed.sp <- function(ids, path = '.') {
  for(sid in ids){ 
    d1<-read.delim(path, sid, '.3.seed') #reads taxonmy by SEED for metagenomes.
    dd<-ddply(.data = d1[1:11],.(query.sequence.id,hit.m5nr.id..md5sum.,alignment.length.,e.value),.fun = extractOTU) #To each splitted data frame (query.sequence.id,hit.m5nr.id..md5sum.,alignment.length.,e.value) applies function "extractOTU" and returns results in a data frame.
    seed.df<-rbind(seed.df,cbind(data.frame(mgrast.id=rep(sid,dim(dd)[1])),dd)) 
    save(seed.df,file = 'SEED.RData')
  }
}





#'Downloading annotations for all samples (metagenomes) 
#' 
#'Gets IDs from file$MG.RAST.ID to put it into myGetMgrastAnnotation function to get annotations and then checks each for download completess.
#'@useage getAllAnnots(file)
#'@param file mdt by default, which should have MGRAST IDs in row MG.RAST.ID
#'@return ind file with indexes of unloaded samples
#'@export
#'@example getAllAnnots(mdt)
getAllAnnots <- function(mdt) {
  ids<-mdt$MG.RAST.ID
  annot<-lapply(ids,myGetMgrastAnnotation,webkey=key)
  ind<-which(!sapply(annot,function(.x)grepl("Download\\s+complete", .x[nrow(.x), 1])))
  ind
}

#'Plots graphs for single specific species
#'
#'After entering biological species name this function will plot mapped KEGG pathway maps and heatmaps for that species.
#'@details 1) gets annoations of reads for single species in "d" file.
#'@details 2) splits rows with several "md5"s  to separate rows with unique "m5" values in "d5" file.
#'@details 3) confers unique rows resulted from merging of "d5" and "d.kres" files at "m5" and "md5" columns respectively. "d.kres" file adds up kegg orthology IDs by each unique "m5" value.
#'@details 4) "adk5" file has enrichments on each KEGG Orthology ID for each metagenome sample.
#'@details 5) "gdk5" has the same information as "adk5", but for 5 samples only.
#'@details 6) outputs pathway maps for each pathway for each metagenome sample.
#'@details 7) generates general heatmap and table of expressions of all KEGG functions for all samples for species entered.
#'@useage extractFigures('Geobacter')
#'@param SpecieName name of biological species for which KEGG pathway maps, heatmap and table of all functions and samples in given specie are going to be plotted.
#'@return mapped KEGG pathway maps
#'@return Heatmap and table of all functions and samples in given specie
#'@example extractFigures('Geobacter')
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


#'Returning appropriate object as a matrix
#'
#'Checks if "obj" is matrix. Then normalizes and/or takes log2 when their corresponding values are TRUE.
#'@useage returnAppropriateObj(obj, norm, log)
#'@param obj should be data of matrix class.
#'@param norm logical value for normalization.
#'@param log logical value for taking logarithm of 2.
#'@return processed "res" file
#'@export
#'@example returnAppropriateObj(obj, norm = FALSE,log = TRUE)
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
#' 
#' @details 
#' @usage plotHeatmap(obj, n, norm, log, fun, ...)
#' @param obj
#' @param n
#' @param norm
#' @param log
#' @param fun
#' @param ...
#' @return 
#' @export
#' @example  plotHeatmap(obj,n,norm=TRUE,log=TRUE,fun=sd,...)
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
#' @details 
#' @usage 
#' @param d3
#' @param sp
#' @param ...
#' @return 
#' @export
#' @example
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
#'
#' @details 
#' @usage 
#' @param obj
#' @param main
#' @param landscape
#' @return 
#' @export
#' @example
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
#' 
#' @details 
#' @usage analyzeMatches(metagenomeID, .path )
#' @param metagenomeID
#' @param .path
#' @return 
#' @export
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

#merges taxonomy with metagenomes' data
taxall<- merge(d.bm,tax.df,all=TRUE,by.x='usp',by.y='strain')[,.(usp,species,genus,family,order,class,phylum,domain,md5,ufun,mgm4714659.3,mgm4714661.3,mgm4714663.3,mgm4714665.3,mgm4714667.3,mgm4714669.3,mgm4714671.3,mgm4714673.3,mgm4714675.3,mgm4714677.3,mgm4714679.3)]
# inputs are: taxlevel and taxName(change in responce to taxlevel)
d.res<-aggregate(.~ufun,as.data.frame(d.res)[,-c(1,3)],FUN = sum)

x.usp <- as.vector(unique(taxall[,"usp"]))

funtree <- read.delim("subsys.txt", header = FALSE, quote = "")
funtree <- funtree[,-5]
colnames(funtree) <- c("FUN4", "FUN3", "FUN2", "FUN1")
funtaxall <- merge(taxall, funtree, by.x = 'ufun', by.y = 'FUN1')[,.(usp,species,genus,family,order,class,phylum,domain,md5,ufun,FUN2,FUN3,FUN4,mgm4714659.3,mgm4714661.3,mgm4714663.3,mgm4714665.3,mgm4714667.3,mgm4714669.3,mgm4714671.3,mgm4714673.3,mgm4714675.3,mgm4714677.3,mgm4714679.3)]

