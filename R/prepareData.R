#'Load metadata of metagenome samples
#'
#'Takes in a file containing metadata and assigns source and origin values depending on MetagenomeID
#'@param file file containing metadata of selected samples ususally exported from MG-RAST
#'@return formatted tab-delimited metadata table called "mdt"
#'@export
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
  ids<-as.character(mdt$MG.RAST.ID)
}

#'If one fragment of DNA matches more than one function or taxon, MG-RAST returns list separated by semicolon. 
#'This function replaces square brackets by space and reads string as a string. 
#'@param .x input list of annotations
#'@return list of formatted annotations as string
#'@export
extractOTU<-function(.x){
  .res<-data.frame(otu=gsub('(\\[|\\])','',
                            unlist(
                              strsplit(
                                as.character(
                                  .x$semicolon.separated.list.of.annotations),
                                ';'))))
}

#' Takes in JSON format as R friendly-formaе and saves in "Biom.RData" file
#' 
#' This function is to make JSON format avalable in R and save in .RData format for further use.
#' Biom additonally gives for each strain its taxonomic lineage.
#' @param file file to be read in.
#' @return Biom.RData file containing R friendly data from JSON file.
#' @export
read.biome <- function(file) {
  json_res<-fromJSON(file)
  biom(json_res)->biom_res
  biom_df<-biom_data(biom_res)
  save(biom_df,biom_res,file = 'Biom.RData')
}
#' Read KO annotation and save it in SEED.RData file. 
#' 
#' @details To each splitted data frame applies function "extractOTU" and returns results in a data frame.
#' @details Takes data-frame arguments and combine by rows.
#' @param file 
#' @return file called "SEED.RData"
#' @export
Read.ko <- function() {
  ko.df<-data.frame(mgrast.id='mgm',query.sequence.id='query.sequence.id',
                    hit.m5nr.id..md5sum.='hit.m5nr.id..md5sum.',
                    alignment.length.=1,
                    e.value=1e-5,
                    otu='out')[FALSE,]
  for(sid in ids){
    d1<-fread('mgm4714675.3.ko')
    dd<-ddply(.data = d1[1:11],.(query.sequence.id,hit.m5nr.id..md5sum.,alignment.length.,e.value),.fun = extractOTU)
    seed.df<-rbind(seed.df,cbind(data.frame(mgrast.id=rep(sid,dim(dd)[1])),dd))
    save(seed.df,file = 'SEED.RData')
  }
}

#' 
#' 
#' @details 
#' @details 
#' @param  
#' @return 
#' @export
d.fseed<-fread('tmp/mgm4714679.3.fseed',sep='\t',header = FALSE,skip = 1L)
names(d.fseed)<-c('sid','md5','pers.id','al.len','num.mis','num.gap','qstart','qend','hstart','hend','e.val','b.score','ann.list')

d.sseed<-fread('tmp/mgm4714679.3.seed',sep='\t',header = FALSE,skip = 1L)
names(d.sseed)<-c('sid','md5','pers.id','al.len','num.mis','num.gap','qstart','qend','hstart','hend','e.val','b.score','ann.list')
keycols<-c('sid','md5')
setkeyv(d.fseed,keycols)
setkeyv(d.sseed,keycols)
d.merge<-merge(d.fseed,d.sseed,all=TRUE,suffixes = c('.fun','.sp'))
num.sp<-sapply(lapply(str_split(d.merge$ann.list.sp,';'),unique),length)
num.fun<-sapply(lapply(str_split(d.merge$ann.list.fun,';'),unique),length)
d.merge.ex<-data.frame(id=rep('',sum(num.fun*num.sp)),ann.fun='',ann.sp='',stringsAsFactors = FALSE)
id1<-which(num.fun==1&num.sp==1)
ki<-0
kk<-sum(num.sp[id1]*num.fun[id1])
d.merge.ex[1:kk+ki,]<-data.table(id=paste(d.merge$sid[id1],d.merge$md5[id1],sep='|md5'),ann.fun=d.merge$ann.list.fun[id1],ann.sp=d.merge$ann.list.sp[id1])
ki<-ki+kk
id1<-which(num.fun==1&num.sp>1)
kk<-sum(num.sp[id1]*num.fun[id1])
d.merge.ex[1:kk+ki,]<-data.table(
  id=rep(paste(d.merge$sid[id1],
               d.merge$md5[id1],sep='|md5'),
         num.sp[id1]),
  ann.fun=rep(sapply(str_split(d.merge$ann.list.fun[id1],';'),unique)
              ,num.sp[id1]),
  ann.sp=unlist(sapply(
    str_split(d.merge$ann.list.sp[id1],';'),
    unique)))
ki<-ki+kk
id1<-which(num.fun>1&num.sp==1)
kk<-sum(num.sp[id1]*num.fun[id1])
d.merge.ex[1:kk+ki,]<-data.table(
  id=rep(paste(d.merge$sid[id1],
               d.merge$md5[id1],sep='|md5'),
         num.fun[id1]),
  ann.fun=unlist(sapply(
    str_split(d.merge$ann.list.fun[id1],';'),
    unique)),
  ann.sp=rep(sapply(str_split(d.merge$ann.list.sp[id1],';'),unique)
             ,num.fun[id1])
)
ki<-ki+kk
id1<-which(num.fun>1&num.sp>1)
kk<-sum(num.sp[id1]*num.fun[id1])
d.merge.ex[1:kk+ki,]<-data.table(
  id=rep(paste(d.merge$sid[id1],
               d.merge$md5[id1],sep='|md5'),
         num.fun[id1]*num.sp[id1]),
  ann.fun=unlist(sapply(id1,
                        function(.x) rep(sapply(
                          str_split(d.merge$ann.list.fun[.x],';'),
                          unique),num.sp[.x]))),
  ann.sp=unlist(sapply(id1,
                       function(.x) rep(sapply(
                         str_split(d.merge$ann.list.sp[.x],';'),
                         unique),num.fun[.x])))
)

#' Read SEED species and save it in "SEED.RData" file. 
#' 
#' @details To each splitted data frame applies function "extractOTU" and returns results in a data frame.
#' @details Takes data-frame arguments and combine by rows.
#' @param file 
#' @return file called "SEED.RData"
#' @export
d.sseed<-data.frame(mgrast.id='mgm',query.sequence.id='query.sequence.id',
                    hit.m5nr.id..md5sum.='hit.m5nr.id..md5sum.',
                    alignment.length.=1,
                    e.value=1e-5,
                    otu='out')[FALSE,]
for(sid in ids){
  d1<-read.delim('tmp/mgm4714675.3.seed')
  dd<-ddply(.data = d1[1:11],.(query.sequence.id,hit.m5nr.id..md5sum.,alignment.length.,e.value),.fun = extractOTU)
  seed.df<-rbind(seed.df,cbind(data.frame(mgrast.id=rep(sid,dim(dd)[1])),dd))
  save(seed.df,file = 'SEED.RData')
}