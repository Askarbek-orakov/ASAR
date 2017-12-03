## Loading libraries for standalone run

if(!require(data.table)){stop('Install library "data.table"')}
if(!require(taxize)){stop('Install library "taxize"')}
if(!require(stringr)){stop('Install library "stringr"')}
if(!require(plyr)){stop('Install library "plyr"')}
#load('../data/funtree.rda')

####--------------------------------------------------####


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
#'@seealso @seealso \code{\link{mergeAnnots}}
#'@export
#'@example expandNamesDT(d.merge)
expandNamesDT<-function(d.merge){
  d.merge[,ab:=.N,by=.(`hit m5nr id (md5sum)`)]
  d<-unique(d.merge[,.(`hit m5nr id (md5sum)`,fun,sp,ab)])
  names(d)<-c('md5sum','fun','sp','ab')
  d.uspfun<-d[,expandFunSP.list(fun,sp),by=.(md5sum,ab)]
  return(d.uspfun)
}

expandFunSP.df<-function(fun,sp){
  if(length(fun)>1){
    return(ldply(1:length(fun),.fun=function(.i)expandFunSP(fun[.i],sp[.i])))
  }
  ufun<-gsub('(\\]|\\[)','',
            unlist( strsplit( fun, "\\]; *\\[" ) ))
  usp<-gsub('(\\]|\\[)','',unlist( strsplit( sp, "\\]; *\\[" ) ))
  res<-unique(expand.grid(ufun=ufun,usp=usp,fun=fun,sp=sp,stringsAsFactors = FALSE))
  return(res)
}

expandFunSP.list<-function(fun,sp){
  return(as.list(expandFunSP.df(fun,sp)))
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

load.data.from.file<-function(path,pattern = '*.3.seed$'){
  path<-sub('/*$','/',path)
  list<-dir(path = path,pattern = pattern)
  cat(paste(list,collapse = '\n'))
  loadMGRAST<-function(fname){
    l <-readLines(textConnection(system(paste0('tail -n -1 ',path,fname),intern=TRUE)))
    if(!grepl('Download complete',tail(l))){
      warning('Loading of the file "',fname,'"was not complete\n')
      return(NULL)
    }
    res<-fread(paste0('ghead -n -1 ', paste0(path,fname)),sep='\t',header = TRUE)
    return(res)
  }
  res<-lapply(list,loadMGRAST)
  idx<-which(!sapply(res,is.null))
  res<-res[idx]
  names(res)<-list[idx]
  return(res)
}

#'Loading files output by MG-RAST with functional and taxonomical analysis by SEED.
#'
#' After entering location of a file (path) and name of the file (pattern), this function reads in information from the file as a list.
#' @usage load.kodata.from.file (path)
#' @param path location of a file that should be input.
#' @return list called "ko"
#' @export
#' @example load.sdata.from.file()
#'@return list "sannot"
#' command ">sannot <- load.sdata.from.file()" should be run
load.sdata.from.file <- function(path = '.') {
#  res<-load.data.from.file(path,pattern = '*.3.seed$')
#  idx<-which(!sapply(res,is.null))
#  sannot<-res[idx]
  sannot<-load.data.from.file(path,pattern = '*.3.seed$')
  return(sannot)
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
  fannot<-load.data.from.file(path,pattern = '*.3.fseed$')
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
  ko<-load.data.from.file(path,pattern = '^m.*.ko$')
  return(ko)
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
our.aggregate <- function(d.res) {
  dcast(setDT(d.res),usp+ufun+md5 ~ mgid,value.var = 'sum',fill = 0)->d.bm

  d.sp<-aggregate(.~usp,as.data.frame(d.bm)[,-c(2,3)],FUN = sum)
  d.fun<-aggregate(.~ufun,as.data.frame(d.bm)[,-c(1,3)],FUN = sum)
  return(d.bm)
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
#command "kres.res <- our.merge()" should be run
our.merge <- function(path='.') {
  res<-list()
  kres<-list()
  fannot <- load.fdata.from.file(path = path)
  sannot <- load.sdata.from.file(path = path)
  ko <- load.kodata.from.file(path = path)
  nms<-gsub('.fseed$','',names(fannot))
  for(i in 1:length(fannot)){
    f<-fannot[[i]]
    f$`query sequence id`<-gsub('\\|SEED$','',f$`query sequence id`)
    s<-sannot[[i]]
    s$`query sequence id`<-gsub('\\|SEED$','',s$`query sequence id`)
    d.k<-unique(ko[[i]][,list(`hit m5nr id (md5sum)`,
                              `semicolon separated list of annotations`)])
    #creates 3rd column with accession number itself only
    d.k1<-d.k[,list(`semicolon separated list of annotations`,
                    ko=unlist(gsub('accession=\\[K([0-9]+)\\].*','K\\1',
                                   unlist(str_split(`semicolon separated list of annotations`,
                                                    ';'))))),
              by=.(`hit m5nr id (md5sum)`)]
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

kres.res <- our.merge()
d.res <- make.d.res(kres.res)
d.kres <- make.d.kres(kres.res)
d.bm <- our.aggregate(d.res)
usp<-unique(d.bm$usp)
tempC <- gnr_resolve(names = usp,best_match_only=TRUE,data_source_ids = ncbi,canonical = TRUE)
taxCLS<-classification(unique(tempC$matched_name2),db='ncbi')

txTab<-ldply(taxCLS,.fun = function(.x){
  if(is.na(.x)){
    res<-data.frame(species=NA,
                    genus=NA,
                    family=NA,
                    order=NA,
                    class=NA,
                    phylum=NA,
                    domain=NA)
  }else{
    res<-data.frame(species=ifelse(any(.x$rank=='species'),.x$name[.x$rank=='species'],NA),
                    genus=ifelse(any(.x$rank=='genus'),.x$name[.x$rank=='genus'],NA),
                    family=ifelse(any(.x$rank=='family'),.x$name[.x$rank=='family'],NA),
                    order=ifelse(any(.x$rank=='order'),.x$name[.x$rank=='order'],NA),
                    class=ifelse(any(.x$rank=='class'),.x$name[.x$rank=='class'],'NA'),
                    phylum=ifelse(any(.x$rank=='phylum'),.x$name[.x$rank=='phylum'],'NA'),
                    domain=ifelse(any(.x$rank=='superkingdom'),.x$name[.x$rank=='superkingdom'],NA))
  }
  return(res)})
taxL<-merge(tempC,txTab,by.x='matched_name2',by.y='.id')
taxall<-merge(d.bm,data.table(taxL),by.y = 'user_supplied_name',by.x='usp')
funtaxallL <- merge(taxall, funtree, by.x = 'ufun', by.y = 'FUN1')
funtaxall<-funtaxallL[,c("usp", "species","genus","family","order","class","phylum","domain","md5","ufun","FUN2","FUN3","FUN4",grep('mgm',names(funtaxall),value = TRUE)), with=FALSE]
