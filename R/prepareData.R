

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
#' @param ids list of IDs to read
#' @return file called "SEED.RData"
#' @export
Read.ko <- function(ids) {
  for(sid in ids){
    d1<-fread(c(sid, '.3.ko'))
    dd<-ddply(.data = d1[1:11],.(query.sequence.id,hit.m5nr.id..md5sum.,alignment.length.,e.value),.fun = extractOTU)
    seed.df<-rbind(seed.df,cbind(data.frame(mgrast.id=rep(sid,dim(dd)[1])),dd))
    save(seed.df,file = 'SEED.RData')
  }
}

#' Read SEED species and save it in "SEED.RData" file. 
#' 
#' @details To each splitted data frame applies function "extractOTU" and returns results in a data frame.
#' @details Takes data-frame arguments and combine by rows.
#' @param file 
#' @return file called "SEED.RData"
read.seed.sp <- function(ids, path = '.') {
for(sid in ids){ 
  d1<-read.delim(path, sid, '.3.seed') #reads taxonmy by SEED for metagenomes.
  dd<-ddply(.data = d1[1:11],.(query.sequence.id,hit.m5nr.id..md5sum.,alignment.length.,e.value),.fun = extractOTU) #To each splitted data frame (query.sequence.id,hit.m5nr.id..md5sum.,alignment.length.,e.value) applies function "extractOTU" and returns results in a data frame.
  seed.df<-rbind(seed.df,cbind(data.frame(mgrast.id=rep(sid,dim(dd)[1])),dd)) #
  save(seed.df,file = 'SEED.RData')
}
}

