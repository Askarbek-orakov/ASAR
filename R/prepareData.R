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

