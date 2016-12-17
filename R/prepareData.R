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
#'
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

#' Takes in JSON format as R friendly-formay and saves in "Biom.RData" file
#' 
#' This function is to make JSON format avalable in R and save in .RData format for further use
#' @param file file to be read in
#' @return Biom.RData file containing R friendly data from JSON file
#' @export
read.biome <- function(file) {
  json_res<-fromJSON(file)
  biom(json_res)->biom_res
  biom_df<-biom_data(biom_res)
  save(biom_df,biom_res,file = 'Biom.RData')
}