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
