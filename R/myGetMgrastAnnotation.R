#'Get MG-RAST annotation
#'
#'@param MetagenomeID Id of your metagenome for which you get annotation from MG-RAST
#'@param evalue = 5 by default
#'@param identity = 60 (%) by default
#'@param length = 15 (bp) by default
#'@param resource has two parameters: source = "KO" by default and type = "ontology" by default
#'@param webkey authentication key to access annotation file from MG-RAST
#'@return annotation file for sample having ID of "MetagenomeID"
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