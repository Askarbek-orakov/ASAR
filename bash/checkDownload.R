#!/usr/local/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
print(args)
webkey<-args[1]
proj.ID<-args[2]
library(data.table)
library(stringr)
library(plyr)
library(biomformat)
mdt<-read.csv(file = paste0(proj.ID,'.meta.csv'))
lines<-c()
i<-1
for(mid in mdt$MG.RAST.ID){
  s<-tryCatch(paste(fread(paste0('tail -n 1 mgm',mid,'.fseed')),collapse = ' '))
  if (inherits(s, "try-error")) {
    lines[i]<-paste0('sbatch getMGRAST.fseed.sh ',webkey,' mgm',mid)
    i<-i+1
  }else if(!grepl('Download complete',s)){
    lines[i]<-paste0('sbatch getMGRAST.fseed.sh ',webkey,' mgm',mid)
    i<-i+1
  }
  s<-tryCatch(paste(fread(paste0('tail -n 1 mgm',mid,'.seed')),collapse = ' '))
  if (inherits(s, "try-error")) {
    lines[i]<-paste0('sbatch getMGRAST.seed.sh ',webkey,' mgm',mid)
    i<-i+1
  }else if(!grepl('Download complete',s)){
    lines[i]<-paste0('sbatch getMGRAST.seed.sh ',webkey,' mgm',mid)
    i<-i+1
  }
  s<-tryCatch(paste(fread(paste0('tail -n 1 mgm',mid,'.ko')),collapse = ' '))
  if (inherits(s, "try-error")) {
    lines[i]<-paste0('sbatch getMGRAST.ko.sh ',webkey,' mgm',mid)
    i<-i+1
  }else if(!grepl('Download complete',s)){
    lines[i]<-paste0('sbatch getMGRAST.ko.sh ',webkey,' mgm',mid)
    i<-i+1
  }
}
if(length(lines)>0){
  fname<-paste0('resubmit',format(Sys.time(), "%Y.%m.%d.%H.%M.%S"),'.sh')
  writeLines(c('#!/bin/bash',lines),fname)
  system(paste0('chmod a+x ',fname))
  cat('The download is incomplete. Run ',fname,'\n')
}else{
  cat('The download is complete. Preparing Rdata file\n')

}
