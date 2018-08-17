#!/bin/bash

echo 'MGRAST download preparation script.\n Usage:\n prepareProject.sh webkey projectID\n'
echo $1
echo $2
R --slave  <<EOR
webkey<-'$1'
prjTMP<-'$2'
cat(webkey,prjTMP)
rmarkdown::render('prepareProject.Rmd','pdf_document')
system(paste0('mkdir project.',proj.ID))
system(paste0('cp submit.sh getMGRAST*.sh prepareProject.pdf ',proj.ID,'.meta.csv project.',proj.ID,'/'))
system(paste0('cp checkDownload.R createPathView.R funtree.rda ../R/keggmappings.Rdata project.',proj.ID,'/'))
system(paste0('cp ',proj.ID,'.biom 	',proj.ID,'/'))
system(paste0('tar cvjf project.',proj.ID,'.tbz2 project.',proj.ID,'/'))
q(save='no')
EOR
