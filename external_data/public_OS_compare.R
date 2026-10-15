suppressPackageStartupMessages(library(jsonlite))
read_json_n <- function(path) {
  z <- readBin(path,'raw',n=file.info(path)$size); s <- rawToChar(z); s <- sub('^\\ufeff','',s); fromJSON(s,simplifyVector=FALSE)
}
collect_cases <- function(path) read_json_n(path)$data$hits
cases <- collect_cases('new/results/E0/public_primary_gdc_clinical.json')
primary <- read.csv('new/results/E0/public_primary_clinical.csv',check.names=FALSE,stringsAsFactors=FALSE)
getv <- function(x) {if(is.null(x)||!length(x)) return(numeric()); suppressWarnings(as.numeric(unlist(x,use.names=FALSE)))}
rows<-lapply(cases,function(ca){
 pid<-ca$submitter_id; p<-primary[primary$patient_id==pid,,drop=FALSE]; if(nrow(p)!=1)return(NULL)
 d<-getv(ca$demographic$days_to_death); diag<-unlist(lapply(ca$diagnoses,function(x)getv(x$days_to_last_follow_up)),use.names=FALSE); fu<-unlist(lapply(ca$follow_ups,function(x)getv(x$days_to_follow_up)),use.names=FALSE)
 stat<-tolower(as.character(ca$demographic$vital_status)); if(!length(stat)) stat<-NA_character_; d<-d[is.finite(d)]; diag<-diag[is.finite(diag)]; fu<-fu[is.finite(fu)]
 expected<-if(stat=='dead') if(length(d)) d[[1]] else NA_real_ else if(length(c(diag,fu))) max(c(diag,fu),na.rm=TRUE) else NA_real_
 data.frame(patient_id=pid,project=ca$project$project_id,status_local=as.numeric(p$OS_status),status_GDC=stat,OS_time_local=as.numeric(p$OS_time),days_to_death=if(length(d)) d[[1]] else NA_real_,max_days_last_followup=max(c(diag,fu),na.rm=TRUE),expected_GDC_days=expected,diff_local_minus_expected=as.numeric(p$OS_time)-expected,exact_local_days=as.numeric(p$OS_time)==expected,stringsAsFactors=FALSE)
})
primary_cmp<-do.call(rbind,rows); primary_cmp$max_days_last_followup[!is.finite(primary_cmp$max_days_last_followup)]<-NA_real_
write.csv(primary_cmp,'new/results/E0/public_primary_OS_match.csv',row.names=FALSE,na='')
cat('PRIMARY n=',nrow(primary_cmp),'status match=',sum((primary_cmp$status_local==1 & primary_cmp$status_GDC=='dead')|(primary_cmp$status_local==0 & primary_cmp$status_GDC=='alive')),'/',nrow(primary_cmp),'exact OS time=',sum(primary_cmp$exact_local_days,na.rm=TRUE),'/',nrow(primary_cmp),'\n')
print(aggregate(cbind(exact_local_days,diff_local_minus_expected)~status_GDC,primary_cmp,function(x)c(n=sum(!is.na(x)),mean=mean(x,na.rm=TRUE),median=median(x,na.rm=TRUE),exact=sum(x==0,na.rm=TRUE))))
# External survival.tsv candidates are evaluated only against public GDC days fields; do not attribute origin from filename.
ext_cases<-collect_cases('new/results/E0/public_external_gdc_clinical.json')
clin<-read.delim('data/livergabage_cancer_data/clinical.tsv',check.names=FALSE,stringsAsFactors=FALSE)
surv<-read.delim('data/livergabage_cancer_data/survival.tsv',check.names=FALSE,stringsAsFactors=FALSE)
clin<-clin[clin$project_id.project=='TCGA-CHOL',]
surv$patient<-substr(surv$sample,1,12)
rows2<-lapply(ext_cases,function(ca){
 pid<-ca$submitter_id; p<-unique(surv[surv$patient==pid,c('OS.time','OS')]); d<-getv(ca$demographic$days_to_death); diag<-unlist(lapply(ca$diagnoses,function(x)getv(x$days_to_last_follow_up)),use.names=FALSE); fu<-unlist(lapply(ca$follow_ups,function(x)getv(x$days_to_follow_up)),use.names=FALSE); stat<-tolower(as.character(ca$demographic$vital_status)); if(!length(stat)) stat<-NA_character_; d<-d[is.finite(d)]; diag<-diag[is.finite(diag)]; fu<-fu[is.finite(fu)]; expected<-if(!is.na(stat) && stat=='dead') if(length(d))d[[1]] else NA_real_ else if(length(c(diag,fu)))max(c(diag,fu),na.rm=TRUE)else NA_real_; if(!nrow(p)) return(NULL)
 data.frame(patient_id=pid,status_GDC=stat,GDC_expected_days=expected,survival_rows=nrow(p),OS_time_values=paste(unique(p$OS.time),collapse='|'),OS_values=paste(unique(p$OS),collapse='|'),diff_days=if(length(unique(p$OS.time))==1)as.numeric(p$OS.time[1])-expected else NA_real_,diff_months=if(length(unique(p$OS.time))==1)as.numeric(p$OS.time[1])-expected/30.4375 else NA_real_,stringsAsFactors=FALSE)
})
extcmp<-do.call(rbind,rows2); write.csv(extcmp,'new/results/E0/public_external_OS_match.csv',row.names=FALSE,na=''); print(summary(extcmp$diff_days)); print(summary(extcmp$diff_months)); cat('external distinct patients=',nrow(extcmp),'survival time rows=',sum(extcmp$survival_rows),'\n')

