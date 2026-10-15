load('data/data.rda')
d<-mrna[mrna$patient_id %in% c('TCGA-XR-A8TC','TCGA-BC-A112'),]
write.csv(d,'new/results/E0/public_LIHC_two_case_input.csv',row.names=FALSE)
