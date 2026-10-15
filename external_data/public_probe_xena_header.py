import gzip,csv
p='new/results/E0/public_TCGA_LIHC_HiSeqV2.gz'
with gzip.open(p,'rt',newline='') as f:
 r=csv.reader(f,delimiter='\t'); h=next(r)
 print(len(h), [(s,h.index(s)) for s in h if s.startswith(('TCGA-XR-A8TC','TCGA-BC-A112'))])
 print('feature first',h[:6])
