import gzip,csv,math,statistics
from pathlib import Path
base=Path('new/results/E0')
targets=['TCGA-XR-A8TC','TCGA-BC-A112']
rda={}
with open(base/'public_LIHC_two_case_input.csv',newline='',encoding='utf-8-sig') as f:
 for row in csv.DictReader(f):
  rda[row['patient_id']]=row
values={s:{} for s in targets}
with gzip.open(base/'public_TCGA_LIHC_HiSeqV2.gz','rt',newline='',encoding='utf-8') as f:
 r=csv.reader(f,delimiter='\t'); header=next(r)
 samplecols={s:header.index(s+'-01') for s in targets}
 for row in r:
  gene=row[0]
  for s in targets:
   try: values[s][gene]=float(row[samplecols[s]])
   except (ValueError,IndexError): pass
summary=[]; out=[]
for s in targets:
 genes=[g for g in values[s] if g in rda[s] and rda[s][g]!='']
 a=[float(rda[s][g]) for g in genes]; x=[values[s][g] for g in genes]
 for g,y,z in zip(genes,a,x): out.append({'case':s,'gene':g,'input':y,'xena_hiseqv2':z,'diff_raw_minus_input':z-y,'diff_log2_1p_minus_input':math.log2(z+1)-y if z>=0 else float('nan'),'diff_2pow_input_minus1_minus_xena':2**y-1-z})
 for trans,xx in [('raw',x),('log2_1p',[math.log2(z+1) if z>=0 else float('nan') for z in x]),('exp2_input_minus1',[2**y-1 for y in a])]:
  pairs=[(float(y),float(z)) for y,z in zip(a,xx) if math.isfinite(y) and math.isfinite(z)]
  aa=[q[0] for q in pairs]; bb=[q[1] for q in pairs]; n=len(aa)
  mean_a=sum(aa)/n;mean_b=sum(bb)/n
  cor=sum((u-mean_a)*(v-mean_b) for u,v in pairs)/math.sqrt(sum((u-mean_a)**2 for u in aa)*sum((v-mean_b)**2 for v in bb)) if n>1 else float('nan')
  rmse=math.sqrt(sum((u-v)**2 for u,v in pairs)/n)
  exact=sum(abs(u-v)<1e-6 for u,v in pairs)
  nonzero=[(u,v) for u,v in pairs if u!=0 and v!=0]
  nzexact=sum(abs(u-v)<1e-6 for u,v in nonzero)
  summary.append({'case':s,'transform':trans,'n':n,'pearson':cor,'rmse':rmse,'within_1e-6':exact,'both_nonzero_n':len(nonzero),'both_nonzero_within_1e-6':nzexact})
with open(base/'public_LIHC_Xena_HiSeqV2_value_match.csv','w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
with open(base/'public_LIHC_Xena_HiSeqV2_match_summary.csv','w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=summary[0].keys());w.writeheader();w.writerows(summary)
print('samplecols',samplecols,'genes',len(genes));print(summary)
