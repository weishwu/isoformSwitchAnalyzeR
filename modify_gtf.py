import pandas as pd
import csv
import sys

gtf=sys.argv[1]
out=sys.argv[2]

gtf=pd.read_csv(gtf,comment='#',sep='\t',dtype=str,header=None)
gtf=gtf.loc[gtf.iloc[:,2]!="gene"].copy()
gtf['gene_id']=gtf.iloc[:,8].apply(lambda x: [y for y in x.split('; ') if y.startswith('gene_id')][0])
gtf['transcript_id']=gtf.iloc[:,8].apply(lambda x: [y for y in x.split('; ') if y.startswith('transcript_id')][0])
gtf['gene_transcript_id']=gtf['gene_id']+'; '+gtf['transcript_id']
gtf.drop([8,'gene_id','transcript_id'],axis=1).to_csv(out,sep='\t',index=False,header=None,quoting=csv.QUOTE_NONE)

