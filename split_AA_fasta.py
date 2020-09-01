import os
import Bio
from Bio import SeqIO

os.mkdir('isoformSwitchAnalyzeR_isoform_AA_split')
os.chdir('./isoformSwitchAnalyzeR_isoform_AA_split')

records=list(SeqIO.parse('../isoformSwitchAnalyzeR_isoform_AA.fasta','fasta'))

for i in records:
    out=open(i.id+'.fa','w')
    out.write('>'+i.id+'\n'+str(i.seq)+'\n')
    out.close()

