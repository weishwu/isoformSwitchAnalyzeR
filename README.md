### 1. Prepare Rsem count data
Create a folder for each sample (e.g. Sample_2/). Put Rsem gene-level and isoform-level count data (e.g. Sample_2.genes.results and Sample_2.isoforms.results) in the sample folders.

### 2. Prepare a sample table. 
The first column is the sample ID and the rest is the conditions.

### 3. Prepare GTF 
    # remove "gene" lines and retain "gene_id" and "transcript_id" fields
    python modify_gtf.py /nfs/turbo/umms-brcfpipeline/references/ENSEMBL_genomes/Homo_sapiens/GRCh38/Homo_sapiens.GRCh38.98.gtf Homo_sapiens.GRCh38.98.mod.gtf`  

### 4. Load data into IsoformSwitchAnalyzeR
    library('IsoformSwitchAnalyzeR')

    # Load RSEM TPM. By default, calculateCountsFromAbundance=T, addIsofomIdAsColumn=T, interLibNormTxPM=T, normalizationMethod=TMM
    data=importIsoformExpression(parentDir = '/nfs/med-bfx-activeprojects/Love_Simulated_DTU/IsoformSwitchAnalyzeR')  
    myDesign=as.data.frame(read.table('smps.txt',header=T))  
    colnames(myDesign)=c('sampleID','condition')  
    gtf='Homo_sapiens.GRCh38.98.mod.gtf'  
    transcript.fa='/nfs/med-bfx-activeprojects/trsaari/sandbox/20200717_Love_sim_comps/analysis_test/alignment_results/04-rsem_star_genome_generate/GRCh38.transcripts.fa'  
    aSwitchList=importRdata(  
    isoformCountMatrix   = data$counts,  
    isoformRepExpression = data$abundance,  
    designMatrix         = myDesign,  
    isoformExonAnnoation = gtf,  
    isoformNtFasta       = transcript.fa,  
    showProgress = FALSE  
    )

### 5. High-level function covering filtering, switch test, and nt & aa sequence extraction
    # default parameters:
    #     isoformSwitchAnalysisPart1(
    #     switchAnalyzeRlist,
    #     alpha = 0.05,
    #     dIFcutoff = 0.1,
    #     switchTestMethod='DEXSeq',
    #     orfMethod = "longest",
    #     genomeObject = NULL,
    #     cds = NULL,
    #     pathToOutput = getwd(),
    #     outputSequences = TRUE,
    #     prepareForWebServers = FALSE,
    #     overwriteORF=FALSE,
    #     quiet=FALSE)
    
    asl_analyzed=isoformSwitchAnalysisPart1(aSwitchList)
    save(asl_analyzed,file='asl_analyzed.rda')
    
    # select significant DTU isoforms
    asl_analyzed_selected=extractTopSwitches(
         asl_analyzed,
         filterForConsequences=FALSE,
         extractGenes=FALSE,
         alpha=0.05,
         dIFcutoff = 0.1,
         n=Inf,
         inEachComparison=FALSE,
         sortByQvals=TRUE
     )
     write.table(asl_analyzed_selected, 'isoforms_DTU.txt', quote=F,sep='\t',row.names=F,col.names=T)


### Sequence analysis
* Docker image containing cpat, pfam, signalp and iupred2a: 
`docker://weishwu/isoform_sequence_analysis:08202020`

#### 1. Coding potential assessment

##### 1.1 CPAT
* Use default parameters and the nucleotide fasta file (_nt.fasta). If the webserver was used, download the tab-delimited result file (available at the bottom of the result page). If a stand-alone version was used, just supply the path to the result file.
* command-line manual: http://rna-cpat.sourceforge.net/#make-hexamer-tab-py
* Download pre-built hexamer files from: https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/
```
source activate cpat
# mouse
cpat.py -g isoform_nucleotide.fasta -d /usr/share/cpat_data/Mouse_logitModel.RData -x /usr/share/cpat_data/Mouse_Hexamer.tsv -o cpat_output.txt
# human
cpat.py -g isoform_nucleotide.fasta -d /usr/share/cpat_data/Human_logitModel.RData -x /usr/share/cpat_data/Human_Hexamer.tsv -o cpat_output.txt
```

##### 1.2 CPC2
* Use default parameters and if required select the most similar species. If the webserver (batch submission) was used, download the tab-delimited result file (via the “Download the result” button). If a stand-alone version was just just supply the path to the result file.

#### 2. Pfam: Prediction of protein domains
* Use default parameters and the amino acid fasta file (_AA.fasta). If the webserver is used you need to copy/paste the result part of the mail you receive into an empty plain text document (notepad, sublimetext, TextEdit or similar (not Word)) and save that to a plain text (txt) file. The path to that file should be supplied. If a stand-alone version was used, just supply the path to the result file. A more detailed walkthrough is found under details in the documentation of the analyzePFAM() function (?analyzePFAM).
* command-line manual: http://avrilomics.blogspot.com/2015/08/pfamscanpl.html
* download HMM files from: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release
* The web-server linked from the IsoformSwitchAnalyzeR manual is actuall hmmscan, whereas the manual states that "pfam ... can be run either locally (using the pfam_scan.pl script as described in the readme found here or via their [webserver](http://lilab.research.bcm.edu/cpat/)". The difference is explained here: https://gist.github.com/olgabot/f65365842e27d2487ad3. The actual pfamscan web-server should be: https://www.ebi.ac.uk/Tools/pfa/pfamscan/
* The results are a little different between command-line version and the web-server. The lines are slightly fewer from the command-line version and the E-values are slightly differet. The web-server uses `-e_seq 10 -e_dom 10` by default but they are not the default parameters in the command-line version (values not shown).
```
source activate pfam
pfam_scan.pl -fasta isoform_aminoacid.fasta -dir /usr/share/pfam_data -outfile pfam_output.txt
```

#### 3. SignalP: Prediction of Signal Peptides — a short N-terminal sequence of a peptide indicating where a protein should be membrane bound or secreted
* Use the amino acid fasta file (_AA.fasta). If using the webserver SignalP should be run with the parameter “Short output (no figures)” under “Output format” and one should select the appropriate “Organism group”. When using a stand-alone version SignalP should be run with the ‘-f summary’ option. If using the webserver the results can be downloaded using the “Downloads” button in the top-right corner where the user should select “Prediction summary” and supply the path to the resulting file to the “pathToSignalPresultFile” argument. If a stand-alone version was just supply the path to the summary result file.
* Download command-line version from: http://www.cbs.dtu.dk/services/SignalP/portable.php
```
signalp -fasta isoform_aminoacid.fasta -org euk -format short -prefix signalp_output
```

#### 4. Prediction of Intrinsically Unstructured Proteins
##### 4.1 IUPred2A 
* Web-server: 1) Go to the webserver 2) Upload the amino acoid file (_AA) created with extractSequence() function. 3) Add your email (you will recieve a notification when the job is done). 4) In the email use the link indicated by “The text file can be found here:” and save the result (right click on a blank space and use “save as” or use the keybord shortcut “ctrl/cmd + s”). 5) Supply a string indicating the path to the downloaded file to the “pathToIUPred2AresultFile” argument. If multiple files are creted (multiple web-server runs) just supply the path to all of them as a string.
* Download command-line version from: https://iupred2a.elte.hu/download_new
* Manual: https://iupred2a.elte.hu/help_new
```
# iupred2a.py takes only one fasta sequence, so the isoform fasta file has to be split to individual isoform files first:
python split_AA_fasta.py

cd isoformSwitchAnalyzeR_isoform_AA_split
count=0
Nproc=20
for f in *.fa
do
python /usr/share/iupred2a/iupred2a.py -a -d /usr/share/iupred2a/ ${f} long >iupred_output_${f}.txt &
let count+=1
[[ $((count%Nproc)) -eq 0 ]] && wait
done
```

##### 4.2 NetSurfP-2 
* Web-server: 1) Go to webserver. 2) Upload the amino acid file (_AA) created with extractSequence() function. 3) Submit your job. 4) Wait till job is finished (if you submit your email you will receive a notification). 5) In the top-right corner of the result site use the “Export All” button to download the results as a CNV file. 6) Supply a string indicating the path to the downloaded csv file directly to the “pathToNetSurfP2resultFile” argument. If you run NetSurfP-2 locally just use the “–csv” argument and provide the resulting csv file to the pathToNetSurfP2resultFile argument.



### Add CPAT analysis
    asl_data = analyzeCPAT(
    switchAnalyzeRlist   = asl_extracted,
    pathToCPATresultFile = "./asl_CPAT.txt",
    codingCutoff         = 0.725, # the coding potential cutoff we suggested for human
    removeNoncodinORFs   = TRUE   # because ORF was predicted de novo
    )

### Add PFAM analysis
    asl_data = analyzePFAM(
    switchAnalyzeRlist   = asl_data,
    pathToPFAMresultFile = "./",
    showProgress=FALSE
    )

### Add SignalP analysis
    asl_data = analyzeSignalP(
    switchAnalyzeRlist       = asl_data,
    pathToSignalPresultFile  = "./asl_SignalP.txt"
    )
#> Added signal peptide information to 17 (10.49%) transcripts

### Add IUPred2A analysis
    asl_data = analyzeIUPred2A(
    switchAnalyzeRlist        = asl_data,
    pathToIUPred2AresultFile = "./asl_IUPred2A.result",
    showProgress = FALSE
    )


### Predicting Alternative Splicing
    asl_data_altspl = analyzeAlternativeSplicing(
    switchAnalyzeRlist = asl_data,
    quiet=TRUE
    )

### overview of number of intron retentions (IR)
    table( asl_data_altspl$AlternativeSplicingAnalysis$IR )


### plot
    ggplot(data=asl_conseq$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
     geom_point(
       aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=1
    ) +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
    geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
    facet_wrap( ~ condition_2) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
    labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
    theme_bw()

    switchPlot(
    asl_conseq,
    gene='ENSMUSG00000030741',
    condition1 = 'Ctrl',
    condition2 = 'Mfn12.Double.Kr',
    localTheme = theme_bw(base_size = 13) # making text sightly larger for vignette
    

    extractSplicingEnrichment(
    asl_conseq,
    returnResult = FALSE # if TRUE returns a data.frame with the results
    )

    extractConsequenceEnrichment(
    asl_conseq,
    consequencesToAnalyze='all',
    analysisOppositeConsequence = TRUE,
    returnResult = FALSE # if TRUE returns a data.frame with the results
    )

    extractConsequenceSummary(
    asl_conseq,
    consequencesToAnalyze='all',
    plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
    asFractionTotal = FALSE      # enables analysis of fraction of significant features
    )

    extractSwitchOverlap(
    asl_conseq,
    filterForConsequences=TRUE
    )

