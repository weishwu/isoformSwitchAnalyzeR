# isoformSwitchAnalyzeR


### prepare GTF
    paste <(awk '{if ($3!="gene") print $0}' /nfs/med-bfx-common/ENSEMBL_references/Mus_musculus/GRCm38/Mus_musculus.GRCm38.98.gtf|grep -v "^#"|cut -f1-8) <(awk '{if ($3!="gene") print $0}' /nfs/med-bfx-common/ENSEMBL_references/Mus_musculus/GRCm38/Mus_musculus.GRCm38.98.gtf|grep -v "^#"|cut -f9|awk '{print $1,$2,$5,$6}') >Mus_musculus.GRCm38.98.mod.gtf  

### in R
    library('IsoformSwitchAnalyzeR')

    # Load RSEM TPM. By default, calculateCountsFromAbundance=T, addIsofomIdAsColumn=T, interLibNormTxPM=T, normalizationMethod=TMM
    data=importIsoformExpression(parentDir = './')  
    myDesign=as.data.frame(read.table('smps.txt',header=T))  
    colnames(myDesign)=c('sampleID','condition')  
    #gtf='/nfs/med-bfx-common/ENSEMBL_references/Mus_musculus/GRCm38/Mus_musculus.GRCm38.98.gtf'  
    gtf='Mus_musculus.GRCm38.98.mod.gtf'  
    transcript.fa='/nfs/med-bfx-activeprojects/Soleimanpour_RS1_weishwu_damki_NovaA-225/analysis_20200318/alignment_results/04-rsem_star_genome_generate/GRCm38.transcripts.fa'  
    aSwitchList=importRdata(  
    isoformCountMatrix   = data$counts,  
    isoformRepExpression = data$abundance,  
    designMatrix         = myDesign,  
    isoformExonAnnoation = gtf,  
    isoformNtFasta       = transcript.fa,  
    showProgress = FALSE  
    )  

### Filtering
### default: 
         # switchAnalyzeRlist,
         # geneExpressionCutoff = 1,
         # isoformExpressionCutoff = 0,
         # IFcutoff=0.01,
         # acceptedGeneBiotype = NULL,
         # acceptedIsoformClassCode = NULL,
         # removeSingleIsoformGenes = TRUE,
         # reduceToSwitchingGenes=FALSE,
         # reduceFurtherToGenesWithConsequencePotential = FALSE,
         # onlySigIsoforms = FALSE,
         # keepIsoformInAllConditions=FALSE,
         # alpha=0.05,
         # dIFcutoff = 0.1,
         # quiet=FALSE

    asl_filtered=preFilter(aSwitchList)

### Testing for Isoform Switches via DEXSeq
### default:
         # switchAnalyzeRlist,
         # alpha = 0.05,
         # dIFcutoff = 0.1,
         # correctForConfoundingFactors=TRUE,
         # overwriteIFvalues=TRUE,
         # reduceToSwitchingGenes = TRUE,
         # reduceFurtherToGenesWithConsequencePotential = TRUE,
         # onlySigIsoforms = FALSE,
         # keepIsoformInAllConditions = TRUE,
         # showProgress = TRUE,
         # quiet = FALSE

    asl_analyzed = isoformSwitchTestDEXSeq(switchAnalyzeRlist = asl_filtered)
    extractSwitchSummary(asl_analyzed)

### Extracting Nucleotide and Amino Acid Sequences
### default:
         # switchAnalyzeRlist,
         # genomeObject  = NULL,
         # onlySwitchingGenes = TRUE,
         # alpha = 0.05,
         # dIFcutoff = 0.1,
         # extractNTseq = TRUE,
         # extractAAseq = TRUE,
         # removeShortAAseq = TRUE,
         # removeLongAAseq  = FALSE,
         # alsoSplitFastaFile = FALSE,
         # removeORFwithStop=TRUE,

    asl_extracted = extractSequence(asl_analyzed, pathToOutput = '/nfs/med-bfx-activeprojects/Soleimanpour_RS1_weishwu_damki_NovaA-225/isoformswitch/rsem')
    save(asl_extracted,file='asl_extracted.RData')

### Sequence analysis

#### 1. Coding potential assessment

##### 1.1 CPAT
* Use default parameters and the nucleotide fasta file (_nt.fasta). If the webserver was used, download the tab-delimited result file (available at the bottom of the result page). If a stand-alone version was used, just supply the path to the result file.
* command-line manual: http://rna-cpat.sourceforge.net/#make-hexamer-tab-py
* Download pre-built hexamer files from: https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/
```
cpat.py -g isoform_nucleotide.fasta -d /usr/share/cpat_data/Mouse_logitModel.RData -x /usr/share/cpat_data/Mouse_Hexamer.tsv -o cpat_output.txt
```

##### 1.2 CPC2
* Use default parameters and if required select the most similar species. If the webserver (batch submission) was used, download the tab-delimited result file (via the “Download the result” button). If a stand-alone version was just just supply the path to the result file.

#### 2. Pfam: Prediction of protein domains
* Use default parameters and the amino acid fasta file (_AA.fasta). If the webserver is used you need to copy/paste the result part of the mail you receive into an empty plain text document (notepad, sublimetext, TextEdit or similar (not Word)) and save that to a plain text (txt) file. The path to that file should be supplied. If a stand-alone version was used, just supply the path to the result file. A more detailed walkthrough is found under details in the documentation of the analyzePFAM() function (?analyzePFAM).
* command-line manual: http://avrilomics.blogspot.com/2015/08/pfamscanpl.html
* download database files from: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release
* Results are different from the web-server results.
```
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
# Input is the AA sequence of one transcript only
python /usr/share/iupred2a/iupred2a.py -a -d /usr/share/iupred2a ENSMUST00000092794_aminoacid.fasta long >iupred_output_ENSMUST00000092794.txt
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

