# isoformSwitchAnalyzeR


### prepare GTF
awk '{if ($3!="gene") print $0}' /nfs/med-bfx-common/ENSEMBL_references/Mus_musculus/GRCm38/Mus_musculus.GRCm38.98.gtf|grep -v "^#"|cut -f1-8 >aa
awk '{if ($3!="gene") print $0}' /nfs/med-bfx-common/ENSEMBL_references/Mus_musculus/GRCm38/Mus_musculus.GRCm38.98.gtf|grep -v "^#"|cut -f9|awk '{print $1,$2,$5,$6}' >bb
paste aa bb >Mus_musculus.GRCm38.98.mod.gtf

### in R
library('IsoformSwitchAnalyzeR')

### Load RSEM TPM. By default, calculateCountsFromAbundance=T, addIsofomIdAsColumn=T, interLibNormTxPM=T, normalizationMethod=TMM
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

# Filtering
# default: 
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

# Testing for Isoform Switches via DEXSeq
# default:
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

# Extracting Nucleotide and Amino Acid Sequences
# default:
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

# Sequence analysis

# CPAT or CPC2: Coding potential assessment
# CPAT : Use default parameters and the nucleotide fasta file (_nt.fasta). If the webserver was used, download the tab-delimited result file (available at the bottom of the result page). If a stand-alone version was used, just supply the path to the result file.
# CPC2 : Use default parameters and if required select the most similar species. If the webserver (batch submission) was used, download the tab-delimited result file (via the “Download the result” button). If a stand-alone version was just just supply the path to the result file.

# Pfam: Prediction of protein domains
# Use default parameters and the amino acid fasta file (_AA.fasta). If the webserver is used you need to copy/paste the result part of the mail you receive into an empty plain text document (notepad, sublimetext, TextEdit or similar (not Word)) and save that to a plain text (txt) file. The path to that file should be supplied. If a stand-alone version was used, just supply the path to the result file. A more detailed walkthrough is found under details in the documentation of the analyzePFAM() function (?analyzePFAM).

# SignalP: Prediction of Signal Peptides — a short N-terminal sequence of a peptide indicating where a protein should be membrane bound or secreted
# Use the amino acid fasta file (_AA.fasta). If using the webserver SignalP should be run with the parameter “Short output (no figures)” under “Output format” and one should select the appropriate “Organism group”. When using a stand-alone version SignalP should be run with the ‘-f summary’ option. If using the webserver the results can be downloaded using the “Downloads” button in the top-right corner where the user should select “Prediction summary” and supply the path to the resulting file to the “pathToSignalPresultFile” argument. If a stand-alone version was just supply the path to the summary result file.

# IUPred2A or NetSurfP-2: Prediction of Intrinsically Unstructured Proteins
# IUPred2A: 1) Go to the webserver 2) Upload the amino acoid file (_AA) created with extractSequence() function. 3) Add your email (you will recieve a notification when the job is done). 4) In the email use the link indicated by “The text file can be found here:” and save the result (right click on a blank space and use “save as” or use the keybord shortcut “ctrl/cmd + s”). 5) Supply a string indicating the path to the downloaded file to the “pathToIUPred2AresultFile” argument. If multiple files are creted (multiple web-server runs) just supply the path to all of them as a string.
# NetSurfP-2: 1) Go to webserver. 2) Upload the amino acid file (_AA) created with extractSequence() function. 3) Submit your job. 4) Wait till job is finished (if you submit your email you will receive a notification). 5) In the top-right corner of the result site use the “Export All” button to download the results as a CNV file. 6) Supply a string indicating the path to the downloaded csv file directly to the “pathToNetSurfP2resultFile” argument. If you run NetSurfP-2 locally just use the “–csv” argument and provide the resulting csv file to the pathToNetSurfP2resultFile argument.

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


# Predicting Alternative Splicing
asl_data_altspl = analyzeAlternativeSplicing(
    switchAnalyzeRlist = asl_data,
    quiet=TRUE
)

### overview of number of intron retentions (IR)
table( asl_data_altspl$AlternativeSplicingAnalysis$IR )

