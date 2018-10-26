###############################################
##    SecretSanta: flexible pipelines for    ##
##    functional secretome prediction        ##
###############################################

#LOAD LIBRARIES AND DATA
library(SecretSanta)

#Verificar Dependencias
c<- manage_paths(in_path = TRUE)

" As in the pipeline the functions are designed to work together by 
producing standardized output as an instance of CBSResult class, so:"

#1.Read FASTA files in an XStringSet object:

Load_fasta<-readAAStringSet("/home/paola/1.Masterproyect/Dataset1/LILI_transcriptome_V5_transdecoder.pep_mod_withoutStopCodonSymbol.fasta")

#2.Create a simple CBSResult object:
"in_fasta: MultiFASTA format in XStringSet object
 out_fasta: File with only positive candidates, i.e those that passed tool filters"

CBSResult_input_LILIGogleva <-CBSResult(in_fasta = Load_fasta,
                                        out_fasta = Load_fasta[1:20])
CBSResult_input_LILIGogleva

#Check that CBSResult instance is valid
validObject(CBSResult_input_LILIGogleva)

#Extracts the corresponding attribute
getInfasta<-getInfasta(CBSResult_input_LILIGogleva)
class(getInfasta)
getOutfasta<- getOutfasta(CBSResult_input_LILIGogleva)
getFastas<-getFastas(CBSResult_input_LILIGogleva)
class<-class(CBSResult_input_LILIGogleva)

################################################
# SIGNALP
################################################

#starter when initiating the secretome pipeline with the above function
#piper for downstream/intermediate steps, so that the function expects an output from another 

CBSResult_input_SignalP <-CBSResult(in_fasta = Load_fasta)


# All the SignalP versions on the same input for comparison
############################################################
# changing value of the version argument and the run_mode

#SignalP2
step1_sp2_Gogleva <- signalp(CBSResult_input_SignalP,
                     version = 2,
                     organism = 'euk',
                     run_mode = "starter",
                     legacy_method = 'hmm')
#Extracts the corresponding attribute for signalP
getInfasta_step1_sp2_Gogleva<-getInfasta(step1_sp2_Gogleva)
class(getInfasta_step1_sp2_Gogleva)
getOutfasta_step1_sp2_Gogleva<- getOutfasta(step1_sp2_Gogleva)
#All and mature_fasta ->Sequences with cleaved N-terminal signal peptides
getFastas_step1_sp2_Gogleva<-getFastas(step1_sp2_Gogleva)
class_step1_sp2_Gogleva<- class(step1_sp2_Gogleva)
#sp_tibble parsed SignalP tabular output for positive candidates
getSPtibble_step1_sp2_Gogleva<-getSPtibble(step1_sp2_Gogleva)

#SignalP3
step1_sp3_Gogleva <- signalp(CBSResult_input_SignalP, 
                     version = 3,
                     organism = 'euk',
                     run_mode = "starter",
                     legacy_method = 'hmm')
#Extracts the corresponding attribute for signalP
getInfasta_step1_sp3_Gogleva<-getInfasta(step1_sp3_Gogleva)
class(getInfastaSignalP)
getOutfasta_step1_sp3_Gogleva<- getOutfasta(step1_sp3_Gogleva)
#All and mature_fasta ->Sequences with cleaved N-terminal signal peptides
getFastas_step1_sp3_Gogleva<-getFastas(step1_sp3_Gogleva)
class_step1_sp3_Gogleva<- class(step1_sp3_Gogleva)
#sp_tibble parsed SignalP tabular output for positive candidates
getSPtibble_step1_sp3_Gogleva<-getSPtibble(step1_sp3_Gogleva)

#SignalP 4.1
#could be not sensitive enough to predict certain classes 
#of secreted oomycete and fungal effectors: sensitive=TRUE
step1_sp4_Gogleva <- signalp(CBSResult_input_SignalP, version = 4,
                            organism = 'euk',
                            run_mode = "starter",
                            sensitive = TRUE,
                            cores = 4)
  
#Extracts the corresponding attribute for signalP
getInfasta_step1_sp4_Gogleva<-getInfasta(step1_sp4_Gogleva)
class(getInfasta_step1_sp4_Gogleva)
getOutfasta_step1_sp4_Gogleva<- getOutfasta(step1_sp4_Gogleva)
#All and mature_fasta ->Sequences with cleaved N-terminal signal peptides
getFastas_step1_sp4_Gogleva<-getFastas(step1_sp4_Gogleva)
class_step1_sp4_Gogleva<- class(step1_sp4_Gogleva)
#sp_tibble parsed SignalP tabular output for positive candidates
getSPtibble_step1_sp4_Gogleva<-getSPtibble(step1_sp4_Gogleva)


#  signalp2 -> signalp3 -> signalp4
############################################################
#switch to the run_mode = 'piper' for the second and other downstream steps.
step2_sp3_Gogleva <- signalp(step1_sp2_Gogleva,
                     version = 3,
                     organism = 'euk',
                     run_mode = 'piper',
                     legacy_method = 'hmm')

step3_sp4_Gogleva.toPipe <- signalp(step2_sp3_Gogleva,
                     version = 4,
                     organism = 'euk',
                     run_mode = "piper",
                     cores = 4)

#Extracts the corresponding attribute for signalP
getInfasta_step3_sp4_Gogleva<-getInfasta(step3_sp4_Gogleva.toPipe)
class(getInfasta_step3_sp4_Gogleva)
getOutfasta_step3_sp4_Gogleva.toPipe<- getOutfasta(step3_sp4_Gogleva.toPipe)
#All and mature_fasta ->Sequences with cleaved N-terminal signal peptides
getFastas_step3_sp4_Gogleva.toPipe<-getFastas(step3_sp4_Gogleva.toPipe)
class_step3_sp4_Gogleva.toPipe<- class(step3_sp4_Gogleva.toPipe)
#sp_tibble parsed SignalP tabular output for positive candidates to pipe
getSPtibble_step3_sp4_Gogleva.toPipe<-getSPtibble(step3_sp4_Gogleva.toPipe)

################################################
# TMHMM
################################################
"Identifies integral membrane proteins based on HMMs.
It is important to exclude proteins with transmembrane 
domains located after signal peptide, as they will
be retained in the membrane "
#with SignalP4 Started
tm__step1_sp4_Gogleva <- tmhmm(step1_sp4_Gogleva, TM = 0)
#with SignalP4 to Pipe
tm__step3_sp4_Gogleva <- tmhmm(step1_sp4_Gogleva, TM = 0)

################################################
# Topcons 
################################################
" "
# extract out_fasta slot with positive candidates and use these sequences to run TOPCONS prediction on a web-server:
writeXStringSet(getOutfasta_step1_sp2_Gogleva, 'sp_out.fasta')
# provide file path/file name with archived results:
p_dir <- system.file("extdata", "heree", package = "SecretSanta")
# integrate TOPCONS predictions:
tpc <- topcons(input_obj = step1_sp2_Gogleva,
               parse_dir = p_dir,
               topcons_mode = "WEB-server",
               TM = 0,
               SP = TRUE)

class(step1_sp2_Gogleva)
?topcons

################################################
# Targetp
################################################
"predicts the subcellular localization of secreted 
eukaryotic proteins based on the presence of signal
peptide (SP), chloroplast transit peptides(cTP) or 
mitochondrial targeting peptides (mTP) in the N-terminus "

# to make sure that a set of selected proteins 
#is not targeted to plastids or mitochondria. i.e > most likely to be secreted. 

targetp_Gogleva_all<- targetp(CBSResult_input_SignalP,  #Number of submitted sequences... 36839
                              network = 'N',
                              run_mode = 'starter')
targetp_Gogleva_all

targetp_Gogleva_pipe<- targetp(step3_sp4_Gogleva.toPipe,  #Number of submitted sequences... 2459
                              network = 'N',
                              run_mode = 'piper')
targetp_Gogleva_pipe

################################################
# WoLF PSORT 
################################################
"predicts the subcellular localization of secreted 
eukaryotic proteins based on the presence of signal
peptide (SP), chloroplast transit peptides(cTP) or 
mitochondrial targeting peptides (mTP) in the N-terminus "

wlf_Gogleva_all <- wolfpsort(CBSResult_input_SignalP, 
                             organism = 'fungi', 
                             run_mode = 'starter')
wlf_Gogleva_all
getWOLFtibble(wlf_Gogleva_all)

wlf_Gogleva_pipe <- wolfpsort(step3_sp4_Gogleva.toPipe, 
                             organism = 'fungi', 
                             run_mode = 'starter')
wlf_Gogleva_pipe
getWOLFtibble(wlf_Gogleva_pipe)

################################################
# Check_khdel 
################################################

# As the protein can also have can have an 
#ER-retention signal in the C-terminal domain that
#prevents protein from being secreted outside the cell.

#pattern = 'prosite'-->ER “[KRHQSA] [DENQ] EL $>” 
#pattern = 'elm'    -->ER “[KRHQSAP] [DENQT] EL $ 
#pattern = 'strict'--> strict matching for N-terminal KDEL or HDEL motifs

ER_result_Gogleva <- check_khdel(step3_sp4_Gogleva.toPipe, pattern = 'prosite')
ER_result_Gogleva
getOutfasta(ER_result_Gogleva)

################################################
# M_slicer
################################################
"takes the input amino acid sequences and 
generates all possible subsequences starting with methionine based.
has 2 running modes:"
  #slice - to simply slice input fasta regardless of it’s origin

slice_Gogleva<- m_slicer(Load_fasta, # a set of amino acid sequences
                         run_mode = 'slice',
                         min_len = 100)  # minimal length of the outputed slices

  #rescue - having output from any other up-stream function it extracts proteins

rescued_Gogleva <- m_slicer(step3_sp4_Gogleva.toPipe,  # signalp2 output
                    run_mode = 'rescue',
                    min_len = 100)




