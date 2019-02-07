###############################################      
##    SecretSanta: flexible pipelines for    ##       
##    functional secretome prediction        ##
###############################################


                                                                ##        ##
                                                                ##    functional secretome prediction        ##
                                                                ###############################################

#LOAD LIBRARIES AND DATA
library(SecretSanta)

################################################
# SIGNALP
################################################

#starter when initiating the secretome pipeline with the above function
#piper for downstream/intermediate steps, so that the function expects an output from another 

Load_fasta<-readAAStringSet("file= data.fasta")
CBSResult_input <-CBSResult(in_fasta = Load_fasta)

#SignalP
step1_sp2 <- signalp(CBSResult_input,
                     version = 2,
                     organism = 'euk',
                     run_mode = "starter",
                     legacy_method = 'hmm')
step2_sp3 <- signalp(step1_sp2,
                     version = 3,
                     organism = 'euk',
                     run_mode = 'starter',
                     legacy_method = 'hmm')
step3_sp4<- signalp(CBSResult_input,
                    version = 4,
                    organism = 'euk',
                    run_mode = "starter", 
                    cores = 4)

################################################
# TMHMM
################################################
"Identifies integral membrane proteins based on HMMs.
It is important to exclude proteins with transmembrane 
domains located after signal peptide, as they will
be retained in the membrane "
#with SignalP4 to pipe
tm_step3_sp4<- tmhmm(step2_sp3, TM = 0)

################################################
# Targetp
################################################
"predicts the subcellular localization of secreted 
eukaryotic proteins based on the presence of signal
peptide (SP), chloroplast transit peptides(cTP) or 
mitochondrial targeting peptides (mTP) in the N-terminus "

# to make sure that a set of selected proteins 
#is not targeted to plastids or mitochondria. i.e > most likely to be secreted. 

targetp_pipe<- targetp(tm_step3_sp4,  #Number of submitted sequences... 2459
                       network = 'N',
                       run_mode = 'piper')

################################################
# Check_khdel 
################################################

# As the protein can also have can have an 
#ER-retention signal in the C-terminal domain that
#prevents protein from being secreted outside the cell.

#pattern = 'prosite'-->ER “[KRHQSA] [DENQ] EL $>” 
#pattern = 'elm'    -->ER “[KRHQSAP] [DENQT] EL $ 
#pattern = 'strict'--> strict matching for N-terminal KDEL or HDEL motifs

ER_result<- check_khdel(targetp_pipe, pattern = 'prosite')

setwd("/home/Effectors/Data/Proteomes_Phytophthora/P_palmivora_(LILI_transcriptome_V5)/initial_secretome/")
writeXStringSet(getOutfasta(ER_result), 'initialSecretome_out_Ppalmivora.fasta')
writeXStringSet(getOutfasta(step1_sp2), 'initialSecretome_out_step1_sp2_Ppalmivora.fasta')
writeXStringSet(getOutfasta(step2_sp3), 'initialSecretome_out_step2_sp3_Ppalmivora.fasta')
writeXStringSet(getOutfasta(step3_sp4), 'initialSecretome_out_step3_sp4_Ppalmivora.fasta')
writeXStringSet(getOutfasta(tm_step3_sp4), 'initialSecretome_out_tm_step3_sp4_Ppalmivora.fasta')
writeXStringSet(getOutfasta(targetp_pipe), 'initialSecretome_out_targetp_pipe_Ppalmivora.fasta')


                        #RESCUED SECRETOME
                        #(initial secretome - initial proteome)
                        #m-slicer with other dependences.

################################################
# M_slicer
################################################
"takes the input amino acid sequences and 
generates all possible subsequences starting with methionine based.
has 2 running modes:"

#rescue - having output from any other up-stream function it extracts proteins


rescued<- m_slicer(step1_sp2,  #(initial secretome - initial proteome)
                          run_mode = 'rescue',
                          min_len = 100)

rescued_input <- CBSResult(in_fasta = rescued)

#SignalP
step1_sp2_rescued <- signalp(rescued_input,
                     version = 2,
                     organism = 'euk',
                     run_mode = "starter",
                     legacy_method = 'hmm')
step2_sp3_rescued <- signalp(step1_sp2_rescued,
                     version = 3,
                     organism = 'euk',
                     run_mode = 'piper',
                     legacy_method = 'hmm')
step3_sp4_rescued<- signalp(step2_sp3_rescued,
                    version = 4,
                    organism = 'euk',
                    run_mode = 'piper',
                    legacy_method = 'hmm',
                     cores = 4)

################################################
# TMHMM
################################################
"Identifies integral membrane proteins based on HMMs.
It is important to exclude proteins with transmembrane 
domains located after signal peptide, as they will
be retained in the membrane "
#with SignalP3 to pipe
tm_step2_sp4_rescued<- tmhmm(step2_sp3_rescued, TM = 0)

################################################
# Targetp
################################################
"predicts the subcellular localization of secreted 
eukaryotic proteins based on the presence of signal
peptide (SP), chloroplast transit peptides(cTP) or 
mitochondrial targeting peptides (mTP) in the N-terminus "

# to make sure that a set of selected proteins 
#is not targeted to plastids or mitochondria. i.e > most likely to be secreted. 

targetp_pipe_rescued<- targetp(tm_step2_sp4_rescued,  #Number of submitted sequences... 2459
                       network = 'N',
                       run_mode = 'piper')

################################################
# Check_khdel 
################################################

# As the protein can also have can have an 
#ER-retention signal in the C-terminal domain that
#prevents protein from being secreted outside the cell.

#pattern = 'prosite'-->ER “[KRHQSA] [DENQ] EL $>” 
#pattern = 'elm'    -->ER “[KRHQSAP] [DENQT] EL $ 
#pattern = 'strict'--> strict matching for N-terminal KDEL or HDEL motifs

ER_result<- check_khdel(targetp_pipe_rescued, pattern = 'prosite')

setwd("file download result")
writeXStringSet(getOutfasta(ER_result), 'file name')
writeXStringSet(getOutfasta(step1_sp2_rescued), 'file name')
writeXStringSet(getOutfasta(step2_sp3_rescued), 'file name')
writeXStringSet(getOutfasta(step3_sp4_rescued), 'file name')
writeXStringSet(getOutfasta(tm_step2_sp4_rescued),'file name')
writeXStringSet(getOutfasta(targetp_pipe_rescued), 'file name')
getwd()
