###############################################
##    SecretSanta: flexible pipelines for    ##
##    functional secretome prediction        ## 
##                                           ##
##          RESCUED SECRETOME                ##
###############################################

library(SecretSanta)
Load_fasta<-readAAStringSet("/home/paola/1.Masterproyect/Dataset1/LILI_transcriptome_V5_transdecoder.pep_mod_withoutStopCodonSymbol.fasta")

#To Cluster
#load_fasta<-readAAStringSet("/hpcfs/home/cp.rojas/MasterProyect_effectors/SecretSanta_ResultGogleva/Archivos_Gogleva
#/LILI_transcriptome_V5_transdecoder.pep_mod_withoutStopCodonSymbol.fasta")

CBSResult_input_SignalP <-CBSResult(in_fasta = Load_fasta)


pipeline_Rescued <- m_slicer(Load_fasta, # a set of amino acid sequences
                                             run_mode = 'slice',
                                             min_len = 100)
  
  
  
  signalp(CBSResult_input_SignalP, version = 4, organism = 'euk', run_mode = 'starter') %>% 
  tmhmm(TM = 0) %>%
  targetp(network = 'N', run_mode = 'piper') %>%
  check_khdel(pattern = 'prosite')
writeXStringSet(getOutfasta(pipeline_all), 'secretome_out.fasta')