###############################################
##    SecretSanta: flexible pipelines for    ##
##    functional secretome prediction        ## 
##                                           ##
##          RESCUED SECRETOME                ##
###############################################

library(SecretSanta)
Load_fasta<-readAAStringSet("/home/paola/1.Masterproyect/04_MAKER/P8084.finalAssembly.all.maker.augustus_masked.proteins.fasta")

#To Cluster
#load_fasta<-readAAStringSet("/hpcfs/home/cp.rojas/MasterProyect_effectors/SecretSanta_ResultGogleva/Archivos_Gogleva
#/LILI_transcriptome_V5_transdecoder.pep_mod_withoutStopCodonSymbol.fasta")

CBSResult_input_SignalP <-CBSResult(in_fasta = Load_fasta)


pipeline_Rescued <- m_slicer(Load_fasta,run_mode = 'slice',min_len = 100) %>%
  signalp(CBSResult_input_SignalP, version = 2, organism = 'euk', run_mode = 'starter',legacy_method = 'hmm' ) %>% 
  signalp(version = 3, organism = 'euk', run_mode = 'piper') %>% tmhmm(TM = 0) %>%
  targetp(network = 'N', run_mode = 'piper') %>%
  check_khdel(pattern = 'prosite')
writeXStringSet(getOutfasta(pipeline_all), 'secretome_RESCUED_out.fasta') 
  
  
  
