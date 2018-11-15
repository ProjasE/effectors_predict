###############################################
##    SecretSanta: flexible pipelines for    ##
##    functional secretome prediction        ## 
##                                           ##
##          INITIAL SECRETOME                ##
###############################################

library(SecretSanta)
Load_fasta<-readAAStringSet("/home/paola/1.Masterproyect/04_MAKER/P8084.finalAssembly.all.maker.augustus_masked.proteins.fasta")

#To Cluster
#load_fasta<-readAAStringSet("/hpcfs/home/cp.rojas/MasterProyect_effectors/SecretSanta_ResultGogleva/Archivos_Gogleva
#/LILI_transcriptome_V5_transdecoder.pep_mod_withoutStopCodonSymbol.fasta")

CBSResult_input_SignalP <-CBSResult(in_fasta = Load_fasta)


pipeline_all <- signalp(CBSResult_input_SignalP, version = 2, organism = 'euk', run_mode = 'starter',legacy_method = 'hmm' ) %>% 
  signalp(version = 3, organism = 'euk', run_mode = 'piper',legacy_method = 'hmm') %>% tmhmm(TM = 0) %>%
  targetp(network = 'N', run_mode = 'piper') %>%
  check_khdel(pattern = 'prosite')
writeXStringSet(getOutfasta(pipeline_all), 'secretome_out.fasta')

