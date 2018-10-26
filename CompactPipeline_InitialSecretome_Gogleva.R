###############################################
##    SecretSanta: flexible pipelines for    ##
##    functional secretome prediction        ## 
##                                           ##
##          INITIAL SECRETOME                ##
###############################################


pipeline_all <- signalp(CBSResult_input_SignalP, version = 4, organism = 'euk', run_mode = 'starter') %>% 
  tmhmm(TM = 0) %>%
  targetp(network = 'N', run_mode = 'piper') %>%
  check_khdel(pattern = 'prosite')
writeXStringSet(getOutfasta(pipeline_all), 'secretome_out.fasta')