#library(parallel)
#cl <- makeCluster(getOption("cl.cores", 4))
system.time(epistasisFn(Data,Y.FA))
#stopCluster(cl)
