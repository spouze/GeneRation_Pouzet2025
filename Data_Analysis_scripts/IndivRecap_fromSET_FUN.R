IndivRecap = function(storeLasts_set = dream_SET, line = 1, clean=T, cleanthresh=0.01){
  
  selected_line = storeLasts_set[line,]
  new_ind = list()
  new_ind$id = selected_line$id
  new_ind$topo = selected_line$topo
  
  
  # Load specific parameters
  L     = 10
  theta = c() ; for (i in 1:L){theta = c(theta, selected_line[[paste("FIT_OPT",  i, sep="")]])}
  s1    = c() ; for (i in 1:L){s1    = c(s1   , selected_line[[paste("FIT_STR",  i, sep="")]])}
  s2    = rep(46000, L)
  
  # Record parameters in ind.
  new_ind$L     = L
  new_ind$theta = theta
  new_ind$s1    = s1
  new_ind$s2    = s2
  
  #get matrix
  newW=c()
  for (i in 1:(L*L)){newW = c(newW, selected_line[[paste("MeanAll", i, sep="")]])} # get interaction values
  for (i in 1:L)    {newW = c(newW, selected_line[[paste("MTrans",  i, sep="")]])} # add coding
  newW=matrix(newW, nrow=L)
  new_ind$mom = newW
  new_ind$dad = newW
  new_ind$ind = newW
  
  # Compute phenotype
  dev = model.M2(newW, steps=20, measure=2, tfreg="notunique!", full=T)
  new_ind$mean = dev$mean
  new_ind$var  = dev$var
  
  new_ind$fitness = addfitness(new_ind, L=L, theta=theta, s1=s1, s2=s2)$fitness
  
  if (clean) {
    #get clean matrix
    cleaned_newW = WClean(newW, threshold = cleanthresh)
    new_ind$indclean = cleaned_newW
    #get corresponding dev
    dev = model.M2(cleaned_newW, steps=20, measure=2, tfreg="notunique!", full=T)
    new_ind$meanclean = dev$mean
    new_ind$varclean  = dev$var
    #to get the fitness, need to create an ind intermediate
    intermediate = list()
    intermediate$mean = dev$mean
    intermediate$var = dev$var
    #Now get fitness thanks to this guy:
    new_ind$fitnessclean = addfitness(intermediate, L=L, theta=theta, s1=s1, s2=s2)$fitness
    }
  
  return(new_ind)
  
}
