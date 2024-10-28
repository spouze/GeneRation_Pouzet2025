## From The last line of the simu[[1]], reconstitute an individual

# WFinale

# L=length(grep(x = colnames(simu[[1]]), pattern = "MPhen.*"))
# simu[[1]][1,] #ligne 1
# gen = nrow(simu[[1]])
# selected_line = simu[[1]][gen,] #last line
# mean_ind = list()
# 
# lastW=c()
# for (i in 1:(L*L)){lastW = c(lastW, selected_line[[paste("MeanAll", i, sep="")]])} # get interaction values
# for (i in 1:L)    {lastW = c(lastW, selected_line[[paste("MTrans",  i, sep="")]])} # add coding
# lastW=matrix(lastW, nrow=L)
# mean_ind$mom = lastW
# mean_ind$dad = lastW
# mean_ind$ind = lastW
# 
# MPhens=c()
# for (i in 1:L) {MPhens = c(MPhens, selected_line[[paste("MPhen",  i, sep="")]])}
# mean_ind$mean = MPhens
# 
# VPhens=c()
# for (i in 1:L) {VPhens = c(VPhens, selected_line[[paste("VPhen",  i, sep="")]])}
# mean_ind$var = VPhens
# 
# mean_ind$fitness = selected_line[["MFit"]]
# 
# 
# ######
# allfits=sapply(simu[[3]], function(i) i$fitness) # extract all fitnesses
# best=which.max(allfits) #ind index w max fitness
# # best=685
# best_ind=simu[[3]][[best]]
# best_ind


#####################################
#########################################
#############################################
##
### TURNING THIS INTO A USABLE FUNCTION:
##
#############################################
#########################################
#####################################


IndivFinale = function(.simu=simu, gen = 1000){
  
  L=length(grep(x = colnames(simu[[1]]), pattern = "MPhen.*"))
  # simu[[1]][1,] #ligne 1
  # gen = nrow(simu[[1]])
  gen = gen
  selected_line = simu[[1]][gen,] #last line
  mean_ind = list()
  
  lastW=c()
  for (i in 1:(L*L)){lastW = c(lastW, selected_line[[paste("MeanAll", i, sep="")]])} # get interaction values
  for (i in 1:L)    {lastW = c(lastW, selected_line[[paste("MTrans",  i, sep="")]])} # add coding
  lastW=matrix(lastW, nrow=L)
  mean_ind$mom = lastW
  mean_ind$dad = lastW
  mean_ind$ind = lastW
  
  MPhens=c()
  for (i in 1:L) {MPhens = c(MPhens, selected_line[[paste("MPhen",  i, sep="")]])}
  mean_ind$mean = MPhens
  
  VPhens=c()
  for (i in 1:L) {VPhens = c(VPhens, selected_line[[paste("VPhen",  i, sep="")]])}
  mean_ind$var = VPhens
  
  mean_ind$fitness = selected_line[["MFit"]]
  
  return(mean_ind)
  
}

# IndivFinale(gen=150)


