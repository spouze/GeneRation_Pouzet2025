setwd("/Users/sylvain/Documents/DOCUMENTS/GeneRation")
#source("GeneRation_Fun.R")
setwd("M2_PAPER/Extracted_data/")
source("../R_SCRIPTS/IndivRecap_fromSET_FUN.R")

dream_SET = read.csv("dream_SIMUSET_Lasts.csv", header = T)
dream_SET = dream_SET[,-1]
# If already loaded
dream_SET = dream_Lasts

######################################
# CREATIon OF THE IndivRecap_fromSET_FUN.R FUNCTION
##########################################################
# Check how MEAN IND FIT compared to the POP MEAN FIT

fitsallinds = c()
for (i in 1:3000){
  fitsallinds = c(fitsallinds, IndivRecap(dream_SET, i, clean = T)$fitnessclean)
  cat(paste0(i, "."))
}
dream_SET$fitsallindsclean = fitsallinds
play()

#topologies = c("RANDOM", "SCALEF", "HIGHCO")
#colrs      = c("green4", "blue3" , "red3" )

plot(NA, NA, xlim=c(0.957,1), ylim=c(0.953,1), xlab = "meaninds fits", ylab="CLEANED meaninds fits")
abline(v=c(0.96, 0.97, 0.98, 0.99, 1), col="grey")
abline(h=c(0.96, 0.97, 0.98, 0.99, 1), col="grey")
lines(c(0,2), c(0,2), lwd=5, col="grey")
points(dream_SET$fitsallinds, dream_SET$fitsallindsclean, pch=20, cex = 0.6,
       col=adjustcolor(colrs[match(dream_SET$topo, topologies)], 0.5))

# NOTE= THE CLEANED UP VERSION WILL ONLY BE USED TO COMPUTE DEGREE ETC. 

######################################################
# ON REPREND CE QUI AVAIT ETE FAIT AVANT pour les degrees etc
######################################################

line = 1003

#INITIALIZATION
Wdegree_df = data.frame(matrix(NA, ncol = 8, nrow = 10))
colnames(Wdegree_df)=c("line", "id", "topo", "gene", 
                       "indegree", "outdregree", "inout_cross", "inout_add")
Wmean = IndivRecap(dream_SET, line = line, clean = T)$indclean
id    = dream_SET[line,]$id
topo  = dream_SET[line,]$topo

for (gene in 1:10) {
  
  #ext pour "extended"
  degree_ext_in    = which(Wmean[gene,1:10]!=0)
  degree_ext_out   = which(Wmean[,gene]!=0)
  degree_ext_inout = unique(c(degree_ext_in, degree_ext_out)) #???
  
  degree_in    = length(degree_ext_in   )
  degree_out   = length(degree_ext_out  )
  degree_inout = length(degree_ext_inout)
  
  Wdegree_df[gene,]$line        = line
  Wdegree_df[gene,]$id          = id
  Wdegree_df[gene,]$topo        = topo
  Wdegree_df[gene,]$gene        = gene
  Wdegree_df[gene,]$indegree    = degree_in
  Wdegree_df[gene,]$outdregree  = degree_out
  Wdegree_df[gene,]$inout_cross = degree_inout
  Wdegree_df[gene,]$inout_add   = degree_in + degree_out
}

Wdegree_df


###########################################################################
###########################################################################
# On loop ca et c'est parti
###########################################################################
###########################################################################

#Create master dataframe
SUPER_Wdegree_df = data.frame(matrix(NA, ncol = 8, nrow = 10*3000))
colnames(SUPER_Wdegree_df)=c("line", "id", "topo", "gene", 
                       "indegree", "outdregree", "inout_cross", "inout_add")



for (line in 1:3000){
  
  #INITIALIZATION
  #Wdegree_df = data.frame(matrix(NA, ncol = 8, nrow = 10))
  #colnames(Wdegree_df)=c("line", "id", "topo", "gene", 
  #                       "indegree", "outdregree", "inout_cross", "inout_add")
  
  #Select matrix
  #Wmean = IndivRecap(dream_SET, line = line, clean = T, cleanthresh = 0.02)$indclean
  #Wmean = WClean(IndivRecap(dream_SET, line = line)$ind, threshold = 0.05)
  # AHH PUREE IL FAUT INTEGRER LE CODING VECTOR.
  Wmean = WConvert(IndivRecap(dream_SET, line = line, clean = T, cleanthresh = 0.01)$indclean)
  
  id    = dream_SET[line,]$id
  topo  = dream_SET[line,]$topo
  
  for (gene in 1:10) {
    
    #ext pour "extended"
    degree_ext_in    = which(Wmean[gene,1:10]!=0)
    degree_ext_out   = which(Wmean[,gene]!=0)
    degree_ext_inout = unique(c(degree_ext_in, degree_ext_out)) #???
    
    degree_in    = length(degree_ext_in   )
    degree_out   = length(degree_ext_out  )
    degree_inout = length(degree_ext_inout)
    
    #Wdegree_df[gene,]$line        = line
    #Wdegree_df[gene,]$id          = id
    #Wdegree_df[gene,]$topo        = topo
    #Wdegree_df[gene,]$gene        = gene
    #Wdegree_df[gene,]$indegree    = degree_in
    #Wdegree_df[gene,]$outdregree  = degree_out
    #Wdegree_df[gene,]$inout_cross = degree_inout
    #Wdegree_df[gene,]$inout_add   = degree_in + degree_out
    
    # Remplit le master dataframe
    SUPER_Wdegree_df[(((line-1)*10)+gene),] = c(line, id, topo, 
                                                gene, 
                                                degree_in, degree_out, degree_inout, 
                                                (degree_in + degree_out))
  }
  cat(paste0(line,".")) #just for notice
}

play()
# write.csv(SUPER_Wdegree_df, file = "Degree_count_dream.csv")
# write.csv(SUPER_Wdegree_df, file = "Degree_count_dreamWClean0.01_integreatedCvec.csv")

# OR LOAD THE DATA 
#SUPER_Wdegree_df = read.csv("Extracted_data/Fig02_Degree/Degree_count_dreamWClean0.01_integreatedCvec.csv")
SUPER_Wdegree_df = read.csv("Extracted_data/Figures_vrac/Fig02_Degree/Degree_count_dreamWClean0.01_integreatedCvec.csv")
View(SUPER_Wdegree_df)

##################################################
############### PLOT DEGREES #####################
##################################################

#topologies = c("RANDOM", "SCALEF", "HIGHCO")
#colrs = c("green4", "blue3",  "red3")

par(mfrow=c(1,3))

for (topo in topologies){
  #topo = "SCALEF"
  #SUPER_Wdegree_df[which(SUPER_Wdegree_df$topo==topo),]$inout_add
  
  #isolate
  table(SUPER_Wdegree_df[which(SUPER_Wdegree_df$topo==topo),]$degree_inout)
  #preprocess
  occurences = as.data.frame(table(SUPER_Wdegree_df[which(SUPER_Wdegree_df$topo==topo),]$inout_add), stringsAsFactors = F)
  occurences$Var1=as.integer(occurences$Var1)
  #add non-existing values
  if (min(occurences$Var1)>0)  {for (i in 0:(min(occurences$Var1)-1)) {occurences = rbind(occurences, c(i,0))}}
  if (max(occurences$Var1)<18) {for (i in (max(occurences$Var1)+1):18){occurences = rbind(occurences, c(i,0))}}
  
  #reorder
  occurences = occurences[order(occurences$Var1),]
  
  barplot(occurences$Freq/sum(occurences$Freq), 
          col=colrs[match(topo, topologies)], main = topo, 
          names.arg = occurences$Var1, ylim=c(0,0.36),
          ylab="Frequency", xlab="Degree (connexions per gene)")
  
}

resetplot()

#################
#############################
###################################################
## On va rajouter les 20 de AUTUMN - script 160
