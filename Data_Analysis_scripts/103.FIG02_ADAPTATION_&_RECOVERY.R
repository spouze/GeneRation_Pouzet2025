# Based on Script 030 Mostly
# Copy de SCRIPT 100 on 20240422 - réduit au minnimum just pour la figure 2. 
# Revoir  SCRIPT 100 pour plus d'options de figures - genre pour ploter la moyenne etc. 
setwd("/Users/sylvain/Documents/DOCUMENTS/GeneRation/M2_PAPER/")
source("../GeneRation_Fun.R")
source("R_SCRIPTS/Comfort_FUN.R")
topologies ; topocolors # quick check

# LOAD DATA
#### DREAM_SET  - Directly load the dream simus for the burnin  #############
store_Fits_BATCH_succ = read.csv("Extracted_data/dreamburnin_SIMUSET_Fits.csv")[,-1]  #enleve la premiere column
nrow(store_Fits_BATCH_succ) ; play()
generations = 2000

###########################
{ ### LOAD THE ADDITIONAL ONES
store_Fits_add_RANDOM = read.csv("Extracted_data/Massive_set_Whole_burnin_extract/storeFits_20240609_1500+1500+3000_burnin_RANDOM",    header = T)[,-1] ; print("RANDOM - Uploaded")
store_Fits_add_SCALEF = read.csv("Extracted_data/Massive_set_Whole_burnin_extract/storeFits_20240609_1500+1500+3000_burnin_SCALEF",    header = T)[,-1] ; print("SCALEF - Uploaded")
store_Fits_add_HIGHCO = read.csv("Extracted_data/Massive_set_Whole_burnin_extract/storeFits_20240609_1500+1500+3000_burnin_HIGHCO",    header = T)[,-1] ; print("HIGHCO - Uploaded")
store_Fits_add_AUTUMN = read.csv("Extracted_data/Massive_set_Whole_burnin_extract/storeFits_20240520_1x1500_burnin_AUTU1_burnin_AUTU", header = T)[,-1] ; print("AUTUMN - Uploaded")
} ; play() # load

{ # FILTER THE SUCCESSFUL ONES AND ADD THEM
thresh1 = 0.95
thresh2 = 0
store_Fits_add_RANDOM_succ = store_Fits_add_RANDOM[store_Fits_add_RANDOM$gen1000>thresh1 & store_Fits_add_RANDOM$gen2000>thresh2,] ; print("RANDOM - success")
store_Fits_add_SCALEF_succ = store_Fits_add_SCALEF[store_Fits_add_SCALEF$gen1000>thresh1 & store_Fits_add_SCALEF$gen2000>thresh2,] ; print("SCALEF - success")
store_Fits_add_HIGHCO_succ = store_Fits_add_HIGHCO[store_Fits_add_HIGHCO$gen1000>thresh1 & store_Fits_add_HIGHCO$gen2000>thresh2,] ; print("HIGHCO - success")
store_Fits_add_AUTUMN_succ = store_Fits_add_AUTUMN[store_Fits_add_AUTUMN$gen1000>thresh1 & store_Fits_add_AUTUMN$gen2000>thresh2,] ; print("AUTUMN - success")

# print(nrow(store_Fits_add_RANDOM_succ))
# print(nrow(store_Fits_add_SCALEF_succ))
# print(nrow(store_Fits_add_HIGHCO_succ))
# print(nrow(store_Fits_add_AUTUMN_succ))

store_Fits_BATCH = rbind(#store_Fits_BATCH_succ, # garder ou degager les dream
                         store_Fits_add_RANDOM_succ, 
                         store_Fits_add_SCALEF_succ,
                         store_Fits_add_HIGHCO_succ,
                         store_Fits_add_AUTUMN_succ)
}# FILTER #######################

# PRE-PROCESSING : 
# SPLIT IN 3  #############
topologies
topologies2 = c(topologies, "AUTUMN")
for (topo in topologies) {
  # SPLIT to get store_Fits_TOPOLOGY
  store_Fits_TOPOLOGY = store_Fits_BATCH[which(store_Fits_BATCH[,"topo"]==topo),]
  # remove the identifiers for now (two first columns)
  store_Fits_TOPOLOGY = store_Fits_TOPOLOGY[,-c(1,2)]
  # Assign
  assign(paste0("store_Fits_", topo), store_Fits_TOPOLOGY)
  print(paste0(topo, " ", nrow(get(paste0("store_Fits_", topo)))))
  }
# # Previously remove the identifiers (two first columns)
# store_Fits_RANDOM # = store_Fits_RANDOM[,-c(1,2)]
# store_Fits_SCALEF # = store_Fits_SCALEF[,-c(1,2)]
# store_Fits_HIGHCO # = store_Fits_HIGHCO[,-c(1,2)]
# store_Fits_AUTUMN # = store_Fits_AUTUMN[,-c(1,2)]


# ############# 
# {###################### # MOUAI CA MARCHE C'est un peu plus beau, lines less shaky.. et ca change pas le result.
# # NOW ON VA FAIRE CA DIFFEREMMENT POUR QUE CE SOIT PLUS BEAU.
# # COMBINE store_Fits_BATCH_Dream avec les 1000 dernieres colonnes de store_Fits_BATCH (chnge with NA)
# store_Fits_BATCH_single = read.csv("Extracted_data/dream_SIMUSET_Fits.csv")[,-1]
# store_Fits_BATCH_single = store_Fits_BATCH_single[,-c(1,2)]
# # on ajoute des lignes de NA pour avoir du 1000 par set
# n_add_total = 1000
# # RANDOM
# n_add = 1000-nrow(store_Fits_RANDOM)
# df_add = data.frame(matrix(nrow = n_add, ncol = 2000)) ; colnames(df_add) = colnames(store_Fits_RANDOM)
# store_Fits_RANDOM = rbind(store_Fits_RANDOM, df_add)
# # SCALEF
# n_add = 1000-nrow(store_Fits_SCALEF)
# df_add = data.frame(matrix(nrow = n_add, ncol = 2000)) ; colnames(df_add) = colnames(store_Fits_SCALEF)
# store_Fits_SCALEF = rbind(store_Fits_SCALEF, df_add)
# # HIGHCO
# n_add = 1000-nrow(store_Fits_HIGHCO)
# df_add = data.frame(matrix(nrow = n_add, ncol = 2000)) ; colnames(df_add) = colnames(store_Fits_HIGHCO)
# store_Fits_HIGHCO = rbind(store_Fits_HIGHCO, df_add)
# # And finally for all of the them we combine:
# store_Fits_RANDOM[,1:1000] = store_Fits_BATCH_single[   1:1000,]
# store_Fits_SCALEF[,1:1000] = store_Fits_BATCH_single[1001:2000,]
# store_Fits_HIGHCO[,1:1000] = store_Fits_BATCH_single[2001:3000,]
# } ; play()
# ##############

# #############
# {###################### # 
#   # ADDDDITION FROM SINGLE DREAM - modif from previous on 202406171025 mais bon
#   store_Fits_BATCH_single = read.csv("Extracted_data/dream_SIMUSET_Fits.csv")[,-1]
#   store_Fits_BATCH_single = store_Fits_BATCH_single[,-c(1,2)]
#   # on ajoute des lignes de NA pour avoir du 1000 par set
#   n_add_total = 1000
#   # RANDOM
#   n_add = n_add_total
#   df_add = data.frame(matrix(nrow = n_add, ncol = 2000)) ; colnames(df_add) = colnames(store_Fits_RANDOM)
#   store_Fits_RANDOM = rbind(df_add, store_Fits_RANDOM)
#   # SCALEF
#   n_add = n_add_total
#   df_add = data.frame(matrix(nrow = n_add, ncol = 2000)) ; colnames(df_add) = colnames(store_Fits_SCALEF)
#   store_Fits_SCALEF = rbind(df_add, store_Fits_SCALEF)
#   # HIGHCO
#   n_add = n_add_total
#   df_add = data.frame(matrix(nrow = n_add, ncol = 2000)) ; colnames(df_add) = colnames(store_Fits_HIGHCO)
#   store_Fits_HIGHCO = rbind(df_add, store_Fits_HIGHCO)
#   # And finally for all of the them we combine:
#   store_Fits_RANDOM[1:1000,1:1000] = store_Fits_BATCH_single[   1:1000,]
#   store_Fits_SCALEF[1:1000,1:1000] = store_Fits_BATCH_single[1001:2000,]
#   store_Fits_HIGHCO[1:1000,1:1000] = store_Fits_BATCH_single[2001:3000,]
# } ; play()
# ##############


###############################################################################
###############################################################################
#####  PLOT                                                                 ###
###############################################################################
###############################################################################

# Chose which to display: 
chosentopologies = topologies#[c(3,1,2)]
colrs = topocolors
#envcols = c("#C29600", "#DEAC00") # brown / yellow
#envcols = c("#7C7C7C", "#989797") # black / grey
#envcols = c("red", "blue") # black / grey
envcols = c("#B89D6E", "#6EA9B8")
onlysuccessfuls = T # ON PREND QUE LES SUCCESSFUL LA.

width = 20 # unit is cm
height =  11 # unit is cm

# width = 10/2 # unit is cm
# height =  11/2 # unit is cm

{ # ALL PDF
pdf(file = paste0("Rplot", format(Sys.time(), "%Y%m%d%H%M%S"), "_", width, "x", height, "_", thresh1, "_", thresh2, ".pdf"), 
    width = width/cmtoin, height = height/cmtoin)

{ # LAUNCH WHOLE FIGURE

par(mgp=c(2.1, 1, 0)) # default c(3, 1, 0)
########################################################
main = ""#"Adaptation from evolved networks"
plot(NA, xlim=c(0,generations), ylim=c(0,1), xlab="Generations", ylab="Population Fitness", 
     main = main, bty ='n') #bty ='n' removes the plot
abline(v=1000, col="grey", lty="dashed") ; 
abline(h=c(0.5, 1), col="grey", lwd=.5)
#abline(h=seq(0,1,.1), col="grey", lwd=.5)

#########################
####### MEDIAN ##########
#########################

# Fond simple
ph = 0.2 # for polygon_height
polygon(c(   0, 0,    1000, 1000), c(0, ph, ph, 0), col=adjustcolor(envcols[1], .2), border = NA) # ENVIRONMENT #1
polygon(c(1000, 1000, 2000, 2000), c(0, ph, ph, 0), col=adjustcolor(envcols[2], .2), border = NA) # ENVIRONMENT #2

# WRITE ENV
if(onlysuccessfuls){textposition=c(500, 1500)}else{textposition=c(700, 1300)}
text(x = textposition[1],  y = 0.1, 'Environment 1', srt = 0, cex = 1.5, col = envcols[1])
text(x = textposition[2],  y = 0.1, 'Environment 2', srt = 0, cex = 1.5, col = envcols[2])


# Q1 and Q3 polygon
for (topo in chosentopologies) {
  store_Fits_XXXXXX = get(paste0("store_Fits_",topo))
  colr = colrs[match(topo, topologies)] ; col2 = "grey"
  q1 = apply(store_Fits_XXXXXX, 2, function(x) quantile(x, probs = 0.25, na.rm = TRUE))
  q2 = apply(store_Fits_XXXXXX, 2, function(x) quantile(x, probs = 0.75, na.rm = TRUE)) # Q3 genre
  #lines(c(1:generations), q1, col=col2)
  #lines(c(1:generations), q2, col=col2)
  polygon(c(1:generations, generations:1), c(q1, rev(q2)), col=adjustcolor(colr, .2)
          , border = NA #adjustcolor(colr, .5)
          )
}

# MEDIANs
for (topo in chosentopologies) {
  store_Fits_XXXXXX = get(paste0("store_Fits_",topo)) ; 
  colr = colrs[match(topo, topologies)]
  # colr = "black"
  qmed = apply(store_Fits_XXXXXX, 2, function(x) quantile(x, probs = 0.5, na.rm = TRUE))
  lines(1:generations, qmed, col=colr, lwd=2)
}

# LEGEND
if(onlysuccessfuls){yposition=0.85}else{yposition=0.45}
text(x = 2000, y = yposition, # "topright", # x=800, y=.4, # 0.85 for successful, 0.45 for all simus
     c("Highly-connected", "\n\n\nRandom", "\n\n\n\n\n\nScale-free"), 
     col=c("red3", "green4", "blue3"), bty="n", pos = 2) # pos 2 = aligned to the right

} # LAUNCH WHOLE FIGURE
dev.off() ; play()# for the end of the PDF
} # ALL PDF
