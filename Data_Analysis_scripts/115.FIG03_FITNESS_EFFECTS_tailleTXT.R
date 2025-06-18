## BASED ON 080.Plot_Tests_02.R

source(here::here("GeneRation_Fun_v1.R"))
source(here::here("Data_Analysis_scripts/Comfort_FUN.R"))
source(here::here("Data_Analysis_scripts/128.New_Mutate_in_new_Env_4.R"))
source(here::here("Data_Analysis_scripts/IndivRecap_fromSET_FUN.R"))

library(RColorBrewer)


# cnames ; cEnames
dream_SET = read.csv("Extracted_data/dream_SIMUSET_Lasts.csv", header = T)[,-1]
###############################################################################

###################
# LOAD TESTS ######
###################
# Load just one quickly:
#allmutations = read.csv("Extracted_data/allmutations_20240417_fig06/allmutations001.csv")[,-1]
allmutations = read.csv("Extracted_data/allmutations_20240314/allmutations00.csv")[,-1]
nrow(allmutations)

{# LOAD ALL ####################################################################

#setwd("Extracted_data/allmutations_20240417_fig06/")
allmut_files = list.files(path=here::here("Extracted_data/allmutations_20240314"), pattern = "allmutations.{2,3}\\.csv", full.names=TRUE) #"allmutations.*\\.csv"
n_batches = 15 #length(allmut_files)
n_batches = length(allmut_files)
all_batches = c()

for (batch_x_file in allmut_files[1:n_batches]){
  batch_x_name = substr(batch_x_file, start = 1, stop = 14)
  all_batches = c(all_batches, batch_x_name)
  assign(batch_x_name, read.csv(batch_x_file)[,-1])
  cat(paste0(sprintf("%02d", match(batch_x_file, allmut_files)), "/", 
             sprintf("%02d", length(allmut_files[1:n_batches])), " - ", batch_x_file, "\n"))
} ; play(1)

allmutations = do.call(rbind, lapply(all_batches, function(x) get(x))) ; play(2) # combine all the allmuations
rm(list = all_batches) #remove individual batches
nrow(allmutations)
play()# LOAD ALL
}

####################################################################
{## PRE-TREATMENT
####################################################################

# Rassemble DUP et DEL for easy ploting
allmutations[allmutations$mut_type=="DEL",]$dupdel_number = -allmutations[allmutations$mut_type=="DEL",]$dupdel_number
allmutations[allmutations$mut_type=="DEL",]$mut_type = "DUP"

# Define different classes for effect on fitness.
effect_class_tresh = rev(c(0.1, 0.01, 0.001, -0.001, -0.01, -0.1))
effect_class= rev(c("BENEF3", "BENEF2", "BENEF1", "NEUTR", "DELET1", "DELET2", "DELET3"))
#effect_class_color = rev(c(colorRampPalette(c("#066400","#989797", "#831D01"))(7)))
effect_class_color = rev(c(colorRampPalette(c("#066400","#BFBFBF", "#831D01"))(7)))
allmutations$fit_effect <- assign_category(allmutations$fit_delta)
} ; play()
####################################################################
## FOR RELATIVE EFFECTS
####################################################################

# ...
# Relative values? 
# (mysubset$mutvalue_final/mysubset$mutvalue_init)[1:15]

####################################################################
## PLOT BARPLOTS EFFECTS ON FITNESS
####################################################################
width = 34 # unit is cm
height =  19 # unit is cm
{ # FULL PDF FIGURE
pdf(file = paste0("Rplot", format(Sys.time(), "%Y%m%d%H%M%S"), "_", width, "x", height, ".pdf"), 
    width = width/cmtoin, height = height/cmtoin)

{ # DRAW WHOLE FIGURE

par(mfrow=c(3,3))
par(oma = c(5.3, 5.5, 3, 0)) # Outer (overall) MArgin
par(mar=c(1,0,1.5,0)+0.1) # plot inside MARgins
draw_box = F #for box("figure", col="forestgreen)

barplot_space_neutr  = 0.14 # default 0.2
barplot_space_others = 0.03# cannot be NA: has to be 0
#barplot_space = 0.2 
barplot_border_color = NA #"white" # default black

plotcounter = 0
for (topology in c("HIGHCO", "RANDOM", "SCALEF")){
  #topology = "RANDOM"
  #Extract data
  pool_topology = allmutations[allmutations$topology==topology,]
  
  for (mut_type in c("REG", "COD", "DUP")){
    #mut_type = "REG"
    #Extract data
    plotcounter = plotcounter+1
    pool_mut_type = pool_topology[pool_topology$mut_type==mut_type,]
    #find the right muteff column
    muteff_var = c("reg_muteff", "cod_muteff", "dupdel_number")[match(mut_type, c("REG", "COD", "DUP"))] # was dup_number before
    #Get the different values tested
    muteff_levels = as.numeric(levels(as.factor(pool_mut_type[,muteff_var])))
    
    # deal with spaces
    barplot_space = rep(barplot_space_others, length(muteff_levels))
    barplot_space[c(which(muteff_levels==0), which(muteff_levels==0)+1)] = c(barplot_space_neutr, barplot_space_neutr)
    
    #Setup empty dataframe for results
    bilan_pool_muteff = data.frame(matrix(ncol = length(muteff_levels), nrow = length(effect_class)))
    colnames(bilan_pool_muteff) = muteff_levels
    rownames(bilan_pool_muteff) = rev(effect_class)
    
    #Compute and fill in the dataframe
    for (used_muteff in colnames(bilan_pool_muteff)){ # 0 0.1 0.2 0.3
      pool = pool_mut_type[pool_mut_type[,muteff_var] == used_muteff,]
      for (fit_effect in rownames(bilan_pool_muteff)){ #"NEUTR" "BENEF" "DELET"
        bilan_pool_muteff[fit_effect, used_muteff] = sum(pool$fit_effect == fit_effect)/nrow(pool)  
      }
    }
    
    cat(paste0(topology, " - ", mut_type))
    bilan_pool_muteff
    
    # Barplots with colored coreders and an extra large black border (2nd barplot)
    coordinates = barplot(as.matrix(bilan_pool_muteff), col=rev(effect_class_color), border = NA, xaxt = "n", yaxt = "n", space = barplot_space)
    barplot(colSums(as.matrix(bilan_pool_muteff)), border = barplot_border_color, add = T, col="NA", xaxt = "n", yaxt = "n", space = barplot_space)
    
    # Axis (with plotcounter)
    
    if (sum(plotcounter==c(1,4,7)) == 1){
      #barplot(colSums(as.matrix(bilan_pool_muteff)), border=NA, col=NA, xaxt = "n", add = T)
      axis(2, at = c(0, 0.25, 0.5, 0.75, 1), tick = T, pos = -0.6,
           labels = c(0, "", 0.5, "", 1), cex.axis = 1.5)
      } # adds y axis
    
    decalage_bas = -0.16
    taille_text_gris = 1.5
    talle_nombres_gris = 1.5
    if (sum(plotcounter==c(7,8,9)) == 1){
      #barplot(colSums(as.matrix(bilan_pool_muteff)), border=NA, col=NA, yaxt = "n", add = T)
      if (mut_type == "REG" | mut_type == "COD"){
        axis(1, at = coordinates[c(1, 6, 11, 16, 21)], tick = T, pos = -0.05, #pos empeche que ce soit collé au barplot
             #labels =c("-1", "","","","", "-0.5", "","","","", "0", "","","","", "+0.5", "","","","", "+1"))}
             labels =c("-1", "-0.5", "0",  "+0.5","+1"), 
             cex.axis = talle_nombres_gris)}
      if (mut_type == "REG"){
        axis(1, at = coordinates[c(6)], tick = F, pos = decalage_bas, 
             labels = c("towards repression"), 
             cex.axis = taille_text_gris, col.axis = "#808080")
        axis(1, at = coordinates[c(16)], tick = F, pos = decalage_bas, 
             labels = c("towards activation"), 
             cex.axis = taille_text_gris, col.axis = "#808080")}
      if (mut_type == "COD"){
        axis(1, at = coordinates[c(6)], tick = F, pos = decalage_bas, 
             labels = c("towards loss of function"), cex.axis = taille_text_gris, col.axis = "#808080")
        axis(1, at = coordinates[c(16)], tick = F, pos = decalage_bas, 
             labels = c("towards gain of function"), cex.axis = taille_text_gris, col.axis = "#808080")
        }
      if (mut_type == "DUP"){
        axis(1, at = coordinates[c(1, 4, 6, 8, 11, 16)], tick = T, pos = -0.05,
             #labels =c("-5", "","","-2","", "0", "","+2","","", "+5", "","","","", "+10"))}
             labels =c("-5","-2", "0", "+2", "+5","+10"), 
             cex.axis = talle_nombres_gris)
        axis(1, at = coordinates[c(4, 11)], tick = F, pos = decalage_bas, labels = c("Deletions", "Duplications"), 
             cex.axis = taille_text_gris, col.axis = "#808080")
        }
      } # adds x axis
    if (plotcounter == 1){
      # Add LEGEND
      text(0.5, 0.85, "Deleterious", adj=0, cex=2, font = 2, col="white")
      text(0.5, 0.34, "Neutral",     adj=0, cex=2, font = 2, col="white")
      text(0.5, 0.11, "Beneficial",  adj=0, cex=2, font = 2, col="white")
    }
    if (draw_box){box("figure", col="forestgreen"); box(which="plot", col="red")}
    mtext(paste0(LETTERS[plotcounter], "."), side = 3, line = 0.3, at = 0 , outer = FALSE, cex = 1.5, font =2)
    cat (" - Completed.\n") ; play(2)
  }# for (mut_type
}# for (topology
play()
}# DRAW WHOLE FIGURE

{# ADD TEXT
mtext('Regulatory mutations'  , side = 3, line = 0.2, at = 0.333/2, outer = TRUE, cex = 1.5, font = 2)
mtext('Coding mutations'      , side = 3, line = 0.2, at = 0.5    , outer = TRUE, cex = 1.5, font = 2)
mtext('Deletion / Duplication', side = 3, line = 0.2, at = 1.666/2, outer = TRUE, cex = 1.5, font = 2)

mtext('Highly Connected', side = 2, line = 3.8, at = 1.666/2, outer = TRUE, col="red3"       , font = 2, cex=1.3)
mtext("Random",           side = 2, line = 3.8, at = 0.5    , outer = TRUE, col="forestgreen", font = 2, cex=1.3)
mtext('Scale-Free',       side = 2, line = 3.8, at = 0.333/2, outer = TRUE, col="blue3"      , font = 2, cex=1.3)

taille_text_noir = 1.2

mtext("Frequency", side = 2, line = 2.2, at = 1.666/2 , outer = TRUE, cex = taille_text_noir) #0.7
mtext("Frequency", side = 2, line = 2.2, at = 0.5     , outer = TRUE, cex = taille_text_noir)
mtext("Frequency", side = 2, line = 2.2, at = 0.333/2 , outer = TRUE, cex = taille_text_noir)

descente_text = 3.8
mtext('Mutation size (Regulation)', side = 1, line = descente_text, at = 0.333/2, outer = TRUE, cex = taille_text_noir) #0.7
mtext('Mutation size (Activity)', side = 1, line = descente_text, at = 0.5    , outer = TRUE, cex = taille_text_noir)
mtext('Number of genes'           , side = 1, line = descente_text, at = 1.666/2, outer = TRUE, cex = taille_text_noir)
}# ADD TEXT

dev.off() ; play()
} # FULL PDF FIGURE
################################################################
####### CHEKING THINGS ##############

# check_mut_type = "COD"
# mysubset = allmutations[allmutations$mut_type==check_mut_type,]
# plot(mysubset$mutvalue_init, pch=20, col=adjustcolor(topocolors[match(mysubset$topology, topologies)], 0.3), main= "INIT")
# plot(mysubset$mutvalue_final, pch=20, col=adjustcolor(topocolors[match(mysubset$topology, topologies)], 0.3), main= "FINAL")
# plot(mysubset$mutvalue_delta, pch=20, col=adjustcolor(topocolors[match(mysubset$topology, topologies)], 0.3), main= "DELTA")
# abline(h=c(0,1, -1))


