
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
allmutations = read.csv("Extracted_data/allmutations_20240417_fig06/allmutations001.csv")[,-1]
nrow(allmutations)

{# LOAD ALL ####################################################################
setwd("/Users/sylvain/Documents/DOCUMENTS/GeneRation/M2_PAPER")
setwd("Extracted_data/allmutations_20240417_fig06/")
allmut_files = list.files(path=here::here("Extracted_data/allmutations_20240417_fig06/"), pattern = "allmutations.{2,3}\\.csv", full.names=TRUE)
#n_batches = 20 #length(allmut_files)
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
play()
} # LOAD ALL

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

# REMOVE NO MUTATIONS # Needed parce qu'on veut pas biaiser nos data (?) - on a pas besoin de la barre à zero
allmutations = allmutations[allmutations$notes != "NO MUTATION",]

# Remove DUP DEL en fait on en a pas besoin
# allmutations = allmutations[allmutations$mut_type != "DUP",]
# allmutations = allmutations[allmutations$mut_type != "DEL",]
} ; play() # pre-treatment

##################
###################################
##################################################
################################################################################
################################################################################
################################################################################
{ # LANCE ET SAVE PDF
width = 23 # unit is cm
height =  20 # unit is cm
pdf(file = paste0("Rplot", format(Sys.time(), "%Y%m%d%H%M%S"), "_", width, "x", height, ".pdf"), 
    width = width/cmtoin, height = height/cmtoin)

{ # DRAW WHOLE FIGURE

# LAYOUT
par(oma = c(3, 5.7, 3, 0)) # Outer (overall) MArgin
mar_plot_pleio = c(0.2, 0.2, 2.5, 0.5)+0.1 # plot inside MARgins
mar_plot_cis   = c(0.1, 0.2, 0.1, 0.5)+0.1   # plot inside MARgins

layout_matrix <- matrix(c(1,3,5,2,4,6,7,9,11,8,10,12,13,15,17,14,16,18), nrow = 6, byrow = TRUE)
layout(layout_matrix, heights = rep(c(4,0.4), 9))  # heights are ratios
#layout.show(18)
draw_box = F #for box("figure", col="forestgreen)

# ESTHETICS
color_cis = adjustcolor("#EDD405", 0.6) #yellow
color_cis_grey = "#DBDBDB"
barplot_space = 0.05 # default 0.2
barplot_border_color = NA #"white" # default black
segment_color = "white"

# PARAMETERS
pleio_levels=0:10
tout_a_100 = F # si on met tout à 100, alors tous les bars sont entre 0 et 100. Sinon, c'est relatif
include_minus_values = T #include nnot only 0.5 but -0.5 too - and not only DUP 1, but also DEL 1
degage_DUP = T # en fait c'est les DEL qu'on degage (si degage_DUP = T alors on dégage les DEL)
focus_genes = "ALL" # "SELECTED" "FREE" # "ALL" # pour la sup figure

######### START PLOT ##################
plotcounter = 0
for (topology in c("HIGHCO", "RANDOM", "SCALEF")){
  #topology = "HIGHCO"
  pool_topology = allmutations[allmutations$topology==topology,]
  
  for (mut_type in c("REG", "COD", "DUP")){
    #mut_type = "REG"
    plotcounter = plotcounter+1
    #plotcounter = 1
    cat(paste0(topology, " - ", mut_type, " - "))
    pool_mut_type = pool_topology[pool_topology$mut_type==mut_type,]
    
    #Find the right muteff column and value
    set_muteffs= c(0.5,          0.5,         1               )
    muteff_var = c("reg_muteff", "cod_muteff", "dupdel_number")[match(mut_type, c("REG", "COD", "DUP"))]
    muteff_val = set_muteffs[match(mut_type, c("REG", "COD", "DUP"))]
    # pool for muteff = 0.5 or dup =1
    if(include_minus_values){imv = -1} else {imv=1}
    if(degage_DUP){if(mut_type=="DUP"){imv=1}} # ICI ON EMPECHE LES DELETIONS
    pool_filter = pool_mut_type[(pool_mut_type[,muteff_var]== muteff_val     |
                                 pool_mut_type[,muteff_var]== muteff_val*imv ),]
    pool_filter$CIS = pool_filter$pleiotropy * pool_filter$Expr_effect_cisness
    
    ####################################
    # Nouveau truc pour séparer 
    
    if (focus_genes=="SELECTED"){focus_genes_list = c(1,2,3,4,5 ); cat("**SELECTED GENES** " )}
    if (focus_genes=="FREE"    ){focus_genes_list = c(6,7,8,9,10); cat("**FREE GENES** ")}
    if (focus_genes=="ALL"){focus_genes_list = c(1,2,3,4,5,6,7,8,9,10)}
    
    pool_filter = pool_filter[(pool_filter[,"i"]            %in% focus_genes_list    |
                                pool_filter[,"dupdel_genes"] %in% focus_genes_list),]
    ####################################
    
    
    #Setup empty dataframes for results
    # PLEIOTROPY
    bilan_pool_pleio = data.frame(matrix(ncol = length(pleio_levels), nrow = length(effect_class)))
    colnames(bilan_pool_pleio) = pleio_levels
    rownames(bilan_pool_pleio) = rev(effect_class)
    #CIS / TRANS
    bilan_pool_cis = data.frame(matrix(ncol = length(pleio_levels), nrow = 1))
    colnames(bilan_pool_cis) = pleio_levels
    rownames(bilan_pool_cis) = "cisness"
    
    
    #Compute and fill in the dataframe
    for (used_pleio in colnames(bilan_pool_pleio)){ # 0:10
      # used_pleio = "0"
      cat(paste0(sprintf("%02d", as.integer(used_pleio)), " "))
      pool = pool_filter[pool_filter$pleiotropy == used_pleio,] # column is affected_genes in old allmutations sets and not pleiotropy (recents sets)
      
      # PLEIOTROPY
      for (fit_effect in rownames(bilan_pool_pleio)){ #"NEUTR" "BENEF" "DELET"
        # fit_effect = "NEUTR"
        if (tout_a_100){
          bilan_pool_pleio[fit_effect, used_pleio] = sum(pool$fit_effect == fit_effect)/nrow(pool) # Everybody at 100
        }else{
          bilan_pool_pleio[fit_effect, used_pleio] = sum(pool$fit_effect == fit_effect)/nrow(pool_filter) # See relative numbers
        } # end if tout à 100
      } # end for fit_effect
      
      # CIS / TRANS
      bilan_pool_cis["cisness", used_pleio] = mean(pool$CIS)
      
    } # end for used_pleio
    
    cat(paste0("- Done (", plotcounter,").\n"))
    
    # PLOT PLEIOTROPY
    par(mar=mar_plot_pleio)
    if (tout_a_100 == F){
      ylim = c(0,0.42)
      if (plotcounter==7){ylim = c(0,0.6)}
        barplot(as.matrix(bilan_pool_pleio), col=rev(effect_class_color), border = NA, xaxt = "n", yaxt = "n", ylim=ylim, space = barplot_space)
        barplot(colSums(as.matrix(bilan_pool_pleio)), add = T, col="NA", yaxt = "n", ylim=ylim, names.arg = NA, space = barplot_space, border=barplot_border_color)
      # y axes
      if (sum(plotcounter==c(1,4)) == 1){
        #barplot(colSums(as.matrix(bilan_pool_pleio)), border=NA, col=NA, xaxt = "n", add = T)
        axis(side=2, at=seq(0, 0.4, by=0.1), labels=c("", "0.1", "0.2", "0.3", "0.4"), cex.axis =1.5)
        }
      if (sum(plotcounter==c(2,5,6,8,9)) == 1){
        axis(side=2, at=seq(0, 0.4, by=0.1), labels=c("", "", "", "", 0.4), cex.axis =1.5)}
      if (plotcounter==3){
        axis(side=2, at=seq(0, 0.4, by=0.1), labels=c("", "", "", "", ""))}
        if (plotcounter==7){
          axis(side=2, at=c(0, 0.1, 0.2, 0.3, 0.4), labels=c("", "", "0.2", "", "0.4"), cex.axis =1.5)}
    }
    
    if (tout_a_100 == T){
      mar_plot_pleio2 = mar_plot_pleio + c(0, 0, 0.4, 0.4)
      mar_plot_cis2   = mar_plot_cis   + c(0, 0, 0.0, 0.4)
      par(mar=mar_plot_pleio2)
      barplot(as.matrix(bilan_pool_pleio), col=rev(effect_class_color), border = rev(effect_class_color), xaxt = "n", yaxt = "n", ylim=c(0,1), space = barplot_space, names.arg = NULL)
      #barplot(colSums(as.matrix(bilan_pool_pleio)), add = T, col="NA", yaxt = "n", ylim=c(0,1), space = barplot_space, border = barplot_border_color, names.arg = NULL)
      # y axes
      if (sum(plotcounter==c(1,4,7)) == 1){      axis(side=2, at=seq(0, 1, by=0.2), labels=c("", "0.2", "0.4", "0.6", "0.8", "1"))}
      if (sum(plotcounter==c(2,3,5,6,8,9)) == 1){axis(side=2, at=seq(0, 1, by=0.2), labels=rep("", 6))}
    }
    mtext(paste0(LETTERS[plotcounter], "."), side = 3, line = 0.3, at = 0 , outer = FALSE, font=2, cex = 1.2)
    if (draw_box){box("figure", col="forestgreen"); box(which="plot", col="red")}
    
    if (plotcounter == 1 & tout_a_100 == F){
      # Add LEGEND
      my_cex = 1.2
      my_font=2
      remov=0.042
      text(0.1, 0.37-remov*0, "Deleterious", adj=0, cex=1.4, font = 2, col=effect_class_color[2])
      text(0.1, 0.37-remov*1, "Neutral",     adj=0, cex=1.4, font = 2, col=effect_class_color[4])
      text(0.1, 0.37-remov*2, "Beneficial",  adj=0, cex=1.4, font = 2, col=effect_class_color[6])
    }
    
    if (plotcounter == 9 & tout_a_100 == F){
      # ADD LEGEND
      par(family = "Courier")
      remov=0.03
      my_cex = 1.2
      my_font=2
      text(midpoints[11]+0.5, 0.4-remov*0, expression("         " * Delta * w * " < -0.1  ") , adj=1, col=effect_class_color[1], cex=my_cex, font=my_font)
      text(midpoints[11]+0.5, 0.4-remov*1, expression("  -0.1 < " * Delta * w * " < -0.01 ") , adj=1, col=effect_class_color[2], cex=my_cex, font=my_font)
      text(midpoints[11]+0.5, 0.4-remov*2, expression(" -0.01 < " * Delta * w * " < -0.001") , adj=1, col=effect_class_color[3], cex=my_cex, font=my_font)
      text(midpoints[11]+0.5, 0.4-remov*3, expression("-0.001 < " * Delta * w * " <  0.001") , adj=1, col=effect_class_color[4], cex=my_cex, font=my_font)
      text(midpoints[11]+0.5, 0.4-remov*4, expression(" 0.001 < " * Delta * w * " <  0.01 ") , adj=1, col=effect_class_color[5], cex=my_cex, font=my_font)
      text(midpoints[11]+0.5, 0.4-remov*5, expression("  0.01 < " * Delta * w * " <  0.1  ") , adj=1, col=effect_class_color[6], cex=my_cex, font=my_font)
      text(midpoints[11]+0.5, 0.4-remov*6, expression("   0.1 < " * Delta * w * "         ") , adj=1, col=effect_class_color[7], cex=my_cex, font=my_font)
      par(family = "")
    }
    
    # PLOT CIS/TRANS
    if (tout_a_100 == T){par(mar=mar_plot_cis2)}else{par(mar=mar_plot_cis)}
    if (mut_type!="DUP"){
      plotable_data_cis = c(1, as.numeric(bilan_pool_cis[1,])[2:11])
      color_cis_array = c(color_cis_grey, rep(color_cis, 10))
      midpoints = barplot(plotable_data_cis, col=color_cis_array, border = NA, axes = F, space = barplot_space, names.arg = NA) #, names.arg = c(0:10)
      barplot(rep(1, 11), col = NA, add=T, names.arg = NA, axes = F, space = barplot_space, border=barplot_border_color)
      segments(midpoints[1]-0.5, 0, midpoints[1]+0.5, 1, col = segment_color)
      if (sum(plotcounter==c(1,4,7)) == 1){axis(side = 2, at = c(0,1), labels = c(0,1))
        } else {axis(side = 2, at = c(0,1), labels = c("", ""))}
    } else { # if DUP
      midpoints = barplot(rep(1, 11), col = color_cis_grey, axes = F, space = barplot_space, border=NA, names.arg = NA)
                  barplot(rep(1, 11), col = NA, axes = F, names.arg = NA, space = barplot_space, border=NA, add=T)
      axis(side = 2, at = c(0,1), labels = c("", ""))
      for (i in c(midpoints)){
        segments(i-0.5, 0, i+0.5, 1, col = segment_color)      
      } # if (mut_type!="DUP")
    } # else if DUP or not
    
    text(midpoints, 0.5, labels = 0:10, font = 2, col="#696969", cex=1.2)
    
    if (draw_box){box("figure", col="forestgreen"); box(which="plot", col="red")}
    
  } # for (mut_type in c("REG", "COD", "DUP"))
} # for (topology in c("HIGHCO", "RANDOM", "SCALEF"))

if (focus_genes=="SELECTED"){cat("ATTENTION - focus_genes is set on SELECTED - genes under selection only - ATTENTION")}
if (focus_genes=="FREE"    ){cat("ATTENTION - focus_genes is set on FREE - genes not-under selection only - ATTENTION")}

} # DRAW WHOLE FIGURE


{# TEXT
  mtext(paste0("Regulatory mutations (+/-", set_muteffs[1], ")"), side = 3, line = 0.7, at = 0.3333*0.5    , outer = TRUE, cex = 1.2, font=2)
  mtext(paste0("Coding mutations (+/-",     set_muteffs[2], ")"), side = 3, line = 0.7, at = 0.3333*1.5    , outer = TRUE, cex = 1.2, font=2)
  mtext(paste0("Duplications (+",           set_muteffs[3], ")"), side = 3, line = 0.7, at = 0.3333*2.5    , outer = TRUE, cex = 1.2, font=2)
  
  mtext('Highly Connected', side = 2, line = 3.9, at = 1.666/2, outer = TRUE, col="red3"       , font = 2, cex=1.2)
  mtext("Random",           side = 2, line = 3.9, at = 0.5    , outer = TRUE, col="forestgreen", font = 2, cex=1.2)
  mtext('Scale-Free',       side = 2, line = 3.9, at = 0.333/2, outer = TRUE, col="blue3"      , font = 2, cex=1.2)
  
  mtext('Pleiotropy', side = 1, line = 0.7, at = 0.3333*0.5    , outer = TRUE, cex = 1.1)
  mtext('Pleiotropy', side = 1, line = 0.7, at = 0.3333*1.5    , outer = TRUE, cex = 1.1)
  mtext('Pleiotropy', side = 1, line = 0.7, at = 0.3333*2.5    , outer = TRUE, cex = 1.1)
    
  if (tout_a_100 == F){
    mtext("Density", side = 2, line = 2.3, at = 1.666/2 , outer = TRUE, cex = 0.9)
    mtext("Density", side = 2, line = 2.3, at = 0.5     , outer = TRUE, cex = 0.9)
    mtext("Density", side = 2, line = 2.3, at = 0.333/2 , outer = TRUE, cex = 0.9)
  }
  if (tout_a_100 == T){
    mtext("Frequency", side = 2, line = 2.2, at = 1.666/2 , outer = TRUE, cex = 0.7)
    mtext("Frequency", side = 2, line = 2.2, at = 0.5     , outer = TRUE, cex = 0.7)
    mtext("Frequency", side = 2, line = 2.2, at = 0.333/2 , outer = TRUE, cex = 0.7)
  }
  
  
  mtext("pba. cis", side = 2, line = 1.5, at = 1.666/2-0.15, outer = TRUE, cex = 0.8, las=1) #las = 1 pour horizontal text
  mtext("pba. cis", side = 2, line = 1.5, at = 0.5    -0.15, outer = TRUE, cex = 0.8, las=1)
  mtext("pba. cis", side = 2, line = 1.5, at = 0.333/2-0.15, outer = TRUE, cex = 0.8, las=1)
  
  mtext("0", side = 2, line = 0.65, at = 1.666/2-0.129, outer = TRUE, cex = 1)#, col="red")
  mtext("0", side = 2, line = 0.65, at = 0.5    -0.129, outer = TRUE, cex = 1)#, col="red")
  mtext("0", side = 2, line = 0.65, at = 0.333/2-0.129, outer = TRUE, cex = 1)#, col="red")
  
}

dev.off() ; play()

} # LANCE ET SAVE PDF

##########

# if (plotcounter == 9){
#   # ADD LEGEND
#   par(family = "Courier")
#   remov=0.03
#   my_cex = 1.2
#   my_font=2
#   text(midpoints[11]+0.5, 0.4-remov*0, "         s < -0.1  " , adj=1, col=effect_class_color[1], cex=my_cex, font=my_font)
#   text(midpoints[11]+0.5, 0.4-remov*1, "  -0.1 < s < -0.01 " , adj=1, col=effect_class_color[2], cex=my_cex, font=my_font)
#   text(midpoints[11]+0.5, 0.4-remov*2, " -0.01 < s < -0.001" , adj=1, col=effect_class_color[3], cex=my_cex, font=my_font)
#   text(midpoints[11]+0.5, 0.4-remov*3, "-0.001 < s <  0.001" , adj=1, col=effect_class_color[4], cex=my_cex, font=my_font)
#   text(midpoints[11]+0.5, 0.4-remov*4, " 0.001 < s <  0.01 " , adj=1, col=effect_class_color[5], cex=my_cex, font=my_font)
#   text(midpoints[11]+0.5, 0.4-remov*5, "  0.01 < s <  0.1  " , adj=1, col=effect_class_color[6], cex=my_cex, font=my_font)
#   text(midpoints[11]+0.5, 0.4-remov*6, "   0.1 < s         " , adj=1, col=effect_class_color[7], cex=my_cex, font=my_font)
#   par(family = "")
# }



