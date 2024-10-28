setwd("/Users/sylvain/Documents/DOCUMENTS/GeneRation/M2_PAPER")
source("../GeneRation_Fun.R")
source("R_SCRIPTS/Comfort_FUN.R") # includes assign
# source("R_SCRIPTS/151.FILM_FUN.R") # includes WNetwork3
# source("R_SCRIPTS/126.New_Mutate_in_new_Env_2.R") # IMPORTANT #includes cnames
# source("R_SCRIPTS/IndivRecap_fromSET_FUN.R")
library(RColorBrewer)
###############################################

####################################################################
## LOAD
####################################################################

# Load just one quickly:
allmutations = read.csv("Extracted_data/allmutations_20240417_fig06/allmutations001.csv")[,-1] ; play(1)
nrow(allmutations)

# LOAD ALL ####################################################################
setwd("/Users/sylvain/Documents/DOCUMENTS/GeneRation/M2_PAPER")
setwd("Extracted_data/allmutations_20240417_fig06/")
allmut_files = list.files(pattern = "allmutations....csv")
#allmut_files = list.files(pattern = "allmutations[0-9]{3}\\.csv")
n_batches = 10 #length(allmut_files)
n_batches = length(allmut_files)
all_batches = c()

for (batch_x_file in allmut_files[1:n_batches]){
  batch_x_name = substr(batch_x_file, start = 1, stop = 14)
  all_batches = c(all_batches, batch_x_name)
  assign(batch_x_name, read.csv(batch_x_file)[,-1])
  cat(paste0(sprintf("%02d", match(batch_x_file, allmut_files)), "/", 
             sprintf("%02d", length(allmut_files[1:n_batches])), " - ", batch_x_file, "\n"))
} ; play()
setwd("/Users/sylvain/Documents/DOCUMENTS/GeneRation/M2_PAPER")

allmutations = do.call(rbind, lapply(all_batches, function(x) get(x))) ; play() # combine all the allmuations
rm(list = all_batches) #remove individual batches
nrow(allmutations)


####################################################################
## PRE-TREATMENT
##################################################

# {# Rassemble DUP et DEL for easy ploting
#   allmutations[allmutations$mut_type=="DEL",]$dupdel_number = -allmutations[allmutations$mut_type=="DEL",]$dupdel_number
#   allmutations[allmutations$mut_type=="DEL",]$mut_type = "DUP"
# }

{# Define different classes for effect on fitness.
  effect_class_tresh = rev(c(0.1, 0.01, 0.001, -0.001, -0.01, -0.1))
  effect_class  = rev(c("BENEF3", "BENEF2", "BENEF1", "NEUTR", "DELET1", "DELET2", "DELET3"))
  effect_class_color = rev(c(colorRampPalette(c("#066400","#989797", "#831D01"))(7)))
  allmutations$fit_effect <- assign_category(allmutations$fit_delta)
}


# REMOVE NO MUTATIONS # Needed parce qu'on veut pas biaiser nos data (?) - on a pas besoin de la barre à zero
allmutations = allmutations[allmutations$notes != "NO MUTATION",]
# Remove DUP DEL en fait on en a pas besoin
allmutations = allmutations[allmutations$mut_type != "DUP",]
allmutations = allmutations[allmutations$mut_type != "DEL",]


################################################################################################
##################################################################################################
######                                                    ##########################################
#####                        PLOT                          ###########################################
######                                                    ##########################################
##################################################################################################
################################################################################################

effect_class
effect_class_color
effect_class_color2 = c(t(outer(effect_class_color, c("",""), paste, sep="")))
mut_types_test = c("REG", "COD") #, "DUP", "DEL")
pleio_levels = 0:10
cex_mtext = 1
effect_class2 = effect_class
cex_mtext = 1.2
effect_class2 = c("---", "--", "-", "0", "+", "++", "+++")

#####################

{ # FULL FIGURE PDF SAVE
#resetplot() # also resets cex_mtext to 1 (default)

chosentopologies = "SCALEF"
#chosentopologies = c("HIGHCO", "RANDOM", "SCALEF")

plot5 = F  ; if(plot5){p5=1}else{p5=0}
plotsize = F ; if(plotsize){psize=1}else{psize=0}
draw_box = F
titles = F

nb_topo = length(chosentopologies)
nb_plot = 3+p5+psize

width = 11*nb_topo # unit is cm
height =  6*nb_plot+1 # unit is cm
pdf(file = paste0("Rplot", format(Sys.time(), "%Y%m%d%H%M%S"), "_", width, "x", height, ".pdf"),
    width = width/cmtoin, height = height/cmtoin)

{ # LAUNCH FIGURE GENERATION ###################################################
  par(mfcol=c(nb_plot,nb_topo))
# par(mfrow=c(3,4))
plotcounter = 0
cex_mtext2 = 0.68
par(oma = c(0,1,2,0))
par(mgp=c(2.1, 1, 0))
par(mar=c(3,4,1.7,1))
old_mar = par("mar")
old_mgp = par("mgp")
line1 = 0.4
##########

for (topology in chosentopologies){ # start for (topology in ...)

pool_topology = allmutations[allmutations$topology==topology,]
#Check
nrow(pool_topology[pool_topology$mut_type=="REG",])
nrow(pool_topology[pool_topology$mut_type=="COD",])

# COLORS
regcod_colors = c("#FFE9B3", "#B3E4FF") #yello et blue
regcod_colors  = c(adjustcolor("#FF9830", 0.3), adjustcolor("grey", 0.5)) # orange et gris
regcod_colors2 = c(adjustcolor("#FF9830", 0.8), adjustcolor("grey", 0.7)) # orange et gris
regcod_colors  = rev(c(adjustcolor("#2B799E", 0.4), adjustcolor("grey", 0.5))) # rev(teal et gris
regcod_colors2 = rev(c(adjustcolor("#267FA8", 0.8), adjustcolor("grey", 0.7))) # rev(teal et gris

##############################################################
##############################################################
# FIRST TAKE OUT DATA # POUR PLOT 1 ET PLOT 5
#Setup empty dataframe for results
bilan_pool_regcod = data.frame(matrix(ncol = length(effect_class)+1, nrow = length(mut_types_test)))
colnames(bilan_pool_regcod) = c("ALL", effect_class)
rownames(bilan_pool_regcod) = mut_types_test

# FILL COLUMN FOR ALL 
effect_cat = "ALL"
for (mut_type in mut_types_test){
  bilan_pool_regcod[mut_type,effect_cat] = nrow(pool_topology[pool_topology$mut_type==mut_type,])
}
# FILL OTHER COLUMNS
for (effect_cat in effect_class){
  for (mut_type in mut_types_test){
    bilan_pool_regcod[mut_type,effect_cat] = nrow(pool_topology[pool_topology$fit_effect==effect_cat & pool_topology$mut_type==mut_type,])
  }
}
# pour PLOT 5
bilan_pool_regcod_distrib = bilan_pool_regcod
# pour avoir un pourcentage:
column_sums = colSums(bilan_pool_regcod)
bilan_pool_regcod = sweep(bilan_pool_regcod, 2, column_sums, "/") # column operation
bilan_pool_regcod


##############################################################
####### PLOT 5 - Density distribution
##############################################################
par(mar = old_mar+c(0.5,0,0,0.2))
par(mgp = old_mgp)

if (plot5){

bilan_pool_regcod_distrib
my_distrib = as.vector(unlist(bilan_pool_regcod_distrib[,-1] , use.names = FALSE))
my_distrib = my_distrib/sum(my_distrib)
# Let's add a fake bar to push things on the side.
my_distrib = c(mean(my_distrib), my_distrib)
my_barplot = barplot(my_distrib, 
                     col=c("white",rep(regcod_colors2,7)),
                     border = F, #c(F,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T),
                     space = c(0, 2.5, rep(c(0.05, 0.5), 10))[1:15],
                     ylab = "Density", axes = F)
axis(2, at = c(0,0.05, 0.1, 0.15, 0.2), labels = c(0, 0.05, 0.1, 0.15, 0.2))
mids = colSums(matrix(my_barplot[-1], nrow = 2, ncol = 7, byrow = F))/2

#axis(1, at = my_barplot[-1], labels = rep(c("reg", "cod"), 7), las=3, tick = F, line = 0-0.5)

fit_labels = c("-0.1>", "-0.01>", "-0.001>", expression(Delta * w), ">0.001", ">0.01", ">0.1")
for (i in 1:8){
  mtext(fit_labels[i] , side=1, line=line1, at = mids[i], font = 2, cex = cex_mtext2, col = effect_class_color[i])
}

plotcounter = plotcounter + 1 ; draw_box_check()
mtext(paste0(LETTERS[plotcounter], "."), side = 3, line = 0.5, at = 0 , outer = FALSE, cex = 1, font =2)

# box()
# abline(v=c((my_barplot[3]+my_barplot[4])/2, 
#            (my_barplot[5]+my_barplot[6])/2,
#            (my_barplot[7]+my_barplot[8])/2,
#            (my_barplot[9]+my_barplot[10])/2,
#            (my_barplot[11]+my_barplot[12])/2,
#            (my_barplot[13]+my_barplot[14])/2))
abline(h=0, col="grey")

text(x = 0, y = 0.20, "Regulatory", col=make_darker_color(regcod_colors2[1]), cex=1.2, adj=0)
text(x = 0, y = 0.175, "Coding",     col=make_darker_color(regcod_colors2[2]), cex=1.2, adj=0)

line2 = 1.8
mtext("Deleterious", side=1, line=line2, at = mids[2], font = 2, cex = 0.9, col = effect_class_color[2])
mtext("Neutral",     side=1, line=line2, at = mids[4], font = 2, cex = 0.9, col = effect_class_color[4])
mtext("Beneficial",  side=1, line=line2, at = mids[6], font = 2, cex = 0.9, col = effect_class_color[6])

} # if plot5

##############################################################
####### PLOT 1 - REG / COD
##############################################################
par(mar = old_mar) ; par(mgp = old_mgp)

if (!plot5){par(mar = old_mar+c(0.5,0,0,0.2))}

mut_types_test

main = ""

spaces = c(0, 0.35 ,0,0,0,0,0,0)
if (titles){main = paste0("Regulatory & Coding: ", topology)}
mids = barplot(as.matrix(bilan_pool_regcod), col = regcod_colors, space = spaces, border = "white", 
        main = main, ylab="Frequency", names.arg = rep("", 8))
abline(h=.5, lty=2, col="#6F6F6F")

plotcounter = plotcounter + 1 ; draw_box_check()
mtext(paste0(LETTERS[plotcounter], "."), side = 3, line = 0.5, at = 0 , outer = FALSE, cex = 1, font =2)


text(x = mids[1], y = .25, "Regulatory", cex = 1, col = make_darker_color(regcod_colors[1], 80), srt = 90)
text(x = mids[1], y = .75, "Coding"    , cex = 1, col = make_darker_color(regcod_colors[2], 80), srt = 90)

mtext("ALL"                , side=1, line=line1, at = mids[1], font = 1, cex = cex_mtext2)
mtext("-0.1>"              , side=1, line=line1, at = mids[2], font = 1, cex = cex_mtext2, col = effect_class_color[1])
mtext("-0.01>"             , side=1, line=line1, at = mids[3], font = 1, cex = cex_mtext2, col = effect_class_color[2])
mtext("-0.001>"            , side=1, line=line1, at = mids[4], font = 1, cex = cex_mtext2, col = effect_class_color[3])
mtext(expression(Delta * w), side=1, line=line1, at = mids[5], font = 1, cex = cex_mtext2, col = effect_class_color[4])
mtext(">0.001"             , side=1, line=line1, at = mids[6], font = 1, cex = cex_mtext2, col = effect_class_color[5])
mtext(">0.01"              , side=1, line=line1, at = mids[7], font = 1, cex = cex_mtext2, col = effect_class_color[6])
mtext(">0.1"               , side=1, line=line1, at = mids[8], font = 1, cex = cex_mtext2, col = effect_class_color[7])

if (!plot5){
  line2 = 1.8
  mtext("Deleterious", side=1, line=line2, at = mids[3], font = 2, cex = 0.9, col = effect_class_color[2])
  mtext("Neutral",     side=1, line=line2, at = mids[5], font = 2, cex = 0.9, col = effect_class_color[4])
  mtext("Beneficial",  side=1, line=line2, at = mids[7], font = 2, cex = 0.9, col = effect_class_color[6])
}

##############################################################
####### PLOT 2 - PLEIO
##############################################################
par(mar = old_mar) ; par(mgp = old_mgp)
pleio_levels
pleio_color_main = "#40015A" # violet un peu pétant
pleio_color_main = "#533774" # violet plus doux
pleio_colors = colorRampPalette(c("white", pleio_color_main))(length(pleio_levels))
#Setup empty dataframe for results
bilan_pool_pleio = data.frame(matrix(ncol = length(effect_class)+1, nrow = length(pleio_levels)))
colnames(bilan_pool_pleio) = c("ALL", effect_class)
rownames(bilan_pool_pleio) = pleio_levels

# FILL COLUMN FOR ALL 
effect_cat = "ALL"
for (pleio in pleio_levels){
  bilan_pool_pleio[pleio+1,effect_cat] = nrow(pool_topology[pool_topology$pleiotropy==pleio,])
}
# pour avoir un pourcentage:
bilan_pool_pleio[,effect_cat] = bilan_pool_pleio[,effect_cat] / sum(bilan_pool_pleio[,effect_cat])


for (effect_cat in effect_class){
  for (pleio in pleio_levels){
    bilan_pool_pleio[pleio+1,effect_cat] = nrow(pool_topology[pool_topology$fit_effect==effect_cat & pool_topology$pleiotropy==pleio,])
  }
  # pour avoir un pourcentage:
  bilan_pool_pleio[,effect_cat] = bilan_pool_pleio[,effect_cat] / sum(bilan_pool_pleio[,effect_cat])
}
bilan_pool_pleio ; "" ; colSums(bilan_pool_pleio)

spaces = c(0, 0.35 ,0,0,0,0,0,0)
if (titles){main = paste0("Pleiotropy: ", topology)}
mids = barplot(as.matrix(bilan_pool_pleio), col = pleio_colors, space = spaces, border = "white", 
        main = main, ylab="Frequency", names.arg = rep("", 8))

plotcounter = plotcounter + 1 ; draw_box_check()
mtext(paste0(LETTERS[plotcounter], "."), side = 3, line = 0.5, at = 0 , outer = FALSE, cex = 1, font =2)

cumsum_vector = cumsum(c(0, bilan_pool_pleio$ALL))
intermediate_values <- cumsum_vector[-length(cumsum_vector)] + diff(cumsum_vector) / 2
pleio_colors_bis = rev(pleio_colors) ; pleio_colors_bis[1:6] = pleio_colors[10]; pleio_colors_bis[7:10] = pleio_colors[1]

for (i in 1:11){
  text(x = mids[1], y = intermediate_values[i], labels = c(0:10)[i], cex = 1, col = pleio_colors_bis[i], srt=0)
}

# mtext("ALL", side=1, line=1, at = (1-0.5), font = 2, cex = cex_mtext2)
# for (i in 1:7){
#   mtext(effect_class2[i], side=1, line=1, at = mids[i+1], col = effect_class_color[i], font = 2, cex = cex_mtext)
# }

# mtext("Deleterious", side=1, line=line2, at = mids[3], font = 2, cex = cex_mtext2, col = effect_class_color[2])
# mtext("Neutral",     side=1, line=line2, at = mids[5], font = 2, cex = cex_mtext2, col = effect_class_color[4])
# mtext("Beneficial",  side=1, line=line2, at = mids[7], font = 2, cex = cex_mtext2, col = effect_class_color[6])

mtext("ALL"                , side=1, line=line1, at = mids[1], font = 1, cex = cex_mtext2)
mtext("-0.1>"              , side=1, line=line1, at = mids[2], font = 1, cex = cex_mtext2, col = effect_class_color[1])
mtext("-0.01>"             , side=1, line=line1, at = mids[3], font = 1, cex = cex_mtext2, col = effect_class_color[2])
mtext("-0.001>"            , side=1, line=line1, at = mids[4], font = 1, cex = cex_mtext2, col = effect_class_color[3])
mtext(expression(Delta * w), side=1, line=line1, at = mids[5], font = 1, cex = cex_mtext2, col = effect_class_color[4])
mtext(">0.001"             , side=1, line=line1, at = mids[6], font = 1, cex = cex_mtext2, col = effect_class_color[5])
mtext(">0.01"              , side=1, line=line1, at = mids[7], font = 1, cex = cex_mtext2, col = effect_class_color[6])
mtext(">0.1"               , side=1, line=line1, at = mids[8], font = 1, cex = cex_mtext2, col = effect_class_color[7])



##############################################################
####### PLOT 3 - cis/trans
##############################################################
par(mar = old_mar) ; par(mgp = old_mgp)
# Make cis/trans computation easier
pool_topology$CIS=pool_topology$pleiotropy*pool_topology$Expr_effect_cisness
cis_color = adjustcolor("#EDD405", 0.6) # en fait faudrait que ce soit le même couleur que dans la figure 4bis.
#Setup empty dataframe for results
bilan_pool_cis = data.frame(matrix(ncol = length(effect_class)+1, nrow = 1))
colnames(bilan_pool_cis) = c("ALL", effect_class)
rownames(bilan_pool_cis) = "pba_CIS"

# Faut retirer les pleio=0 ici
pool_cistrans = pool_topology[pool_topology$pleiotropy != 0,]

# FILL COLUMN FOR ALL 
effect_cat = "ALL"
bilan_pool_cis[1, effect_cat] = mean(pool_cistrans[,]$CIS)


for (effect_cat in effect_class){
  bilan_pool_cis[1, effect_cat] = mean(pool_cistrans[pool_cistrans$fit_effect==effect_cat,]$CIS)
}
bilan_pool_cis

spaces = c(0, 0.35 ,0,0,0,0,0,0)
if (titles){main = paste0("Max. cis-effect: ", topology)}
mids = barplot(as.matrix(bilan_pool_cis), col = cis_color, space = spaces, border = "white", ylim=c(0,1), 
        main = main, ylab = "Probability", names.arg = rep("", 8))

abline(h=bilan_pool_cis$ALL, lty = 2, col="grey")
plotcounter = plotcounter + 1 ; draw_box_check()
mtext(paste0(LETTERS[plotcounter], "."), side = 3, line = 0.5, at = 0 , outer = FALSE, cex = 1, font =2)

# text(x = mids[1], y = bilan_pool_cis$ALL/2+0.04, labels = "probability", cex = 1, col = "black")
# text(x = mids[1], y = bilan_pool_cis$ALL/2-0.04, labels = "cis-effect" , cex = 1, col = "black")

text(x = mids[1], y = bilan_pool_cis$ALL/2, labels = "probability cis-effect" , cex = 1, 
     col = make_darker_color(cis_color, 100), srt=90)


# mtext("ALL", side=1, line=1, at = (1-0.5), font = 2, cex = cex_mtext2)
# for (i in 1:7){
#   mtext(effect_class2[i], side=1, line=1, at = mids[i+1], col = effect_class_color[i], font = 2, cex = cex_mtext)
# }

# mtext("Deleterious", side=1, line=line2, at = mids[3], font = 2, cex = cex_mtext2, col = effect_class_color[2])
# mtext("Neutral",     side=1, line=line2, at = mids[5], font = 2, cex = cex_mtext2, col = effect_class_color[4])
# mtext("Beneficial",  side=1, line=line2, at = mids[7], font = 2, cex = cex_mtext2, col = effect_class_color[6])

mtext("ALL"                , side=1, line=line1, at = mids[1], font = 1, cex = cex_mtext2)
mtext("-0.1>"              , side=1, line=line1, at = mids[2], font = 1, cex = cex_mtext2, col = effect_class_color[1])
mtext("-0.01>"             , side=1, line=line1, at = mids[3], font = 1, cex = cex_mtext2, col = effect_class_color[2])
mtext("-0.001>"            , side=1, line=line1, at = mids[4], font = 1, cex = cex_mtext2, col = effect_class_color[3])
mtext(expression(Delta * w), side=1, line=line1, at = mids[5], font = 1, cex = cex_mtext2, col = effect_class_color[4])
mtext(">0.001"             , side=1, line=line1, at = mids[6], font = 1, cex = cex_mtext2, col = effect_class_color[5])
mtext(">0.01"              , side=1, line=line1, at = mids[7], font = 1, cex = cex_mtext2, col = effect_class_color[6])
mtext(">0.1"               , side=1, line=line1, at = mids[8], font = 1, cex = cex_mtext2, col = effect_class_color[7])



##############################################################
####### PLOT 4 - Mutation size
##############################################################
par(mar = old_mar) ; par(mgp = old_mgp)
effect_class
mut_types_test
new_mar = old_mar+c(0,1,0,1)
par(mar = new_mar)
new_mgp = old_mgp+c(1,1,1)
par(mgp = new_mgp)

if(plotsize){
  
library(vioplot) #install.packages("vioplot", type = "binary")
# vioplot(rnorm(100, mean = 5), rnorm(100, mean = 10), names = c("Group 1", "Group 2"), col = c("gold", "skyblue"), main = "Violin Plot of Two Groups")
#Setup empty dataframe for results
bilan_pool_size = data.frame(matrix(ncol = (length(effect_class)+1)*2, nrow = nrow(pool_topology)/2))
colnames(bilan_pool_size) = c(t(outer(c("ALL", effect_class), mut_types_test, paste, sep="-")))

# FILL COLUMN FOR ALL 
effect_cat = "ALL"
for (mut_type in mut_types_test){
  the_values = pool_topology[pool_topology$mut_type==mut_type,]$mutvalue_delta
  bilan_pool_size[[paste0(effect_cat, "-", mut_type)]] = c(the_values, rep(NA, nrow(pool_topology)/2 - length(the_values)))
}

for (effect_cat in effect_class){
  for (mut_type in mut_types_test){
    the_values = pool_topology[pool_topology$fit_effect==effect_cat & pool_topology$mut_type==mut_type,]$mutvalue_delta
    bilan_pool_size[[paste0(effect_cat, "-", mut_type)]] = c(the_values, rep(NA, nrow(pool_topology)/2 - length(the_values)))
  }
}
#boxplot(abs(bilan_pool_size), ylim=c(0,1))
#vioplot(abs(bilan_pool_size), ylim=c(0,1))

###### PARAMETERS FOR PLOT: ###
absolute_values = T
condensed = T #if TRUE, launch 4B.
central_point = "mean" # or "median
space_val = 0.4

if (absolute_values){ ylim = c(-0.02, 1.02) ; bilan_pool_size = abs(bilan_pool_size)
} else { ylim = c(-1.02, 1.02) }

increments = c(0.3, 0, 0.4+0.5, rep(c(0, space_val), 8))[1:16]
at_specified = 1:16 + cumsum(increments)
xlim = c(1,max(at_specified))

if (titles){main = paste0("Mutation size (REG/COD): ", topology)}
plot(NA, NA, ylim=ylim, xlim=xlim, xaxt = "n", xlab="", 
     ylab="Mutation size (absolute value)", main=main, bty="n")
box(col = "#8F8F8F", lwd=.5) #grey
abline(h=.55, col='grey', lty = 2)
if (!absolute_values){abline(h=c(0, -0.5), col='grey', lty = 2)}

###### 4.A - (NOT CONDENSED was here before)
###### 4.B - 
par(bty = 'n')
# vioplot(bilan_pool_size, ylim=c(-1,1), add=T, at = at_specified, border=NA,
#         col=adjustcolor("white", 0), colMed = adjustcolor(effect_class_color2, 0))
boxplot(bilan_pool_size, add=T, at = at_specified, col=regcod_colors2, 
        names = rep("", 16), axes=F, lty=1, staplelty = 1, boxwex=.3, 
        whiskcol="white", staplecol = "white")
par(bty = 'o') # restore default


# # Point central
# reduce_size=0.5
# if (central_point == "median") {central_point_values = apply(bilan_pool_size, 2, median, na.rm = TRUE)}
# if (central_point == "mean")   {central_point_values = apply(bilan_pool_size, 2, mean,   na.rm = TRUE)}
# points(at_specified, central_point_values, col = "white", pch = 20, cex=4)
# points(at_specified, central_point_values, cex=3-reduce_size)
# #points(at_specified, central_point_values, col = effect_class_color2, pch = 20, cex=3)
# points(at_specified, central_point_values, col = regcod_colors2, pch = 20, cex=3-reduce_size)
# points(at_specified, central_point_values, cex=2-reduce_size)

###### AXES etc.
#axis(1, at = at_specified, labels = rep(c("reg", "cod"), 8), las = 2)
text(x = at_specified+0.1, y = 0-0.02, labels = c("reg", "cod"), srt=90, pos=3)
at_midpoints <- (at_specified[seq(1, length(at_specified) - 1, by = 2)] + at_specified[seq(2, length(at_specified), by = 2)]) / 2
#axis(1, at = at_midpoints, labels = paste0("-", c("ALL",effect_class), "-"), tick = F, line = 1, col = c("black",effect_class_color))


at_lines <- (at_specified[seq(2, length(at_specified) - 1, by = 2)] + at_specified[seq(3, length(at_specified), by = 2)]) / 2
abline(v=at_lines[2:7], col="grey")
abline(v=at_lines[1]-0.2, col="grey")
abline(v=at_lines[1]+0.2, col="grey")

plotcounter = plotcounter + 1 ; draw_box_check()
mtext(paste0(LETTERS[plotcounter], "."), side = 3, line = 0.5, at = 0 , outer = FALSE, cex = 1, font =2)

# mtext("ALL", side=1, line=1, at = at_midpoints[1], font = 2, cex = cex_mtext2)
# for (i in 1:7){
#   mtext(effect_class2[i], side=1, line=0.75, at = at_midpoints[i+1], col = effect_class_color[i], font = 2, cex = cex_mtext)
#   # line était 2 avant
#   #mtext("BENEF", side=1, line=1, at = 0.5, col="red")
# }

mids = at_midpoints

mtext("ALL"                , side=1, line=line1, at = mids[1], font = 1, cex = cex_mtext2)
mtext("-0.1>"              , side=1, line=line1, at = mids[2], font = 1, cex = cex_mtext2, col = effect_class_color[1])
mtext("-0.01>"             , side=1, line=line1, at = mids[3], font = 1, cex = cex_mtext2, col = effect_class_color[2])
mtext("-0.001>"            , side=1, line=line1, at = mids[4], font = 1, cex = cex_mtext2, col = effect_class_color[3])
mtext(expression(Delta * w), side=1, line=line1, at = mids[5], font = 1, cex = cex_mtext2, col = effect_class_color[4])
mtext(">0.001"             , side=1, line=line1, at = mids[6], font = 1, cex = cex_mtext2, col = effect_class_color[5])
mtext(">0.01"              , side=1, line=line1, at = mids[7], font = 1, cex = cex_mtext2, col = effect_class_color[6])
mtext(">0.1"               , side=1, line=line1, at = mids[8], font = 1, cex = cex_mtext2, col = effect_class_color[7])

} # end if plotsize

##############################################################
####### FINAL STUFF
##############################################################

# Y-titles
sep = 0.25/2
midpoints = seq(0,1, length.out = nb_plot+1)[1:(nb_plot)]
midpoints = midpoints+((midpoints[2]-midpoints[1])/2)
i = 1
if(plot5){   mtext("Distribution",   side = 2, line = -0.65, at = rev(midpoints)[i]+0.009, outer = T, font = 2) ; i=i+1}
             mtext("Mutation type",  side = 2, line = -0.65, at = rev(midpoints)[i]+0.009, outer = T, font = 2) ; i=i+1
             mtext("Pleiotropy",     side = 2, line = -0.65, at = rev(midpoints)[i]+0.009, outer = T, font = 2) ; i=i+1
             mtext("Cis-effect",     side = 2, line = -0.65, at = rev(midpoints)[i]+0.009, outer = T, font = 2) ; i=i+1
if(plotsize){mtext("Mutation size",  side = 2, line = -0.65, at = rev(midpoints)[i]+0.009, outer = T, font = 2) ; i=i+1}

# X-titles (topologies)
chosentopologies
topologies_fullnames = c("Highly Connected", "Random", "Scale-Free") 
tot = length(chosentopologies)
evenly_spacing <- function(n){return((2 *(1:n)-1)/(2*n))}
at_coord = evenly_spacing(tot)
for (i in 1:tot){
  if (tot==1){addnet=" Networks"}else{addnet=""}
  mtext(paste0(topologies_fullnames[match(chosentopologies[i], topologies)], addnet), 
        side = 3, line = 0.7, at = at_coord[i], outer = T, font = 2, 
        col = topocolors[match(chosentopologies[i], topologies)]) 
}


##############

} # for (topology in ...)
} # ALL FIGURE GENERATION

 dev.off() ; play()

} # FULL FIGURE PDF SAVE
