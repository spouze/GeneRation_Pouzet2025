# COPY SCRIPT 090 pour mettre le truc en mode vertical
{# Import basic tools
source(here::here("GeneRation_Fun_v1.R"))
source(here::here("Data_Analysis_scripts/Comfort_FUN.R"))
source(here::here("Data_Analysis_scripts/IndivFinale_FUN.R"))
source(here::here("Data_Analysis_scripts/IndivRecap_fromSET_FUN.R"))
source(here::here("Data_Analysis_scripts/128.New_Mutate_in_new_Env_4.R"))
# source(here::here("Data_Analysis_scripts/097.WTab2_WTab_modified.R"))
library(MASS)
}
#####################################################################
{######################### LOAD DATA #################################
#####################################################################

### LOAD FITS
dream_Fits = read.csv("Extracted_data/dream_SIMUSET_Fits.csv")[,-1] ; play(1) #takes 20 seconds
batch_november_Fits = read.csv("Extracted_data/2.storeFits/storeFits_november_100+100+100_burnin")[,-1][,1:1002]
# batch_november_Fits = read.csv("Extracted_data/Batch_november_Fits.csv")[,-1] # ou SCRIPT 032 lol
generations = ncol(dream_Fits)-2
# SPLIT IN 3
for (topo in topologies) {
  assign(paste0("dream_Fits_", topo), dream_Fits[which(dream_Fits[,"topo"]==topo),])
}
nrow(dream_Fits_RANDOM) ; nrow(dream_Fits_SCALEF) ; nrow(dream_Fits_HIGHCO)

### LOAD LASTS
dream_Lasts = read.csv("Extracted_data/dream_SIMUSET_Lasts.csv") ; play(1) #takes 3 seconds
dream_Lasts = dream_Lasts[,-1] # remove first column which is just numbers
# View(dream_Lasts)

### Load degree info
#SUPER_Wdegree_df = read.csv("Extracted_data/Figures_vrac/Fig02_Degree/Degree_count_dreamWClean0.01_integreatedCvec.csv")
SUPER_Wdegree_df = read.csv("Figures/Figures_vrac/Fig02_Degree/Degree_count_dreamWClean0.01_integreatedCvec.csv")

} ; play() # load data

# PARAMETERS # Go at it sequentially:
batch="november"
simus = sample(1:1000, 100)
candidates = c("0612_215808_HIGHCO_012_november",
               "0613_220243_RANDOM_099_november", 
               "0614_083747_SCALEF_077_november")
# "0613_091244_RANDOM_002_november"
# "0614_052534_SCALEF_054_november

#####################################################################
######################### FULL FIGURE ###############################
#####################################################################
# resetplot()
# Check for Unclosed Devices
# FIX Error in mtext(paste0(LETTERS[letter_order[plotcounter]], "."), side = 3, :plot.new has not been called yet
if (dev.cur() != 1) dev.off()

{ # LAUNCH WHOLE FIGURE

#Save the figure in 20 x 11 (inches) avant
size_factor = 1.3
width = 17*size_factor # unit is cm
height = 20*size_factor # unit is cm
pdf(file = paste0("Rplot", format(Sys.time(), "%Y%m%d%H%M%S"), "_", width, "x", height, ".pdf"), 
    width = width/cmtoin, height = height/cmtoin)

{# LAYOUT VERTICAL
plotcounter = 0
par(mfrow=c(5,3)) # 3 rows, 3 columns
par(mfrow=c(5,3), oma = c(3, 3, 4, 0)) # 5 rows, 3 columns
layout_order = c(1,5,9, 2,6,10, 3,7,11, 4,8,12, 13,14,15)
letter_order = c(1:4, 6:9, 11:14, 5,10,15)
layout_matrix = matrix(layout_order, nrow = 5, byrow = TRUE)
layout_heights = c(1.75, 1, 1.75, 1.75, 1)
layout(layout_matrix, heights = layout_heights)  # heights are ratios
#layout.show(15)
draw_box = F #for box("figure", col="forestgreen)
} # Layout  

# Obligé d'y aller à la main sequentiellement, sinon ca saute des étapes. 
topo = "HIGHCO" ; play(1)
topo = "RANDOM" ; play(1)
topo = "SCALEF" ; play(1)

for (topo in topologies){

{## 0. SELECT SIMU & 0b. LOAD SIMU
simu_id = candidates[match(topo, topologies)]
#simu = LoadSimu(paste0(list.files(pattern = batch)[1], "/SIMUS/", simu_id, "/"))
simu = LoadSimu(paste0("SIMULATIONS/", list.files(pattern = batch, path = "SIMULATIONS/")[1], "/SIMUS/", simu_id, "/"))
print("Completed Step 0")

## 1. INITIAL MATRIX
par(mar=c(0.5, 4, 3, 1) + 0.1)
WTab3(WFinale(gen = 1))
axis(3, at = 10, labels = "10", pos = 10.5+0.15, las = 0, tck = 0)
if (topo=="HIGHCO"){
  mtext(text = "from", side = 3, line = 1+0.18, at = -0.1, cex = 1/size_factor)
  mtext(text = "on",   side = 2, line = 2, at = 10+0.04, cex = 1/size_factor, las=1)}
} ; print("Completed Step 1")  ; play(1)
{#mtext("from", side = 3, line=1.3, at=0.2, cex=0.7) # A AJUSTER
plotcounter = plotcounter +1 ; draw_box_check()
mtext(paste0(LETTERS[letter_order[plotcounter]], "."), side = 3, line = 2.2, at = 0 , outer = FALSE, cex = 1, font =2)
} ; play(2)

{## 2. EVOLUTION HISTORY
par(mar=c(2.5,4,1.7,1)+0.1, mgp=c(1.5, 1, 0))
if (topo == "HIGHCO"){ylab="Pop. fitness"}else{ylab=""}
plot(NA, xlim=c(0,generations), ylim=c(0,1), xlab="Generations", ylab=ylab, main="", 
     bty="n", xaxt="n", yaxt="n")
x_labels = c("0", "", "", "", "", "1000")
axis(1, at = seq(0,1000,200), labels = x_labels)
axis(2, at = seq(0,1,0.2), labels = c("0", "", "", "", "", "1"), las=1)
# abline(v=1000, col="grey", lty="dashed") ; abline(h=.5, col="grey", lwd=.5)
colr = topocolors[match(topo, topologies)]
# SUCCESSFUL SIMUS
#dream_Fits_XXXXXX = get(paste0("dream_Fits_",topo))
#for (i in simus){lines(1:generations,dream_Fits_XXXXXX[i,c(3:1002)], col=adjustcolor(colr,.4))}
# RANDOM SIMUS
n_simus = 30 # number of simus to display
simuset = batch_november_Fits[batch_november_Fits$topo==topo,]
relous = which(simuset$gen1>0.1)
#print(paste0("Removed ", length(relous), " relous."))
simuset = simuset[-relous,] # on enleve les relous qui démarrent haut
for (i in 1:n_simus){lines(1:generations,simuset[i,c(3:1002)], col=adjustcolor(colr,.4))}
# MEAN RANDOM (?)
#lines(1:generations, colMeans(dream_Fits_XXXXXX[,c(3:1002)]), col="white", lwd=5)
#lines(1:generations, colMeans(dream_Fits_XXXXXX[,c(3:1002)]), col="black", lwd=3)
restore_mgp()
} ; print("Completed Step 2")  ; play(1)
{plotcounter = plotcounter +1 ; draw_box_check()
mtext(paste0(LETTERS[letter_order[plotcounter]], "."), side = 3, line = 0.3, at = -60 , outer = FALSE, cex = 1, font =2)
}  ; play(2)

{## 3. FINAL MATRIX
par(mar=c(0.5, 4, 3, 1) + 0.1)
WTab3(WFinale(gen = 1000))
axis(3, at = 10, labels = "10", pos = 10.5+0.15, las = 0, tck = 0)
} ; print("Completed Step 3")  ; play(1)
{#mtext("from", side = 3, line=1.3, at=0.2, cex=0.7) # A AJUSTER
plotcounter = plotcounter +1 ; draw_box_check()
mtext(paste0(LETTERS[letter_order[plotcounter]], "."), side = 3, line = 2.2, at = 0 , outer = FALSE, cex = 1, font =2)
} ; play(2)

{## 4. FINAL NETWORK
#source("../GeneRation_Fun.R")
par(mar=c(1, 2, 1.2, 0) + 0.1)
WNetwork2(WFinale(), vcex = 1.2, dist = 2.3, lab_col = "#979695")
}; print("Completed Step 4")  ; play(1)
{plotcounter = plotcounter +1 ; draw_box_check()
  mtext(paste0(LETTERS[letter_order[plotcounter]], "."), side = 3, line = -1, at = -1.3 , outer = FALSE, cex = 1, font =2)
} ; play(2) ; play(2)

} # for topo in topologies

############################################
{## 5. DRAW TOPOLOGY / DEGREE (single-step)
par(mar=c(2,4,0.3,1)+0.1, mgp=c(2.2, 0.5, 0))
black_line_lwd = 3

# --------------------------------------------------------------------------------------------
topo = "HIGHCO"
x = SUPER_Wdegree_df[which(SUPER_Wdegree_df$topo==topo),]$inout_add
# sum(x==15) # Check
x = c(x, 0:17) # pour que tous aient jusque 15
myhist = hist(x, prob = TRUE, 
              ylim=c(0,0.36), main = "",
              col=c("grey", rep(adjustcolor(topocolors[match(topo, topologies)], 0.7), 20)), 
              xaxt="n", yaxt="n", xlab="", ylab="",
              right=F) # needed otherwise aggregates 0 and 1 values
# axis(1, at=myhist$mids, labels = rep("", 17))
# axis(1, at=myhist$mids[8 ], labels = "7", tick = F)
# axis(1, at=myhist$mids[9 ], labels = "8", tick = F)
# axis(1, at=myhist$mids[10], labels = "9", tick = F)
axis(1, at=myhist$mids, labels = rep("", 17))
axis(1, at=myhist$mids[c(1,6,11,16)], labels = c(0,5,10,15), tick = F)

axis(2, at=c(0,0.1,0.2, 0.3), labels = c("0", "", "", "0.3"), las=1)
points(myhist$mids, myhist$density, pch=20)

fit <- fitdistr(x, "normal")
para <- fit$estimate
BIC(fit) # the model with the lowest BIC is preferred # 45 832.08
xx = c(seq(0,16, 0.05))
lines(xx+0.5, dnorm(xx, para[1], para[2]), 
      col = adjustcolor("black", 0.7), lwd = black_line_lwd, pch=20, type = "l")
plotcounter = plotcounter +1 ; draw_box_check()
mtext(paste0(LETTERS[letter_order[plotcounter]], "."), side = 3, line = 0.3, at = -1.3 , outer = FALSE, cex = 1, font =2)
mtext("Degree",  side = 1, line = 1.3, cex = 1/size_factor)
mtext("Density", side = 2, line = 1.3, cex = 1/size_factor)


# --------------------------------------------------------------------------------------------
topo = "RANDOM"
x = SUPER_Wdegree_df[which(SUPER_Wdegree_df$topo==topo),]$inout_add
x = c(x, 0:17) # pour que tous aient jusque 15
myhist = hist(x, prob = TRUE, 
              ylim=c(0,0.36), main = "",
              col=c("grey", rep(adjustcolor(topocolors[match(topo, topologies)], 0.7), 20)), 
              xaxt="n", yaxt="n", ylab="", xlab="",
              right=F) # needed otherwise aggregates 0 and 1 values
# axis(1, at=myhist$mids, labels = rep("", 17))
# axis(1, at=myhist$mids[3 ], labels = "2", tick = F)
# axis(1, at=myhist$mids[4 ], labels = "3", tick = F)
axis(2, at=c(0,0.1,0.2, 0.3), labels = c("0", "", "", "0.3"), las=1)
axis(1, at=myhist$mids, labels = rep("", 17))
axis(1, at=myhist$mids[c(1,6,11,16)], labels = c(0,5,10,15), tick = F)
points(myhist$mids, myhist$density, pch=20)
### ALL GOOD TILL THERE

fit <- fitdistr(x, "poisson")
para <- fit$estimate
xx = c(0:16)
lines(spline(xx+0.5, dpois(xx, para[1]), n = 300), # Need a smooth because poisson is discrete 
     col = adjustcolor("black", 0.7), lwd = black_line_lwd)
# lines(xx+0.5, dpois(c(0:16), para[1]),
#       col = "black", lwd = 2, pch=20, type = "o")
plotcounter = plotcounter +1 ; draw_box_check()
mtext(paste0(LETTERS[letter_order[plotcounter]], "."), side = 3, line = 0.3, at = -1.3 , outer = FALSE, cex = 1, font =2)
mtext("Degree", side = 1, line = 1.3, cex = 1/size_factor)


# --------------------------------------------------------------------------------------------
topo = "SCALEF"
x = SUPER_Wdegree_df[which(SUPER_Wdegree_df$topo==topo),]$inout_add
x = c(x, 0:17) # pour que tous aient jusque 15
myhist = hist(x, prob = TRUE, 
              ylim=c(0,0.36), main = "",
              col=c("grey", rep(adjustcolor(topocolors[match(topo, topologies)], 0.7), 20)), 
              xaxt="n", yaxt="n", xlab="", ylab="",
              right=F) # needed otherwise aggregates 0 and 1 values
axis(1, at=myhist$mids, labels = rep("", 17))
# par(col.axis="grey")
#   axis(1, at=myhist$mids[1], labels = "0", tick = F)
# par(col.axis="black")
#for (i in 1:8){axis(1, at=myhist$mids[i+1], labels = i, tick = F)}
axis(1, at=myhist$mids[c(1,6,11,16)], labels = c(0,5,10,15), tick = F)

axis(2, at=c(0,0.1,0.2, 0.3), labels = c("0", "", "", "0.3"), las=1)
points(myhist$mids, myhist$density, pch=20)

# PREVIOUS EXPONENTIAL FIT
# fit <- fitdistr(x, "exponential")
# para <- fit$estimate
# xx = c(seq(0,16, 0.05))
# lines(xx+0.5, dexp(xx, rate = para), 
#       col = adjustcolor("black", 0.7), lwd = black_line_lwd, pch=20, type = "l")
# BIC(fit) # the model with the lowest BIC is preferred # 32908.34
# polygon(x = c(-1,3,3,-1), y = c(0.375, 0.375, 0.5, 0.5), border = "white", col = "white") # degager le trop haut
# plotcounter = plotcounter +1 ; draw_box_check()
# mtext(paste0(LETTERS[letter_order[plotcounter]], "."), side = 3, line = 0.3, at = -1.3 , outer = FALSE, cex = 1, font =2)
# mtext("Degree", side = 1, line = 1.3, cex = 1/size_factor)

# NEW POWERLAW FIT
############# POWER LAW
xx = c(1:16)[1:7]
yy = my_hist$density[2:17][1:7]
fit = lm(log(yy) ~ log(xx))
b = fit$coef[2]
a = exp(fit$coef[1])
x_plot = seq(0.1,17, .1)
y_powerlaw = a * x_plot^b
hist_decalage = 0.5
lines(x_plot+hist_decalage, y_powerlaw,
      col = adjustcolor("black", 0.7), lwd = black_line_lwd)
plotcounter = plotcounter +1 ; draw_box_check()
mtext(paste0(LETTERS[letter_order[plotcounter]], "."), side = 3, line = 0.3, at = -1.3 , outer = FALSE, cex = 1, font =2)
mtext("Degree", side = 1, line = 1.3, cex = 1/size_factor)



} ; play(2) ## 5. DRAW TOPOLOGY / DEGREE (single-step)

############################################

{##6. ## AND ADD THE COLUMN NAMES ETC. 
repere = cumsum(rev(layout_heights))/sum(layout_heights)
midpoints = c(repere[1]/2, (repere[-1] + repere[-length(repere)])/2)

# ROW TITLES IN 2 PARTS
mtext('Initial',      side = 2, line = 0.7, at = midpoints[5], outer = T, font = 1)
mtext('Evolution',    side = 2, line = 0.7, at = midpoints[4], outer = T, font = 1)
mtext('Evolved',      side = 2, line = 0.7, at = midpoints[3], outer = T, font = 1)
mtext('Evolved',      side = 2, line = 0.7, at = midpoints[2], outer = T, font = 1)
mtext('Total Degree', side = 2, line = 0.7, at = midpoints[1], outer = T, font = 1)

mtext('Interaction Matrix', side = 2, line = -0.6, at = midpoints[5], outer = T, font = 1)
mtext('Pattern',            side = 2, line = -0.6, at = midpoints[4], outer = T, font = 1)
mtext('Interaction Matrix', side = 2, line = -0.6, at = midpoints[3], outer = T, font = 1)
mtext('Network',            side = 2, line = -0.6, at = midpoints[2], outer = T, font = 1)
mtext('Distribution',       side = 2, line = -0.6, at = midpoints[1], outer = T, font = 1)

#mtext('A', side = 1, line = 0.5, at = 1, outer = TRUE, col="red3", font = 2)

add_val = 0.03
mtext("\"full-active\"",    side = 3, line = 0.5, at = 0.333/2   +add_val, outer = TRUE, font = 2)
mtext("\"empty-active\"",   side = 3, line = 0.5, at = 0.5       +add_val, outer = TRUE, font = 2)
mtext("\"empty-inactive\"", side = 3, line = 0.5, at = 0.5+0.333 +add_val, outer = TRUE, font = 2)

mtext('Highly Connected', side = 1, line = 1, at = 0.333/2   +add_val, outer = TRUE, col="red3",        font = 2)
mtext("Random",           side = 1, line = 1, at = 0.5       +add_val, outer = TRUE, col="forestgreen", font = 2)
mtext('Scale-Free',       side = 1, line = 1, at = 0.5+0.333 +add_val, outer = TRUE, col="blue3",       font = 2)

#mtext('Initial conditions',  side = 3, line = 0.5, at = 0.5, outer = TRUE,       font = 2)

#mtext('Generations', side = 1, line = 0.5, at = 0.315, outer = TRUE, cex = 0.7, col="black")
#mtext('Edge number', side = 1, line = 0.5, at = 0.915, outer = TRUE, cex = 0.7, col="black")

} ; play(2) # Add text

#mtext("2",    side = 3, line = 0.5, at = 1/3/2+0.05,   outer = TRUE, font = 2)



dev.off() ; play()# for the end of the PDF


} # LAUNCH WHOLE FIGURE

