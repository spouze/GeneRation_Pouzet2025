# BASED ON SCRIPT 099

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Set working directory where the file is.
# getwd()
source("src/GeneRation_Fun_v1.R")
source("src/GeneRation_DataAnalysis_FUN.R")
require(MASS)

#####################################################################
{######################### LOAD DATA #################################
#####################################################################

### LOADING DATA FITS

# 1. PAPER DATA
# dream_Fits = read.csv("USED_DATA/Extracted_from_successful_simulations/dream_SIMUSET_Fits.csv")[,-1] ; play(1) #can take a while
# batch_november_Fits = read.csv("USED_DATA/Extracted_from_successful_simulations/storeFits_november_100+100+100_burnin")[,-1][,1:1002]
# dream_Lasts = read.csv("USED_DATA/Extracted_from_successful_simulations/dream_SIMUSET_Lasts.csv") ; play(1) #takes 3 seconds
# # THOSE ARE THE ONES THAT WILL BE PLOTED
# candidates = c("USED_DATA/SIMULATIONS_small_set/0612_215808_HIGHCO_012_november/",
#                "USED_DATA/SIMULATIONS_small_set/0613_220243_RANDOM_099_november/", 
#                "USED_DATA/SIMULATIONS_small_set/0614_083747_SCALEF_077_november/")
  
# 2. FEW SIMULATIONS
# dream_Fits = read.csv("USED_DATA/Extracted_data_small_set/storeFits_november")[,-1] #can take a while
# batch_november_Fits = dream_Fits
# dream_Lasts = read.csv("USED_DATA/Extracted_data_small_set/storeLasts_november") 
# candidates = c("USED_DATA/SIMULATIONS_small_set/0612_215808_HIGHCO_012_november/",
#                "0USED_DATA/SIMULATIONS_small_set/613_220243_RANDOM_099_november/", 
#                "USED_DATA/SIMULATIONS_small_set/0614_083747_SCALEF_077_november/")


# 3. TUTO_SIMUS
dream_Fits = read.csv("Extracted_data/storeFits_TUTO_SIMUS")[,-1] #can take a while
batch_november_Fits = dream_Fits
dream_Lasts = read.csv("Extracted_data/storeLasts_TUTO_SIMUS")
## HERE YOU NEED TO PUT THE NAMES OF 3 SIMULATIONS FROM THE SIMULATION FOLDER (HANDPICKED)
# THOSE ARE LOCATIONS, include the "/" at the end.
candidates = c("SIMULATIONS/0625_114103_HIGHCO_001_TUTO_SIMUS/", # HIGHCO
               "SIMULATIONS/0625_114117_RANDOM_002_TUTO_SIMUS/",  # RANDOM
               "SIMULATIONS/0625_114128_SCALEF_003_TUTO_SIMUS/") # SCALEF

################################################################
################################################################

generations = ncol(dream_Fits)-2
# SPLIT IN 3
for (topo in topologies) {
  assign(paste0("dream_Fits_", topo), dream_Fits[which(dream_Fits[,"topo"]==topo),])
}
nrow(dream_Fits_RANDOM) ; nrow(dream_Fits_SCALEF) ; nrow(dream_Fits_HIGHCO)
dream_Lasts = dream_Lasts[,-1] # remove first column which is just numbers
# View(dream_Lasts)
simus = c(1:100)#sample(1:1000, 100)


################################################################
################################################################
################################################################
## PLOT FIGURE 2 (E,J,O) - Network distributions
################################################################
################################################################
################################################################
# BASED ON SCRIPT 060

dream_SET = dream_Lasts

#Create master dataframe
SUPER_Wdegree_df = data.frame(matrix(NA, ncol = 8, nrow = 10*nrow(dream_SET)))
colnames(SUPER_Wdegree_df)=c("line", "id", "topo", "gene", 
                             "indegree", "outdregree", "inout_cross", "inout_add")
for (line in 1:nrow(dream_SET)){
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
    # Remplit le master dataframe
    SUPER_Wdegree_df[(((line-1)*10)+gene),] = c(line, id, topo, 
                                                gene, 
                                                degree_in, degree_out, degree_inout, 
                                                (degree_in + degree_out))
  } ; cat(paste0(line,".")) #just for notice
}


} ; play() # load data




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
  pdf(file = paste0("figures/Rplot", format(Sys.time(), "%Y%m%d%H%M%S"), "_fig2_", width, "x", height, ".pdf"), 
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
      simu = LoadSimu(paste0(simu_id))
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
      if (length(relous) > 0) {simuset = simuset[-relous,]} # on enleve les relous qui démarrent haut
      for (i in 1:n_simus){lines(1:generations,simuset[i,c(3:(generations+2))], col=adjustcolor(colr,.4))}
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
      WTab3(WFinale(gen = generations))
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
    require(MASS)
    par(mar=c(2,4,0.3,1)+0.1, mgp=c(2.2, 0.5, 0))
    black_line_lwd = 3
    
    # --------------------------------------------------------------------------------------------
    topo = "HIGHCO"
    x = SUPER_Wdegree_df[which(SUPER_Wdegree_df$topo==topo),]$inout_add
    # sum(x==15) # Check
    x = as.numeric(c(x, 0:17)) # pour que tous aient jusque 15
    myhist = hist(x, prob = TRUE, 
                  ylim=c(0,0.36), main = "",
                  col=c("grey", rep(adjustcolor(topocolors[match(topo, topologies)], 0.7), 20)), 
                  xaxt="n", yaxt="n", xlab="", ylab="",
                  right=F) # needed otherwise aggregates 0 and 1 values
    # axis(1, at=myhist$mids, labels = rep("", 17))
    # axis(1, at=myhist$mids[8 ], labels = "7", tick = F)
    # axis(1, at=myhist$mids[9 ], labels = "8", tick = F)
    # axis(1, at=myhist$mids[10], labels = "9", tick = F)
    axis(1, at=myhist$mids, labels = rep("", length(myhist$mids)))
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
    x = as.numeric(c(x, 0:17)) # pour que tous aient jusque 15
    myhist = hist(x, prob = TRUE, 
                  ylim=c(0,0.36), main = "",
                  col=c("grey", rep(adjustcolor(topocolors[match(topo, topologies)], 0.7), 20)), 
                  xaxt="n", yaxt="n", ylab="", xlab="",
                  right=F) # needed otherwise aggregates 0 and 1 values
    # axis(1, at=myhist$mids, labels = rep("", 17))
    # axis(1, at=myhist$mids[3 ], labels = "2", tick = F)
    # axis(1, at=myhist$mids[4 ], labels = "3", tick = F)
    axis(2, at=c(0,0.1,0.2, 0.3), labels = c("0", "", "", "0.3"), las=1)
    axis(1, at=myhist$mids, labels = rep("", length(myhist$mids)))
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
    x = as.numeric(c(x, 0:17)) # pour que tous aient jusque 15
    myhist = hist(x, prob = TRUE, 
                  ylim=c(0,0.36), main = "",
                  col=c("grey", rep(adjustcolor(topocolors[match(topo, topologies)], 0.7), 20)), 
                  xaxt="n", yaxt="n", xlab="", ylab="",
                  right=F) # needed otherwise aggregates 0 and 1 values
    axis(1, at=myhist$mids, labels = rep("", length(myhist$mids)))
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
    yy = myhist$density[2:17][1:7]
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
  cat("FIG 2 COMPLETED.\n")
  
} # LAUNCH WHOLE FIGURE

