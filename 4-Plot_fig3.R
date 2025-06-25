# BASED ON SCRIPT 103

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Set working directory where the file is.
# getwd()
source("src/GeneRation_Fun_v1.R")
source("src/GeneRation_DataAnalysis_FUN.R")

##########################
# LOAD DATA
##########################

# 1. PAPER DATA
# store_Fits_BATCH_succ = read.csv("USED_DATA//Extracted_from_successful_simulations/dreamburnin_SIMUSET_Fits.csv")[,-1] # takes 10 seconds

# 2. FEW SIMULATIONS
#store_Fits_BATCH_succ = read.csv("USED_DATA/Extracted_data_small_set/storeFits_november")[,-1]

# 3. TUTO_SIMUS
store_Fits_BATCH_succ = read.csv("Extracted_data/storeFits_TUTO_SIMUS")[,-1]


##########################
# PRE-PROCESSING : 
##########################

nrow(store_Fits_BATCH_succ) ; play()
generations = ncol(store_Fits_BATCH_succ)-2
store_Fits_BATCH = store_Fits_BATCH_succ


# SPLIT IN 3  #############
topologies
for (topo in topologies) {
  # SPLIT to get store_Fits_TOPOLOGY
  store_Fits_TOPOLOGY = store_Fits_BATCH[which(store_Fits_BATCH[,"topo"]==topo),]
  # remove the identifiers for now (two first columns)
  store_Fits_TOPOLOGY = store_Fits_TOPOLOGY[,-c(1,2)]
  # Assign
  assign(paste0("store_Fits_", topo), store_Fits_TOPOLOGY)
  print(paste0(topo, " ", nrow(get(paste0("store_Fits_", topo)))))
}


###############################################################################
###############################################################################
#####  PLOT                                                                 ###
###############################################################################
###############################################################################

# Chose which to display: 
chosentopologies = topologies#[c(3,1,2)]
colrs = topocolors
envcols = c("#B89D6E", "#6EA9B8")
onlysuccessfuls = T # ON PREND QUE LES SUCCESSFUL LA.

width = 20 # unit is cm
height =  11 # unit is cm
{ # ALL PDF
  pdf(file = paste0("figures/Rplot", format(Sys.time(), "%Y%m%d%H%M%S"), "_fig3_", width, "x", height,".pdf"), 
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
  cat("FIG 3 COMPLETED.\n")
} # ALL PDF
