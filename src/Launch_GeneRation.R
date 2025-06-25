#!/usr/bin/env Rscript

#####################################################
#
# Simul_Prog.R
#
# R clone of the C++ Simul_Prog simulation program
# So far, many functions and parameters are missing
#
######################################################
system(paste("echo ---", sep=""))

cmd <- commandArgs(trailingOnly=FALSE)
my.path <- cmd[grep(cmd, pattern="--file=")]
my.path <- dirname(strsplit(my.path, split='=')[[1]][2])
date=format(Sys.time(), "%m%d_%H%M%S_")

# -p parameterfile
which.p <- which(cmd=="-p")
param.file <- cmd[which.p+1]
if (length(which.p)!=1){
	system("echo 'SIMUL_GENER\t100\nINIT_PSIZE\t100\nGENET_NBLOC\t5\nINIT_BASAL \t0.2\nINIT_CLONAL\tnotclonal\nINIT_CONDIAG\t0\nINIT_ALLELES\t0.0\t0.00002\nINIT_TRANSALLELES\t1.0 0.00002\nTYPE_ALLELES\tzero\nGENET_MUTRATES\t0.01\nGENET_TRANSMUTRATES \t0.01\nGENET_MUTSD\t0.5 sqrt\nFITNESS_OPTIMUM\trandom\nFITNESS_STRENGTH\t10 10 0 0 0\nFITNESS_STABSTR\t46000\nSIMUL_OUTPUT\t1\nDEV_TIMESTEPS\t20\nDEV_CALCSTEPS\t2\nTF_REG\tboth\nGENET_MUTTYPE	\tindividual\nGENET_RECRATES\t0.5\nGENET_SELFING\t0.0\nGENET_CLONAL\t0.0\nGENET_PLOIDY\t2\nGENET_EPIGENET\t0.0\nFITNESS_TYPE\tgaussian\nINIT_CONNECT\t1\nTYPE_SO\tbasal\nFITNESS_STAB\texponential_stab\nOUT_UNSTAB\tyes\nOUT_GENO\tyes\nOUT_CANAL_TESTS\t0\nOUT_CANAL_MUTSD\t0.5\nOUT_HERIT_TESTS\t0\nOUT_DIREPI_TESTS\t0\nINIT_RECURRENCE\t0\nFITNESS_FLUCT\tno_fluctuation\nSIMUL_MAXGEN\t5000\nTYPE_ARCHI\tm2\n' > defaultparam.txt")
	system("echo Default parameter file generated.")
	param.file="defaultparam.txt"
}else{
	stopifnot(length(which.p)==1, length(cmd) > which.p)
	param.file <- cmd[which.p+1]
}

# -ipop initialpop
which.ipop <- which(cmd=="-ipop") ##load pop initiale
init.pop <- cmd[which.ipop+1]			##load pop initiale

# -o outputname
output.file <- ""
which.o <- which(cmd=="-o")
if (length(which.o > 0) && length(cmd) > which.o) {
	output.file <- cmd[which.o+1]
}
if (output.file==""){output.file="simulation"}


# -fun : function file to be selected
functions_file <- "GeneRation_Fun_v1.R" ######### OR MODIFY BY HAND DEFAULT HERE
which.fun <- which(cmd=="-fun")
if (length(which.fun > 0) && length(cmd) > which.fun) {
  functions_file <- cmd[which.fun+1]
  system(paste0("echo Using GeneRation Functions file: ", functions_file, "."))
}


#Create a folder in which the files will be moved, and copy params
system(paste("mkdir ", date, output.file, sep=""))
system(paste("cp ", param.file, " ",  date, output.file, ".param", sep=""))
system(paste("mv ",
							date, output.file, ".param ",
							date, output.file,
							sep=""))
param.file=paste(date, output.file, "/", date, output.file, ".param", sep="")
#################################################################################################################################""

source(paste(my.path, functions_file, sep="/"))
library(rlist)
cat("Running simulation ",date, output.file,": \n",sep="")
#Rprof("prof.prof")
#Here we verify and use if needed the initial population
if (length(init.pop)==0){sim <- launchprogevol(paramfile = param.file, initialpop=NULL) #if no pop initial loadee
									 }else{sim <- launchprogevol(paramfile = param.file, initialpop=list.load(init.pop))}
#Rprof(NULL)

#OUTPUTS
write.table(sim[[1]], file=paste(date, output.file, ".table", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
list.save(sim[[2]], file=paste(date, output.file, ".initialpop.rds", sep=""))
list.save(sim[[3]], file=paste(date, output.file, ".finalpop.rds", sep=""))

#move files to folder
system(paste("mv ",
							date, output.file, ".table ",
							date, output.file, ".initialpop.rds ",
							date, output.file, ".finalpop.rds ",
							date, output.file,
							sep=""))

# simu=LoadSimu(paste(date, output.file, "/", sep=""))
# png(paste (date, output.file, "/", simu[[4]], ".00.graphall", sep=""), width=1000, height=600)
# 	graphall(simu[[1]]) ; invisible(dev.off())
# 	relinit=Quantifycoding(simu, param = simu[[5]] ,reference ="INIT", write=F)$relevant
# png(paste (date, output.file, "/", simu[[4]], ".04.relevant_toInit.heatmap", sep=""), width=700, height=700)
#   WTab(relinit, cellnote = TRUE, main="Relevant matrix compared to Init", devnew=F) ; invisible(dev.off())
# relinit=Quantifycoding(simu, param = param ,reference ="INIT", write=F)$relevant
# mat=WConvert(relinit) ; mat[abs(mat)<0.01]=0
# if(sum(mat)!=0){
# 	png(paste (date, output.file, "/", simu[[4]], ".04.relevant_toInit.network", sep=""), width=700, height=500)
# 	options(warn=-1) ; WNetwork(relinit, tkplot = FALSE) ;  invisible(dev.off()) ; options(warn=0)
# }else{cat("NO INTERACTION DETECTED. No plot generated.\n")}

simu=LoadSimu(paste(date, output.file, "/", sep=""))
png(paste (date, output.file, "/", simu[[4]], ".00.graphall.png", sep=""), width=1000, height=600)
	graphall(simu[[1]])
	invisible(dev.off())
	relinit=Quantifycoding(simu, param = simu[[5]] ,reference ="INIT", write=F)$relevant
png(paste (date, output.file, "/", simu[[4]], ".01.WInit.png", sep=""), width=1000, height=450)
  	par(mfrow=c(1,2))
  	W = WFinale(output=simu[[1]], gen=1)
  	WTab3(W, cellnote = T, from = T, round_nb = 2, on = T)
  	options(warn=-1) ; WNetwork4(W) ; invisible(dev.off()) ; options(warn=0)
png(paste (date, output.file, "/", simu[[4]], ".02.WFinal.png", sep=""), width=1000, height=450)
  	par(mfrow=c(1,2))
  	W = WFinale(output=simu[[1]])
  	WTab3(W, cellnote = T, from = T, round_nb = 2, on = T)
  	options(warn=-1) ; WNetwork4(W) ; invisible(dev.off()) ; options(warn=0)
png(paste (date, output.file, "/", simu[[4]], ".03.WFinal_Kinetics.png", sep=""), width=750, height=450)
  W = WFinale(output=simu[[1]])	
  WKinetics3(W, end = 19, optima = gv("FITNESS_OPTIMUM", simu[[5]])[which(gv("FITNESS_STRENGTH", simu[[5]])!=0)])


#Message de courtoisie pour la sortie
system(paste("echo Final fitness of pop: ", round(sim[[1]]$MFit[nrow(sim[[1]])],3), sep=""))
system(paste("echo JOB COMPLETED", sep=""))
system(paste("echo ---", sep=""))

#Petit bip de notification
#if(file.exists("~/Music/snap.mp3")){system("play ~/Music/snap.mp3 >/dev/null 2>&1")} #hide the output heheeee!!
#if(file.exists("~/Music/snap.mp3")){system("afplay /Users/sylvain/Music/snap.mp3")}
