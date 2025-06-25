##Modified version thats supposed to work without sed, via R commands

launchprogevol = function(paramfile="param.txt", initialpop=NULL) { #launch the program WITH A PARAMETER FILE  ""
  #read parameters
  parameters=read.param(paramfile)
  #check vector length (if too small, copy the first term as many times as needed)
  L=      gv("GENET_NBLOC", parameters)
  theta = gv("FITNESS_OPTIMUM", parameters)
  if (length(theta)==1){
    if (theta=="random"){
      theta=runif(L)
      paramms  <- readLines(paramfile)
      tx2  <- gsub(pattern = "random", replace = paste(theta, collapse = ' '), x = paramms)
      writeLines(tx2, con=paramfile)
      system("echo Fitness optima randomly determined: the param file was modified.")
    }
  }
    if (length(theta)<L){ for (i in (length(theta)+1):L){theta[length(theta)+1]=theta[1]}}
  s1=     gv("FITNESS_STRENGTH", parameters) ; if (length(s1)<L){ for (i in (length(s1)+1):L){s1[length(s1)+1]=s1[1]}}
  s2=     gv("FITNESS_STABSTR", parameters)  ; if (length(s2)<L){ for (i in (length(s2)+1):L){s2[length(s2)+1]=s2[1]}}
  muteff=gv("GENET_MUTSD", parameters)
  if(length(muteff)==1){muteff=c(muteff, muteff)}
  if (muteff[2]=="sqrt"){muteff=c(as.numeric(muteff[1]), sqrt(1/L)*as.numeric(muteff[1]))}
  #cat(muteff, "\n")

  #Launch simulation
  runevo=runevolution(N=gv("INIT_PSIZE", parameters),
                      gen=gv("SIMUL_GENER", parameters),
                      L=L,
                      mucis=gv("GENET_MUTRATES", parameters),
                      mutrans=gv("GENET_TRANSMUTRATES", parameters),
                      muteff=muteff,
                      a=gv("INIT_BASAL", parameters),
                      theta=theta,
                      s1=s1,
                      s2=s2,
                      pas=gv("SIMUL_OUTPUT", parameters),
                      steps=gv("DEV_TIMESTEPS", parameters),
                      measure=gv("DEV_CALCSTEPS", parameters),
                      diag=gv("INIT_CONDIAG", parameters),
                      clon=gv("INIT_CLONAL", parameters),
                      initall=gv("INIT_ALLELES", parameters),
                      inittrans=gv("INIT_TRANSALLELES", parameters),
                      tfreg=gv("TF_REG", parameters),
                      initialpop)
  return(runevo)
} #end launchprogevol

read.param = function(paramfile) { #from Estelle: turn un the param file into a list
	filterpar = function(line) {
		line = line[-1] # Remove the name tag
		# Two possibilities: the line contains only numbers -> returns a numeric vector ; the line contains at least a non-number -> returns a character vector
		ans = suppressWarnings(as.numeric(line))
		if(any(is.na(ans))) ans = line
		return(ans)
	}
	stopifnot(is.character(paramfile), file.exists(paramfile))
	ss  = scan(file=paramfile, what=character(), sep="\n", quiet=TRUE)
	ss  = strsplit(ss, split="\\s+")
	nam = sapply(ss, "[", 1)
	ans = lapply(ss, filterpar)
	names(ans) = nam
	return(ans)
} #end read.param

gv=function(name, param){ # gv: get value: extracts the variable from parameters
  stopifnot(is.character(name), length(name)==1) ; stopifnot(is.list(param))
  if (! name %in% names(param)) stop("Tag ", name, " not found in the parameter file.")
  return(param[[name]])
} #end gv

#============================================================================================================================================================
#============================================================================================================================================================

runevolution = function(N=100, gen=100, L=5, mucis=0.5, mutrans=0.5, a=0.2, pas=1, steps=20, measure=2, diag=0, muteff=c(0.5, 0.5), clon="notclonal", initall=c(0,0.01), inittrans=c(1,0.01), initialpop=NULL, tfreg="unique",...){

  #CREATE INITIAL POPULATION - from file if asked
  if(length(initialpop)!=0){ #i.e., if I asked to use an initial pop
    pop=initialpop ; #load an initial pop if asked for, check that the matrices correspond
    stopifnot(nrow(pop[[1]]$ind)==L)
    #stopifnot(length(initialpop)==N) #un jour il faudra degager N et L du model, cad qu'ils soient mesurEs
    cat("Loaded: initial population of",length(initialpop), "individuals.\n")
    if(length(pop)!=N){cat("New population of", N, "individuals created.\n")}
    }else{
      #or create a new initial population
      if (clon=="clonal"){    #~~mettre en init_clonal yesss ##hard
            pop = development(populclonale(N, L,diag, initall, inittrans, tfreg, ...), N, a, L, steps, measure, tfreg, ...)
      }else{pop = development(populrandom(N, L, diag, initall, inittrans, tfreg, ...), N, a, L, steps, measure, tfreg, ...)}
  } #end iniitalpop test
  initpop=pop

  #create and name columns of the output dataframe
  output = data.frame(matrix(ncol = 1+L*6+2+L*L*2, nrow = 0))
  names=c("Gen")
  for (i in c("MPhen","VPhen","MUnstab","VUnstab", "MTrans", "VTrans")){
    for (j in 1:L){names[length(names)+1]=paste(i, j, sep="")}}
  names[length(names)+1] = "MFit"
  names[length(names)+1] = "VFit"
  for (i in c("MeanAll","VarAll")){
    for (j in 1:squared(L)){names[length(names)+1]=paste(i, j, sep="")}}
  colnames(output) = names

  myseq=dna(50) ; count=c()
  cat(myseq, "\n") #Affichage
  #HERE IT STARTS:
  for (gg in 1:gen){ #for each generation
    #FIRST, THE INFO FROM THE POPULATION IS RECORDED
    if(gg==1 || gg == gen || gg %% pas==0) { #definition du pas pour sortir les datas, pour avoir la premiere, la derniere et sorties intermediaires
      output[nrow(output)+1,] = NA #add empty line to dataframe
      output[nrow(output), "Gen"] = gg
      output[nrow(output), "MFit"] = mean(sapply(1:length(pop), function(i) pop[[i]]$fitness))
      output[nrow(output), "VFit"] = var(sapply(1:length(pop), function(i) pop[[i]]$fitness))
      for (i in 1:L){
        N2=length(pop)
        output[nrow(output), paste("MPhen",   i, sep="")] = mean(sapply(1:N2, function(j) pop[[j]]$mean[[i]]))
        output[nrow(output), paste("VPhen",   i, sep="")] =  var(sapply(1:N2, function(j) pop[[j]]$mean[[i]]))
        output[nrow(output), paste("MUnstab", i, sep="")] = mean(sapply(1:N2, function(j) pop[[j]]$var[[i]]))
        output[nrow(output), paste("VUnstab", i, sep="")] =  var(sapply(1:N2, function(j) pop[[j]]$var[[i]]))
        output[nrow(output), paste("MTrans",  i, sep="")] = mean(sapply(1:N2, function(j) pop[[j]]$ind[[i,L+1]]))
        output[nrow(output), paste("VTrans",  i, sep="")] =  var(sapply(1:N2, function(j) pop[[j]]$ind[[i,L+1]]))
      }
      for (i in 1:squared(L)){
        output[nrow(output), paste("MeanAll", i, sep="")] = mean(sapply(1:length(pop), function(j) pop[[j]]$ind[i]))
        output[nrow(output), paste("VarAll",  i, sep="")] =  var(sapply(1:length(pop), function(j) pop[[j]]$ind[i]))
      }
    } #end of recording

    #THEN, A NEW POPULATION IS CREATED
    #if in order to get the last generation in pop, not the last+1.
    if (gg != gen) {pop = development(newpopul(pop, N, tfreg, L, mucis, mutrans, diag, muteff), N, a, L, steps, measure, tfreg, ...)}
    #Affichage
    if (gg==1 | (gg %% ceiling(gen/50))==0 | gg==gen){count=c(1,count) ; cat(substr(comp(myseq),length(count),length(count)))} ; if (gg == gen){cat("\n")} #Affichage
  } #end of a generation
  simu=list(output,initpop,pop)
  return(simu) #return simu
} #end runevolution




        #######    ##     ##    ####      ###   #F#####   #F#########  ###    #######     ####      ###    #######
      #FF#######  ###     ###  ######     ###  #FF######  F##########  ###  ###########  ######     ###  ###########
      #FF    ###  ###     ###  #######    ###  #FF    ##      #FF      ###  ###     ###  #######    ###  ###     ###
      #FF         ###     ###  ###   ##   ###  #FF            #FF      ###  ###     ###  ###   ##   ###  ###
      #FF######   ###     ###  ###    ### ###  #FF            #FF      ###  ###     ###  ###    ### ###  #####
      #FF######   ###     ###  ###     ######  #FF            #FF      ###  ###     ###  ###     ######    #######
      #FF         ###     ###  ###      #####  #FF            #FF      ###  ###     ###  ###      #####        #####
      #FF         ###     ###  ###       ####  #FF            #FF      ###  ###     ###  ###       ####          ###
      #FF         ###     ###  ###        ###  #FF     ##     #FF      ###  ###     ###  ###        ###  ###     ###
      #FF         ###########  ###        ###  #FF#######     #FF      ###  ###########  ###        ###  ###########
      #FF           #######    ###        ###   #F######      #FF      ###    #######    ###        ###    #######



#POPULATIONS===========================================================================================================

indiv = function(L=5, diag=0, initall=c(0,0.001), inittrans=c(1,0.001), tfreg="unique", homoz=FALSE, ...){
  ind=list()                  #INIT_ALLELES
  if (tfreg=="unique"){
    ind$mom = matrix(rnorm(L*L, initall[1], initall[2]), nrow=L)
    coding=c(); for (i in 1:L){coding=c(coding, rnorm(1, sample(c(inittrans[1], -inittrans[1]),1), inittrans[2]))}
    ind$mom=cbind(ind$mom,coding)
    if (homoz==TRUE){ind$dad = ind$mom
               }else{ind$dad = matrix(rnorm(L*L, initall[1], initall[2]), nrow=L)
                     coding=c(); for (i in 1:L){coding=c(coding, rnorm(1, sample(c(inittrans[1], -inittrans[1]),1), inittrans[2]))}
                     ind$dad=cbind(ind$dad,coding)}
  }else{
    ind$mom = matrix(rnorm(L*L, initall[1], initall[2]), nrow=L)
    coding=rnorm(L, inittrans[1], inittrans[2])
    ind$mom=cbind(ind$mom,coding)
    if (homoz==TRUE){ind$dad = ind$mom
               }else{ind$dad = matrix(rnorm(L*L, initall[1], initall[2]), nrow=L)
                 coding=rnorm(L, inittrans[1], inittrans[2])
                 ind$dad=cbind(ind$dad,coding)}
  }
  if (diag==0){for (i in 1:L){ind$mom[i,i]=0 ; ind$dad[i,i]=0}}
  ind$ind = (ind$mom+ind$dad)/2
  return(ind)
}

populrandom = function(N=10, L=5, diag=0, initall=c(0,0.001), inittrans=c(1,0.001), tfreg, developed=FALSE, ...){
  pop=list()
  pop = lapply(1:N, function(i) pop[[i]]=indiv(L, diag, initall, inittrans, tfreg, ...))
  if (developed==TRUE){
    pop=development(pop, tfreg)
    #for (i in 1:10){pop[[i]]$fitness=runif(1)}
  }
  return(pop)
}

populclonale = function(N=10, L=5, diag=0, initall=c(0,0.001), inittrans=c(1,0.001), tfreg, ...){
  pop = list()
  Lelu = indiv(L, diag, initall, inittrans, homoz=TRUE, tfreg, ...) #car cest lui l'elu
  pop = lapply(1:N, function(i) pop[[i]]=Lelu)
  return(pop)
}

populdevfromW=function(W, N=10, tfreg="unique"){ #generates a developed clonal population from a single matrix
  ind=list()
  ind$mom=W; ind$dad=W; ind$ind=W
  pop=list()
  pop = lapply(1:N, function(i) pop[[i]]=ind)
  pop = development(pop, tfreg)
  return(pop)
}

#DEVELOPMENT===========================================================================================================

development = function(pop, N=10, a=0.2, L=5, steps=40, measure=10, tfreg, ...){
  pop = lapply(pop, function(i){ #~mclapply
    dev = model.M2(i$ind, a, steps, measure, tfreg=tfreg)
    i$mean = dev$mean
    i$var  = dev$var
    i = addfitness(i, L, ...)
    })#, mc.cores=3)
  return(pop)
}

sigma.M2p = function(x, aam1, l1am1) { 1. / (1. + exp((-x/aam1)+l1am1)) }
model.M2 = function(W, a=0.2, steps=40, measure=10, init= rep(a,nrow(W)), S0=init, full=FALSE, varFUN=function(x) mean((x-mean(x))^2), tfreg="unique") {
#note# steps:pas de temps du developement  #measure:on fait les mesure sur les x derniers steps du development
#note# init: runif(nrow(W), min=0, max=1), ##ici toutes les valuers initiales sont aleatoires : entre individus et entre generaztions...attention
#note# S0:#minimum,maximum,median,random_binary,random,basal
  aam1 = a*(1-a) ; l1am1 = log(1/a-1)
  sto = matrix(NA, nrow=length(S0), ncol=steps+1)
  sto[,1] = S0
  L=nrow(W)
  LL=ncol(W)
  if(L==LL){W2=W}else{
    transeffect=W[,LL] ; W=W[,1:L]  #on separe les arguments qui vont passer dans l'equation
    transeffect[transeffect>2]=2
    if (tfreg=="unique"){
      W[W<0]=0
      transeffect[transeffect<-2]=-2
    } #Tous ce qui est en desous de zero correspond à des vrais zeros. sinon -x-=+
    W2 = t(t(W)*transeffect)
  }
  for (i in 1:steps) {
    S0 = sigma.M2p(W2 %*% S0, aam1, l1am1) 	#not S0 %*% W
    sto[,i+1] = S0
  }
  ans = list()
  ans$mean = apply(sto[,(steps+1-measure):(steps+1)], 1, mean)
  ans$var  = apply(sto[,(steps+1-measure):(steps+1)], 1, varFUN)
  if (full) ans$full = sto
  return(ans)
}

addfitness = function(ind, L=5, theta=rep(1,L), s1=rep(10,L), s2=rep(4600,L), ...) {
#optimum pour les 3 genes a la fois et #force de selection
  fitmean = 0 ; fitvar = 0
  for (i in 1:L) {
    fitmean = fitmean + (-s1[i]*(squared(ind$mean[i]-theta[i])))
    fitvar  = fitvar  + (-s2[i]*(ind$var[i]))
    }
  ind$fitness = exp(fitvar)*exp(fitmean)
  return(ind)
}

#NEWPOP=============================================================================================================

newpopul = function(pop, N=10, tfreg, ...){
  newpop = list()
  newpop = lapply(1:N, function(i) newpop[[i]]=reproduction(pop, N, tfreg, ...))
  return(newpop)
}

reproduction = function(pop, N, tfreg, ...){
#create new individu from pop
  ind=list()
  #xx=TRUE ; while (xx==TRUE) { #Control autofecondation
    mom = selectforfitness(pop, N)
    dad = selectforfitness(pop, N)
  #  if (mom != dad){xx=FALSE}
  #}
  ind$mom = gametogenesis(pop[[mom]], tfreg=tfreg, ...)
  ind$dad = gametogenesis(pop[[dad]], tfreg=tfreg, ...)
  ind$ind = (ind$mom+ind$dad)/2
  return(ind)
}

selectforfitness = function(pop, N=10) {
	fitnesses <- sapply(pop, "[", "fitness")
  if (sum(fitnesses==0)==length(pop)){fitnesses=rep(1,length(pop))} #to avoid Error in sample.int(...) : too few positive probabilities
	return(sample.int(length(pop), 1, prob=fitnesses))
}

gametogenesis = function(ind, L=5, mucis=0.002, mutrans=0.001, diag=1, muteff, tfreg, ...){
  gamete=matrix(NA, L,(L+1))
  mutdad=ind$dad ; mutmom=ind$mom
  mutdad=cismutation(mutdad, L, mucis, diag, muteff)
  mutmom=cismutation(mutmom, L, mucis, diag, muteff)
  mutdad=transmutation(mutdad, L, mutrans, diag, muteff, tfreg)
  mutmom=transmutation(mutmom, L, mutrans, diag, muteff, tfreg)
  for (i in 1:L){
    r=runif(1)
    if (r<0.5){gamete[i,]=mutdad[i,]}else{gamete[i,]=mutmom[i,]}   ##TAUX DE RECOMBINAISON GENET_RECRATES
  } #~ taux de recombinaison marche pas pour different de 0.5
  return(gamete)
}

cismutation = function(W, L=5, mucis=0.001, diag=1, muteff, r=rpois(1,mucis)){
  transeffect=W[,L+1] ; W=W[,1:L]
  if (r>0){ # on tire un nombre de mutation sur une loi de poisson - on utilise si egal a 1 ou si superieur a 1???
    for (i in 1:r){ #car il peut tirer r=2 donc 2 mutations sur l'individu
      gene=sample(c(1:L),1) ; bindsite=sample(c(1:L),1)#tire quel gene va etre mutE, quelle ligne
      if (diag==0){while (gene==bindsite){ bindsite=sample(c(1:L),1) }} #on garde les diag à zero
      combien=rnorm(1,0,muteff[1]) #on tire de combien on va muter
      W[bindsite,gene]=W[bindsite,gene]+combien
    }
  }
  W=cbind(W,transeffect)
  return(W)
}

transmutation = function(W, L=5, mutrans=0.01, diag=1, muteff, tfreg="unique",r=rpois(1,mutrans)){
  transeffect=W[,L+1] ; W=W[,1:L]
  if (r>0){
    for (i in 1:r){
      gene=sample(c(1:L),1)
      combien=rnorm(1,0,muteff[2])
      if (tfreg=="unique"){
        if(transeffect[gene]<0){sign=-1}else{sign=1} #on prend le sign de la region codante (activator or repressor)
        transeffect[gene]=abs(transeffect[gene])+combien #on lui ajoute la mutation
        if(transeffect[gene]<0){transeffect[gene]=0.000001} #ici on empeche le switch activator/repressor, on pourrait mettre un rebound instead of zero - on met pas vraiment ZERO car on veut lui réattribuer un sign
        transeffect[gene] = transeffect[gene]*sign #on conserve le sign (un activator reste un activator, un represseur reste un represseur, mm si inactif)
        }
        else{
          transeffect[gene]=transeffect[gene]+combien
          if (transeffect[gene]<0){transeffect[gene]=0}
        } #important de garder ca sinon on va plus rien comprendre en regardant la matrice. Donc pas de robustness qd on reach zero
      }
  } #on pourrait le faire aller a -2 ssi les cismut, elles, ne pouvaient etre negatives; car -x-=+
  W=cbind(W,transeffect)
  return(W)
}

########################################################################################################################################
##### LOAD SIMULATIONS ETC #############################################################################################################
########################################################################################################################################

LoadSimu=function(name){
  library(rlist)
  rank=unlist(gregexpr(pattern ='/',name)) #chech where "/" are in the address

  if (file.exists(paste(name, ".param", sep=""))==TRUE){ #ca veut dire qu'il a loadé l'interieur
    if(length(rank)==1 & rank[1]==-1){sim=""; realname=name} #loadé depuis l'interieur
    if(length(rank)==1 & rank[1]==1){sim=""; realname=substr(name,1,rank[length(rank)]-1)} #loadé depuis just above
    if(length(rank)>1){sim=""; realname=substr(name,rank[length(rank)-1]+1,rank[length(rank)]-1)} #loadé d'aileurs
  }else{
    if (length(rank)==1){sim=substr(name,1,rank[length(rank)]-1) ;realname=sim}#IF JUST ABOVE THE SIMULATION FOLDER
    if (length(rank)>1){sim=substr(name,rank[length(rank)-1]+1,rank[length(rank)]-1) ; realname=sim}#IF SOMEWHERE ELSE
  }
  simu=list(read.table(paste(name,sim,".table", sep=""), header=TRUE),
              list.load(paste(name,sim, ".initialpop.rds", sep="")),
              list.load(paste(name,sim, ".finalpop.rds", sep="")),
              realname, #stores the name
              read.param(paste(name,sim, ".param", sep=""))) #stores the name
  return(simu)
}

balance=function(simu=simu, ref="INIT", compare=FALSE){ #pour avoir un resumE du truc
  rel=Quantifycoding(simu, reference=ref)$relevant
  Quantifycistrans(simu)
  WNetwork(rel)
  if (compare==TRUE){WCompareNet(mat=rel, ref=WFromPop(simu[[2]]), integrate=TRUE)}
  Quantifycoding(simu, reference=ref)
  Quantifycistrans(simu)
  dev.new() ; graphall(simu[[1]])
  WTab(rel, cellnote = TRUE, main = simu[[4]])
  }#end balance




      #F#####      #######     #######      ######    ###     ###    #######
     #FF######   #FF#######  ###########  #FF#######  ###     ###  ###########
     #FF    ##   #FF    ###  ###     ###  #FF    ###  ###     ###  ###     ###
     #FF         #FF    ###  ###     ###  #FF    ###  ###     ###  ###
     #FF         #FF######   ###########  #FF#######  ###########  #####
     #FF  ####   #FF###      ###########  #FF#####    ###########    #######
     #FF  #####  #FF####     ###     ###  #FF         ###     ###        #####
     #FF  ## ##  #FF  ###    ###     ###  #FF         ###     ###          ###
     #FF     ##  #FF   ###   ###     ###  #FF         ###     ###  ###     ###
     #FF#######  #FF   ###   ###     ###  #FF         ###     ###  ###########
      #F######   #FF   ####  ###     ###  #FF         ###     ###    #######    from simu


#UN PEU PLUS VIF:
bgcol=function(){bg='white'; return(bg)} #899caa ##e7ebf4
mycolors=function(){mycolors=c('#b53122','#892daf','#2c85c2','#25bd66','#e09b2d', '#45b39d', '#dd5c25','#6f6666','#cd30cc') ; return(mycolors)}
varcol=function(){varcol='#dbdbed'; return(varcol)} #8092a6

##DOESNT HURT MY EYES:
#bgcol=function(){bg='#899caa'; return(bg)}
#mycolors=function(){mycolors=c('#78281f','#4a235a','#1b4f72','#186a3b','#b9770e','#b9410e','#0b0b0a','#9e3c9d','#3c9e7e') ; return(mycolors)}
#varcol=function(){varcol='#8092a6'; return(varcol)}

# Here "output" represent the output dataframe used before

graphfit = function(output=simu[[1]]){
  L=length(grep(x = colnames(output), pattern = "MPhen.*"))
  par(bg=bgcol())
  plot(NA, xlim=c(0,output$Gen[nrow(output)]) ,ylim=c(0,1), main="Fitness of the population", xlab="Generations", ylab="Fitness of pop")
  #polygon(c(0,output$Gen,output$Gen[nrow(output)]), c(0,output$MFit,0), col=adj('#9fb1b6',2), border=NA)
  polygon(c(output$Gen, rev(output$Gen)), c(output$MFit+output$VFit, rev(output$MFit-output$VFit)), col=varcol(), border=NA)
  lines(output$Gen, output$MFit, type = 'l', col='#001f3f')
}

graphexp=function(output=simu[[1]]){
  L=length(grep(x = colnames(output), pattern = "MPhen.*"))
  par(bg=bgcol())
  plot(NA, xlim=c(0,output$Gen[nrow(output)]) ,ylim=c(0,1), main="Phenotype - Gene expressions (mean)", xlab="Generations", ylab="Expression")
  mycolrs=mycolors() ; if (length(mycolrs)<L) {for (i in (length(mycolrs)+1):L){mycolrs[i]=mycolrs[i-9]}}
  for (i in 1:L){ polygon(c(output$Gen, rev(output$Gen)), c(output[[paste("MPhen", i, sep="")]]+output[[paste("VPhen", i, sep="")]], rev(output[[paste("MPhen", i, sep="")]]-output[[paste("VPhen", i, sep="")]])), col=varcol(), border=NA)}
  for (i in 1:L){ lines(output$Gen, output[[paste("MPhen", i, sep="")]], type = 'l', col=mycolrs[i])}
  legend("topleft", lty=1, col=mycolrs, legend=c(1:L), bty = "n", cex = 0.75, lwd=2)
}

graphmtrxvalues = function(output=simu[[1]]){ #GRAPH les valeurs de la matrice
  par(bg=bgcol())
  L=length(grep(x = colnames(output), pattern = "MPhen.*"))
  tt=output[,grep(x = colnames(output), pattern = "MeanAll.*")] #extracts values of the MeanAll
  plot(NA, xlim=c(0,output$Gen[nrow(output)]) ,ylim=c(min(tt)-0.1,max(tt)+0.1), main="Regulatory values in population (mean)", xlab="Generations", ylab="Strength of interaction")
  mycolrs=mycolors() ; if (length(mycolrs)<L) {for (i in (length(mycolrs)+1):L){mycolrs[i]=mycolrs[i-9]}}
  newcolors=c()
  #colors par colonnes (par effet DU gene sur les autres):
  for (j in 1:L){for (i in 1:L){newcolors[length(newcolors)+1]=mycolrs[j]}} ; mtext("Effect OF the colored gene ON others", 3, line=.3)
  #colors par lignes (par effet SUR gene des autres):
  #for (j in 1:L){for (i in 1:L){newcolors[length(newcolors)+1]=mycolrs[i]}} ; mtext("Effect ON the colored gene BY others", 3, line=.3)
  for (i in 1:(L*L)){ polygon(c(output$Gen, rev(output$Gen)), c(output[[paste("MeanAll", i, sep="")]]+output[[paste("VarAll", i, sep="")]], rev(output[[paste("MeanAll", i, sep="")]]-output[[paste("VarAll", i, sep="")]])), col=varcol(), border=NA)}
  for (i in 1:(L*L)){ lines(output$Gen, output[[paste("MeanAll", i, sep="")]], type = 'l', col=newcolors[i])}
  legend("topleft", lty=1, col=mycolrs, legend=c(1:L), bty = "n", cex = 0.75, lwd=2)
}

graphtransvector = function(output=simu[[1]]){
  par(bg=bgcol())
  L=length(grep(x = colnames(output), pattern = "MPhen.*"))
  tt=output[,grep(x = colnames(output), pattern = "MTrans.*")]
  plot(NA, xlim=c(0,output$Gen[nrow(output)]) ,ylim=c(-2,2), main="Activity (coding) values in pop (mean)", xlab="Generations", ylab="Strength of interaction / activation force") #ylim=c(min(tt)-0.1,max(tt)+0.1)
  mycolrs=mycolors() ; if (length(mycolrs)<L) {for (i in (length(mycolrs)+1):L){mycolrs[i]=mycolrs[i-9]}}
  for (i in 1:L){ polygon(c(output$Gen, rev(output$Gen)), c(output[[paste("MTrans", i, sep="")]]+output[[paste("VTrans", i, sep="")]], rev(output[[paste("MTrans", i, sep="")]]-output[[paste("VTrans", i, sep="")]])), col=varcol(), border=NA)}
  for (i in 1:L){ lines(output$Gen, output[[paste("MTrans", i, sep="")]], type = 'l', col=mycolrs[i])}
  legend('topleft',  col=mycolrs, legend=c(1:L), bty = "n", cex = 0.75, lwd=2)
}

graphall=function(output=simu[[1]]){
  par(mfrow=c(2,2), bg=bgcol())
  graphfit(output)
  graphmtrxvalues(output)
  graphexp(output)
  graphtransvector(output)
}


            ##      ##      ##
           ###     ####     ###
           ###     ####     ###
           ###     ####     ###
           ###     ####     ###
           ###     ####     ###
           ###     ####     ###   #####
           ###     ####     ###  ##
           ###     ####     ###   ####
           ##################       ##
            ######    ######    #####      Extraction, visualisation et traitement des matrices d'interaction W.



WFinale = function(output=simu[[1]], gen=nrow(output), noCOD=FALSE, tfreg="unique"){
#extraire la matrice moyenne de la gen (dernière par default) generation du dataframe de sortie
  L=length(grep(x = colnames(output), pattern = "MPhen.*"))
  lastW=c()
  for (i in 1:(L*L)){lastW = c(lastW, output[[paste("MeanAll", i, sep="")]][gen])} #get interaction values
  if (length(output[["MTrans1"]])>0){
    for (i in 1:L){ lastW = c(lastW, output[[paste("MTrans",  i, sep="")]][gen]) }
    lastW=matrix(lastW, nrow=L)
    if (noCOD==TRUE){
      transeffect=lastW[,(L+1)] ; lastW=lastW[,1:L]
      if (tfreg=="unique"){lastW[lastW<0]=0}
      lastW=t(t(lastW)*transeffect) #Merged Coding vector in order to only have interactions
    }
  }else{lastW=matrix(lastW, nrow=L)}
  return(lastW)
}

WFromPop = function(pop=simu[[2]], median=FALSE){
#extraire la matrice moyenne/mediane d'une population
  L=nrow(pop[[1]]$ind)
  LL=ncol(pop[[1]]$ind)
  matmean=c() ; matmedian=c()
  for (k in 1:L){
    for (j in 1:(LL)){
      matmean[length(matmean)+1]=mean(sapply(pop, function(i) i$ind[k,j]))
      matmedian[length(matmedian)+1]=median(sapply(pop, function(i) i$ind[k,j]))
    }
  }
  matmean=matrix(matmean, nrow=L, ncol=LL, byrow = TRUE)
  matmedian=matrix(matmedian, nrow=L, ncol=LL, byrow = TRUE)
  if (median==TRUE){return(matmedian)}else{return(matmean)}
}

WConvert=function(W, tfreg="unique"){ #integrate the coding vector / merge / convert / transform
  L=nrow(W)
  if(ncol(W)==L){stop("Matrice is square, no convertion possible.")}
  transeffect=W[,(L+1)] ; W=W[,1:L]
  if (tfreg=="unique"){W[W<0]=0}
  W=t(t(W)*transeffect)
  return(W)
}

WClean=function(W, threshold=0.01){ #permets d'enlever toutes les petites interactions parasites qui ont une force en dessous de 0.01.
 L=nrow(W)
 LL=ncol(W)
 if (LL==L){W[abs(W)<0.01]=0}else{#clean matrice carree
   for (j in 1:L){#clean matrix non carree
     for (i in 1:L){
       if (abs(W[i,j]*W[j,(L+1)])<threshold){W[i,j]=0}
     }
     if (mean(W[,j])==0){W[j,(L+1)]=0}
   }
 }
 return(W)
} # W[abs(W)<threshold]=0

#########################################
#  VISUALIZE                            #
#########################################

WKinetics=function(W, start=1, end=10, main=paste("Kinetics of ",deparse(substitute(W))), tfreg="unique"){
  L=nrow(W)
  clr=c('#b53122','#892daf','#2c85c2','#25bd66','#e09b2d', '#45b39d', '#dd5c25','#6f6666','#cd30cc')
  if (length(clr)<L) {for (i in (length(clr)+1):L){clr[i]=clr[i-9]}}
  par(mar=c(5,4,2,6))
  plot(NA, xlim=c(start,end), ylim=c(0,1), xlab="Timesteps", ylab="Gene expression", las=1, main=main)
  for (i in seq(0,1,0.1)){abline(h=i, col=adj(mygrey,3))}
  for (i in seq(start,end,1)){abline(v=i, col=adj(mygrey,3))}
  for (i in c(1:L)){lines(model.M2(W, full = T, steps=end, tfreg=tfreg)$full[i,start:end], col=clr[i], lwd=2)} #tfreg="unique" ERROR WAS HERE @202403061258
  legend("topleft", legend=c(1:L), lty=1, lwd=2, inset=c(1,0), xpd=TRUE, bty="n", col=clr[1:L])
}#end WKinetics

WKinetics3 = function(W, start=0, end=19, 
                      main="", #paste("Kinetics of ",deparse(substitute(W))), 
                      tfreg="NOTunique", 
                      optima = simu[[5]]$FITNESS_OPTIMUM[1:5], 
                      the_points=""){
  
  # from SCRIPT 131
  # ajouter un truc pour gérer le fait que optima est parfois vide.
  # -> JUst mettre optima = NA
  
  # parce que ca commence en fait à 1 on va faire un shift
  start = start+1
  end = end+1
  
  L = nrow(W)
  clr=c('#b53122','#892daf','#2c85c2','#25bd66','#e09b2d', '#45b39d', '#dd5c25','#6f6666','#cd30cc')
  mygrey = "#6F636399"
  if (length(clr)<L) {for (i in (length(clr)+1):L){clr[i]=clr[i-9]}}
  
  par(mar=c(5,4,2,6)) ; par(mgp=c(2.2, 1, 0))
  
  plot(NA, xlim=c(start,end+.5), ylim=c(0,1), 
       xlab="", ylab="Normalized gene expression", las=1, main=main)
  
  polygon(x=c(end-2.5, end-2.5, end-.5, end-.5)+1, y=c(-1, 2, 2, -1), col=adjustcolor("grey", 0.4), border = NA)
  for (i in seq(0,1,0.1)){abline(h=i, col=adj(mygrey,30))}
  for (i in seq(start,end,1)){abline(v=i, col=adj(mygrey,30))}
  abline(h=optima, col = adj(clr, 20), lwd = 10)
  
  for (i in c(1:L)){
    the_dev = model.M2(W, full = T, steps=end, tfreg=tfreg)$full[i,start:end]
    lines(the_dev, col=clr[i], lwd=2)
    if (the_points == "all"){points(the_dev, col=clr[i], pch=20)}
    if (the_points == "end"){points(c(end-1,end), the_dev[c(end-1,end)], col=clr[i], pch=20)}
  }
  
  legend("topleft", legend=c(1:L), lty=1, lwd=2, inset=c(1,0), xpd=TRUE, bty="n", col=clr[1:L])
  title(xlab="Developmental timesteps")
  #return default
  par(mar = c(5, 4, 4, 2)+0.1) ; par(mgp=c(3, 1, 0))
} #end Wkinetics 3 


WPlot=function(W, main=deparse(substitute(W))){
  L=nrow(W) ; LL=ncol(W)
  if(LL==L+1){cols=c(rep("#c0392b",L), "#42110b", rep("#9b59b6",L), "#2d0e39", rep("#2980b9",L), "#0c2839", rep("#27ae60",L),"#0a361c", rep("#f39c12",L), "#47310c", rep("#45b39d",L), "#043028", rep("#dd5c25",L), "#5b2007", rep("#6f6666",L), "#000000", rep("#c727c0",L), "#460a50")}
  if(LL==L){cols=c(rep("#c0392b",L), rep("#9b59b6",L), rep("#2980b9",L), rep("#27ae60",L), rep("#f39c12",L), rep("#45b39d",L), rep("#dd5c25",L), rep("#6f6666",L), rep("#c727c0",L))}
  plot(as.vector(t(W)), col="white",ylim=c(-1.5,2.2), xaxt="n", ylab="", xlab="", main=main)
  for (i in 0:(LL-1)){abline(v = i*(LL)+(LL)+0.5, col="grey", lty=5)}
  abline(h=1, col="#17202a", lty=3) ; abline(h=0, col="#2c3e50", lty=3)
  points(as.vector(t(W)), col=cols, pch=19,ylim=c(-1.5,2.2), xaxt="n")
}#end WPlot

#WTab=function(matrx, cellnote=TRUE, main=deparse(substitute(matrx)), devnew=TRUE){
#  if (devnew==TRUE){dev.new()}
#  library(gplots,warn.conflicts = FALSE, quietly=TRUE) ; library(RColorBrewer)
#  my_palette = colorRampPalette(c("#dc2424", "#ededed", "#3db221"))(n = 299)
#  #library(grid) ; library(gridGraphics)
#  if (cellnote==TRUE){matrix2=matrx ; matrix2[matrix2==0]=NA ; cellnote=round(matrix2,2)}else{matrix2=matrx ; matrix2[matrix2!=277]=NA ; cellnote=matrix2}
#  heatmap.2(matrx,
#            cellnote = cellnote,  # same data set for cell labels #drole de truc avec les signes, donc je prends la valeur absolue et les couleurs sont inversees.
#            main = main,          # heat map title
#            notecol="black",      # change font color of cell labels to black
#            density.info="none",  # turns off density plot inside color legend
#            trace="none",         # turns off trace lines inside the heat map
#            margins =c(12,9),     # widens margins around plot
#            col=my_palette,       # use on color palette defined earlier
#            dendrogram="none",    # only draw a row dendrogram
#            Rowv = FALSE,
#            Colv = FALSE,
#            sepcolor="white",
#            sepwidth=c(0.01,0.01, 0.2),
#            colsep=1:ncol(matrx),
#            rowsep=1:nrow(matrx),
#            symm=F,symkey=F,symbreaks=T, scale="none")#to set the color scale la
#}#end WTab function

#New function to avoid the use of the "gplots" car cest dela merde. 
WTab=function(W, cellnote=TRUE, main="", devnew=TRUE){ #main=deparse(substitute(W)) #????
  if (devnew==TRUE){dev.new()}
  L=nrow(W)
  require("RColorBrewer", quietly = T)
  my_palette = colorRampPalette(c("#dc2424", "#ededed", "#3db221"))(n = 299)
  if(ncol(W)==L+1){colnames(W) = c(paste("from ", 1:L, sep=""), "Coding")}
  if(ncol(W)==L  ){colnames(W) = c(paste("from ", 1:L, sep=""))} # so that it works for squared matrices @202402141351
  rownames(W) <- paste("on ", 1:L, sep="")
  zlim=c(-2,2)
  # Need to adjust the colors if they are values exceeding the palette zlim @202402141204
  adjustedW = W
  adjustedW[adjustedW >  2] =  2
  adjustedW[adjustedW < -2] = -2
  image(1:ncol(W), 1:nrow(W), t(apply(adjustedW,2,rev)), col = my_palette, axes = F, ann=F, zlim=zlim)
  title(main=main)
axis(3, 1:ncol(W), colnames(W), las=1)
  axis(2, 1:nrow(W), rev(rownames(W)), las=1)
  if (cellnote){
    for (x in 1:ncol(W))
      for (y in 1:nrow(W))
        text(x, y, round(apply(W,2,rev)[y,x], 2))
  }
  invisible(dev.set()) #@202402141351
}#end WTab function

WTab3=function(W, cellnote=F, main="", devnew=F, round_nb = 1, space_cod = T, from=F, delimitation=T, on=F){
  ## FROM SCRIPT 097
  L=nrow(W)
  if (ncol(W)==L){square_mat = T #is a square met
  } else if (ncol(W)==L+1){square_mat = F #is NOT a square mat -> our classical matrices
  } else {stop("ERROR MATRIX FORMAT!")}
  require("RColorBrewer", quietly = T)
  my_palette = colorRampPalette(c("#dc2424", "#ededed", "#3db221"))(n = 299)
  zlim=c(-2,2)
  if (devnew==TRUE){dev.new()}
  # LABELS
  if (on) {rownames(W) = c("on 1", 2:L)}else{rownames(W) = c(1:L)}#paste("on ", 1:L, sep="")
  if (from){
    if(!square_mat){colnames(W) = c("from 1", 2:L ,"Coding")}
    if( square_mat){colnames(W) = c("from 1", 2:L          )} # so that it works for squared matrices @202402141351
  } else {
    if(!square_mat){colnames(W) = c(1:L ,"Coding")}
    if( square_mat){colnames(W) = c(1:L          )}
  }
  # Adjust values for colors if value > |zlim|
  adjustedW = W
  adjustedW[adjustedW >  2] =  2
  adjustedW[adjustedW < -2] = -2
  at_adjust = 0
  if (space_cod & !square_mat){
    W         = cbind(        W[1:L,1:L], rep(NA, L),         W[,L+1])
    adjustedW = cbind(adjustedW[1:L,1:L], rep(NA, L), adjustedW[,L+1])
    at_adjust = 1
  }
  # PLOT
  image(1:ncol(W), 1:nrow(W), t(apply(adjustedW,2,rev)), 
        col = my_palette, axes = F, ann=F, zlim=zlim)
  title(main=main)
  axis(3, at = 1:L, labels = colnames(W)[1:L],       las=1, pos = L+.5+0.15) #default pos = 10.5
  axis(3, at = L+1+at_adjust, labels = "Coding", tick = F,      las=1, pos = L+.5-0.2)
  axis(2, at = 1:nrow(W), labels = rev(rownames(W)), las=1, pos = 0.5-0.15)  #default pos = 0.5
  if (delimitation){
    for (i in 1:nrow(W)) {
      for (j in 1:ncol(W)) { # lwd for the width of the border
        rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5, border = "white", lwd=1)
      }
    }
  }
  if (cellnote){
    W_text = as.data.frame(round(apply(W,2,rev), round_nb))
    if (space_cod){W_text[,L+1]=rep("__", L)}
    for (x in 1:ncol(W))
      for (y in 1:nrow(W))
        text(x, y, W_text[y,x])
  }else{ # No cell notes but still the "--" separation
    W_text = W
    W_text[,] = NA
    if (space_cod){W_text[,L+1]=rep("__", L)}
    for (x in 1:ncol(W))
      for (y in 1:nrow(W))
        text(x, y, W_text[y,x], col="grey")
  }
  invisible(dev.set()) #@202402141351
}#end WTab3 function

WNetwork=function(mat=WFinale(simu[[1]]), tkplot=TRUE, tfreg="unique"){
  library(igraph,warn.conflicts = FALSE, quietly=TRUE)
  #browser()
  L=nrow(mat)
  if (ncol(mat)!=L){ #ya la colonne trans et on veut une matrix carree
    cod=mat[,(L+1)] ; W=mat[,1:L]
    if (tfreg=="unique"){W[W<0]=0}
    mat=t(t(W)*cod)
    cat("The trans vector has been integrated to generate a square matrix.\n")
  }
  mat[abs(mat)<0.01]=0
  weight=as.vector(mat)
  weight=weight[!weight %in% 0] #weight contains the forces of interactions (quantitative info)

  arrows=c() ;
  #if (length(weight)==0){weight=rep(0,L);arrows=rep(0,L)}else{
  if (length(weight)==0){cat("NO INTERACTION in the matrix. No plot generated.");opt <- options(show.error.messages=FALSE);on.exit(options(opt));stop()}
  #just to avoid the stop() error message #stfu shut
    for (i in 1:length(weight)){if (weight[i]<=0){arrows[i]=0}else{arrows[i]=2}}
    weight=abs(weight)
    #}################################################################

  mat[mat!=0]=1 #mat contains the interactions (qualitative info)
  g = graph.adjacency(t(mat), mode="directed", weighted=NULL) # For directed networks

  mycolors=mycolors() ; if (length(mycolors)<L) {for (i in (length(mycolors)+1):L){mycolors[i]=mycolors[i-9]}}
  edge.start = ends(g, es=E(g), names=F)[,1] #get the starting node for each edge
  edgecolrs=c() ; for (i in 1:length(edge.start)){edgecolrs[i]=mycolors[edge.start[i]]}

  E(g)$arrow.size = .7
  E(g)$arrow.mode = arrows
  E(g)$width = weight*5
  E(g)$color = edgecolrs
  E(g)$curved = .2
  V(g)$frame.color = "white"
  V(g)$label.color = "black"
  V(g)$color=mycolors
  if (tkplot==TRUE){tkplot(g)}else{plot.igraph(g)}
  #plot(g)
  #tkplot(g)
  #http://kateto.net/netscix2016
  detach("package:igraph")
}#end WNetwork function


WNetwork2 = function(W=W, tkplot=F, vcex=1, dist = 0, angles="", lab_col = "black"){ #function on2 02402151600
  library(igraph)
  # W should be WFinal()
  tfreg="NOTunique"
  
  # Integrate trans vector
  L=nrow(W)
  if (ncol(W)!=L){ #ya la colonne trans et on veut une matrix carree donc on va integrer
    cod=W[,(L+1)] ; reg=W[,1:L]
    if (tfreg=="unique"){reg[reg<0]=0} 
    mat=t(t(reg)*cod)
    cat("The trans vector has been integrated to generate a square matrix.\n")
  }
  
  # get rid of too small interactions
  mat[abs(mat)<0.01]=0
  # Check if there are any interaction left in the matrix: # A TESTER !!!!!!!!!!!!!!!!!!!!
  CHECK_INTERACTIONS = T
  if (sum(mat)==0){
    cat("NO INTERACTION in the matrix. No plot generated.");
    opt <- options(show.error.messages=FALSE);
    on.exit(options(opt));
    CHECK_INTERACTIONS = F
    # stop()
    }
  
  if (CHECK_INTERACTIONS){
    
    # DEFINE THE INTERACTIONS:
    interactions = mat
    interactions[interactions!=0] = 1 # Detected interaction -> interactions contains the interactions (qualitative info)
    g = graph.adjacency(t(interactions), mode="directed", weighted=NULL) # For directed networks
    # plot(g)
    # We have created the graph object.
    
    # Create a dataframe to store all the infos for the graph
    alledges = as.data.frame(ends(g, es=E(g)))
    colnames(alledges)=c("start.E", "target.E")
    intnb = nrow(alledges)
    allvertices = data.frame(vertices = cbind(V(g)))
    
    # Add weights
    alledges$interac = mat[mat!=0]
    alledges$weights = abs(mat[mat!=0])
    if (ncol(W)!=L){ #if we have some coding
      allvertices$size = cod
    }else{allvertices$size = 1}
    
    # Activation or Repression (arrows)
    alledges$arrow.mode <- ifelse(alledges$interac <= 0, 0, 2)
    alledges$arrow.lty  <- ifelse(alledges$interac <= 0, 2, 1) # 1 “solid”, 2 (“dashed”), 3 (“dotted”), 4 (“dotdash”), 5 (“longdash”), 6 (“twodash”)
    
    # define colors - the good amount and for interactions
    mycolrs=mycolors()
    allvertices$color = rep(mycolrs, length.out = nrow(allvertices))
    alledges$color = allvertices$color[alledges$start.E] # Assign color according to start.E 
    
    # Define the angles
    if (angles == "" && dist != 0){
      #angles = c(0, -0.2, -0.4, -0.6, -0.8, 1, 0.8, 0.6, 0.4, 0.2)*pi # for the layout_in_circle
      angles = c(-0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 9, -0.9, -0.7, -0.5)*pi # for the clockwise_layout
    } # works only for 10 genes here and a circular representation
    
    # VERTICES 
    V(g)$color        = allvertices$color
    V(g)$frame.color  = "white"
    V(g)$label.color  = "black"
    V(g)$size         = allvertices$size*15
    V(g)$label.cex    = vcex
    V(g)$label.family = "sans"
    V(g)$label.dist   = dist #distance from center of the node
    V(g)$label.degree = angles # angle from center of the node
    V(g)$label.color  = lab_col
    
    # EDGES
    # Trait
    E(g)$width      = alledges$weights*5
    E(g)$color      = alledges$color
    E(g)$curved     = .1
    E(g)$lty        = alledges$arrow.lty
    # Arrow
    E(g)$arrow.mode  = alledges$arrow.mode
    E(g)$arrow.size  = 0 # seq(from=0.2, to=2, length.out=21) # alledges$weights*10 # 0.5 # default=1 #### AHA ! peutetre que c'est just ca le pb.
    E(g)$arrow.width = 1 # alledges$weights*2*3 # 0.5*3 # default=1
    
    # DEFINE A CURSTOM LAYOUT: 1st node on top (instead of on the right), and clockwize distrib.
    layout <- layout_in_circle(make_ring(L)) # Generate the circular layout
    angle <- pi/2  # Find the angle for rotating the layout # 90 degrees to place node 1 at the top
    rotation_matrix <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), ncol = 2) # Create a rotation matrix
    rotated_layout <- layout %*% rotation_matrix # Apply the rotation matrix to the layout
    clockwise_layout <- rotated_layout[order(-seq_len(nrow(rotated_layout))), ] # Reverse the order of the layout to make it clockwise
    
    #g$layout <- layout_in_circle
    g$layout <- clockwise_layout
    
    if (tkplot) {
      tkplot(g)
    }else{
      plot(g)}
    
  }else{
    plot(1, 1, type = "p", axes = FALSE, xlab = "", ylab = "", main = "", pch = "/")
  }
  
  detach("package:igraph")
} # end of Wnetwork2

## COPY OF WNetwork2 => FOR FILM
# Goal is to have the network oriented correctly and numbers out
## SCRIPT 151
WNetwork4 = function(W=W, tkplot=F, vcex="vertices_size", dist = 0.01, angles="", lab_col = "white", speak=F){ #function on2 02402151600
  suppressPackageStartupMessages(library(igraph))
  # W should be WFinal()
  tfreg="NOTunique"
  
  # Integrate trans vector
  L=nrow(W)
  empty_cod=F
  if (ncol(W)!=L){ #ya la colonne trans et on veut une matrix carree donc on va integrer
    cod=W[,(L+1)] ; reg=W[,1:L]
    if(sum(cod)<0.01){empty_cod=T}
    if (tfreg=="unique"){reg[reg<0]=0} 
    mat=t(t(reg)*cod)
    #cat("The trans vector has been integrated to generate a square matrix.\n")
  }
  
  # get rid of too small interactions
  mat[abs(mat)<0.01]=0
  # Check if there are any interaction left in the matrix: # A TESTER !!!!!!!!!!!!!!!!!!!!
  
  # If there is no interaction in the matrix, stop everything:
  # if (sum(mat)==0){
  #   cat("NO INTERACTION in the matrix. No plot generated.");
  #   opt <- options(show.error.messages=FALSE);
  #   on.exit(options(opt));
  #   stop()}
  
  # If there is no interaction in the matrix, USES A FAKE MATRIX AND just change all col to white.
  if (sum(mat)==0){
    if(speak){cat("NO INTERACTION in the matrix.")}
    #mat = matrix(c(0,rep(c(rep(0.5,10),0), 9)), nrow = 10, ncol = 10)
    mat = matrix(c(0,rep(c(rep(0,L),0), (L-1) )), nrow = L, ncol = L)
    mat[L,(L/2)]=0.5
    empty_reg = T
  }else{empty_reg = F}
  
  # DEFINE THE INTERACTIONS:
  interactions = mat
  interactions[interactions!=0] = 1 # Detected interaction -> interactions contains the interactions (qualitative info)
  g = graph.adjacency(t(interactions), mode="directed", weighted=NULL) # For directed networks
  # plot(g)
  # We have created the graph object.
  
  # Create a dataframe to store all the infos for the graph
  alledges = as.data.frame(ends(g, es=E(g)))
  colnames(alledges)=c("start.E", "target.E")
  intnb = nrow(alledges)
  allvertices = data.frame(vertices = cbind(V(g)))
  
  # Add weights
  alledges$interac = mat[mat!=0]
  alledges$weights = abs(mat[mat!=0])
  if (ncol(W)!=L){ #if we have some coding
    allvertices$size = cod
  }else{allvertices$size = 1}
  
  # Activation or Repression (arrows)
  alledges$arrow.mode <- ifelse(alledges$interac <= 0, 0, 2)
  alledges$arrow.lty  <- ifelse(alledges$interac <= 0, 2, 1) # 1 “solid”, 2 (“dashed”), 3 (“dotted”), 4 (“dotdash”), 5 (“longdash”), 6 (“twodash”)
  
  # define colors - the good amount and for interactions
  mycolrs=mycolors()
  allvertices$color = rep(mycolrs, length.out = nrow(allvertices))
  alledges$color = allvertices$color[alledges$start.E] # Assign color according to start.E 
  
  # Define the angles
  if (angles == "" && dist != 0){
    #angles = c(0, -0.2, -0.4, -0.6, -0.8, 1, 0.8, 0.6, 0.4, 0.2)*pi # for the layout_in_circle
    angles = c(c(-0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 9, -0.9, -0.7, -0.5)*pi)[1:L] # for the clockwise_layout
  } # works only for 10 genes here and a circular representation
  
  # VERTICES / Nodes
  V(g)$size         = abs(allvertices$size*15) #*15 because the default is 15 #### C'est ce ABS qui est necessaire et qui m'a causé du soucis on 202407091719
  V(g)$color        = allvertices$color
  V(g)$frame.color  = "white"
  
  V(g)$label.family = "sans"
  V(g)$label.dist   = dist #distance from center of the node
  V(g)$label.degree = angles # angle from center of the node
  if(empty_cod){lab_col=rep("grey", L)} #needed when tout à zero
  V(g)$label.color  = lab_col
  if (vcex=="vertices_size"){
    V(g)$label.cex  = allvertices$size
  }else{
    V(g)$label.cex  = vcex 
  }
  
  # EDGES / traits
  E(g)$width      = alledges$weights*5
  E(g)$color      = alledges$color
  if(empty_reg){E(g)$color = "white"}
  E(g)$curved     = .1
  E(g)$lty        = alledges$arrow.lty
  
  # Arrow
  E(g)$arrow.mode  = alledges$arrow.mode
  E(g)$arrow.size  = 0 # seq(from=0.2, to=2, length.out=21) # alledges$weights*10 # 0.5 # default=1 #### AHA ! peutetre que c'est just ca le pb.
  E(g)$arrow.width = 1 # alledges$weights*2*3 # 0.5*3 # default=1
  
  # DEFINE A CURSTOM LAYOUT: 1st node on top (instead of on the right), and clockwize distrib.
  layout <- layout_in_circle(make_ring(L)) # Generate the circular layout
  angle <- pi/2  # Find the angle for rotating the layout # 90 degrees to place node 1 at the top
  rotation_matrix <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), ncol = 2) # Create a rotation matrix
  rotated_layout <- layout %*% rotation_matrix # Apply the rotation matrix to the layout
  clockwise_layout <- rotated_layout[order(-seq_len(nrow(rotated_layout))), ] # Reverse the order of the layout to make it clockwise
  
  #g$layout <- layout_in_circle
  g$layout <- clockwise_layout
  
  #Now let's add a background graph in slight grey
  g0 = g
  size_factor=0.9
  V(g0)$size         = 15*size_factor
  V(g0)$color        = "white" #c(rep("grey",5), rep("white", 5))
  V(g0)$frame.color  = "grey"
  V(g0)$label.family = "sans"
  V(g0)$label.color  = "grey" #c(rep("grey",5), rep("white", 5))
  V(g0)$label.cex    = 1*size_factor
  E(g0)$color        = "white"
  
  if (tkplot) {
    tkplot(g)
  }else{
    plot(g0)
    plot(g, add = T)}
  
  detach("package:igraph")
} #end WNetwork 4


#########################################
#  TREAT AND COMPARE                    #
#########################################

WCompareNet=function(mat=rel, ref=WFromPop(simu[[2]]), tkplot=TRUE, integrate=FALSE, tfreg="unique"){ #needed: full matrices, ie, CIS+TRANS
  # Will compare two matrices, and draw the difference between them as a Network with colors denoting: red=new link created, orange=modification
  # Cool for network that needs to adapt to a novel environment, we can compare the initial matrix to the new one.
  library(igraph)
  L=nrow(mat)
  matcarree=mat[,1:L] ; matTrans=mat[,(L+1)] ; if (tfreg=="unique"){matcarree[matcarree<0]=0}
  refcarree=ref[,1:L] ; refTrans=ref[,(L+1)] ; if (tfreg=="unique"){refcarree[refcarree<0]=0}
  #cat("The trans vector has NOT been integrated to generate a square matrix.\nHere we identify the place of the mutation, not the final effect.\nThe edges correpond to CIS, the nodes to TRANS.\n")
  if (integrate==TRUE){
    matcarree=t(t(matcarree)*matTrans)
    refcarree=t(t(refcarree)*refTrans)
    cat("The trans vector has been integrated to generate a square matrix.\n")
    }

    matcarree[abs(matcarree)<0.01]=0 ; weight=as.vector(matcarree)
    refcarree[abs(refcarree)<0.01]=0 ; ref=as.vector(refcarree)

  colors=c("#586063", "#c87929", "#da1818", "#1883da")
  edgecolors=c()
  for (i in 1:length(weight)){
    if (ref[i]==0 && weight[i]==0){
    }else if (weight[i]==ref[i]){edgecolors[length(edgecolors)+1]=colors[1]
    }else if (ref[i]==0){edgecolors[length(edgecolors)+1]=colors[3]
    }else if (ref[i]!=0){edgecolors[length(edgecolors)+1]=colors[2]
    }else{edgecolors[length(edgecolors)+1]=colors[4]}
  }
  nodecolors=c()
  for (i in 1:length(matTrans)){
    if (matTrans[i]==refTrans[i]){nodecolors[i]=colors[1]
    }else if (refTrans[i]==0){nodecolors[i]=colors[3]
    }else if (refTrans[i]!=0){nodecolors[i]=colors[2]
    }else{nodecolors[i]=colors[4]}
  }

  weight=weight[!weight %in% 0] #weight contains the forces of interactions (quantitative info)
  arrows=c() ; for (i in 1:length(weight)){if (weight[i]<=0){arrows[i]=0}else{arrows[i]=2}}
  weight=abs(weight)
  matcarree[matcarree!=0]=1 #mat contains the interactions (qualitative info)
  g = graph.adjacency(t(matcarree), mode="directed", weighted=NULL) # For directed networks
  edge.start = ends(g, es=E(g), names=F)[,1] #get the starting node for each edge

  E(g)$arrow.size = .7
  E(g)$arrow.mode = arrows
  E(g)$width = weight*5
  E(g)$color = edgecolors
  E(g)$curved = .2
  V(g)$frame.color = "white"
  V(g)$label.color = "black"
  V(g)$color= nodecolors
  if (tkplot==TRUE){tkplot(g)}else{plot.igraph(g)}
  #plot(g)
  #tkplot(g)
  #http://kateto.net/netscix2016
  detach("package:igraph")
}#end WCompareNet


WBoxplot=function(pop, file, stripchart=FALSE, alpha=0.1, cex=1, main="Distribution of alleles in population"){ #can be launched with either a file containing a pop, or directly a pop
  library(rlist)  #Plots all matrices of all ind of pops #shows ditribution of alleles in the pop
  if (missing(file)){
    test=pop
    N=length(test) ; L=nrow(test[[1]]$ind)
    title=paste("fitness of pop", mean(sapply(1:N, function(i) test[[i]]$fitness)))
  }else{
    test=list.load(file = file)
    N=length(test) ; L=nrow(test[[1]]$ind)
    title=paste(file, "- fitness of pop", mean(sapply(1:N, function(i) test[[i]]$fitness)))
  }
  dat=c()
  for (j in 1:L){
    for (i in 1:(L+1)){
      val=c(); pval=c()
      pval=sapply(1:N, function(k) val[length(val)+1]=test[[k]]$ind[j,i])
      dat=cbind(dat,pval)
    }
  }
  cols=c(rep("#c0392b",L), "#42110b", rep("#9b59b6",L), "#2d0e39", rep("#2980b9",L), "#0c2839", rep("#27ae60",L),"#0a361c", rep("#f39c12",L), "#47310c", rep("#45b39d",L), "#043028", rep("#dd5c25",L), "#5b2007", rep("#6f6666",L), "#000000", rep("#c727c0",L), "#460a50")
  if (stripchart==TRUE){stripchart(as.data.frame(dat), vertical = TRUE, col="white", ylim=c(-1.5,2.2),xaxt="n", main=main)
  }else{
  boxplot(dat[,c(1:(L*(L+1)))], ylim=c(-1.5,2.2),xaxt="n", medcol="white", whiskcol="white",staplecol="white",boxcol="white",outcol="white",outcex=0.5,outpch=1, main=main)}
  for (i in 0:(L-2)){abline(v = i*(L+1)+(L+1)+0.5, col="grey", lty=5)}
  abline(h=1, col="#17202a", lty=3) ; abline(h=0, col="#2c3e50", lty=3)
  if (stripchart==TRUE){stripchart(as.data.frame(dat), add=TRUE, vertical = TRUE, main=main, method = "jitter", cex=cex, pch=20, col=adjustcolor(cols, alpha), ylim=c(-1.5,2.2))}else{
  boxplot(dat[,c(1:(L*(L+1)))], ylim=c(-1.5,2.2), add=TRUE, xaxt="n", main=main, medcol=cols, whiskcol=cols,staplecol=cols,boxcol=cols,outcol="grey",outcex=0.5,outpch=1)}
  #fitness dans la pop
  a=paste("fitness of pop:", round(mean(sapply(1:N, function(i) test[[i]]$fitness)),2))
  #matrix moyenne in pop
  b=round(Reduce('+', lapply(1:N, function(i) test[[i]]$ind))/N, 2)
  return(list(a,b))
} # end WBoxplot


WDegree=function(W, plot=FALSE){
  L=nrow(W)
  W=WClean(W)
  newmat=W[,1:L]
  #sum(newmat!=0) #48 connexions
  gene=c(); connex=c() ;
  for (i in 1:L){
    gene[i]=i ; connex[i]=sum(sum(newmat[,i]!=0), sum(newmat[i,]!=0))
  }
  deg=cbind(gene, connex)[order(connex, decreasing = TRUE),]
  if (plot==TRUE){plot(table(deg[,2]), type = "p", pch=16, col="grey65", xlab="Nombre de connexions", ylab="Nomnre de genes")}
  return(deg)
}#end WDegree function



            #########      #######    #F#########    #######          #F#########    #######     #######      #######    #F#########
            ###########  ###########  F##########  ###########        F##########  #FF#######  ###########  ###########  F##########
            ###     ###  ###     ###      #FF      ###     ###            #FF      #FF    ###  ###          ###     ###      #FF
            ###     ###  ###     ###      #FF      ###     ###            #FF      #FF    ###  ###          ###     ###      #FF
            ###     ###  ###########      #FF      ###########  #######   #FF      #FF######   ###########  ###########      #FF
            ###     ###  ###########      #FF      ###########  #######   #FF      #FF###      ###########  ###########      #FF
            ###     ###  ###     ###      #FF      ###     ###            #FF      #FF####     ###          ###     ###      #FF
            ###     ###  ###     ###      #FF      ###     ###            #FF      #FF  ###    ###          ###     ###      #FF
            ###     ###  ###     ###      #FF      ###     ###            #FF      #FF   ###   ###          ###     ###      #FF
            ###########  ###     ###      #FF      ###     ###            #FF      #FF   ###   ###########  ###     ###      #FF
            #########    ###     ###      #FF      ###     ###            #FF      #FF   ####    #######    ###     ###      #FF  -ment



Quantifycoding=function(simulation=simu, gen=nrow(simulation[[1]]) ,param=paste(simulation[[4]], "/", simulation[[4]], ".param", sep=""), reference="INIT", write=TRUE, method=1, tfreg="unique"){
  stopifnot(reference=="ZERO" || reference=="INIT" || reference=="PARAM") #we can choose here, which reference we want to take to compute the  relevant mutations: from ZERO, from the initial individual (luca - last universal common ancester) or whats written in the param (first value).
#load parameters used for simulation to check mutations
  parameters=simulation[[5]]
  L=nrow(simulation[[2]][[1]]$ind)
  theta = gv("FITNESS_OPTIMUM", parameters)  ; if (length(theta)<L){ for (i in (length(theta)+1):L){theta[length(theta)+1]=theta[1]}} # on verifie que theta a le bon nombre de parametres, sinon on lui rajoute le 1er parametre le nombre de fois qu'il faut
  s1=gv("FITNESS_STRENGTH", parameters) ; if (length(s1)<L){ for (i in (length(s1)+1):L){s1[length(s1)+1]=s1[1]}}
  s2=gv("FITNESS_STABSTR", parameters)  ; if (length(s2)<L){ for (i in (length(s2)+1):L){s2[length(s2)+1]=s2[1]}}
  initall=gv("INIT_ALLELES", parameters)
  inittrans=gv("INIT_TRANSALLELES", parameters)
#prepare data
        if (reference=="ZERO") { ref=matrix(rep(0, L*(L+1)), nrow=L, ncol=(L+1))
  }else if (reference=="PARAM"){ ref=cbind(matrix(rep(initall[1],L*L), nrow = L), c(rep(inittrans[1],L))) #matrix with initial values
  }else if (reference=="INIT") { ref=WFinale(simulation[[1]], 1) #simulation[[2]][[1]]$ind #~a revoir, ce serait bien que ce soit une moyenne de la pop, ici cest le premier ind de initpop
  }
  mat=WFinale(simulation[[1]], gen) #matrix from last gen - or required gen
  moyfit=addfitness(model.M2(mat, tfreg=tfreg), L=L, theta=theta, s1=s1, s2=s2) #ind moyen fitness
#result becomes the matrix containing the difference of fitness w/ or w/o the mutation
  difffit=matrix(NA, nrow=L, ncol =(L+1))
  for(i in 1:L){ for(j in 1:(L+1)){
    mat2=mat
    mat2[i,j]=ref[i,j]
    moyfit2=addfitness(model.M2(mat2, tfreg=tfreg), L=L, theta=theta, s1=s1, s2=s2)
    difffit[i,j]=moyfit$fitness-moyfit2$fitness}}
  difffit[difffit < 0.01] = 0 #get rid of smaller effects
  round(difffit,2)
#reconstruct consensus matrx for testing ie matrix with only relevant mutations determined from diffinit
  relevant=ref
  for(i in 1:L){ for(j in 1:(L+1)){ if(difffit[i,j]!=0){relevant[i,j]=mat[i,j]}}}
  consfit=addfitness(model.M2(relevant, tfreg=tfreg), L=L, theta=theta, s1=s1, s2=s2)$fitness #fitness de l'indiv moyen
  #exit the found results
  if (write==TRUE){cat(paste("Compared to ", reference, ":\n", sum(difffit!=0)," apparent relevant mutations: ", sum(difffit[,1:L]!=0), " CIS-reg and ", sum(difffit[,(L+1)]!=0), " CODING.\n", sep=""))}
  if (sum(simulation[[1]]$MFit>0.9)==0){if (write==TRUE){cat(paste("The max population's fitness is ", round(max(simulation[[1]]$MFit),4), " at generation ", match(max(simulation[[1]]$MFit), simulation[[1]]$MFit), ".", sep=""))}}else{#mean of the individual's fitnesses = of inds of pop
  if (write==TRUE){cat(paste("The population's fitness reaches 0.90 at generation ", simulation[[1]][min(which(simulation[[1]]$MFit>0.9)),]$Gen, ".", sep=""))}}#mean of the individual's fitnesses = of inds of pop
  relevant
  best=sapply(simulation[[3]], function(i) i$fitness) ; best=match(max(best), best) #extract fitnesses; #numero de l'individu avec la max fitness
  if (write==TRUE){
    cat(paste("\nFitness: of initial pop: ", round(simulation[[2]][[1]]$fitness,4), #mean of the individual's fitnesses
            "\n         of inds of pop: ", round(mean(sapply(1:length(simulation[[3]]), function(i) simulation[[3]][[i]]$fitness)),4),#mean of the individual's fitnesses
            "\n        of the mean ind: ", round(moyfit$fitness,4), #fitness of the mean individual, ie, of the mean matrix of the individuals
            "\n        of the best ind: ", round(max(sapply(simulation[[3]], function(i) i$fitness)),4), " [[", best, "]]", #indicate the which individual is the best  [[ ... ]]
            "\n          from relevant: ", round(consfit,4), sep=""))
  }
  output=list()
  output$meanindmat=mat
  output$difffit=difffit
  output$relevant=relevant
  if(method==1){output$CIS=sum(difffit[,1:L]!=0) #Here we differenciate counting cisreg mutations by cases (maximum number of mutations CIS=100, TRANS=10)
    }else if(method==2){output$CIS=sum(rowSums(difffit[,1:L])!=0)} #versus counting along loci (max mut CIS=10, TRANS=10).
  output$TRANS=sum(difffit[,(L+1)]!=0)
  output$fit_MeanINdsPop_IndMoy_Relevant=c(round(mean(sapply(1:length(simulation[[3]]), function(i) simulation[[3]][[i]]$fitness)),4),
                                           round(moyfit$fitness,4),
                                           round(consfit,4))
  return(output)
} # end Quantifycoding

Quantifycistrans = function(simulation=simu, threshold=0.01, total=TRUE){ #--> eQTL Study
#Allows for determination of the effect of mutations for a variation in expression
#What genes are differentially expressed in those two populations?
  L=nrow(simulation[[3]][[1]]$ind)
  a=0 ; genes=c() ; impact=c()
  for (j in 1:L){ #check differences in mean gene expressions between pop initial and finale
    k=abs(mean(sapply(simulation[[2]], function(i) i$mean[j]))-mean(sapply(simulation[[3]], function(i) i$mean[j])))
    impact[length(impact)+1]=k
    if (k>threshold){ #threshold above which it's considered a change in expression
    a=a+1
    genes[length(genes)+1]=j  # a is the counter, genes get the gene numbers differentially expressed
    }
  }
  cat("Changes in expression in ", a, " genes: ", paste(as.character(genes), collapse = " "), ".\n\n", sep="")

  #Now: are this differences in expression due to cis or trans effects?
  init=round(WClean(WFromPop(simulation[[2]])),2)
  final=round(WClean(WFromPop(simulation[[3]])),2)
  initeffect=WConvert(init)
  finaleffect=WConvert(final)
  totalinit=rowSums(abs(initeffect))
  totalfinal=rowSums(abs(finaleffect))

  results=data.frame(matrix(ncol=7, nrow=0)); colnames(results)=c("onGene","fromGene","type","CISvalue", "TRANSvalue", "effectvalue", "impact")
  #CIS/TRANSvalue is relative to the effect on the matrix/of the regulation (positive or negative), final-init.
  #impact: on gene expression (positive) of the mutation cis/trans/cistrans: divisE (selon les effets) si plusieurs mutations impliquees pour le changement d'expression
  #effectvalue: total cis x trans, case de la matrice carree (final-init)
  for (i in genes){
    #print(i)
    for (j in 1:L){
            if (final[i,j]==0 && init[i,j]==0){
          #cat("00000\n")
      }else if (final[i,j]==init[i,j] && final[j,(L+1)]!=init[j,(L+1)]){
          #cat("trans- effect", finaleffect[i,j]/totalfinal[i]*impact[i], "\n")
          xrow=nrow(results)+1 ; results[xrow,] = NA #adds a new line to dataframe #xrow est la last row
          results$onGene[xrow]=i ; results$fromGene[xrow]=j ; results$type[xrow]="trans"
          results$impact[xrow]=abs(finaleffect[i,j])/totalfinal[i]*impact[i]
          results$effectvalue[xrow]=finaleffect[i,j]-initeffect[i,j]
          results$TRANSvalue[xrow]=final[j,(L+1)]-init[j,(L+1)] ; results$CISvalue[xrow]=0
      }else if (final[i,j]!=init[i,j] && final[j,(L+1)]==init[j,(L+1)]){
          #cat("cis - effect", finaleffect[i,j]/totalfinal[i]*impact[i], "\n")
          xrow=nrow(results)+1 ; results[xrow,] = NA
          results$onGene[xrow]=i ; results$fromGene[xrow]=j ; results$type[xrow]="cis"
          results$impact[xrow]=abs(finaleffect[i,j])/totalfinal[i]*impact[i]
          results$effectvalue[xrow]=finaleffect[i,j]-initeffect[i,j]
          results$TRANSvalue[xrow]=0 ; results$CISvalue[xrow]=final[i,j]-init[i,j]
      }else if (final[i,j]!=init[i,j] && final[j,(L+1)]!=init[j,(L+1)]){
          #cat("cistrans- effect", finaleffect[i,j]/totalfinal[i]*impact[i], "\n")
          xrow=nrow(results)+1 ; results[xrow,] = NA
          results$onGene[xrow]=i ; results$fromGene[xrow]=j ; results$type[xrow]="cistrans"
          results$impact[xrow]=abs(finaleffect[i,j])/totalfinal[i]*impact[i]
          results$effectvalue[xrow]=finaleffect[i,j]-initeffect[i,j] #OK
          results$TRANSvalue[xrow]=final[j,(L+1)]-init[j,(L+1)] ; results$CISvalue[xrow]=final[i,j]-init[i,j]
      }
    }
  }

  if (total==FALSE){return(results)}else{
    total=data.frame(matrix(ncol=5, nrow=0))
    for (i in genes){
      isolated=results[results$onGene==i,]
      carry=c(i, colSums(cbind(isolated$CISvalue, isolated$TRANSvalue, isolated$effectvalue, isolated$impact))) #total de lignes de matrices
      total=rbind(total, carry)
    }
    colnames(total)=c("onGene","CISvalue", "TRANSvalue", "effectvalue", "impact")
    return(total)
  }
} #end Quantifycistrans


QuantifyExpressionChange=function(meanind.=meanind, mut.=mut, threshold=0.01, write=TRUE, tfreg="unique"){
  #threshold above which it's considered a change in expression
  #Need to be in an environment where params have been imported
  diff=meanind-mut ; diff2=as.vector(which(diff!=0, arr.ind=TRUE)) #localize the mutation
  if (nrow(which(diff!=0, arr.ind=TRUE))>1){cat("Error: more than one mutation.\n")}

  meanindmeans = addfitness(model.M2(meanind, tfreg=tfreg), L=L, theta=theta, s1=s1, s2=s2)$mean
  mutmeans = addfitness(model.M2(mut, tfreg=tfreg), L=L, theta=theta, s1=s1, s2=s2)$mean
  a=0 ; genes=c() ; impact=c()
  for (i in 1:L){ #check differences in mean gene expressions between pop initial and finale
    k=abs(mutmeans[i]-meanindmeans[i])
    impact[length(impact)+1]=k
    if (k>threshold){a=a+1 ; genes[length(genes)+1]=i} #a is the counter, genes get the gene numbers differentially expressed
  }
  if (write==TRUE){cat("Changes in expression in ", a, " genes: ", paste(as.character(genes), collapse = " "), ".\n\n", sep="")}

  results=data.frame(matrix(ncol=5, nrow=0)); colnames(results)=c("onGene","fromGene","Location", "Effect", "ExpressionChange")
  for (i in genes){
    results[nrow(results)+1,] = NA #new empty row
    lastrow=nrow(results)
    results[lastrow, "onGene"] = i
    results[lastrow, "fromGene"] = diff2[1]
    if (diff2[2]==(L+1)){results[lastrow, "Location"]="coding"}else{results[lastrow, "Location"]="cisreg"}
    if (i == diff2[1]){results[lastrow, "Effect"]="CIS"}else{results[lastrow, "Effect"]="TRANS"}
    results[lastrow, "ExpressionChange"] = mutmeans[i]-meanindmeans[i]
  }
  return(results)
}






#Shortcuts & co.
squared=function(n){q=n*n;return(q)}
dna=function(n=2500){paste(sample(c("A", "T", "G", "C"), n, replace=TRUE), collapse = '')}
comp=function(seq=dna(10)){seqcomp=gsub("D","G", gsub("G","C", gsub("C","D", gsub("B","T", gsub("T","A", gsub("A","B", seq)))))) ; return(seqcomp)}
rt=function(){par(mfrow=c(1,1))}#just a shortcut pour remettre narmol figures
do=function(){invisible(dev.off())}#just a shortcut for dev.off()
dn=function(){dev.new()}
#http://patorjk.com/software/taag/#p=display&f=Banner3&t=coucou
adj=function(col="", alpha=50){return(adjustcolor(col, alpha.f = alpha*0.01))} #color transparent
myred=adjustcolor("red4", alpha.f = 0.6)
mygrey=adjustcolor("#6f6363", alpha.f = 0.6)
