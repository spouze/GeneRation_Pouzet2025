##### IndivRecap_fromSET_FUN.R
IndivRecap = function(storeLasts_set = dream_SET, line = 1, clean=T, cleanthresh=0.01){
  
  selected_line = storeLasts_set[line,]
  new_ind = list()
  new_ind$id = selected_line$id
  new_ind$topo = selected_line$topo
  
  
  # Load specific parameters
  L     = 10
  theta = c() ; for (i in 1:L){theta = c(theta, selected_line[[paste("FIT_OPT",  i, sep="")]])}
  s1    = c() ; for (i in 1:L){s1    = c(s1   , selected_line[[paste("FIT_STR",  i, sep="")]])}
  s2    = rep(46000, L)
  
  # Record parameters in ind.
  new_ind$L     = L
  new_ind$theta = theta
  new_ind$s1    = s1
  new_ind$s2    = s2
  
  #get matrix
  newW=c()
  for (i in 1:(L*L)){newW = c(newW, selected_line[[paste("MeanAll", i, sep="")]])} # get interaction values
  for (i in 1:L)    {newW = c(newW, selected_line[[paste("MTrans",  i, sep="")]])} # add coding
  newW=matrix(newW, nrow=L)
  new_ind$mom = newW
  new_ind$dad = newW
  new_ind$ind = newW
  
  # Compute phenotype
  dev = model.M2(newW, steps=20, measure=2, tfreg="notunique!", full=T)
  new_ind$mean = dev$mean
  new_ind$var  = dev$var
  
  new_ind$fitness = addfitness(new_ind, L=L, theta=theta, s1=s1, s2=s2)$fitness
  
  if (clean) {
    #get clean matrix
    cleaned_newW = WClean(newW, threshold = cleanthresh)
    new_ind$indclean = cleaned_newW
    #get corresponding dev
    dev = model.M2(cleaned_newW, steps=20, measure=2, tfreg="notunique!", full=T)
    new_ind$meanclean = dev$mean
    new_ind$varclean  = dev$var
    #to get the fitness, need to create an ind intermediate
    intermediate = list()
    intermediate$mean = dev$mean
    intermediate$var = dev$var
    #Now get fitness thanks to this guy:
    new_ind$fitnessclean = addfitness(intermediate, L=L, theta=theta, s1=s1, s2=s2)$fitness
    }
  
  return(new_ind)
  
}

########### Comfort_FUN.R

play = function(file=3){
  if (file==3){system("afplay /Users/sylvain/Music/snap.mp3")}
  if (file==1){system("afplay /Users/sylvain/Music/demisnap.mp3")}
  if (file==2){system("afplay /Users/sylvain/Music/demisnap2.mp3")}
  if (file=="mario"){system("afplay /Users/sylvain/Music/mario_end.mp3")}
  }

topologies = c("HIGHCO", "RANDOM", "SCALEF")
topocolors = c("red3", "green4", "blue4")

colrs = c("green4", "blue4", "red3")
# colrs[match("SCALEF", topologies)]

resetplot = function(){
  dev.off()
  dev.new()
  cex_mtext <<- 1
  par(mgp=c(3, 1, 0))
}

assign_category <- function(value, thresholds=effect_class_tresh, names=effect_class) {
  index <- findInterval(value, thresholds, rightmost.closed = T)
  return(names[index + 1])  # Adding 1 since R indexing starts at 1
}

check_folder_exists = function(folder_to_check){ # Check if the folder exists
  if (!dir.exists(folder_to_check)) { 
    dir.create(folder_to_check, recursive = TRUE)
    cat("Folder created at:", folder_to_check)
  } else{
    cat("Folder exists.")
  }
}

draw_box_check = function(){
  if (draw_box){
    box("figure", col="forestgreen")
    box(which="plot", col="red")
    box("outer", col="blue")
    }
}

restore_mgp = function(){
  par(mgp=c(3, 1, 0))
}

cmtoin=2.54 # 1 in = 2.54 cm

make_darker_color = function(color="#B3E4FF", amount=30){
  return(rgb(pmax(col2rgb(color) - amount, 0)[1,], pmax(col2rgb(color) - amount, 0)[2,], pmax(col2rgb(color) - amount, 0)[3,], maxColorValue = 255))
}


################## 128.New_Mutate_in_new_Env_4.R

#copy of 126 # ici on va ajouter les tests pour ce qui est de cis/trans
# DUPLICATED 128 (safeguard) to handle the cis/trans for DUP DEL (ignore cis/trans in those conditions actualyl)

# IF USED AS STANDALONE
#setwd("/Users/sylvain/Documents/DOCUMENTS/GeneRation/M2_PAPER")
#source("../GeneRation_Fun.R")
#source("R_SCRIPTS/Comfort_FUN.R")
#source("R_SCRIPTS/IndivRecap_fromSET_FUN.R")


#dream_SET = read.csv("Extracted_data/dream_SIMUSET_Lasts.csv", header = T)[,-1]

cnames = c("id", "name", "topology", 
           "mut_type", "i", "j", 
           "mutvalue_init", "mutvalue_final", "mutvalue_delta",
           "dupdel_number", "dupdel_genes",
           "fit_pop", "fit_indE1", "fit_initE2", "fit_mutE2", "fit_delta", "fit_effect",
           "affected_genes", "affected_genes_tresh", 
           "reg_muteff", "cod_muteff",
           "notes", 
           "gene_rand", "dExpr_rand", "Expr_effect_rand", 
           "gene_min" , "dExpr_min" , "Expr_effect_min" ,
           "gene_max" , "dExpr_max" , "Expr_effect_max" , 
           "dExpr_mean", "Expr_effect_cisness")

# allmutations = data.frame(matrix(ncol = length(cnames), nrow = 0))
# colnames(allmutations)=cnames

# cEnames=c("id", "onGene", "fromGene", "Location",
#           "Effect", "ExpressionChange", "fit_effect")

# addfitness2 = function(ind, theta, s1, s2) {
#   fitmean = 0 ; fitvar = 0 ; L = nrow(ind$ind)
#   for (i in 1:L) {
#     fitmean = fitmean + (-s1[i]*(squared(ind$mean[i]-theta[i])))
#     fitvar  = fitvar  + (-s2[i]*(ind$var[i]))}
#   fitness = exp(fitvar)*exp(fitmean)
#   return(fitness)
# }

####################################################################
####################################################################
# for TEST
# dream_SET = read.csv("Extracted_data/dream_SIMUSET_Lasts.csv")[,-1]
# ind = IndivRecap()
####################################################################

Mutate_in_new_environment=function(ind, 
                                   mut_type, #"REG" "COD" "DUP" "DEL"
                                   cod_muteff = NA,
                                   reg_muteff = NA,
                                   dupdel_number = NA, #number of genes to be duplicated or deleted
                                   fit_thresh = 0.001,
                                   affected_genes_tresh = 0.00037, 
                                   id = NA){
  
  # check arguments
  if (is.na(match(mut_type, c("REG", "COD", "DUP", "DEL")))){
    print(paste0("mut_type should be REG COD DUP or DEL. You wrote: ", mut_type))
    stop()}

  
  ### 00 Load Parameters
  ##############################################
  name     = ind$id
  topology = ind$topo
  fit_pop  = ind$fitness
  notes    = ""
  
  L       = ind$L
  theta   = ind$theta
  s1      = ind$s1
  s2      = ind$s2
  diag    = 0 # important, such that if 0, i and j cannot be equal for reg mut. 
  
  ### 04 Change its environment # inspired from the burnin script.
  ##############################################
  newtheta=theta
  n = 2 #We change the optima of 2 of the 5 genes
  selected_genes = sample(c(1:5),n) # pour pas galerer avec moins de genes sous selection, mettre le 5 à 3, mais bon. 
  for (i in selected_genes){
    optimum=theta[i]
    s=s1[i]
    d=sqrt(1/(n*s)*log(2))
    newoptimum = optimum + d*sample(c(1, -1), 1)    # soit on ajoute, soit on retire
    if (newoptimum < 0){newoptimum = newoptimum + 2*d} # if results<0, on inverse
    if (newoptimum > 1){newoptimum = newoptimum - 2*d} # if results>1, on inverse
    newtheta[i] = newoptimum
  }
  fit_indE1  = addfitness(ind, L=L, theta=theta,    s1=s1, s2=s2)$fitness
  fit_initE2 = addfitness(ind, L=L, theta=newtheta, s1=s1, s2=s2)$fitness
  # n ; fit_indE1 ; fit_initE2 ; newtheta==theta # to check
  
  ### 05 Do the mutation
  ##############################################
  clone_ind=list(ind = ind$ind)
  
  
  if(mut_type == "REG"){ # ############################################################## 1. REGULATORY mutation
    if (is.na(reg_muteff)){print("REG mutation: reg_muteff MUST be specified.") ; stop()} #check
    if (reg_muteff==0){notes="NO MUTATION"}
    if(diag==0){
      ij = sample(c(1:L), 2) # prevents i = j (so we keep diags at zero)
      i  = ij[1] # gene i touché
      j  = ij[2] # at binding site of gene j
    }else{
      i = sample(c(1:L), 1)
      j = sample(c(1:L), 1)
    }
    clone_ind$ind[i,j] = clone_ind$ind[i,j]+reg_muteff
    # clone_ind$ind == ind$ind # check
    mutvalue_init  = ind$ind[i,j]
    mutvalue_final = clone_ind$ind[i,j]
    mutvalue_delta = mutvalue_final-mutvalue_init
    cod_muteff = NA
    dupdel_number = NA
    dupdel_genes = NA
  }# end REG mutation
  
  if(mut_type == "COD"){ # ############################################################## 2. CODING mutation
    if (is.na(cod_muteff)){print("COD mutation: cod_muteff MUST be specified.") ; stop()} #check
    if (cod_muteff==0){notes="NO MUTATION"}
    i = sample(c(1:L), 1) # gene i touché
    j = 11 # coding value
    clone_ind$ind[i,j] = clone_ind$ind[i,j]+cod_muteff
    mutvalue_init  = ind$ind[i,j]
    mutvalue_final = clone_ind$ind[i,j]
    mutvalue_delta = mutvalue_final-mutvalue_init
    reg_muteff = NA
    dupdel_number = NA
    dupdel_genes = NA
  }# end COD mutation
  
  if(mut_type == "DUP"){ # ############################################################## 3. DUPLICATION
    if (is.na(dupdel_number)){print("DUPDEL mutation: dupdel_number MUST be specified.") ; stop()} #check
    if (dupdel_number==0){notes="NO MUTATION"}
    dupdel_genes = sample(c(1:L), dupdel_number)
    for (gene in dupdel_genes){
      starting_L = nrow(clone_ind$ind)
      clone_ind$ind = cbind(clone_ind$ind[,1:starting_L], clone_ind$ind[,gene], clone_ind$ind[,starting_L+1]) # add column
      clone_ind$ind = rbind(clone_ind$ind, clone_ind$ind[gene,]) # add row
      newtheta = c(newtheta,0)
      s1 = c(s1, 0)
      s2 = c(s2, 0)
    }
    dupdel_genes = paste(as.character(dupdel_genes), collapse = "")
    i = NA
    j = NA
    mutvalue_init  = NA
    mutvalue_final = NA
    mutvalue_delta = NA
    reg_muteff = NA
    cod_muteff = NA
  }# end DUP mutation
  
  if(mut_type == "DEL"){ # ############################################################## 4. DELETION
    if (is.na(dupdel_number)){print("DUPDEL mutation: dupdel_number MUST be specified.") ; stop()} #check
    if (dupdel_number>5){print("DUPDEL > 5, We cannot remove genes that are under selection.") ; stop()} #check
    if (dupdel_number==0){dupdel_genes = 0 ; notes="NO MUTATION"
    }else{
      dupdel_genes = sample(c(6:10), dupdel_number) #only genes not under selection
      clone_ind$ind = clone_ind$ind[-dupdel_genes, -dupdel_genes]
      newtheta = newtheta[-dupdel_genes]
      s1 = s1[-dupdel_genes]
      s2 = s2[-dupdel_genes]
      dupdel_genes = paste(as.character(dupdel_genes), collapse = "")
    }
    i = NA
    j = NA
    mutvalue_init  = NA
    mutvalue_final = NA
    mutvalue_delta = NA
    reg_muteff = NA
    cod_muteff = NA
  }# end DEL mutation
  
  ### 06 Consequences of mutation
  ##############################################
  dev = model.M2(clone_ind$ind, steps=20, measure=2, tfreg="notunique!", full=T)
  clone_ind$mean    = dev$mean
  clone_ind$var     = dev$var
  clone_ind$fitness = addfitness(clone_ind, L=nrow(clone_ind$ind), theta=newtheta, s1=s1, s2=s2)$fitness
  
  fit_mutE2 = clone_ind$fitness
  fit_delta = fit_mutE2 - fit_initE2
  
  if(abs(fit_delta) > fit_thresh) {
    if(fit_delta > 0){fit_effect = "BENEF"}
    if(fit_delta < 0){fit_effect = "DELET"}
  } else {            fit_effect = "NEUTR"}
  
  # Nombre de genes affectés
  
  affected_genes = which(abs(ind$mean - clone_ind$mean[1:L]) > affected_genes_tresh)
  pleiotropy = length(affected_genes)
  
  # # Contribution to instability? 
  # fitmean = 0 ; fitvar = 0
  # for (gene in 1:L) {
  #   fitmean = fitmean + (-s1[gene]*(squared(clone_ind$mean[gene]-newtheta[gene])))
  #   fitvar  = fitvar  + (-s2[gene]*(clone_ind$var[gene]))
  # }
  # fitmean = exp(fitmean)
  # fitvar  = exp(fitvar)
  # fitfinal = fitvar*fitmean
  
  ## INFOS for CIS / TRANS Effects
  # fromGene = i
  # onGene   = gene
  
  if (pleiotropy == 0 || mut_type == "DUP" || mut_type == "DEL") {
    gene_rand = dExpr_rand = Expr_effect_rand = NA 
    gene_min  = dExpr_min  = Expr_effect_min  = NA
    gene_max  = dExpr_max  = Expr_effect_max  = NA 
    dExpr_mean = Expr_effect_cisness = NA
  } else {
    #Random
    gene_rand = sample(affected_genes, 1)
    dExpr_rand = ind$mean[gene_rand] - clone_ind$mean[gene_rand]
    if ( gene_rand == i){Expr_effect_rand = "cis"}else{Expr_effect_rand = "trans"}
    #min
    gene_min = affected_genes[which.min(abs(ind$mean[affected_genes] - clone_ind$mean[affected_genes]))]
    dExpr_min = ind$mean[gene_min] - clone_ind$mean[gene_min]
    if ( gene_min == i){Expr_effect_min = "cis"}else{Expr_effect_min = "trans"}
    #max
    gene_max = affected_genes[which.max(abs(ind$mean[affected_genes] - clone_ind$mean[affected_genes]))]
    dExpr_max = ind$mean[gene_max] - clone_ind$mean[gene_max]
    if ( gene_max == i){Expr_effect_max = "cis"}else{Expr_effect_max = "trans"}
    #mean
    dExpr_mean = mean(abs(ind$mean[affected_genes] - clone_ind$mean[affected_genes])) # ABSOLUTE VALUE
    if (sum(affected_genes==i)==1){ # si on a un cis parmi les effets
      Expr_effect_cisness = 1/pleiotropy
    }else{Expr_effect_cisness = 0}
    # c'est donc un score de 0 à 1 ici. ...? 
  }
  
  return(
    data.frame(id, name, topology, 
               mut_type, i, j, 
               mutvalue_init, mutvalue_final, mutvalue_delta,
               dupdel_number, dupdel_genes,
               fit_pop, fit_indE1, fit_initE2, fit_mutE2, 
               fit_delta, fit_effect, fit_thresh,
               pleiotropy, affected_genes_tresh, 
               reg_muteff, cod_muteff,
               notes,
               #fitmean, fitvar, fitfinal
               gene_rand, dExpr_rand, Expr_effect_rand, 
               gene_min , dExpr_min , Expr_effect_min ,
               gene_max , dExpr_max , Expr_effect_max , 
               dExpr_mean, Expr_effect_cisness
               ) 
  )#end return
  
}# end Mutate_in_new_environment

#####################
#############################
#############################################
# 
# selected_simu = 1001 #index indream_SET
# id = dream_SET[selected_simu,]$id
# ind = IndivRecap(storeLasts_set = dream_SET, line = selected_simu, clean = F)
# # ind
# 
# allmutations = data.frame(matrix(ncol = length(cnames), nrow = 0))
# colnames(allmutations)=cnames
# 
# allmutations = rbind(allmutations,
#                      Mutate_in_new_environment(ind = ind, mut_type = "REG", reg_muteff = 0, id =1),
#                      Mutate_in_new_environment(ind = ind, mut_type = "COD", cod_muteff = 0, id =1),
#                      Mutate_in_new_environment(ind = ind, mut_type = "DUP", dupdel_number = 0, id =1),
#                      Mutate_in_new_environment(ind = ind, mut_type = "DEL", dupdel_number = 0, id =1))
# 
# Mutate_in_new_environment(ind = ind, mut_type = "REG", reg_muteff = 0.5, id =1)
# Mutate_in_new_environment(ind = ind, mut_type = "COD", cod_muteff = 0.5, id =1)
# Mutate_in_new_environment(ind = ind, mut_type = "DUP", dupdel_number = 2, id =1)
# Mutate_in_new_environment(ind = ind, mut_type = "DEL", dupdel_number = 4, id =1)
# 
# for (check in 1:20){
#   allmutations = rbind(allmutations, Mutate_in_new_environment(ind = ind, mut_type = "REG", reg_muteff = 0.5, id=check))
# }
# 
# allmutations = read.csv("Extracted_data/dream_allmutations.csv")[,-1]



# Assuming pleiotropy and mut_type are predefined variables
# and 'do that' is a placeholder for the actual operation you want to perform



