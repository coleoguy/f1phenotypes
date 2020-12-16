

simulate <- function(N, loci, effect.size, afreq, gsize,
                     iter, s.size, epipair, epi.type, hset, mating, verbose){


  #Make a phenotyping function that inputs a variable, and outputs the
  # sum + (effect size * h)
  # x is a matrix 2 by number of loci
  #need to decide which sites to be epistatic
  cv <- function(x){
    sd(x)/mean(x)
  }
  
  phenotyper <- function(x, cur.loci, h, esize, epipair, epi.type){

    y <- x[, cur.loci]
    
    if(epipair == 0){
      pheno <- 0
      for(i in 1:length(cur.loci)){
        switch((sum(y[, i]) + 1),
               pheno <- pheno,
               pheno <- pheno + esize[i] * h[i],
               pheno <- pheno + esize[i])
      }
    }
    if(epipair > 0){
      epi.loci <- sample(1:ncol(y), size=2*epipair)
      if(epi.type=="addbyadd"){
        epi.impact <- matrix(c(1,0,-1,0,0,0,-1,0,1),9)
      }
      if(epi.type=="addbydom"){
        epi.impact <- matrix(c(-1,1,-1,0,0,0,1,-1,1), 9)
      }
      if(epi.type=="dombyadd"){
        epi.impact <- matrix(c(-1,0,1,1,0,-1,-1,0,1), 9)
      }
      if(epi.type=="dombydom"){
        epi.impact <- matrix(c(-1,1,-1,1,-1,1,-1,1,-1), 9)
      }
      pheno <- 0
      for(i in 1:length(cur.loci)){
        if(! i %in% epi.loci){
          switch((sum(y[, i]) + 1),
                 pheno <- pheno,
                 pheno <- pheno + esize[i] * h[i],
                 pheno <- pheno + esize[i])
        }
      }
      counter <- 0
      for(i in seq(from=1,by=2,length.out=epipair)){
        counter <- counter +1
        cgeno <- y[,c(epi.loci[i],epi.loci[i+1])]
        switch(paste(as.character(colSums(cgeno)),collapse=""),
               "22" = pheno <- pheno + esize[i]*epi.impact[1],
               "21" = pheno <- pheno + esize[i]*epi.impact[2],
               "20" = pheno <- pheno + esize[i]*epi.impact[3],
               "12" = pheno <- pheno + esize[i]*epi.impact[4],
               "11" = pheno <- pheno + esize[i]*epi.impact[5],
               "10" = pheno <- pheno + esize[i]*epi.impact[6],
               "02" = pheno <- pheno + esize[i]*epi.impact[7],
               "01" = pheno <- pheno + esize[i]*epi.impact[8],
               "00" = pheno <- pheno + esize[i]*epi.impact[9])
      }
    }
    return(pheno)
  }

  htracker <- vector(length=iter)
  esize.tracker <- vector(length=iter)
  results <- vector(length=iter, mode="list")
  for(k in 1:iter){
    h <- NULL
    #set a dominance coefficient for each locus. Randomly assign h values,
    #generated from a uniform distribution between 0 and 1, for every locus
    if(hset=="runif"){
      h <- runif(min=0, max=1, loci)
    }
    if(hset=="all_dom"){
      h <- rep(1, loci)
    }
    if(hset=="all_add"){
      h <- rep(0.5, loci)
    }
    if(hset=="half_dom"){
      h <- rep(c(0.5, 1), loci/2) 
    }
    htracker[k] <- paste(h, collapse="_")
    esize <- NULL
    if(mode(effect.size)!= "numeric"){
      if(effect.size=="allequal_stephen"){
        esize <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
      }
      if(effect.size=="onelarge_stephen"){
        esize <- c(0.91, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)
      }
      if(effect.size=="runif"){
        esize <- runif(min=0, max=1, loci)
      }
      if(effect.size=="neg_binom"){
        esizedistro <- rnbinom(n = loci, size = 1, prob = 0.01)
        esize <- esizedistro/sum(esizedistro)
      }
    }
    if(mode(effect.size)=="numeric"){
      esize <- effect.size
    }
    esize.tracker[k] <- paste(esize, collapse="_")
    #Make a list object named "SpeciesA"
    SpeciesA <- vector(length=N, mode="list")
    #Simulate "N" diploid individuals, each individual has 2 chromosomes (rows) of
    #"gsize" loci (columns). All of the values are 0 for this initial matrix.
    for(i in 1:N){
      SpeciesA[[i]] <- matrix(data = 0, nrow = 2, ncol = gsize)
    }

    #Make a list object named "SpeciesB"
    SpeciesB <- vector(length=N, mode="list")
    #Simulate "N" diploid individuals, each individual has 2 chromosomes (rows) of
    # "gsize" loci (columns). All of the values are 0 for this initial matrix.
    for(i in 1:N){
      SpeciesB[[i]] <- matrix(data =0, nrow = 2, ncol = gsize)
    }

    #Take a random sample of loci from the genome and save to object named
    #"cur.loci". These loci are the impacted loci.
    cur.loci <- sample(1:gsize, loci)
    #set major allele frequency for population 1
    p <- afreq[1]
    #make hardy-weinberg genotype frequencies for population 1
    genosA <- c(p^2, p*(1-p), p*(1-p), (1-p)^2)
    #Label the genotypes by allele combinations
    names(genosA) <-  c("11","10","01","00")

    #set major allele frequency for population 2
    p <- afreq[2]
    #make hardy-weinberg genotype frequencies for population 2
    genosB <- c(p^2, p*(1-p), p*(1-p), (1-p)^2)
    #Label the genotypes by allele combinations
    names(genosB) <- c("11","10","01","00")

    #loop through all individuals in a population
    for(i in 1:N){
      #loops through all impacted loci and assigns genotypes
      for(j in 1:loci){

        x <- sample(1:4, 1, prob=genosA)
        switch(x,
               SpeciesA[[i]][,cur.loci[j]] <- c(1,1),
               SpeciesA[[i]][,cur.loci[j]] <- c(1,0),
               SpeciesA[[i]][,cur.loci[j]] <- c(0,1),
               SpeciesA[[i]][,cur.loci[j]] <- c(0,0))
      }
      for(j in 1:loci){
        x <- sample(1:4, 1, prob=genosB)
        switch(x,
               SpeciesB[[i]][,cur.loci[j]] <- c(1,1),
               SpeciesB[[i]][,cur.loci[j]] <- c(1,0),
               SpeciesB[[i]][,cur.loci[j]] <- c(0,1),
               SpeciesB[[i]][,cur.loci[j]] <- c(0,0))
      }
    }


    
    #Simulate a SpeciesA gamete with no recombination. This code first linearizes the chromosomes (c(SpeciesA[[i]][1,],SpeciesA[[i]][2,])). Then the code creates a new vector of length 20, and randomly adds either a 0 or a 20 to each position. This new "position vector" is then used to determine if the gamete will have an element from chromosome 1 or 2.
    # gam1 <- c(SpeciesA[[i]][1,],SpeciesA[[i]][2,])[1:20 + sample(c(0, 20), 20, replace = T)]
    
    #Make a for loop to simulate N hybrid individuals that were produced by gametogenesis with free recombination and random mating (which includes poisson distributed family size)
    #Simulate parent population gametes with free recombination and random mating of these gametes (poisson distribution)
    if(mating=="random"){
    SpeciesHyb <- list()
    for(i in 1:N){
      #randomly choose N individuals for gametogenesis and random mating (with replacement so that family size is poisson distributed)
      x <- sample(1:N, 1, replace = T)
      #Unite gametes (produced by free recombination, or no recombination??????) from Species A and Species B
      SpeciesHyb[[i]] <- matrix(c(c(SpeciesA[[x]][1, ], SpeciesA[[x]][2, ])[1:gsize + sample(c(0, gsize), gsize, replace = T)],
                                  c(SpeciesB[[x]][1, ], SpeciesB[[x]][2, ])[1:gsize + sample(c(0, gsize), gsize, replace = T)]),
                                2, gsize, byrow = T)
    }
}
  if(mating=="assortative"){
    mpheno <- matrix(nrow = N, ncol = 2)
    colnames(mpheno) <- c("SpeciesA", "SpeciesB")
    for(i in 1:N){
      mpheno[i,1] <- phenotyper(SpeciesA[[i]], cur.loci, h, esize, epipair)
      mpheno[i,2] <- phenotyper(SpeciesB[[i]], cur.loci, h, esize, epipair)
    }
    distmat <- matrix(NA, nrow=nrow(mpheno), ncol=nrow(mpheno))
    for(i in 1:nrow(mpheno)){
      for(j in 1:nrow(mpheno)){
        distmat[i,j] <- abs(mpheno[i,1] - mpheno[j,2])
      }
    }
    ass.distmat <- round(1-(distmat/max(distmat)), digits = 3)
    chosen.parents <- matrix(NA,0,2)
    colnames(chosen.parents) <- c("SpeciesA", "SpeciesB")
    parents <- 5
    if(type=="assort"){
      par.index <- sample(1:length(ass.distmat), size = parents,
                          prob = as.vector(ass.distmat))
      counter <- 1
      for(i in 1:ncol(ass.distmat)){
        for(j in 1:nrow(ass.distmat)){
          if(counter %in% par.index){
            chosen.parents <- rbind(chosen.parents, c(j,i))
          }
          counter <- counter + 1
        }
      }
    }
    for(i in 1:N){
    SpeciesHyb[[i]] <- matrix(c(c(SpeciesA[chosen.parents][1, ],
                                  SpeciesA[chosen.parents][2, ])[1:gsize + sample(c(0, gsize), gsize, replace = T)],
                                c(SpeciesB[chosen.parents][1, ],
                                  SpeciesB[chosen.parents][2, ])[1:gsize + sample(c(0, gsize), gsize, replace = T)]),
                              2, gsize, byrow = T)
    }
    
    
    
      par.index <- sample(1:length(ass.distmat), size = parents,
                        prob = as.vector(ass.distmat))
      counter <- 1
      for(w in length(ass.distmat)){
        for(i in 1:ncol(ass.distmat)){
          for(j in 1:nrow(ass.distmat)){
           if(counter %in% par.index){
            chosen.parents <- rbind(chosen.parents, c(j,i))
          }
        }
      }
      counter <- counter + 1
    }
    counter <- 1
    for(i in 1:ncol(ass.distmat)){
      for(j in 1:nrow(ass.distmat)){
        if(counter %in% par.index){
          chosen.parents <- rbind(chosen.parents, c(j, i))
        }
        
    
    
    for(i in 1:N){
      x <- sample(distmat, size = 1, replace = T, prob=distmat)
      SpeciesHyb[[i]] <- matrix(c(c(SpeciesA[[x]][1, ],SpeciesA[[x]][2, ])[1:gsize + sample(c(0, gsize), gsize, replace = T)],
                                  c(SpeciesB[[x]][1, ], SpeciesB[[x]][2, ])[1:gsize + sample(c(0, gsize), gsize, replace = T)]),
                                2, gsize, byrow = T)
    }
  }
    }}    
    
    
    #Make a matrix that has N rows and 3 columns
    vals <- matrix(,s.size,3)
    #Name the columns "SpeciesA" and "SpeciesB" and "SpeciesHyb"
    colnames(vals) <- c("SpeciesA","SpeciesB", "SpeciesHyb")
    #Get the phenotype values for N individuals of each species
    counter <- 1



    for(i in sample(1:N, s.size)){
      vals[counter,1] <- phenotyper(SpeciesA[[i]], cur.loci, h, esize, epipair, epi.type)
      vals[counter,2] <- phenotyper(SpeciesB[[i]], cur.loci, h, esize, epipair, epi.type)
      vals[counter,3] <- phenotyper(SpeciesHyb[[i]], cur.loci, h, esize, epipair, epi.type)
      counter <- counter + 1
    }
    vals <- as.data.frame(vals)
    results[[k]] <- vals
    #Also want to print columns for number of loci and average dominance???
    if(verbose==T) print(k)
  }
  if(verbose==T) print(paste("recording results of last",iter, "simulations"))
  # get stuff to plot
  # mean phenotype from each iteration
  # and the variance from each iteration
  dat.plot <- as.data.frame(matrix(NA,1,10))
  colnames(dat.plot) <- c("value","pop","stat","loci","esize","afreq","h",
                          "epipair","s.size", "epi.type")
  counter <- 1
  for(i in 1:length(results)){
    dat.plot[counter, 1] <- mean(results[[i]]$SpeciesA)
    dat.plot[counter, 2] <- "popA"
    dat.plot[counter, 3] <- "mean"
    counter <- counter +1
    dat.plot[counter, 1] <- var(results[[i]]$SpeciesA)
    dat.plot[counter, 2] <- "popA"
    dat.plot[counter, 3] <- "var"
    counter <- counter +1
    dat.plot[counter, 1] <- mean(results[[i]]$SpeciesB)
    dat.plot[counter, 2] <- "popB"
    dat.plot[counter, 3] <- "mean"
    counter <- counter +1
    dat.plot[counter, 1] <- var(results[[i]]$SpeciesB)
    dat.plot[counter, 2] <- "popB"
    dat.plot[counter, 3] <- "var"
    counter <- counter +1
    dat.plot[counter, 1] <- mean(results[[i]]$SpeciesHyb)
    dat.plot[counter, 2] <- "popH"
    dat.plot[counter, 3] <- "mean"
    counter <- counter +1
    dat.plot[counter, 1] <- var(results[[i]]$SpeciesHyb)
    dat.plot[counter, 2] <- "popH"
    dat.plot[counter, 3] <- "var"
    counter <- counter +1
    dat.plot[counter, 1] <- (mad(results[[i]]$SpeciesHyb))/(mad(results[[i]]$SpeciesA))
    dat.plot[counter, 2] <- "allpops"
    dat.plot[counter, 3] <- "MADratio_hybtoA"
    counter <- counter +1
    dat.plot[counter, 1] <- (mad(results[[i]]$SpeciesHyb))/(mad(results[[i]]$SpeciesB))
    dat.plot[counter, 2] <- "allpops"
    dat.plot[counter, 3] <- "MADratio_hybtoB"
    counter <- counter +1
    dat.plot[counter, 1] <- if(!is.na(((mad(results[[i]]$SpeciesHyb))/(mad(results[[i]]$SpeciesA))<1)&((mad(results[[i]]$SpeciesHyb))/(mad(results[[i]]$SpeciesB))<1))){
      print("SmallerVariance")
    }else{print("LargerVariance")}
    dat.plot[counter, 2] <- "allpops"
    dat.plot[counter, 3] <- "MADratio_Smallerthanbothparents"
    counter <- counter +1
  }
  dat.plot$loci <- loci
  dat.plot$esize <- rep(esize.tracker, each=9)
  dat.plot$afreq <- paste(afreq, collapse="_")
  dat.plot$s.size <- s.size
  dat.plot$h <- rep(htracker, each=9)
  dat.plot$epipair <- epipair
  dat.plot$epi.type <- epi.type
  return(dat.plot)
}




