
simulate <- function(N, loci, esize, afreq, gsize,
                     iter, s.size, epipair, epitype, hset, verbose){
  
  
  #Make a phenotyping function that inputs a variable, and outputs the 
  # sum + (effect size * h)
  # x is a matrix 2 by number of loci
  #need to decide which sites to be epistatic
  phenotyper <- function(x, cur.loci, h, esize,epipair, epitype){
    y <- x[, cur.loci]
    pheno <- 0
    for(i in 1:length(cur.loci)){
      switch((sum(y[, i]) + 1),
             pheno <- pheno,
             pheno <- pheno + esize * h[i],
             pheno <- pheno + esize)
    }
    pheno
  }
  
  
  
  
  htracker <- vector(length=iter)
  results <- vector(length=iter, mode="list")
  for(k in 1:iter){
    h <- NULL
    #set a dominance coefficient for each locus. Randomly assign h values,
    #generated from a uniform distribution between 0 and 1, for every locus
    if(hset=="randunif_0_1"){
      h <- runif(min=0, max=1, loci)
    }
    if(hset=="all_dom"){
      h <- rep(1, loci)
    }
    if(hset=="all_add"){
      h <- rep(0.5, loci)
    }
    htracker[k] <- paste(h, collapse="_")
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
    SpeciesHyb <- list()
    for(i in 1:N){
      #randomly choose N individuals for gametogenesis and random mating (with replacement so that family size is poisson distributed)
      x <- sample(1:N, 1, replace = T)
      #Unite gametes (produced by free recombination) from Species A and Species B
      SpeciesHyb[[i]] <- matrix(c(c(SpeciesA[[x]][1, ], SpeciesA[[x]][2, ])[1:gsize + sample(c(0, gsize), gsize, replace = T)],
                                  c(SpeciesB[[x]][1, ], SpeciesB[[x]][2, ])[1:gsize + sample(c(0, gsize), gsize, replace = T)]),
                                2, gsize, byrow = T)
    }
    
    
    
    #Make a matrix that has N rows and 3 columns
    vals <- matrix(,s.size,3)
    #Name the columns "SpeciesA" and "SpeciesB" and "SpeciesHyb"
    colnames(vals) <- c("SpeciesA","SpeciesB", "SpeciesHyb")
    #Get the phenotype values for N individuals of each species
    counter <- 1
    for(i in sample(1:N, s.size)){
      vals[counter,1] <- phenotyper(SpeciesA[[i]], cur.loci, h, esize)
      vals[counter,2] <- phenotyper(SpeciesB[[i]], cur.loci, h, esize)
      vals[counter,3] <- phenotyper(SpeciesHyb[[i]], cur.loci, h, esize)
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
  dat.plot <- as.data.frame(matrix(NA,1,8))
  colnames(dat.plot) <- c("value","pop","stat","loci","esize","afreq","h","s.size")
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
  }
  dat.plot$loci <- loci
  dat.plot$esize <- esize
  dat.plot$afreq <- paste(afreq, collapse="_")
  dat.plot$s.size <- s.size
  dat.plot$h <- rep(htracker, each=6)
  return(dat.plot)
}
