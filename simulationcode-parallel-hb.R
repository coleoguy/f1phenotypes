setwd("~/Documents/GitHub/f1phenotypes")
source("functions.R")
library(doMC)
registerDoMC(cores = 6)


N <- 50
loci <- 10
gsize <- 20
s.size <- 50
iter <- 10
verbose <- F
mating <- "random"


#create lists of parameter values
afreq <- list(c(0.5,0.5), c(0.6, 0.4), c(0.7, 0.3), c(0.8, 0.2), c(0.9,0.1), c(1,0))
hset <- c("all_add", "all_dom", "half_dom")
esize <- list(c(0.91, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
              rep(0.1, 10))
epi.type <- c("noepi", "addbyadd", "addbydom", "dombyadd", "dombydom")
epipair <- c(0,1,2,5)

#Here i am setting up my output structure

x <- foreach(i=1:length(afreq), .combine="c") %dopar% { #this loops through the afreqs, and parallelizes by afreq
#for(i in 1:length(afreq)){ #this loops through allele frequencies
  output <- list()
  counter <- 1
  for(j in 1:length(hset)){ #this loops through different hsets
    for(k in 1:length(esize)){ #this loops through different esizes
      print(counter)
      for(m in 1:length(epipair)){ #this loops through different epipairs
        if(epipair[m]==0){ #if epipair is 0, then set an arbitrary epi.type and record the output
          output[[counter]] <- simulate(N = N, loci = loci, effect.size = esize[[k]],
                                        afreq = afreq[[i]], gsize = gsize,
                                        iter = iter, s.size = s.size, epipair = epipair[m],
                                        epi.type = epi.type[1], hset = hset[j], mating = mating,
                                        verbose = verbose)
          names(output)[counter] <- paste("freqs=", paste(afreq[[i]],collapse="_"), 
                                          "h=", hset[[j]], 
                                          "esize=", paste(esize[[k]], collapse = "_"), 
                                          "epipair=", epipair[[m]], sep="")
          counter <- counter +1
        }
        if(epipair[m]>0){ #if epipair is greater than 0, then record the output
          for(n in 2:length(epi.type)){ #this loops through different epi.types
            output[[counter]] <- simulate(N = N, loci = loci, effect.size = esize[[k]],
                                          afreq = afreq[[i]], gsize = gsize,
                                          iter = iter, s.size = s.size, epipair = epipair[m],
                                          epi.type = epi.type[n], hset = hset[j], mating = mating,
                                          verbose = verbose)
            names(output)[counter] <- paste("freqs=", paste(afreq[[i]],collapse="_"), 
                                            "h=", hset[[j]], 
                                            "esize=", paste(esize[[k]], collapse = "_"), 
                                            "epipair=", epipair[[m]],
                                            "epitype=", epi.type[n], sep="")
            counter <- counter +1
          }
        }
      }
    }
  }
  output
}

newoutput <- do.call(rbind, x)

scatter1 <- ggplot(subset(newoutput, stat %in% c("var")), aes(x = afreq, y = value))
scatter1 + geom_path(aes(color = esize)) + facet_wrap( ~ pop)

                   
#Jamie Notes Need to DO
#Jamie needs to figure out plotting
#Calculate a variance difference variable & coefficient of variance
#Figure out epistasis addbyadd, maybe add another variable?
#Regression tree + random forest to see what minimizes variance of F1s
#



