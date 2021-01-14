setwd("~/Documents/GitHub/f1phenotypes")
source("functions.R")
library(ggplot2)
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
#For some reason I'm getting lot's of NaN's, especially for Hybrids.
is.na(newoutput$stat)

min(newoutput$stat=="var")

cvratiohist <- ggplot(subset(newoutput, stat %in% c("CVratio")), aes(x=value)) + geom_histogram()
cvratiohist

newoutput[is.na(newoutput)] <- 0

point5point5 <- newoutput[newoutput$afreq=="0.5_0.5", ]
point6point4 <- newoutput[newoutput$afreq=="0.6_0.4", ]
point7point3 <- newoutput[newoutput$afreq=="0.7_0.3", ]
point8point2 <- newoutput[newoutput$afreq=="0.8_0.2", ]
point9point1 <- newoutput[newoutput$afreq=="0.9_0.1", ]
fixedandlost <- newoutput[newoutput$afreq=="1_0", ]


scatter1 <- ggplot(subset(point5point5, stat %in% c("CV")), aes(x = esize, y = value))
scatter1 + geom_boxplot(aes(color = pop)) + facet_wrap( ~ epi.type) + ggtitle(label = "point5point5")

scatter2 <- ggplot(subset(point6point4, stat %in% c("CV")), aes(x = esize, y = value))
scatter2 + geom_boxplot(aes(color = pop)) + facet_wrap( ~ epi.type) + ggtitle(label = "point6point4")

scatter3 <- ggplot(subset(point7point3, stat %in% c("CV")), aes(x = esize, y = value))
scatter3 + geom_boxplot(aes(color = pop)) + facet_wrap( ~ epi.type) + ggtitle(label = "point7point3")

scatter4 <- ggplot(subset(point8point2, stat %in% c("CV")), aes(x = esize, y = value))
scatter4 + geom_boxplot(aes(color = pop)) + facet_wrap( ~ epi.type) + ggtitle(label = "point8point2")

scatter5 <- ggplot(subset(point9point1, stat %in% c("CV")), aes(x = esize, y = value))
scatter5 + geom_boxplot(aes(color = pop)) + facet_wrap( ~ epi.type) + ggtitle(label = "point9point1")

scatter6 <- ggplot(subset(fixedandlost, stat %in% c("CV")), aes(x = esize, y = value))
scatter6 + geom_boxplot(aes(color = pop)) + facet_wrap( ~ epi.type) + ggtitle(label = "fixed_lost")

                   
#Jamie Notes Need to DO
#Jamie needs to figure out plotting
#Calculate a variance difference variable
#Regression tree + random forest to see what minimizes variance of F1s
#



