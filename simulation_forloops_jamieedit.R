
source("functions.R")

#set variables
N <- 50
gsize <- 20
s.size <- 50
iter <- 1
verbose <- F
mating <- "random"
baselevelpheno <- 10
afreq <- list(c(0.5,0.5), c(0.6, 0.4), c(0.7, 0.3), c(0.8, 0.2), c(0.9,0.1), c(1,0))
hset <- c("all_add", "all_dom", "half_dom")
epi.type <- c("noepi", "addbyadd", "addbydom", "dombyadd", "dombydom")
loci <- matrix(c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
total_epi <- 5
epipair <- list()
esize <- list()
esize_test <- list()
for(i in 1:length(loci)){
  epipair[[i]] <- round(seq.int(from = 0, to = loci[i]/2, length.out = total_epi), digits = 0)
  for(k in 1:10){
    vec <- sample(x = 10000, size = loci[i], replace = TRUE)
    newvec <- vec/sum(vec)
    esize_test[[k]] <- newvec
  }
  esize[[i]] <- esize_test
}
output <- list()

for(i in 1:length(loci)){
  for(j in 1:length(esize)){
    for(k in 1:length(epipair)){
      for(m in 1:length(hset)){
        for(n in 1:length(afreq)){
          for(o in 1:length(epi.type)){
            output[[i]] <- simulate(N = N, loci = loci[i], effect.size = esize[[k]],
                                    afreq = afreq[n], gsize = gsize, iter = iter,
                                    epi.type = epi.type[o], hset = hset[m],
                                    mating = mating, baselevelpheno = baselevelpheno,
                                    verbose = verbose)
          }
        }
      }
    }
  }
}



