getwd()
setwd("~/Documents/GitHub/f1phenotypes")
source("functions.R")
# set model parameters
#set population size
N <- 50
#set number of impacted loci (QTL)
loci <- 3
# set effect size this is the difference in phenotype between the two
# alternative homozygous genotypes
esize <-runif(n=loci, min=0, max=1)
# set allele frequency of major allele in 2 parent species. First number is
# SpeciesA, second number is SpeciesB
afreq <- c(.9,.1)
#set genome size
gsize <- 10
#set the number of simulations
iter <- 100

# sample size equal for parents and hybrid pops
s.size <- 50
# number of epistatic gene pairs
epipair <- 2
hset <- "randunif_0_1"
hset <- "all_dom"
hset <- "all_add"

test.results <- list()
for(i in 2:3){
  # set effect size this is the difference in phenotype between the two
  # alternative homozygous genotypes
  esize <-runif(n=loci, min=0, max=1)
  episize <-runif(n=epipair, min=20, max=20)
  test.results[[i]] <- simulate(N, loci, esize,
                                afreq, gsize,
                                iter, s.size, epipair, episize,
                                hset, verbose=T)
}


library(ggplot2)
foo <- test.results[[10]][test.results[[10]]$stat=="var",]
a <- ggplot(foo, aes(x =value, fill=pop))
a + geom_density()
t.test(foo$value[foo$pop=="popA"])
t.test(foo$value[foo$pop=="popB"])
t.test(foo$value[foo$pop=="popH"])
fit <- glm(foo$value~foo$pop)
summary(fit)
foo <-test.results[[10]][test.results[[10]]$stat=="mean",]
a <- ggplot(foo, aes(x =value, fill=pop))
a + geom_density()













# ####SCRATCH
# #Below is just scratch notes
#
#
# bw <- .1
# plot(density(vals[,1],bw=bw),xlim=c(0,3),ylim=c(0,2))
# lines(density(vals[,2],bw=bw))
# lines(density(vals[,3],bw=bw))
#
# polygon(density(vals[,1],bw=bw), col=rgb(1,0,0,.2))
# polygon(density(vals[,2],bw=bw), col=rgb(0,1,0,.2))
# polygon(density(vals[,3],bw=bw), col=rgb(0,0,1,.2))
#
# #Calculate variance of phenotype values
# phen.var <- matrix(,iter,3)
# colnames(phen.var) <- c("SpeciesA","SpeciesB", "SpeciesHyb")
#
# for(i in 1:iter){
#   phen.var[i,1] <- var(results[[i]]$SpeciesA)
#   phen.var[i,2] <- var(results[[i]]$SpeciesB)
#   phen.var[i,3] <- var(results[[i]]$SpeciesHyb)
# }
# boxplot(phen.var)
# #Need a stat to see if two variances are equal - Levene's test? yes
#
#
# #Make a function for epistasis (interaction among alleles between loci)
# #Which type of epistasis? Maybe randomly choose a group of loci, and a type of
# #epistasis?
#
#

# #This is the recombination code of 5% recombination rate.
# #sample(c(0,1), size=20, replace=T, prob=c(.95,.05))
#
#
