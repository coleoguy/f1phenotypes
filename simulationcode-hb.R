getwd()
setwd("~/Documents/GitHub/f1phenotypes")
source("functions.R")
# set model parameters
#set population size
N <- 50
#set number of impacted loci (QTL)
loci <- 10
# set effect size this is the difference in phenotype between the two
# alternative homozygous genotypes
# set allele frequency of major allele in 2 parent species. First number is
# SpeciesA, second number is SpeciesB
afreq <- c(0.5,0.5)
#set genome size
gsize <- 30
#set the number of simulations
iter <- 100
# sample size equal for parents and hybrid pops
s.size <- 50
# number of epistatic gene pairs
epipair <- 0
epi.type <- ""
#hset <- "all_dom"
hset <- "all_add"
#hset <- "runif"
#hset <- "halfdom"
#effect.size <- "neg_binom"
#effect.size <- "runif"
effect.size <- "allequal_stephen"
#effect.size <- "onelarge_stephen"
verbose = F
#verbose = T
mating = "random"
#mating = "assortative"
#simulate <- function(N, loci, effect.size, afreq, gsize,
#                    iter, s.size, epipair, hset, verbose)



output <- list()
for(i in 1:10){
  # set effect size this is the difference in phenotype between the two
  # alternative homozygous genotypes
  output[[i]] <- simulate(N = N, loci = loci, effect.size = effect.size,
                          afreq = afreq, gsize = gsize,
                          iter = iter, s.size = s.size, epipair = epipair,
                          epi.type = epi.type, hset = hset, mating = mating,
                          verbose = verbose)
}


install.packages("ggraptR")
library(ggraptR)
ggraptR(x[[1]])

#plotting
library(ggplot2)
output <- do.call(rbind, output)
png("10loci_differentallelefreqs_domrunif_epi5_randommating_mean.png", width = 700, height = 700)
meandistribution <- ggplot(subset(output, stat %in% c("mean")), aes(x = value, fill=pop))
meandistribution <- meandistribution + geom_density(alpha = 0.2)
meandistribution <- meandistribution + xlab("Mean Phenotype") + ylab("Density")
meandistribution <- meandistribution + ggtitle("10loci_differentallelefreqs_domrunif_epi5_randommating")
meandistribution <- meandistribution + theme(
  axis.title.x=element_text(size=20, colour = "black"),
  axis.title.y=element_text(size=20, colour = "black"),
  axis.text.x=element_text(size=18, colour = "black"),
  axis.text.y=element_text(size=18, colour = "black"),
  legend.title = element_text(size = 18),
  legend.text = element_text(size = 18)
)
meandistribution
dev.off()

png("10loci_differentallelefreqs_domrunif_epi5_randommating_variance.png", width = 700, height = 700)
varianceboxplot <- ggplot(subset(output, stat %in% c("var")), aes(x = pop, y = value))
varianceboxplot <- varianceboxplot + geom_boxplot(aes(fill = pop), alpha = 0.6) +
  geom_jitter(color = "black", alpha = 0.2) +
  xlab("Population") +
  ylab("Variance of Phenotype") +
  ggtitle("10loci_differentallelefreqs_domrunif_epi5_randommating") +
  theme(
    axis.title.x=element_text(size=20, colour = "black"),
    axis.title.y=element_text(size=20, colour = "black"),
    axis.text.x=element_text(size=18, colour = "black"),
    axis.text.y=element_text(size=18, colour = "black"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    panel.background = element_rect(fill = "transparent", color = "black") 
  )
varianceboxplot
dev.off()






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
