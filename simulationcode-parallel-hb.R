setwd("~/Documents/GitHub/f1phenotypes")
source("functions.R")
library(ggplot2)
library(doMC)
registerDoMC(cores = 6)


N <- 50
loci <- 10
gsize <- 20
s.size <- 50
iter <- 1
verbose <- F
mating <- "random"
baselevelpheno <- 10


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
                                        baselevelpheno = baselevelpheno, verbose = verbose)
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
                                          baselevelpheno = baselevelpheno, verbose = verbose)
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

#Convert list to dataframe
newoutput <- do.call(rbind, x)
#Make new variable called allelefreqdif. This will be the absolute value of the
#difference of allele frequencies between Pop A and Pop B

library(stringr)
foo <- data.frame(do.call('rbind', strsplit(as.character(newoutput$afreq), '_', fixed = TRUE)))
foo$X1 <- as.numeric(foo$X1)
foo$X2 <- as.numeric(foo$X2)
newoutput <- cbind(newoutput, foo)
newoutput$allelefreqdif <- abs(newoutput$X1 - newoutput$X2)

#Make new variable called propepi. This will be the proportion of loci that are
#involved in epistasis
newoutput$propepi <- newoutput$epipair/newoutput$loci




#Classification tree instructions - from https://www.guru99.com/r-decision-trees.html
#From newoutput, need to extract only the rows where stat=="cvcomparison"
onlypopcomparisons <- newoutput[newoutput[, "stat"] == "cvcomparison", ]
#I should shuffle the data around, because it is sorted
shuffle_index <- sample(1:nrow(onlypopcomparisons))
onlypopcomparisons <- onlypopcomparisons[shuffle_index, ]
#Install and attach packages
install.packages("rpart")
library(rpart)
install.packages("rpart.plot")
library(rpart.plot)
#Create training/testing datasets
create_train_test <- function(data, size = 0.8, train = TRUE) {
  n_row = nrow(data)
  total_row = size * n_row
  train_sample <- 1:total_row
  if (train == TRUE) {
    return (data[train_sample, ])
  } else {
    return (data[-train_sample, ])
  }
}
data_train <- create_train_test(data = onlypopcomparisons, size = 0.8, train = TRUE)
data_test <- create_train_test(data = onlypopcomparisons, size = 0.8, train = FALSE)
#Build the model
fit <- rpart(value ~ esize + afreq + h + epipair + epi.type, data = data_train, method="class")
summary(fit)
#Plot the model
rpart.plot(fit, extra = 106, fallen.leaves = T)
#Make Predictions
predict_unseen <- predict(fit, data_test, type = "class")
#Test hybrids that had smaller CV and larger CV than both parents
table_mat <- table(data_test$value, predict_unseen)
table_mat
#Measure performance with a confusion matrix
accuracy_Test <- sum(diag(table_mat))/sum(table_mat)
print(paste('Accuracy for test', accuracy_Test))

#The below code draws esizes from a normal distribution of their CVs. User needs
#to input the number of loci and the number of effect sizes they would like
loci <- 10
numberofeffectsize <- 2
effectsizelist <- list()
for(i in 1:1000){
  vec <- sample(x = 10000, size = loci, replace = TRUE)
  newvec <- vec/sum(vec)
  effectsizelist[[i]] <- newvec
}
effectsizelist <- as.data.frame(matrix(unlist(effectsizelist), ncol = 10, byrow = TRUE))
effectsizelist <- transform(effectsizelist, SD=apply(effectsizelist, 1, sd, na.rm = TRUE))
effectsizelist <- transform(effectsizelist, MEAN=apply(effectsizelist, 1, mean, na.rm = TRUE))
effectsizelist$CV <- effectsizelist$SD/effectsizelist$MEAN
hist(effectsizelist$CV)
effectsizelist_stepsize <- 1000/numberofeffectsize
effectsizelist <- effectsizelist[order(effectsizelist$CV),]
effectsizelist_new <- round(effectsizelist[seq(1, nrow(effectsizelist), effectsizelist_stepsize), ], digits = 4)
newesize <- subset(effectsizelist_new, select = -c(SD, MEAN, CV))
esize <- as.list(as.data.frame(t(newesize)))

###The below code makes a list of epipairs, depending on the number of loci and the amount of epistases
#you want to analyze. I still need to edit the function code so I can input epipair as a list.
loci <- matrix(c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
loci
total_epi <- 5

epipairlist <- list()
nrow(loci)
for(i in 1:length(loci)){
  epipairlist[[i]] <- round(seq.int(from = 0, to = loci[i]/2, length.out = total_epi), digits = 0)
}


##Trying below to make new esize code compatible with loci matrix/vector. I'm 
#not going to sort or anything by CV, and just take a random sample.
loci <- matrix(c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
esize_test <- list()
esize <- list()
epipair <- list()
total_epi <- 5

for(i in 1:length(loci)){
  epipair[[i]] <- round(seq.int(from = 0, to = loci[i]/2, length.out = total_epi), digits = 0)
  for(k in 1:10){
    vec <- sample(x = 10000, size = loci[i], replace = TRUE)
    newvec <- vec/sum(vec)
    esize_test[[k]] <- newvec
  }
  esize[[i]] <- esize_test
}



