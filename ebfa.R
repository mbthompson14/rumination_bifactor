library(psych)
library(GPArotation)

load("rumitems.RData")
rumitems.cor <- cor(rumitems)

# ======= Exploratory oblique bi-factor analysis =========

# the rotation function I use below does an orthogonal bifactor rotation
# is there an oblique equivalent? e.g. replace GPForth() with GPFoblq()?

# method described here: https://escholarship.org/uc/item/3w14k56d#article_main
# I'm following this example: https://www.tandfonline.com/doi/full/10.1080/10705511.2019.1622421?casa_token=CB6g7O-U4wkAAAAA%3AqkJZVL7Bkf-rS859WFyzEjPzOvTFHd8vEMyzAKlcz2VS_KHzaicdvdpQmO_YvhGZe49Fb7aaEDSm
# and corresponding syntax: https://osf.io/5tdr2
# FindBifactor.orth() from: http://sites.baylor.edu/lvm5

# parallel analysis--principal axis
rum.pa.mean <- fa.parallel(
  rumitems.cor,
  n.obs = 86,
  fa = "fa",
  fm = "pa"
  )

# extract 5 factors
f5.pa.smc <- fa(
  rumitems.cor,
  nfactors = 5,
  n.obs = 86,
  rotate = "none",
  fm = "pa",
  max.iter = 100,
  SMC = smc(rum_sub2.cor)
  )

# extract 5 factors
f5.pa.smc.biquart <- fa(
  rumitems.cor,
  nfactors = 5,
  n.obs = 86,
  rotate = "biquartimin",
  fm = "pa",
  max.iter = 100,
  SMC = smc(rum_sub2.cor)
  )


# function to complete EBFA using arbitrary number of replications taken from:
# Loehlin, J. C., & Beaujean, A. A. (2016). Syntax companion for Latent Variable Models: An Introduction yo Factor, Path, And Structural Equation Analysis (5th Ed.) (5th ed.). Waco, TX: Baylor Psychometric Laboratory.

# rotation options are bifactor (for biquartimin) and geomin (for bigeomin)
FindBifactor.orth <- function(A,reps=10,rotation="bifactor",solutions=1,round=8,maxit=1000, seed=NA)
{
  # A <- A[[1]] # this is only needed if come from psych package
  require(GPArotation)
  require(parallel)
  m <- dim(A)[2]
  seed <- ifelse(!is.na(seed),round(seed),ceiling(runif(1, 0, 10^9)))
  set.seed(seed)
  ran.mat <- replicate(reps,Random.Start(m),simplify=FALSE)
  y <- mclapply(ran.mat, function(z) GPForth(A,method=rotation,Tmat =z,maxit=maxit))  # GPForth
  y <- y[lapply(y, function(z) z$convergence) ==TRUE]
  criterion <- sapply(y, function(z) min(z$Table[,2]))
  results <- lapply(y, function(z) z$loadings)
  criterion <- round(criterion,round)
  results <- lapply(results, function(x) round(x,2))
  index.val <- 1:length(y)
  criterion.index <- data.frame(criterion,index.val)
  criterion.index <- criterion.index[order(criterion),]
  criterion.index.u <- criterion.index[match(unique(criterion.index$criterion), criterion.index$criterion),]
  index.keep <- criterion.index.u$index.val[1:solutions]
  output <- list(criterion=criterion[index.keep], loadings=results[index.keep])
  return(output)
}

# extract factor loadings
f5.pa.smc.load <- f5.pa.smc$loadings

# bi-quartimin rotation
f5.biquart <- FindBifactor.orth(
  f5.pa.smc.load,
  reps = 1000,
  rotation = "bifactor",
  solutions = 10
  )

# bi-geomin rotation
f5.bigeomin <- FindBifactor.orth(
  f5.pa.smc.load,
  reps = 1000,
  rotation = "geomin",
  solutions = 10
  )

# psych biquartimin: implements the oblique bifactor rotation introduced by Jennrich and Bentler (2011)
# https://personality-project.org/r/psych/help/Promax.html
f5.psych.biquart <- psych::biquartimin(f5.pa.smc.load, eps=1e-5, maxit=1000)

f5.psych.biquart$loadings %>%
  as.data.frame() %>%
  apply(2, function(x) ifelse(abs(x) < 0.3, NA, x)) %>%
  apply(2, function(x) round(x,digits = 2)) %>%
  write.csv(paste("f5.psych.biquart",format(Sys.time(),"%Y%m%d%H%M%S"),"csv",sep="."))

f5.psych.biquart.scores <- factor.scores(rumitems, f5.psych.biquart$loadings, method = "Bartlett")
cor(f5.psych.biquart.scores$scores,use = "pairwise")
