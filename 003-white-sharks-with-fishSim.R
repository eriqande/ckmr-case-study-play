

#devtools::install_github(repo = "SMBaylis/fishSim")

library(fishSim)

repro_age <-  14
max_age <- 55
batchsize <- 9   # mean number of offspring per mature female per breeding attempt, poisson distributed (set below)
max_litter <- 12
age_mort <- c(rnorm(repro_age, (1-0.7), 0.001),   # mortality = 1-survival
              rep(rnorm(max_age-repro_age, (1-0.95), 0.001 )))
sim_years <- 20  # if this runs too long as-is, error comes up

indiv <- makeFounders(pop = 2000,
                      osr = c(0.5, 0.5),
                      stocks = c(1), # only one stock here
                      maxAge = max_age,
)


for(y in 1:sim_years){
  indiv <- altMate(indiv,
                   firstBreed = repro_age,
                   batchSize = batchsize,
                   fecundityDist = "poisson",
                   maxClutch = max_litter,
                   type = "flat",
                   year = y
  )
  # doesn't remove animals from indiv, just gives them a death year
  indiv <- mort(indiv,
                maxAge = max_age,
                type = "age",
                ageMort = age_mort
  )
  indiv <- birthdays(indiv)

  ## lethal sampling of 5 animals, either sex, not sure if this is set up correctly

}

indiv <- capture(indiv,
                 n=6,
                 year = y,
                 fatal = TRUE,
                 age = 5:12)  # age of sampled/killed animals
tail(indiv)

check_growthrate(mateType = "flat",
                 mortType = "age",
                 batchSize = batchsize,
                 firstBreed = repro_age,
                 ageMort = age_mort # I think this input is correct
)

# 0.38 indicates annual growth rate is -0.62. Yikes


PoNG(mateType = "flat",
     mortType = "age",
     ageMort = age_mort,
     batchSize = batchsize,
     maxClutch = max_litter,
     firstBreed = repro_age
)

# $root indicates rate at which we would have null growth with these parameters




## mark n individuals alive at the end of the simulation as captured (and killed)

for(year in 1:sim_years){  # I don't know if this is set up correctly to sample in each of the years we have
  indiv <- capture(indiv,
                   n = 5,
                   year = year,
                   fatal = TRUE)
}
## look up each animal's ancestors and look for shared ancestors between each pair of sampled animals:

pairs <- findRelativesPar(indiv = indiv,
                          sampled = TRUE,
                          nCores = 2)

POPs <- pairs[pairs$OneTwo == 1,]   ## Parent-Offspring pairs
GGPs <- pairs[pairs$OneThree == 1,] ## Grandparent-Grandoffspring pairs
HSPs <- pairs[pairs$TwoTwo == 1,]   ## Half-sibling pairs
FSPs <- pairs[pairs$TwoTwo == 2,]   ## Full Sibling pairs (self-comparisons are automatically excluded)
FCPs <- pairs[pairs$ThreeThree == 2 & pairs$TwoTwo != 1,] ## Full Cousin pairs

## look at the number of shared ancestors at each ancestral
## generation, for one of the half-sibling pairs.
lookAtPair(HSPs[1,])
