library(tidyverse)
library(plotly)


N_adults <- function(N0, t, lambda) {
  N0 * exp(lambda * t)
}


A_traj <- tibble(
  year = 1:60
) %>%
  mutate(
    num_adults = N_adults(800, year, 0.01)
  )


# simulate sampling 10 sharks a year from the years
# 20 to 40, with ages ranging from 3 to 8, So let's
# store this in a tibble
set.seed(5)
samples <- tibble(
  samp_year = rep(20:40, each=10)
) %>%
  mutate(
    age = sample(3:8, size = n(), replace = TRUE, prob = c(8:3)),
    born_year = samp_year - age
  ) %>%
  left_join(A_traj, by = c("born_year" = "year"))

n <- nrow(samples)
s1 <- samples[rep(1:n, each = n),]
s2 <- samples[rep(1:n, times = n),]
names(s2) <- paste0(names(s2), ".old")

pairs <- bind_cols(s1, s2) %>%
  filter(born_year.old < born_year)

# that tibble is sufficient to calculate our HSP probabilities, HSP(i, j)

# here is the adult survival:
phiA <- 0.94


pairs2 <- pairs %>%
  mutate(
    HSP_prob = (4/num_adults) * (phiA ^ (born_year - born_year.old))
  )

# make a plot of these things!
ggplot(pairs2, aes(x = born_year, y = HSP_prob, fill = factor(born_year - born_year.old))) +
  geom_point(shape = 21)


# now, what about simulating whether each pairs is a HSP or not
pairs3 <- pairs2 %>%
  mutate(
    isHSP = rbernoulli(n(), p = HSP_prob),
    age_diff = born_year - born_year.old
  )

pair_counts <- pairs3 %>%
  count(born_year, age_diff, HSP_prob, isHSP) %>%
  pivot_wider(names_from = isHSP, values_from = n, values_fill = 0L) %>%
  rename(n_UP = `FALSE`, n_HSP = `TRUE`)



#' @param pars a  vector (Ninit, lambda, phiA)
#' @param X a tibble like pair counts
hsp_nll <- function(pars, X) {
  N0 <- pars[1]
  L <- pars[2]
  P <- pars[3]

  NLL <- X %>%
    mutate(
      N = N0 * exp(L * born_year),
      hspp = (4 / N) * P ^ age_diff,
      logl = log(hspp) * n_HSP + log(1 - hspp) * n_UP
    )

  -sum(NLL$logl)
}

# # starting conditions
# pars <- c(Ninit = 4000, lambda = 0.0, phiA = 0.9)
#
# optim(
#   par = pars,
#   fn = hsp_nll,
#   X = X,
#   method = "L-BFGS-B",
#   lower = c(Ninit = 200, lambda = 0.0001, phiA = 0.6),
#   upper = c(Ninit = 5000, lambda = -0.0001, phiA = 1.0)
#   )

# let's make some three-D plots
# the true values were:

# N_init <- 800
# phiA 0.94
# lambda 0.01


# First assume the correct growth rate. lambda = 0.01
# and calculate the likelihood on a grid of possible
# values of N_init and phiA.
LL_tib_lambda_0.01 <- expand_grid(
  N_init = seq(200, 2000, by = 20),
  phiA = seq(0.85, 1, by = 0.01)
) %>%
  mutate(
    nll = map2_dbl(
      .x = N_init,
      .y = phiA,
      .f = function(x, y) {
        pars <- c(x, 0.01, y)
        -hsp_nll(pars, pair_counts)
      }
    )
  )

Nivals <- unique(LL_tib_lambda_0.01$N_init)
phivals <- unique(LL_tib_lambda_0.01$phiA)
nll <- matrix(LL_tib_lambda_0.01$nll, nrow = length(phivals), ncol = length(Nivals))

# here is a three-d plot
plot_ly(x = ~Nivals, y = ~phivals, z = ~nll) %>%
  add_surface()

# it might be easier to see as a contour plot
plot_ly(
  x = Nivals,
  y = phivals,
  z = nll,
  type = "contour")

# but, if we are going to do that, we may as well just use ggplot, right?
ggplot(LL_tib_lambda_0.01, aes(x = N_init, y = phiA, z = nll)) +
  geom_contour_filled(binwidth = 1) +
  geom_contour(binwidth = 1, colour = "black") +
  theme(legend.position = "none")

# OK, so, what I see from this is that it is going to be best to just
# use ggplot for all of this, and we are going to be able to also facet
# it by assumed lambda values.  Cool.

# First assume the correct growth rate. lambda = 0.01
# and calculate the likelihood on a grid of possible
# values of N_init and phiA.
LL_tib2 <- expand_grid(
  N_init = seq(200, 2000, by = 20),
  lambda = seq(-0.04, 0.04, by = 0.01),
  phiA = seq(0.85, 1, by = 0.01)
) %>%
  mutate(
    nll = pmap_dbl(
      .l = list(N0 = N_init, Lam = lambda, phiA = phiA),
      .f = function(N0, Lam, phiA) {
        pars <- c(N0, Lam, phiA)
        -hsp_nll(pars, pair_counts)
      }
    )
  )


ggplot(LL_tib2, aes(x = N_init, y = phiA, z = nll)) +
  #geom_contour_filled(binwidth = 1) +
  geom_contour(binwidth = 1, colour = "gray") +
  theme(legend.position = "none") +
  facet_wrap(~lambda) +
  theme_bw()

