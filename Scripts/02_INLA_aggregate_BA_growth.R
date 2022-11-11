# Buonanduci et al. 2022. Fine-scale spatial heterogeneity shapes 
#   compensatory responses of a subalpine forest to severe bark 
#   beetle outbreak. Landscape Ecology.

# Aggregate BA growth INLA analysis

require(INLA)
require(tidyverse)

# Import and standardize data ---------------------------------------------

# Import data
dat <- read_csv("Data/input_aggregate_raster.csv")

# Standardizing the predictors by subtracting their means and 
# dividing them by two times their standard deviations
stand <- function(x){
  (x-mean(x,na.rm=T))/(2*sd(x,na.rm=T))
}

nostand <- c("X", "Y", "Plot", 
             "grow_BA_89_04", "grow_BA_04_18", "grow_BA_dif",
             "ingrow_stems_89_04", "ingrow_stems_04_18", "ingrow_stems_dif") # variables NOT to standardize
nostand <- which(names(dat) %in% nostand) # convert to column numbers

# Drop non-host/late-seral proportion NAs
dat_pre <- drop_na(dat, surv_BA_LS_89_prop) 
dat_post <- drop_na(dat, surv_BA_LS_04_prop)

# Standardize all data
dat_std <- cbind(dat[,c(nostand)], apply(dat[,-c(nostand)],2,stand))
dat_pre_std <- cbind(dat_pre[,c(nostand)], apply(dat_pre[,-c(nostand)],2,stand))
dat_post_std <- cbind(dat_post[,c(nostand)], apply(dat_post[,-c(nostand)],2,stand))

# Build mesh --------------------------------------------------------------

# Extract coordinates
coords <- dplyr::select(dat, X, Y) %>% as.matrix()

# Boundary
bound <- inla.nonconvex.hull(coords, convex = 40)

# Build the mesh
mesh <- inla.mesh.2d(boundary = bound, 
                     cutoff = 1, 
                     max.edge = c(20, 150),
                     min.angle=c(30, 21),
                     max.n=c(48000, 16000), ## Safeguard against large meshes.
                     max.n.strict=c(128000, 128000)) ## Don't build a huge mesh!

# Once the mesh has been defined, the projector matrix is computed
A <- inla.spde.make.A(mesh, loc=as.matrix(dplyr::select(dat_std, X, Y)))
A_pre <- inla.spde.make.A(mesh, loc=as.matrix(dplyr::select(dat_pre_std, X, Y)))
A_post <- inla.spde.make.A(mesh, loc=as.matrix(dplyr::select(dat_post_std, X, Y)))

# Model formulation -------------------------------------------------------

# The SPDE model 
spde <- inla.spde2.pcmatern(
  # Mesh and smoothness parameter
  mesh = mesh, alpha = 2,
  # P(practic.range < 30m) = 0.5
  prior.range = c(30, 0.5),
  # P(sigma > 1) = 0.01
  prior.sigma = c(0.1, 0.01)) ##penalized complexity -- prior shifts sigma towards zero


# Create data stacks

# Pre-outbreak growth
stk_pre_G <- inla.stack(
  data = list(y = dat_pre_std$grow_BA_89_04),
  A = list(A_pre, 1), 
  effects = list(list(s = 1:spde$n.spde), #spatial random effects
                 data.frame(Intercept = 1, 
                            BA_surv = dat_pre_std$surv_BA_89,
                            BA_mort_prop = dat_pre_std$mort_BA_89_04_prop,
                            BA_LS_prop_surv = dat_pre_std$surv_BA_LS_89_prop,
                            Stems_surv = dat_pre_std$surv_stems_89,
                            Plot = dat_pre_std$Plot,
                            TPI = dat_pre_std$TPI,
                            Slope = dat_pre_std$Slope)), #fixed effects
  tag = 'est')


# Post-outbreak growth
stk_post_G <- inla.stack(
  data = list(y = dat_post_std$grow_BA_04_18),
  A = list(A_post, 1), 
  effects = list(list(s = 1:spde$n.spde), #spatial random effects
                 data.frame(Intercept = 1, 
                            BA_surv = dat_post_std$surv_BA_04,
                            BA_mort_prop = dat_post_std$mort_BA_04_18_prop,
                            BA_LS_prop_surv = dat_post_std$surv_BA_LS_04_prop,
                            Stems_surv = dat_post_std$surv_stems_04,
                            Plot = dat_post_std$Plot,
                            TPI = dat_post_std$TPI,
                            Slope = dat_post_std$Slope)), #fixed effects
  tag = 'est')

# Growth release
stk_dif_G <- inla.stack(
  data = list(y = dat_post_std$grow_BA_dif),
  A = list(A_post, 1), 
  effects = list(list(s = 1:spde$n.spde), #spatial random effects
                 data.frame(Intercept = 1, 
                            BA_surv = dat_post_std$surv_BA_04,
                            BA_mort_prop = dat_post_std$mort_BA_04_18_prop,
                            BA_LS_prop_surv = dat_post_std$surv_BA_LS_04_prop,
                            Stems_surv = dat_post_std$surv_stems_04,
                            Plot = dat_post_std$Plot,
                            TPI = dat_post_std$TPI,
                            Slope = dat_post_std$Slope)), #fixed effects
  tag = 'est')


# Model fitting ------------------------------------------------------------

form_null <- y ~ f(Plot, model="iid") + f(s, model=spde)

form <- y ~ BA_mort_prop + BA_LS_prop_surv + BA_surv + 
  Stems_surv + I(Stems_surv^2) + Slope + TPI + f(Plot, model="iid") + f(s, model=spde)


# Evaluate model fit -------

# The sum of the log CPO (Pettit 1990) can be computed for each model using 
# the following function:
slcpo <- function(m, na.rm = TRUE) {
  sum(log(m$cpo$cpo), na.rm = na.rm)
}


## Null model for growth release ---------
res_dif_G_null <- inla(form_null, 
                       family = "gaussian",
                       data = inla.stack.data(stk_dif_G), 
                       quantiles=c(0.025, 0.5, 0.975),
                       control.compute = list(cpo = TRUE, waic=TRUE),
                       control.predictor = list(A = inla.stack.A(stk_dif_G), 
                                                compute = TRUE))

## Pre-outbreak model ---------
res_pre_G <- inla(form, 
                  family = "gaussian",
                  data = inla.stack.data(stk_pre_G), 
                  quantiles=c(0.025, 0.5, 0.975),
                  control.compute = list(cpo = TRUE, waic=TRUE),
                  control.predictor = list(A = inla.stack.A(stk_pre_G), 
                                           compute = TRUE))

## Post-outbreak model ---------
res_post_G <- inla(form, 
                   family = "gaussian",
                   data = inla.stack.data(stk_post_G), 
                   quantiles=c(0.025, 0.5, 0.975),
                   control.compute = list(cpo = TRUE, waic=TRUE),
                   control.predictor = list(A = inla.stack.A(stk_post_G), 
                                            compute = TRUE))

## Growth release model ---------
res_dif_G <- inla(form, 
                  family = "gaussian",
                  data = inla.stack.data(stk_dif_G), 
                  quantiles=c(0.025, 0.5, 0.975),
                  control.compute = list(cpo = TRUE, waic=TRUE),
                  control.predictor = list(A = inla.stack.A(stk_dif_G), 
                                           compute = TRUE))


