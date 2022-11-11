# Buonanduci et al. 2022. Fine-scale spatial heterogeneity shapes 
#   compensatory responses of a subalpine forest to severe bark 
#   beetle outbreak. Landscape Ecology.

# Individual diameter growth INLA analysis

require(INLA)
require(tidyverse)

# Import and standardize data ---------------------------------------------

# Import data
dat <- read_csv("Data/input_individual_trees.csv") %>%
  mutate(Species = factor(Species, levels=c("ABLA", "PICO", "PIEN"))) %>%
  mutate(DBH_89_unstd = DBH_89) %>%
  mutate(DBH_04_unstd = DBH_04) %>%
  mutate(DBH_18_unstd = DBH_18)

dat_pre <- dat %>% drop_na(Dia_89_04)
dat_post <- dat %>% drop_na(Dia_04_18)
dat_dif <- dat %>% drop_na(Dia_dif)

# Standardizing the predictors by subtracting their means and 
# dividing them by two times their standard deviations
stand <- function(x){
  (x-mean(x,na.rm=T))/(2*sd(x,na.rm=T))
}

nostand <- c("Plot", "TreeID", "Tag", "Species", "X", "Y", 
             "Dia_89_04", "Dia_04_18", "Dia_dif",
             "DBH_89_unstd", "DBH_04_unstd", "DBH_18_unstd") # variables NOT to standardize
nostand <- which(names(dat) %in% nostand) # convert to column numbers

stand_df <- function(df){
  cbind(df[,c(nostand)], apply(df[,-c(nostand)],2,stand))
}

# Separate lodgepole, spruce, fir
dat_L_pre <- filter(dat_pre, Species=="PICO")
dat_L_post <- filter(dat_post, Species=="PICO")
dat_L_dif <- filter(dat_dif, Species=="PICO")

dat_S_pre <- filter(dat_pre, Species=="PIEN")
dat_S_post <- filter(dat_post, Species=="PIEN")
dat_S_dif <- filter(dat_dif, Species=="PIEN")

dat_F_pre <- filter(dat_pre, Species=="ABLA")
dat_F_post <- filter(dat_post, Species=="ABLA")
dat_F_dif <- filter(dat_dif, Species=="ABLA")

# Standardize all data
dat_pre_std <- stand_df(dat_pre)
dat_post_std <- stand_df(dat_post)
dat_dif_std <- stand_df(dat_dif)

dat_L_pre_std <- stand_df(dat_L_pre)
dat_L_post_std <- stand_df(dat_L_post)
dat_L_dif_std <- stand_df(dat_L_dif)

dat_S_pre_std <- stand_df(dat_S_pre)
dat_S_post_std <- stand_df(dat_S_post)
dat_S_dif_std <- stand_df(dat_S_dif)

dat_F_pre_std <- stand_df(dat_F_pre)
dat_F_post_std <- stand_df(dat_F_post)
dat_F_dif_std <- stand_df(dat_F_dif)




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
A_pre <- inla.spde.make.A(mesh, loc=as.matrix(dplyr::select(dat_pre, X, Y)))
A_post <- inla.spde.make.A(mesh, loc=as.matrix(dplyr::select(dat_post, X, Y)))
A_dif <- inla.spde.make.A(mesh, loc=as.matrix(dplyr::select(dat_dif, X, Y)))

A_L_pre <- inla.spde.make.A(mesh, loc=as.matrix(dplyr::select(dat_L_pre, X, Y)))
A_L_post <- inla.spde.make.A(mesh, loc=as.matrix(dplyr::select(dat_L_post, X, Y)))
A_L_dif <- inla.spde.make.A(mesh, loc=as.matrix(dplyr::select(dat_L_dif, X, Y)))

A_S_pre <- inla.spde.make.A(mesh, loc=as.matrix(dplyr::select(dat_S_pre, X, Y)))
A_S_post <- inla.spde.make.A(mesh, loc=as.matrix(dplyr::select(dat_S_post, X, Y)))
A_S_dif <- inla.spde.make.A(mesh, loc=as.matrix(dplyr::select(dat_S_dif, X, Y)))

A_F_pre <- inla.spde.make.A(mesh, loc=as.matrix(dplyr::select(dat_F_pre, X, Y)))
A_F_post <- inla.spde.make.A(mesh, loc=as.matrix(dplyr::select(dat_F_post, X, Y)))
A_F_dif <- inla.spde.make.A(mesh, loc=as.matrix(dplyr::select(dat_F_dif, X, Y)))

# Model formulation -------------------------------------------------------

# The SPDE model 
spde <- inla.spde2.pcmatern(
  # Mesh and smoothness parameter
  mesh = mesh, alpha = 2,
  # P(practic.range < 30m) = 0.5
  prior.range = c(30, 0.5),
  # P(sigma > 0.02) = 0.01
  prior.sigma = c(0.02, 0.01)) ##penalized complexity -- prior shifts sigma towards zero


# Create data stacks

# Pre-outbreak growth data stacks
inla.stack_pre <- function(dataframe, response, Amatrix, spdeobject) {
  
  inla.stack(
    data = list(y = eval(substitute(response), dataframe)),
    A = list(Amatrix, 1),
    effects = list(list(s = 1:spdeobject$n.spde), #spatial random effects
                   data.frame(Intercept = 1,
                              Species = dataframe$Species,
                              DBH = dataframe$DBH_89,
                              DBHrw = dataframe$DBH_89_unstd, # for fitting random walk to DBH
                              Dens_live = dataframe$NDI_89_live,
                              Dens_mort_prop = dataframe$NDI_89_propmort,
                              Dens_consp_prop = dataframe$NDI_89_propconsp,
                              Plot = dataframe$Plot,
                              TPI = dataframe$TP,
                              Slope = dataframe$Slope)), #fixed effects
    tag = 'est')
  
}

stk_pre <-   inla.stack_pre(dat_pre_std,   Dia_89_04, A_pre,   spde)
stk_L_pre <- inla.stack_pre(dat_L_pre_std, Dia_89_04, A_L_pre, spde)
stk_S_pre <- inla.stack_pre(dat_S_pre_std, Dia_89_04, A_S_pre, spde)
stk_F_pre <- inla.stack_pre(dat_F_pre_std, Dia_89_04, A_F_pre, spde)

# Post-outbreak growth data stacks
inla.stack_post <- function(dataframe, response, Amatrix, spdeobject) {
  
  inla.stack(
    data = list(y = eval(substitute(response), dataframe)),
    A = list(Amatrix, 1),
    effects = list(list(s = 1:spdeobject$n.spde), #spatial random effects
                   data.frame(Intercept = 1,
                              Species = dataframe$Species,
                              DBH = dataframe$DBH_04,
                              DBHrw = dataframe$DBH_04_unstd, # for fitting random walk to DBH
                              Dens_live = dataframe$NDI_04_live,
                              Dens_mort_prop = dataframe$NDI_04_propmort,
                              Dens_consp_prop = dataframe$NDI_04_propconsp,
                              Plot = dataframe$Plot,
                              TPI = dataframe$TPI,
                              Slope = dataframe$Slope)), #fixed effects
    tag = 'est')
  
}

stk_post <- inla.stack_post(dat_post_std, Dia_04_18, A_post, spde)
stk_dif <- inla.stack_post(dat_dif_std, Dia_dif, A_dif, spde)

stk_L_post <- inla.stack_post(dat_L_post_std, Dia_04_18, A_L_post, spde)
stk_L_dif <- inla.stack_post(dat_L_dif_std, Dia_dif, A_L_dif, spde)

stk_S_post <- inla.stack_post(dat_S_post_std, Dia_04_18, A_S_post, spde)
stk_S_dif <- inla.stack_post(dat_S_dif_std, Dia_dif, A_S_dif, spde)

stk_F_post <- inla.stack_post(dat_F_post_std, Dia_04_18, A_F_post, spde)
stk_F_dif <- inla.stack_post(dat_F_dif_std, Dia_dif, A_F_dif, spde)

# Model fitting ------------------------------------------------------------

# The sum of the log CPO (Pettit 1990) can be computed for each model 
# using the following function:
slcpo <- function(m, na.rm = TRUE) {
  sum(log(m$cpo$cpo), na.rm = na.rm)
}


# Model fitting function ---------
inla_fit <- function(form, data_stack){
  inla(form,
       family = "gaussian",
       data = inla.stack.data(data_stack),
       quantiles=c(0.025, 0.5, 0.975),
       control.compute = list(cpo = TRUE, waic=TRUE),
       control.predictor = list(A = inla.stack.A(data_stack),
                                compute = TRUE))
}


# Null models ----------

# Intercept only, all species combined
form_null <- y ~ Species + f(Plot, model="iid") + f(s, model=spde)

res_pre_null <- inla_fit(form_null, stk_pre)
res_post_null <- inla_fit(form_null, stk_post)
res_dif_null <- inla_fit(form_null, stk_dif)

# Intercept only, species-specific models for growth release
form_null_species <- y ~ f(Plot, model="iid") + f(s, model=spde)

res_L_dif_null <- inla_fit(form_null_species, stk_L_dif)
res_S_dif_null <- inla_fit(form_null_species, stk_S_dif)
res_F_dif_null <- inla_fit(form_null_species, stk_F_dif)


# Pre-outbreak lodgepole growth model -------------
form_L_pre <- y ~ DBH*Dens_live + DBH*Dens_mort_prop + Dens_consp_prop + Slope + TPI + 
  f(Plot, model="iid") + f(s, model=spde)
res_L_pre <- inla_fit(form_L_pre, stk_L_pre)

# Pre-outbreak spruce growth model -------------
form_S_pre <- y ~ Dens_live + DBH*Dens_mort_prop + Dens_consp_prop + Slope + TPI + 
  f(Plot, model="iid") + f(s, model=spde)
res_S_pre <- inla_fit(form_S_pre, stk_S_pre)

# Pre-outbreak fir growth model -------------
form_F_pre <- y ~ Dens_live + DBH*Dens_mort_prop + Dens_consp_prop + Slope + TPI + 
  f(Plot, model="iid") + f(s, model=spde)
res_F_pre <- inla_fit(form_F_pre, stk_F_pre)

# Post-outbreak lodgepole growth model -------------
form_L_post <- y ~ 0 + Intercept + f(DBHrw, model="rw1", scale.model=TRUE) + 
  Dens_live + Dens_mort_prop + Dens_consp_prop + Slope + TPI + 
  DBH:Dens_live + DBH:Dens_mort_prop +
  f(Plot, model="iid") + f(s, model=spde)
res_L_post <- inla_fit(form_L_post, stk_L_post)

# Post-outbreak spruce growth model -------------
form_S_post <- y ~ 0 + Intercept + f(DBHrw, model="rw1", scale.model=TRUE) + 
  Dens_live + Dens_mort_prop + Dens_consp_prop + Slope + TPI + 
  DBH:Dens_mort_prop + 
  f(Plot, model="iid") + f(s, model=spde)
res_S_post <- inla_fit(form_S_post, stk_S_post)

# Post-outbreak fir growth model -------------
form_F_post <- y ~ 0 + Intercept + f(DBHrw, model="rw1", scale.model=TRUE) + 
  Dens_live + Dens_mort_prop + Dens_consp_prop + Slope + TPI + 
  DBH:Dens_live + DBH:Dens_consp_prop + DBH:Slope +
  f(Plot, model="iid") + f(s, model=spde)
res_F_post <- inla_fit(form_F_post, stk_F_post)

# Post-outbreak lodgepole growth release model -------------
form_L_dif <- y ~ DBH*Dens_live + DBH*Dens_mort_prop + Dens_consp_prop + Slope + TPI + 
  f(Plot, model="iid") + f(s, model=spde)
res_L_dif <- inla_fit(form_L_dif, stk_L_dif)

# Post-outbreak spruce growth release model -------------
form_S_dif <- y ~ DBH + Dens_live + DBH*Dens_mort_prop + Dens_consp_prop + Slope + TPI + 
  f(Plot, model="iid") + f(s, model=spde)
res_S_dif <- inla_fit(form_S_dif, stk_S_dif)

# Post-outbreak fir growth release model -------------
# Remove Dens_consp_prop due to collinearity with Dens_mort_prop
form_F_dif <- y ~ DBH + Dens_live + Dens_mort_prop + Slope + TPI + 
  f(Plot, model="iid") + f(s, model=spde)
res_F_dif <- inla_fit(form_F_dif, stk_F_dif)


