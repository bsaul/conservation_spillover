# ---------------------------------------------------------------------------- #
# Create plasmode simulators
# 
# ---------------------------------------------------------------------------- #

library(dplyr)
library(mgcv)
library(conservationSpillover)

# Basis dataset for plasmode simulations for effect on 2012 outcomes ####
covars <- rlang::syms(list("ha", "slope", "p_wet", "river_frontage",
                           "lake_frontage", "travel", "township",
                           "coast_2500", "pop_dens_tract_1990"))

basisdt <-
  adt %>%
  tidyr::unnest(
    cols = data
  ) %>%
  select(
    pid, A, A_tilde,
    !!! covars,
    geometry, neighbors, boundary_length, neighbor_boundary_lengths,
    year, forested, undeveloped, 
    starts_with("undeveloped"),
    starts_with("forest")
  ) %>%
  filter(year %in% c(2006, 2012)) %>%
  tidyr::pivot_wider(
    names_from = c("year"),
    values_from = c("forested", "undeveloped")
  ) %>%
  {
    hold  <- .
    
    dplyr::mutate(
      hold,
      # Identify neighbor ids
      neighbor_pids = purrr::map(
        .x = neighbors,
        .f = ~ hold[["pid"]][.x]
      ),
      neighbor_p_wet = purrr::map_dbl(
        .x = neighbors,
        .f = ~ {
          ngbha <- hold[["ha"]][.x]
          ngbpw <- hold[["p_wet"]][.x]
          `if`(length(.x) == 0,
               0,
               sum(ngbha * ngbpw)/sum(ngbha))
          
        }
      )
    )
  } %>%
  dplyr::mutate(
    log_ha = log(ha),
    p_wet_0 = (p_wet == 0),
    neighbor_p_wet_0 = (neighbor_p_wet == 0),
    log1p_coast_2500 = log1p(coast_2500),
    d_forested    = forested_2012 - forested_2006,
    d_undeveloped = undeveloped_2012 - undeveloped_2006,
    forest_not_decrease  = d_forested >= 0,
    any_forest_change    = d_forested != 0,
    develop_not_increase = d_undeveloped <= 0,
    any_develop_change   = d_undeveloped != 0
  )

# Prediction models for simulation ####
expm <- gam(
  formula = A ~ log(ha)*slope + 
    I(p_wet == 0) + p_wet +
    I(neighbor_p_wet == 0) + neighbor_p_wet +
    log1p(coast_2500) + travel +
    s(log(ha)) +
    s(slope) +
    s(p_wet) +
    s(travel),
  data   = basisdt,
  family = binomial
)

expm <- glm(
  formula = A ~ log(ha)*slope + 
    I(slope == 0) + 
    I(p_wet == 0) + p_wet + I(p_wet^2) + I(p_wet^3)  +
    I(neighbor_p_wet == 0) + neighbor_p_wet + I(neighbor_p_wet^2) + I(neighbor_p_wet^3)  +
    log1p(coast_2500)*travel +
    I(log(ha)^2) + I(log(ha)^3) +
    I(slope^2) + I(slope^3),
  data   = basisdt,
  family = binomial
)

summary(expm)
summary(expm2)

exposure_simulator <- 
  create_model_simulator(
    model = expm,
    post   = function(x) { rbinom(length(x), size = 1, prob = x) }
  )

# p <- predict(expm, data = basisdt, type = "response")
# w <- mean(basisdt$A)/dbinom(basisdt$A, size = 1, prob = p)

cmodel2 <- glm(
  formula = I(d_forested * 100) ~ A*A_tilde + 
    log(ha)*slope + 
    I(slope == 0) + 
    I(p_wet == 0) + p_wet + I(p_wet^2) + I(p_wet^3)  +
    I(neighbor_p_wet == 0) + neighbor_p_wet + I(neighbor_p_wet^2) + I(neighbor_p_wet^3)  +
    log1p(coast_2500)*travel +
    I(log(ha)^2) + I(log(ha)^3) +
    I(slope^2) + I(slope^3),
  data   = basisdt %>% filter(undeveloped_2006 > 0.5),
  family = gaussian
)

coutcome_simulator <- 
  create_model_simulator(
    model = cmodel,
    post   = function(x) { drop(x) + rnorm(length(x), sd = 3.5) } 
  )

bmodel <- 
  glm(
    formula = forest_not_decrease ~ A*A_tilde +
      log(ha)*slope + 
      I(slope == 0) + 
      I(p_wet == 0) + p_wet   +
      I(neighbor_p_wet == 0) + neighbor_p_wet  +
      log1p(coast_2500)*travel +
      I(log(ha)^2) + I(log(ha)^3) +
      I(slope^2) + I(slope^3),
    data    = basisdt,
    family  = binomial,
  )

outcome_simulator <- 
  create_model_simulator(
    model = bmodel,
    post   = function(x) { rbinom(length(x), size = 1, prob = x) }
  )

simulator0 <- create_simulator(exposure_simulator, coutcome_simulator)
simulator1 <- create_simulator(exposure_simulator, outcome_simulator)

saveRDS(basisdt, file = "data/simulation_basisdt.rds")
saveRDS(simulator0, file = "data/plasmode_simulator_0.rds")
saveRDS(simulator1, file = "data/plasmode_simulator_1.rds")

rm(p, w)
