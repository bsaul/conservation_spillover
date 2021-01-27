# ---------------------------------------------------------------------------- #
# Create plasmode simulators
# 
# ---------------------------------------------------------------------------- #

library(dplyr)
library(mgcv)

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
    d_forested    = forested_2012 - forested_2006,
    d_undeveloped = undeveloped_2012 - undeveloped_2006
  )

# Functions ####

# Prediction models for simulation
exposure_simulator <- 
  create_model_simulator(
    formula = A ~ log(ha)*slope + 
      I(p_wet == 0) + p_wet + 
      I(neighbor_p_wet == 0) + neighbor_p_wet + 
      log1p(coast_2500) + travel +
      s(log(ha)) +
      s(slope) + 
      s(p_wet) + 
      s(travel),
    data   = basisdt,
    family = binomial,
    post   = function(x) { rbinom(length(x), size = 1, prob = x) }
  )

outcome_simulator <- create_model_simulator(
  formula = d_forested ~
    A*A_tilde + 
    log(ha)*slope +
    I(p_wet == 0) + p_wet +
    I(neighbor_p_wet == 0) + neighbor_p_wet + 
    log1p(coast_2500) + travel +
    s(log(ha)) +
    s(slope) +
    s(p_wet) +
    s(travel),
  data   = basisdt,
  family = gaussian,
  post   = drop
)

simulator1 <- create_simulator(exposure_simulator, outcome_simulator)

saveRDS(basisdt, file = "data/simulation_basisdt.rds")
saveRDS(simulator1, file = "data/plasmode_simulator_1.rds")
