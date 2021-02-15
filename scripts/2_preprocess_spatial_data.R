# ---------------------------------------------------------------------------- #
# Preprocesses spatial data in preparation for analysis
# 
# * pc_MA spatial data downloaded from https://doi.org/10.7910/DVN/8RCE84
# * township data downloaded from:
#   https://docs.digital.mass.gov/dataset/massgis-data-community-boundaries-towns-survey-points
# ---------------------------------------------------------------------------- #

library(dplyr)
library(sf)

source("R/prep_data.R")

outFile <- "data/preprocessed_spatial.rds" # location for saving file
nrows <- NULL # set to 100 (e.g.) for testing pipeline

# Load parcel data ####
sp_parcels <- 
  st_read(
    dsn = here::here("data", "extdata", "pc_MA"), 
    layer = "pc_MA",
    stringsAsFactors = FALSE
  ) %>%
  st_transform(crs = 5070) 

# Load township data ####
sp_townships <- 
  sf::st_read(
    dsn = here::here("data", "extdata", "townssurvey_shp"),
    layer = "TOWNSSURVEY_POLY",
    stringsAsFactors = FALSE
  ) %>%
  st_transform(crs = 5070) %>%
  group_by(TOWN, POP2000) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  select(town = TOWN, pop2000 = POP2000, geometry)

# Add spatial variables needed for analytic dataset ####
sp_parcels %>% 
  filter({ if(!is.null(nrows)) (row_number() <= nrows) else TRUE } ) %>% 
  mutate(
    centroid  = st_centroid(geometry, of_largest_polygon = TRUE),
    boundary  = st_cast(geometry, "MULTILINESTRING", warn = FALSE),
    
    # Length of boundary in meters
    boundary_length = st_length(boundary),
    
    # Indices of the other parcels touching each parcel 
    neighbors = purrr::map(st_touches(geometry), ~ .x),

    # Name of township in which the centroid (buffered by 2km) of a parcel lies.
    township  = sf::st_within(centroid, sp_townships) %>%
      # Handle case that centroid is not in any township and
      # --just to be sure-- the case it is in multiple townships
      purrr::map_int(.f = ~ if(length(.x) == 0) { NA_integer_ } else { .x[1] } ) %>%
      as.data.frame(sp_townships)[., "town"]
  ) %>%
  # { # TODO: this approach is a bit computationally intensive
  #   df <- .
  #   mutate(
  #     df,
  #     
  #     # Indices of parcels contained by a 10k buffer of each parcel, excluding
  #     #   the target parcel itself. This can be used to simplify later spatial
  #     #   operations by quickly subsetting to parcels within 10k.
  #     indices_parcels_within_10k = purrr::imap(
  #       .x = geometry,
  #       .f = ~ {
  #         st_buffer(.x, dist = 10000) %>%
  #           st_contains(df[-.y, ])
  #       }
  #     )
  #   )
  # } %>%
  # Measure shared boundary lengths
  st_set_geometry("boundary") %>%
  {
    hold <- .
    hold$neighbor_boundary_lengths <-
      purrr::imap(
        .x = hold$neighbors,
        .f = ~ {
          shared_boundaries_length(hold, .y, .x)
        })
    hold
  } %>%
  st_set_geometry("geometry") ->
  out

# Checks ###
# Check that all parcels have a township
stopifnot(!anyNA(out$township))

# Check that the sum of neighbor boundaries is within 1% of boundary length
stopifnot(all(purrr::map2_lgl(
  .x = out$boundary_length,
  .y = out$neighbor_boundary_lengths,
  .f = ~ sum(.y) < (.x * 1.01)
)))


saveRDS(out, file = outFile)

