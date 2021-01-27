# ---------------------------------------------------------------------------- #
# Create analysis dataset
# 
# ---------------------------------------------------------------------------- #

library(dplyr)
library(sf)
config <-
  list(
    include_any_fee = TRUE,
    prot_protected_min = 80,
    prot_protected_max = 20,
    tr_year_min = 1985,
    tr_year_max = 2006
  )

spdt <- readRDS(file = "data/preprocessed_spatial.rds")
dt   <- readRDS(file = "data/preprocessed_data.rds")

adt <- 
  dt %>% 
  left_join(spdt, by = "pid") %>%
  mutate(
    
    protmech = case_when(
      (e_p >= f_p) & (e_p > 0) ~ "e",
      f_p > 0 & (f_fee_prot==1) * config$include_any_fee ~ "f",
      TRUE ~ NA_character_
    ),
    
    tr_type = case_when(
      protmech == "e" & e_orgtype != '' ~ paste0('e_', e_orgtype),
      protmech == "f" & f_orgtype != '' ~ paste0('f_', f_orgtype),
      TRUE ~ NA_character_
    ),
    
    tr_year = case_when(
      protmech == "e" ~ e_year, 
      protmech == "f" ~ f_year,
      TRUE ~ NA_real_), 
    
    p_prot_protmech = case_when(
      protmech == "e" ~ e_p, 
      protmech == 'f' ~ f_p, 
      TRUE ~ NA_real_),
    
    tr = (p_prot_protmech > config$prot_protected_min & 
            !is.na(tr_year) & 
            tr_year >= config$tr_year_min &
            tr_year <= config$tr_year_max) * 1L
  ) %>%
  {
    hold  <- .
    
    mutate(
      hold,
      # Identify treatment of neighbors
      tr_neighbors = purrr::map(
        .x = neighbors,
        .f = ~ hold[["tr"]][.x]
      )
    )
  } %>%
  mutate(
    A = tr,
    
    # Proportion of boundary shared with treated units
    A_tilde = purrr::pmap_dbl(
      .l = list(
        x = boundary_length,
        a = tr_neighbors,
        b = neighbor_boundary_lengths
      ),
      .f = function(x, a, b){
        `if`(length(b) == 0, 0, sum(a * b)/x)
      })
  ) 

saveRDS(file = "data/analysis_data.rds")
rm(spdt, dt)
  