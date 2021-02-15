# ---------------------------------------------------------------------------- #
# Create analysis dataset
# 
# ---------------------------------------------------------------------------- #

library(dplyr)
library(sf)
config <-
  list(
    include_any_fee = TRUE,
    prot_protected_min = 0.8,
    prot_protected_max = 0.2,
    tr_year_min = 1985,
    tr_year_max = 2006
  )

spdt <- readRDS(file = "data/preprocessed_spatial.rds") %>%
  select(-fips, -gisjoin, -geometry, -centroid)
dt   <- readRDS(file = "data/preprocessed_data.rds")

adt0 <- 
  dt %>% 
  mutate(
    
    protmech = case_when(
      (e_p >= f_p) & (e_p > 0) ~ "e",
      (e_p > 0) & is.na(f_p)   ~ "e",
      f_p > 0   & (f_fee_prot==1) * config$include_any_fee ~ "f",
      TRUE ~ NA_character_
    ),
    
    e_year = case_when(
      e_year == -99  ~ 0,
      !is.na(e_year) ~ e_year,
      !is.na(e_p)    ~ 0,
      TRUE ~ NA_real_
    ),
    
    f_year = case_when(
      f_year == -99  ~ 0,
      !is.na(f_year) ~ f_year,
      !is.na(f_p)    ~ 0,
      TRUE ~ NA_real_
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
      protmech == "e" ~ e_p/100, 
      protmech == 'f' ~ f_p/100, 
      TRUE ~ NA_real_),
    
  ) %>%
  mutate(
    data = purrr::pmap(
      .l = list(d = data, ty = tr_year, trp = p_prot_protmech),
      .f = function(d, ty, trp) {
        d %>%
          mutate(
            has_tr = if (is.na(ty)) 0L else (as.integer(ty) <= as.integer(year)) * 1L,
            Ac     = if_else(has_tr == 1L, trp, 0),
            A = (Ac > config$prot_protected_min) * 1L
          )
      }
    )
  ) %>%
  left_join(spdt, by = "pid")

neighbor_info <- 
  adt0 %>%
  select(pid, neighbors, neighbor_boundary_lengths, boundary_length) 

adj <- Matrix::sparseMatrix(
  i = unlist(purrr::imap(neighbor_info$neighbors, ~ rep.int(.y, length(.x)))),
  j = unlist(neighbor_info$neighbors)
)

# Add neighbor average covariates ####
# * set to 0 if no neighbors
adt0 <- 
  adt0 %>%
  mutate(
    log_ha              = log(ha),
    neighbors_avg_ha   = as.numeric((adj %*% (ha)) * (adj %*% rep.int(1, nrow(adt0)))^(-1)),
    neighbors_avg_ha   = if_else(is.nan(neighbors_avg_ha), 0, neighbors_avg_ha),
    
    neighbors_avg_p_wet = as.numeric((adj %*% (p_wet * ha)) * (adj %*% ha)^(-1)),
    neighbors_avg_p_wet = if_else(is.nan(neighbors_avg_p_wet), 0, neighbors_avg_p_wet),
    
    neighbors_avg_slope = as.numeric((adj %*% (slope)) * (adj %*% rep.int(1, nrow(adt0)))^(-1)),
    neighbors_avg_slope = if_else(is.nan(neighbors_avg_slope), 0, neighbors_avg_slope),
  )

# Impute missing demographics ####
# - if parcel has neighbors use average of neighbors
# - else average of parcel's township
impute_missing_demo <- function(x){
  hold <- as.numeric((adj %*% x) * (adj %*% rep.int(1, length(x)))^(-1))
  case_when(
    !is.na(x) ~ x,
    is.na(x) & !is.nan(hold) ~ hold,
    is.na(x) & is.nan(hold)  ~ NA_real_
  )
}

demo_vars <- 
  adt0 %>%
  select(pid, township, starts_with("pop_dens"), starts_with("hh_inc")) %>%
  mutate(
    across(
      .cols = starts_with("pop_dens") | starts_with("hh_inc"),
      .fns =  ~ impute_missing_demo(.x)
    )
  ) %>%
  group_by(township) %>%
  mutate(
    across(
      .cols = starts_with("pop_dens") | starts_with("hh_inc"),
      .fns =  ~ if_else(is.na(.x), mean(.x, na.rm = TRUE), .x)
    )
  ) %>%
  ungroup() %>%
  select(-township) %>%
  tidyr::pivot_longer(
    cols = starts_with("pop_dens") | starts_with("hh_inc")
  ) %>%
  mutate(
    variable = stringr::str_extract(name, "(pop_dens|hh_inc_med)"),
    scale    = stringr::str_extract(name, "(bg|tract)"),
    year     = stringr::str_extract(name, "(?<=_)[0-9]{4}(?=($|\\.))")
  ) %>%
  select(-name) %>%
  group_nest(pid, .key = "demo")

neighbor_bndrys <- Matrix::sparseMatrix(
  i = unlist(purrr::imap(neighbor_info$neighbors, ~ rep.int(.y, length(.x)))),
  j = unlist(neighbor_info$neighbors),
  x = unlist(neighbor_info$neighbor_boundary_lengths)
)

adt1 <- 
  adt0 %>%
  select(pid, data) %>%
  tidyr::unnest(cols = "data") %>%
  group_nest(year) %>%
  mutate(
    data = purrr::map(
      .x = data,
      .f = ~ {
        .x %>%
          mutate(
            # Proportion of boundary shared with treated units
            A_tilde = as.numeric(neighbor_bndrys %*% A) /
                       as.numeric(neighbor_info$boundary_length)
          )
      }
    )
  ) %>%
  tidyr::unnest(cols = data) %>%
  group_nest(pid)

adt <- 
  list(
    select(adt0, -data, -starts_with("pop_dens"), -starts_with("hh_inc")),
    adt1,
    demo_vars
  ) %>%
  purrr::reduce(.f = ~ left_join(.x, .y, by = "pid"))

  
saveRDS(adt, file = "data/analysis_data.rds")
saveRDS(neighbor_bndrys, file = "data/neighbor_boundary_matrix.rds")

rm(spdt, dt, adt0, adt1)
  