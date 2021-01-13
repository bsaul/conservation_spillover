# ---------------------------------------------------------------------------- #
# Preprocesses data in preparation for analysis
# 
# * data.csv downloaded from https://doi.org/10.7910/DVN/8RCE84
# 
# ---------------------------------------------------------------------------- #

library(dplyr)
library(tidyr)

outFile <- "data/preprocessed_data.rds" # location for saving file
nrows <- NULL # set to 100 (e.g.) for testing pipeline

# Source data ####
dt <- read.csv(
  file = "data/data.csv",
  stringsAsFactors = FALSE
)
 
longerSpec <- dt %>%
  build_longer_spec(
  cols = matches("p_(f|u|prot)_"),
  names_to = c("variable", "year", "buffer"),
  names_pattern = "p_(f|u|prot)_([0-9]{4})_*($|[0-9]{2,5}$)",
  names_ptypes = list(outcome = "character",
                      year    = "integer",
                      buffer    = "integer"),
  values_to = "value") %>% 
  mutate(
    variable = case_when(
      variable == "f" ~ "forested",
      variable == "u" ~ "undeveloped",
      variable == "prot" ~ "neighbor_protected_area_"
    )
  )


dt %>% 
  filter({ if(!is.null(nrows)) (row_number() <= nrows) else TRUE } ) %>% 
  mutate_at(
    .vars = vars(p_wet, slope, river_frontage, lake_frontage, lake_ha, travel,
                 coast_2500    # TODO: is this appropriate for coast_2500?
    ),
    .f    = ~ {ifelse(is.na(.), 0, .)}
  ) %>% 
      {
        df <- .
        ntdat <- select(df, pid, !matches("p_(f|u|prot)_"))

        tdat <- 
          select(df, pid, matches("p_(f|u|prot)_")) %>%
          pivot_longer_spec(longerSpec) %>%
          mutate(value = value/100) %>% # rescale to proportion for my own sanity
          pivot_wider(
            id_cols = c(pid, year), 
            names_from = c(variable, buffer), 
            values_from = value, 
            names_sep = ''
          ) %>%
          group_by(pid) %>%
          nest() 
        
        left_join(ntdat, tdat, by ="pid")

    } -> out

saveRDS(out, file = outFile)
  
