# ---------------------------------------------------------------------------- #
# Sample to create spatial dataset for simulations
# 
# ---------------------------------------------------------------------------- #

frame <- sf::as_Spatial(sf::st_cast(basisdt$geometry, "POLYGON"))
smpl  <- SDraw::sdraw(frame, n = 3000, type = "BAS")
ids   <- basisdt$pid[as.integer(gsub("ID", "", smpl$geometryID))]
saveRDS(ids, file = "data/bas_3000_sample.rds")
