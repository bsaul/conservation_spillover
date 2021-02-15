
#' Create a simulator function of a GAM
#' 
#' @param model 
#' @param post a `function(x) {...}` applied to the output of the resulting 
#'    function
#' @return a closure of `newdata` and `newparms` that returns a vector of simulated
#'    data from the fitted GAM model
#' @importFrom mgcv gam
#' @export
create_model_simulator <- function(model, post = identity){
  fm   <- rlang::expr_text(formula(model))
  beta <- coef(model)
  fam  <- model[["family"]]
  inv  <- fam[["linkinv"]]

  on.exit({ rm(model) })
  
  function(newdata, parameter_update = function(b, x) { b }){
    X <- mgcv::gam(as.formula(fm), data = newdata, family = fam, fit = FALSE)[["X"]]
    beta <- parameter_update(b = beta, x = X)
    post(inv(X %*% beta))
  }
  

}

#' Create a simulator
#' 
#' @param exposure_sim the result of `create_model_simulator` for the exposure 
#'     of interest
#' @param outcome_sim the result of `create_model_simulator` for the outcome 
#'     of interest
#' @export
create_simulator <- function(exposure_sim, outcome_sim){
  force(exposure_sim)
  force(outcome_sim)
  function(newdata, exposure_newparms = NULL, outcome_newparms = NULL){
    df <- newdata
    df[["A"]] <- exposure_sim(df, exposure_newparms)
    df %>%
      {
        hold  <- .
        
        dplyr::mutate(
          hold,
          # Update treatment of neighbors
          tr_neighbors = purrr::map(
            .x = neighbor_pids,
            .f = ~ hold[["A"]][hold[["pid"]] %in% .x]
          )
        )
      } %>%
      dplyr::mutate(
        # Proportion of boundary shared with treated units
        A_tilde = purrr::map2_dbl(
          .x = tr_neighbors,
          .y = neighbor_boundary_lengths,
          .f = function(x, y){
            `if`(length(x) == 0 || sum(as.numeric(y)) == 0, 
                 0, sum(x * y)/sum(y))
          })
      ) ->
      df
    df[["Y"]] <- outcome_sim(df, outcome_newparms)
    
    df
  }
}

#' Creates a function which can produce the basis for simulation datasets
#' 
#' @param basisdt a dataset containing `pid`, `geometry`, `neighbor_pids`,
#'  `neighbor_boundary_lengths`
#' @export
create_sim_basis_maker <- function(basisdt){
  function(seed_id, km_buffer){
    
    loc <- sf::st_buffer(basisdt$geometry[basisdt$pid == seed_id, "geometry"], 
                         units::as_units(km_buffer, "km"))
    prd <- 
      drop(sf::st_intersects(basisdt$geometry, loc, sparse = FALSE)) |
      drop(sf::st_within(basisdt$geometry, loc, sparse = FALSE)) 
    
    basisdt[prd, ] %>% 
      {
        hold <- .
        hold %>%
          dplyr::mutate(
            neighbors_pids = purrr::map(
              .x = neighbor_pids,
              .f = ~ which(hold$pid %in% .x)
            ),
            neighbors = purrr::map(
              .x = neighbor_pids,
              .f = ~ which(hold$pid %in% .x)
            ),
            neighbor_boundary_lengths = purrr::map2(
              .x = neighbor_boundary_lengths,
              .y = neighbor_pids,
              .f = ~ .x[.y %in% hold[["pid"]]]
            )
          )
      } 
  }
}

#' Create a simulation basis dataset
#' @export
make_sim_basis_data <- function(basis_maker, id, buffer, 
                                checker    = function(hash) { FALSE },
                                redirector = function(out, hash) { out }){
  sha <- digest::sha1(list(basis_maker, id, buffer))
  
  if (checker(sha)) { return(invisible(NULL)) }
  
  dt <- 
    basis_maker(seed_id = id, km_buffer = buffer) %>%
    dplyr::select(-geometry)
  
  out <-
    list(
      data = dt,
      sha = sha,
      id = id,
      buffer = buffer,
      n  = nrow(dt)
    )
  
  redirector(out, sha)
}


#' Create a simulation dataset
#' @export
make_sim <- function(basedt, simulator, 
                     parms = list(exposure_newparms = identity, outcome_newparms = identity)){
  simdt <- do.call(simulator, args = c(list(newdata = basedt), parms))
  
  list(
    data   = simdt,
    mean_A = mean(simdt$A),
    mean_A_tilde = mean(simdt$A_tilde),
    mean_Y = mean(simdt$Y),
    parms  = parms
  )
}