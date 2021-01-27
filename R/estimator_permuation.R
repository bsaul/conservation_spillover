#' Permutes the treatment vector and recomputes the proportion of neighbors 
#' treated.
#' @export
permute_A <- function(dt){
  
  permutation <- sample.int(nrow(dt), replace = FALSE)
  
  dt %>%
    dplyr::mutate(
      A = A[permutation]
    ) %>%
    {
      hold  <- .
      
      dplyr::mutate(
        hold,
        # Identify treatment of neighbors
        tr_neighbors = purrr::map(
          .x = neighbors,
          .f = ~ hold[["A"]][.x]
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
    )
}


#' Stack the observed data with the permuted data
#' @export
permute_and_stack <- function(dt){
  dplyr::bind_rows(
    dplyr::mutate(dt,  C = 0),
    dplyr::mutate(permute_A(dt), C = 1)
  ) 
}

#' Estimate the density ratio by modeling whether an observation is from 
#' the permuted dataset or original dataset
.pw <- function(dt, rhs_formula, fitter){
  m  <- fitter(formula = update(C ~ 1, rhs_formula), 
               data    = permute_and_stack(dt))
  w  <- predict(object = m, newdata = dt, type = "response") 
  w/(1 - w)
}

#' Estimate permutation weights B times and average the results
#' @export
get_permutation_weights <- function(dt, B, rhs_formula, fitter){
  replicate(
    n    = B,
    expr = .pw(dt, rhs_formula = rhs_formula, fitter  = fitter)
  ) %>% 
    apply(1, mean)
}

#' Estimate weights by GLM
#' @export
glm_fitter <- function(formula, data){
  glm(formula = formula, data = data, family = binomial)
}

#' Estimate weights by GAM
#' @export
gam_fitter <- function(formula, data){
  mgcv::gam(formula = formula, data, family = binomial)
}

#' Create a permutation weighted estimator for the marginal structural model
#' @export
make_pw_estimator <- function(fitter, rhs_formula, B, msm_formula = Y ~ A*A_tilde){
  function(data){
    w <- get_permutation_weights(
      dt = data, B = B, rhs_formula = rhs_formula, fitter  = fitter
    )
    
    lm(Y ~ A*A_tilde, data = data, weights = w)
  }
}
