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
permute_and_stack <- function(dt, permuter){
  dplyr::bind_rows(
    dplyr::mutate(dt,  C = 0),
    dplyr::mutate(permuter(dt), C = 1)
  ) 
}

#' Estimate the density ratio by modeling whether an observation is from 
#' the permuted dataset or original dataset
.pw <- function(dt, rhs_formula, fitter, modify_C, predictor, fitter_args, permuter){
  pdt <- permute_and_stack(dt, permuter)
  pdt$C <- modify_C(pdt$C)
  # browser()
  m  <- do.call(
    fitter,
    args = c(list(formula = update(C ~ 1, rhs_formula), data = pdt), fitter_args))
  
  w  <- predictor(object = m, newdata = dt) 
  w/(1 - w)
}

#' Estimate permutation weights B times and average the results
#' @export
get_permutation_weights <- function(dt, B, rhs_formula, fitter, modify_C = identity, 
                                    predictor, fitter_args = list(),
                                    permuter){
  replicate(
    n    = B,
    expr = .pw(dt, rhs_formula = rhs_formula, fitter = fitter, modify_C = modify_C,
               predictor = predictor, fitter_args = fitter_args,
               permuter = permuter)
  ) %>% 
    apply(1, mean)
}

#' Estimate weights by GLM
#' @export
glm_fitter <- function(formula, data, ...){
  glm(formula = formula, data = data, family = binomial)
}

#' Estimate weights by GAM
#' @export
gam_fitter <- function(formula, data, ...){
  mgcv::gam(formula = formula, data, family = binomial, ...)
}

#' Estimate weights by XGBoost
#' @export
xg_fitter <- function(formula, data, ...){
  vars <- attr(terms(formula), "term.labels")
  X <- data.matrix(data[ , vars])
  y <- data[[rlang::expr_text(rlang::f_lhs(formula))]]
  
  xgboost::xgboost(
    data = X,
    label = y,
    ...
  )
}

#' Estimate weights by randomForest
#' @export
rf_fitter <- function(formula, data, ...){
  randomForest::randomForest(formula = formula, data = data)
}

#' Create a permutation weighted estimator for the marginal structural model
#' @export
make_pw_estimator <- function(fitter, rhs_formula, B, 
                              modify_C = identity,
                              predictor = function(object, newdata) { predict(object = object, newdata = newdata, type = "response") },
                              fitter_args = list(),
                              permuter){
  function(data){
    
    w <- get_permutation_weights(
      dt = data, B = B, rhs_formula = rhs_formula, fitter  = fitter,
      modify_C = modify_C, predictor = predictor, fitter_args = fitter_args,
      permuter = permuter
    )
    
    w
  }
}

#' Make an ipw estimator
#' @export
make_ipw_estimator <- function(fitter, rhs_formula,
                              predictor = function(object, newdata) { predict(object = object, newdata = newdata, type = "response") },
                              fitter_args = list()){
  function(data){
    m  <- do.call(
      fitter,
      args = c(list(formula = update(A ~ 1, rhs_formula), data = data), fitter_args))
    
    p <- predictor(object = m, newdata = data) 
    w <- dbinom(data$A, size = 1, prob = p)
    w <- mean(data$A)/w
    w 
  }
}
