cbps_ee <- function( data, object, ...){

  X  <- stats::model.matrix(object$formula, data = data, ...)
  Y  <- as.numeric(stats::model.frame(grab_response_formula(object), data = data)[[1]])
  n  <- length(Y)
  p  <- length(stats::coef(object))
  w  <- object$weights
  family <- object$family$family
  link   <- object$family$link
  invlnk <- object$family$linkinv
  varfun <- object$family$variance
  
  function(theta){
    # browser()
    lp <- drop(X %*% theta[1:p]) # linear predictor
    f  <- invlnk(lp)  # fitted values
    r  <- (Y - f)     # residuals
    V  <- diag(varfun(f), nrow = n, ncol = n)
    xlp <- stats::dlogis(lp)
    w <- (Y - f)/(f*(1 - f))
    D <- apply(X, 2, function(x) x * xlp)
    
    ipw <- 1/dbinom(Y, 1, prob = f)
    
    D2 <- apply(X[ , -1, drop = FALSE], 2, function(x) x * ipw)
    D3 <- apply(X[ , -1, drop = FALSE], 2, function(x) x * w)
    if (!is.matrix(D))  { D <- matrix(D, nrow=1) } ## apply will undesireably coerce to vector

    c(
      t(D) %*% solve(V) %*% (r),
      t(D3)
      # , # score eqns,
      # ((cbw*(Y == 1)) - (cbw*(Y == 0)))
    )
  }
}

dt <- tibble(
  a = rbinom(100, 1, prob = 0.5),
  x = rnorm(100)
)

m <- glm(a ~ x, data = dt, family = binomial)

library(geex)

summary(m)

m_estimate(
  estFUN = cbps_ee,
  data   = dt,
  outer_args = list(object = m),
  root_control = setup_root_control(start = c(beta0 = 0, beta1 = 0, mean = .4))
)

fit <- CBPS(treat ~ age + educ, data = LaLonde, ATT = 0)
summary(fit)
hist(fit$fitted.values)

smd(LaLonde[c("age", "educ")], g = LaLonde$treat, w = fit$weights)
  
m <- glm(treat ~ age + educ, data = LaLonde, family = binomial)
summary(m)


mm <- m_estimate(
  estFUN = cbps_ee,
  data   = LaLonde,
  outer_args = list(object = m),
  root_control = setup_root_control(start = c(beta0 = 3.61, beta1 = -0.11,
                                              beta2 = -0.23, beta1 = -0.11,
                                              beta2 = -0.23), 
                                    # atol = c(0.1, 0.1, 0.1, 0.1, 0.1),
                                    # ctol = 0.1,
                                    verbose = TRUE
  )
)
summary(mm)
