# ---------------------------------------------------------------------------- #
# Toying around with analysis
# 
# ---------------------------------------------------------------------------- #

library(smd)
library(mgcv)
library(dplyr)
library(CBPS)
# temp

covars <- rlang::syms(list("ha", "slope", "p_wet", "river_frontage",
                           "lake_frontage", "travel", "township",
                           "coast_2500", "pop_dens_tract_1990"))
temp <-
  adt %>%
  mutate(
    pop_dens_tract_1990 = if_else(is.na(pop_dens_tract_1990), 0, pop_dens_tract_1990),
    township = case_when(
      is.na(township) ~ sample(unique(township), 1),
      TRUE ~ township
    )
  ) %>%
  mutate(
    forest_1990 = purrr::map_dbl(data, ~ .x$forested[.x$year == "1990"]),
    undev_1990  = purrr::map_dbl(data, ~ .x$undeveloped[.x$year == "1990"]),
    forest_2012 = purrr::map_dbl(data, ~ .x$forested[.x$year == "2012"]),
    undev_2012  = purrr::map_dbl(data, ~ .x$undeveloped[.x$year == "2012"]),
    avg_chg_forest = purrr::map_dbl(data, ~ mean(diff(.x$forested[as.integer(.x$year) > 2006]))),
    avg_chg_undev = purrr::map_dbl(data, ~ mean(diff(.x$undeveloped[as.integer(.x$year) > 2006])))
  ) %>%
  dplyr::select(A, !!! covars, forest_1990, undev_1990, forest_2012, undev_2012,
                starts_with("avg_chg"))



fitters <- list(
  A = function(formula, data){
    glm(formula, data = data, family = binomial(link = "logit"))
  },
  Atilde0 = function(formula, data){
    glm(formula, data = data, family = binomial(link = "logit"))
  },
  Atilde1 = function(formula, data){
    glm(formula, data = data, family = Gamma, subset = A_tilde > 0)
  }
)

formulas <- list(
  A       = rlang::expr(A ~ !! purrr::reduce(covars, ~ rlang::expr(!! .x + !! .y))),
  Atilde0 = rlang::expr(I(A_tilde > 0) ~ !! purrr::reduce(covars, ~ rlang::expr(!! .x + !! .y))),
  Atilde1 = rlang::expr(A_tilde ~ !! purrr::reduce(covars, ~ rlang::expr(!! .x + !! .y)))
)


fit_models <- function(fitters, formulas, data){
  purrr::map2(
    .x = fitters,
    .y = formulas,
    .f = ~ { .x(formula = .y, data = data)  }
  )
  
}

models <- fit_models(
  fitters, formulas, adt
)


predictions <- 
  purrr::map(
    .x = models,
    .f = ~ predict(.x, type = "response")
  )

library(smd)

m1 <- glm(
  A ~
    log(ha)*slope*travel +
    I(forest_1990 == 0) +
    forest_1990 +
    I(forest_1990 == 1) +
    undev_1990 +
    I(undev_1990 == 1) +
    I(travel < 20) +
    I(p_wet == 0) + 
    p_wet +
    I(log1p(river_frontage) == 0) + 
    log1p(river_frontage) +
    I(log1p(lake_frontage) == 0) + 
    log1p(lake_frontage) +
    I(coast_2500 == 0) +
    travel*log1p(coast_2500) +
    forest_1990*p_wet + 
    I(undev_1990 == 1)*p_wet + 
    I(undev_1990 == 1)*log1p(pop_dens_tract_1990) + 
    undev_1990*p_wet + 
    undev_1990*log1p(pop_dens_tract_1990) +
    township,
  data = temp,
  family = binomial(link = "logit"),
)

summary(m1)

ps <- predict(m1, type = "response")
ip <- 1/dbinom(temp$A, 1, prob = as.numeric(ps))
hjk <- tapply(ip, temp$A, function(x) sum(1/x))/tapply(ip, temp$A, length)

ip2 <- ifelse(
  temp$A == 1,
  ip * hjk[2],
  ip * hjk[1]
)

hist(ip)
hist(ip2)
ip2[ip2 > 2000]

smd(dplyr::select(temp, !!! covars, forest_1990, undev_1990), pull(temp, A))
smd(dplyr::select(temp, !!! covars, -township, forest_1990, undev_1990), 
    pull(temp, A),
    w = ip2)

temp %>%
  mutate(
    yhat = forest_2012 * ip2
  ) %>%
  group_by(A) %>%
  summarise(
    y_naive = mean(forest_2012),
    y = mean(yhat)
  )


summary(lm(avg_chg_forest*100  ~ A, data = temp))
summary(lm(avg_chg_forest*100 ~ A + (1/ip2) + township, weights = temp$ha, data = temp))
summary(gam(avg_chg_forest*100~ A + te(I(1/ip2)) + township,  data = temp))

summary(lm(avg_chg_undev*100 ~ A + (1/ip2) + township, weights = temp$ha, data = temp))
