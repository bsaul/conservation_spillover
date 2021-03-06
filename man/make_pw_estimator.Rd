% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimator_permuation.R
\name{make_pw_estimator}
\alias{make_pw_estimator}
\title{Create a permutation weighted estimator for the marginal structural model}
\usage{
make_pw_estimator(
  fitter,
  rhs_formula,
  B,
  modify_C = identity,
  predictor = function(object, newdata) {     predict(object = object, newdata =
    newdata, type = "response") },
  fitter_args = list()
)
}
\description{
Create a permutation weighted estimator for the marginal structural model
}
