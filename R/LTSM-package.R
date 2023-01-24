#' Estimation of latent time shift models
#' 
#' Functions for the estimation of latent time shift models (LTSM).
#' Continuous outcomes are handled for the longitudinal part,
#' whereas the survival part considers multiple competing events.
#' The likelihood is computed using Monte-Carlo integration. Estimation is achieved
#' by maximizing the log-likelihood using a robust iterative algorithm.
#'
#' Please report to the LTSM-team any question or suggestion regarding the package
#' via github only (https://github.com/TiphaineSaulnier/LTSM/issues).
#' 
#' @name LTSM-package
#' @docType package
#' @author Tiphaine Saulnier, Cecile Proust-Lima, Viviane Philipps
#' 
#' @references
#' Kühnel, Raket, Åström, Berger, Hansen, Krismer, ... & EMSA‐SG Natural History Study Investigators. (2022).
#' Disease Progression in Multiple System Atrophy—Novel Modeling Framework and Predictive Factors. Movement Disorders.
#' 
#' Philipps, Hejblum, Prague, Commenges, Proust-Lima (2021).
#' Robust and efficient optimization using a Marquardt-Levenberg algorithm with 
#' R package marqLevAlg, The R Journal 13:2.
#'
#' @keywords package
#' @importFrom graphics axis hist lines matlines matplot mtext par plot points segments polygon
#' @importFrom grDevices rainbow rgb col2rgb n2mfrow
#' @importFrom stats as.formula formula get_all_vars integrate median model.frame model.matrix na.fail na.omit na.pass pchisq pnorm qnorm quantile rnorm sd terms residuals vcov fitted coef update
#' @importFrom survival Surv untangle.specials
#' @importFrom randtoolbox sobol
#' @importFrom stringr str_detect
#' @importFrom parallel clusterEvalQ clusterExport clusterSetRNGStream makeCluster parApply stopCluster
#' @importFrom marqLevAlg mla
#' @useDynLib LTSM, .registration=TRUE, .fixes="C_"
NULL









