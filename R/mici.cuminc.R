#' Function to estimate the cumulative incidence from multiply imputed data
#' 
#' This function estimates the cumulative incidence from multiply imputed data
#' following the approach described in Chase et al. (2024+). Specifically, the 
#' function takes as input a list of M datasets in which all censoring has been
#' replaced with imputed event times and types. Then estimating the cumulative 
#' incidence at time t is just counting the number of events from the cause-of-interest
#' by time t and dividing by the sample size n.  
#' 
#' @param imp_obj A list containing M imputed datasets free of censoring with columns
#' ftime and ftype containing the imputed event times and types, respectively. This 
#' is outputted from mici.impute(). 
#' @param event_interest The event indicator of the event of interest. The default is
#' 1.
#' @param times A length l vector of times at which cumulative incidence estimates are 
#' requested. If none are provided, the cumulative incidence at every unique event time 
#' in the data will be returned.
#' @param int.type The type of uncertainty interval to use. Options include "Wald",
#' "Wilson", "Bayes", or "Wald_untrans". The default (and our recommendation for most 
#' settings) is "Wilson". Briefly, "Wald" is the log-log transformed Wald interval 
#' around the cumulative incidence. "Wilson" is a Wilson score interval around the 
#' cumulative incidence. "Bayes" is a Bayesian beta-binomial interval around the 
#' cumulative incidence. "Wald_untrans" is the untransformed Wald interval around the 
#' cumulative incidence (note that this is not constrained to fall in [0,1]). 
#' @param conf.width The function will return a conf.width*100% uncertainty interval. 
#' The default is 0.95. 
#' @param bayes_alpha If using the "Bayes" option for int.type, this is the alpha
#' shape parameter for the Beta prior on the cumulative incidence at all timepoints. 
#' The default is 0.8.
#' @param bayes_beta If using the "Bayes" option for int.type, this is the beta
#' shape parameter for the Beta prior on the cumulative incidence at all timepoints. 
#' The default is 1.2.
#' @param bayes_samps If using the "Bayes" option for int.type, this is the number 
#' of psoterior samples to draw for the beta-binomial interval on each imputation.
#' The default is 1000.
#' 
#' @return The function returns a list with elements:
#' \describe{
#'  \item{n}{The sample size.}
#'  \item{time}{The length l vector of times.}
#'  \item{cuminc}{The cumulative incidence estimate at each timepoint.}
#'  \item{int.type}{The type of uncertainty interval being returned.}
#'  \item{conf.width}{The uncertainty interval width being returned.}
#'  \item{lower}{The lower bound of the uncertainty interval estimate at each timepoint.}
#'  \item{upper}{The upper bound of the uncertainty interval estimate at each timepoint.}
#' }
#' 
#' @examples
#' ftime <- sample(c(1:10), size = 100, replace = TRUE)
#' ftype <- sample(c(0:2), size = 100, replace = TRUE)
#' 
#' myimps <- mici.impute(ftime, ftype, M = 10, scheme = "KMI")
#' mycuminc <- mici.cuminc(myimps, int.type = "Wilson")
#' 
#' @examples
#' mydata <- data.frame("time" = sample(c(1:10), size = 100, replace = TRUE),
#' "type" = sample(c(0:2), size = 100, replace = TRUE))
#' 
#' myimps <- mici.impute(time, type, data = mydata, scheme = "RSI")
#' mycuminc <- mici.cuminc(myimps, int.type = "Wald")
#' 
#' @importFrom stats qnorm quantile rbeta var
#' @import dplyr
#' @import binom
#' @export
mici.cuminc <- function(imp_obj = NULL,
                        event_interest = 1,
                        times = NULL,
                        int.type = "Wilson",
                        conf.width = 0.95,
                        bayes_alpha = 0.8,
                        bayes_beta = 1.2,
                        bayes_samps = 1000){
  
  if (is.null(imp_obj)){stop("You must input the output from mici.impute as imp_obj!")}
  
  event_types <- unique(imp_obj[[1]]$ftype)
  if (0 %in% event_types){warning("These data still contain censoring: this is only sensible if the remaining
                                  censoring is after the last observed event time or if the original data only
                                  contained censored observations.")}
  
  if (event_interest %notin% event_types){warning("The event of interest indicator (event_interest) does not
                                               appear in the observed data.")}
  
  if (conf.width >= 1 | conf.width <= 0){stop("conf.width must be a number between 0 and 1!")}
  
  crit_value <- qnorm((1 + conf.width)/2)
  
  if (int.type == "Bayes" & bayes_alpha <= 0){stop("bayes_alpha must be greater than 0.")}
  if (int.type == "Bayes" & bayes_beta <= 0){stop("bayes_beta must be greater than 0.")}
  if (int.type == "Bayes" & bayes_samps <= 0){stop("bayes_samps must be greater than 0.")}
  if (int.type == "Bayes" & !is.integer(bayes_samps)){stop("bayes_samps must be an integer.")}
  
  if (is.null(times)){
    times <- sort(unique(imp_obj[[1]]$time[imp_obj[[1]]$type != 0]))
  }
  
  n <- nrow(imp_obj[[1]])
  M <- length(imp_obj)
  l <- length(times)
  
  cuminc <- matrix(NA, nrow = l, ncol = M)
  if (int.type=="Bayes"){b_beta <- NULL}
  
  for (i in 1:M){
    subinc <- NULL
    disdat <- imp_obj[[i]]
    if (int.type=="Bayes"){bayes_cuminc <- matrix(nrow = bayes_samps, ncol = l)}
    for (j in 1:l) {
      subinc[j] <- length(which(disdat$ftime <= times[j] & disdat$ftype == event_interest)) / n
      if (int.type=="Bayes"){
        num_events <- length(which(disdat$ftime <= times[j] & disdat$ftype == event_interest))
        bayes_cuminc[,j] <- rbeta(bayes_samps, 
                                  bayes_alpha + num_events, 
                                  bayes_beta + n - num_events)
      }
    }
    cuminc[,i] <- subinc
    if (int.type=="Bayes"){b_beta <- rbind(b_beta, bayes_cuminc)}
  }
  
  if (int.type=="Wald_untrans"){
    v <- cuminc * (1 - cuminc) / n
    est.cuminc <- apply(cuminc, 1, mean)
    within_var <- apply(v, 1, mean)
    between_var <- apply(cuminc, 1, var)
    overall_var <- within_var + (1 + 1 / M) * between_var
    
    lower <- est.cuminc - crit_value*sqrt(overall_var)
    upper <- est.cuminc + crit_value*sqrt(overall_var)
    
  } else if (int.type=="Wald"){
    log_v <- (1-cuminc)/(n*cuminc*log(cuminc)^2)
    log_cuminc <- log(-log(cuminc))
    
    est.cuminc <- apply(cuminc, 1, mean)
    
    within_log_var <- apply(log_v, 1, mean)
    between_log_var <- apply(log_cuminc, 1, var)
    overall_log_var <- within_log_var + (1 + 1 / M) * between_log_var
    
    lower <- est.cuminc^(exp(crit_value*sqrt(overall_log_var)))
    upper <- est.cuminc^(exp(-1*crit_value*sqrt(overall_log_var)))
    
  } else if (int.type=="Wilson"){
    if (M > 1){
      est.cuminc <- apply(cuminc, 1, mean)
      v <- cuminc * (1 - cuminc) / n
      
      within_var <- apply(v, 1, mean)
      between_var <- apply(cuminc, 1, var)
      
      r <- (1 + 1/M)*between_var/within_var
      r[is.nan(r)] <- 0
      t_n <- (crit_value^2)/n
      tr_n <- (crit_value^2)*r/n
      q_sum <- 2*est.cuminc + t_n + tr_n
      one_sum <- 2*(1 + t_n + tr_n)
      
      centerpoint <- q_sum/one_sum
      summand <- sqrt((q_sum^2/one_sum^2)-(2*(est.cuminc^2)/one_sum))
    } else{
      est.cuminc <- cuminc
      centerpoint <- (cuminc + ((crit_value^2)/(2*n)))/(1 + ((crit_value^2)/n))
      summand <- (crit_value*sqrt((cuminc*(1-cuminc) + ((crit_value^2)/(4*n)))/n))/(1 + ((crit_value^2)/n))
    }
    
    lower <- centerpoint - summand
    upper <- centerpoint + summand
    
  } else if (int.type=="Bayes"){
    est.cuminc <- apply(cuminc, 1, mean)
    lower <- apply(b_beta, MARGIN = 2, quantile, (1-conf.width)/2)
    upper <- apply(b_beta, MARGIN = 2, quantile, (1-(1-conf.width)/2))
  }
  
  results <- list(
    "n" = n,
    "time" = times,
    "cuminc" = est.cuminc,
    "int.type" = int.type,
    "conf.width" = conf.width,
    "lower" = lower,
    "upper" = upper)
  
  return(results)
}