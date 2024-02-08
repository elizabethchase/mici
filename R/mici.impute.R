#' Function to generate multiple imputations
#' 
#' This function generates multiple imputations according to the nonparametric
#' multiple imputation scheme described in Chase et al. (2024+). In particular,
#' for censored individuals it generates an imputed event time and type by sampling
#' from the Kaplan-Meier distribution of overall event times and then randomly
#' sampling an event type among those observed at the imputed event time. The
#' function outputs a length M list of imputed datasets free of censoring. 
#' 
#' @param ftime Either the name of the variable in the data that contains the event
#' times, or a length n vector of event times.
#' @param ftype Either the name of the variable in the data that contains the event
#' types, or a length n vector of event types. Note that event = 0 corresponds to
#' censoring, event = 1 corresponds to the event of interest, and any event indicator
#' greater than 1 corresponds to competing event(s). 
#' @param data Optional: if the ftime and ftype inputs are variable names (rather than 
#' numeric vectors), this is the dataframe in which those variables are stored. 
#' @param M The number of imputations desired; the default is 200. 
#' @param scheme The imputation approach to be used. Options are "RSI" for risk 
#' set imputation and "KMI" for Kaplan-Meier imputation. KMI is the default.
#' 
#' @return The function returns a length M list, where each object of the list 
#' is a dataframe with columns ftime, ftype containing the imputed event times and
#' event types. If data was inputted to the original function, all original columns
#' of the dataset are kept in each of the M dataframes of the list. If only vectors
#' were inputted to the original function, the columns time and type in each of the 
#' M dataframes store the original inputs (without imputation). 
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
#' @import dplyr
#' @import survival
#' @export
mici.impute <- function(ftime = NULL,
                        ftype = NULL,
                        data = NULL,
                        M = 200,
                        scheme = "KMI"){
  if (!is.null(data)){
    ftime_name <- substitute(ftime)
    if(is.symbol(ftime_name)) ftime_name <- deparse(ftime_name)
    ftime <- data[[ftime_name]]
    ftype_name <- substitute(ftype)
    if(is.symbol(ftype_name)) ftype_name <- deparse(ftype_name)
    ftype <- data[[ftype_name]]
    
    if (!is.vector(ftime)){stop("ftime must contain numeric event times.")}
    if (!is.vector(ftype)){stop("ftype must contain numeric event types.")}
    
    u <- data
    u$ftime <- ftime
    u$ftype <- ftype
    u$id <- c(1:length(ftime))
    
  } else{
    if (!is.vector(ftime)){stop("ftime must contain numeric event times.")}
    if (!is.vector(ftype)){stop("ftype must contain numeric event types.")}
    
    u <- data.frame("time" = ftime,
                    "type" = ftype,
                    "ftime" = ftime,
                    "ftype" = ftype,
                    "id" = c(1:length(ftime)))
  }
  
  dc <- u[u$ftype != 0, ] ## complete cases
  dd <- u[u$ftype == 0, ] ##cases that need imputation
  
  if (nrow(dd)==nrow(u)){
    stop("All observations in this dataset are censored; imputation
         cannot be performed.")
  }
  if (nrow(dc)==nrow(u)){
    warning("There was no censoring in this dataset. 
            Imputation was not needed.")
    ipd <- select(u, -id)
    myimps <- list(ipd)
    
    return(myimps)
  }
  
  xt <- dd$ftime
  n <- nrow(u)
  myimps <- list()
  
  if (scheme=="KMI"){
    g <- summary(survfit(Surv(ftime, ftype != 0) ~ 1, data = u))
    gm <- g$surv[length(g$surv)]
    w <- g$time
    wp <- -diff(c(1, g$surv))
    
    if (gm > 0) {
      wp <- c(wp, gm)
      w <- c(w, max(u$ftime) + 1)
    }
    
    for (j in 1:M) {
      cts <- NULL
      cevent <- NULL
      for (jj in 1:length(xt)) {
        sub = w > xt[jj]
        if (gm > 0) {
          if (length(w[sub]) == 1) {
            cts[jj] <- w[sub]
            cevent[jj] <- resample(c(1:2), size = 1)
          } else {
            cts[jj] <- resample(w[sub], 1, replace = TRUE, prob = wp[sub])
            if (cts[jj] <= max(dc$ftime)) {
              cevent[jj] <-
                resample(dc$ftype[near(dc$ftime, cts[jj])], size = 1)
            } else{
              cevent[jj] <- resample(c(1:2), size = 1)
            }
          }
        } else{
          if (length(w[sub]) == 0) {
            cts[jj] <- max(u$ftime) + 1
            cevent[jj] <- resample(c(1:2), size = 1)
          } else if (length(w[sub]) == 1) {
            cts[jj] <- w[sub]
            cevent[jj] <- resample(dc$ftype[near(dc$ftime, cts[jj])], size = 1)
          } else {
            cts[jj] <- resample(w[sub], 1, replace = TRUE, prob = wp[sub])
            cevent[jj] <- resample(dc$ftype[near(dc$ftime, cts[jj])], size = 1)
          }
        }
      }
      
      dd$ftime = cts
      dd$ftype = cevent
      ipd = bind_rows(dd, dc) %>% arrange(id) %>% select(-id)
      
      myimps <- append(myimps, list(ipd))
    }
  } else if (scheme=="RSI"){
    xt <- sort(unique(dd$ftime))
    for (j in 1:M) {
      ipd <- u
      for (jj in 1:length(xt)) {
        sub = ipd$ftime > xt[jj]
        if (length(which(sub)) > 1){
          inds <- resample(which(sub), size = length(which(near(ipd$ftime, xt[jj]) & 
                                                             ipd$ftype==0)), replace = TRUE)
          ipd[near(ipd$ftime, xt[jj]) & ipd$ftype==0,c("ftime", "ftype")] <- ipd[inds,c("ftime", "ftype")]
        } 
      }
      last_event <- max(u$ftime[u$ftype==1])
      last_cens <- max(ipd$ftime[ipd$ftype==0])
      if (length(last_cens)>0){
        if (last_cens <= last_event){
          ipd$ftype[ipd$ftype==0] <- resample(u$ftype[u$ftime==last_event & u$ftype!=0], 
                                              size = length(which(ipd$ftype==0)), replace = TRUE)
        } 
      }
    
      ipd <- ipd %>% arrange(id) %>% select(-id)
      myimps <- append(myimps, list(ipd))
    }
  }
  
  return(myimps)
}
