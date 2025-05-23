#' Stratified nowcasts of incomplete, time-stamped reporting data.
#'
#' Produces nowcasts stratified by a single variable of interest, e.g. by geographic unit (province/state/region) or by age group.
#'
#' @param data A time series of reporting data in line list format (one row per case), with a column \code{onset_date} indicating date of case onset, and a column \code{report_date} indicating date of case report.
#' @param now An object of datatype \code{Date} indicating the date at which to perform the nowcast.
#' @param units Time scale of reporting. Options: "1 day", "1 week".
#' @param onset_date In quotations, the name of the column of datatype \code{Date} designating the date of case onset. e.g. "onset_week"
#' @param report_date In quotations, the name of the column of datatype \code{Date} designating the date of case report. e.g. "report_week"
#' @param strata In quotations, the name of the column indicating the stratifying variable.
#' @param moving_window Size of moving window for estimation of cases (numeric). The moving window size should be specified in the same date units as the reporting data (i.e. specify 7 to indicate 7 days, 7 weeks, etc). Default: NULL, i.e. takes all historical dates into consideration.
#' @param max_D Maximum possible delay observed or considered for estimation of the delay distribution (numeric). Default: (length of unique dates in time series)-1 ; or, if a moving window is specified, (size of moving window)-1
#' @param cutoff_D Consider only delays d<=\code{max_D}? Default: TRUE. If \code{cutoff_D=TRUE}, delays beyond \code{max_D} are ignored. If \code{cutoff_D=FALSE}, \code{max_D} is interpreted as delays>=\code{max_D} but within the moving window given by \code{moving_window}.
#' @param add_dow_cov Whether or not to add day-of-week covariates to the model
#' @param proportion_reported A decimal greater than 0 and less than or equal to 1 representing the proportion of all cases expected to be reported. Default: 1, e.g. 100 percent of all cases will eventually be reported. For asymptomatic diseases where not all cases will ever be reported, or for outbreaks in which severe under-reporting is expected, change this to less than 1.
#' @param quiet Suppress all output and progress bars from the JAGS process. Default: TRUE.
#' @param specs A list with arguments specifying the Bayesian model used: \code{dist} (Default: "Poisson"), \code{beta.priors} (Default: 0.1 for each delay d), \code{nSamp} (Default: 10000), \code{nBurnin} (Default: 1000), \code{nAdapt} (Default: 1000), \code{nChains} (Default: 1), \code{nThin} (Default: 1), \code{alphat.shape.prior} (Default: 0.001), \code{alphat.rate.prior} (Default: 0.001), \code{alpha1.mean.prior} (Default: 0), \code{alpha1.prec.prior} (Default: 0.001), \code{gamma.mean.prior} (Default: 0 for each day of the week (Monday-Saturday) - i.e. assuming initially no difference from Sunday), \code{gamma.prec.prior} (Default: 0.25 for each day of the week), \code{dispersion.prior} (Default: NULL, i.e. no dispersion. Otherwise, enter c(shape,rate) for a Gamma distribution.), \code{conf} (Default: 0.95), \code{quantiles} (Default: 5 quantiles for median, 50\% PI and 95\% PI), \code{param_names} (Default: NULL, i.e. output for all parameters is provided: c("lambda","alpha","beta.logged","tau2.alpha"). See McGough et al. 2019 (https://www.biorxiv.org/content/10.1101/663823v1) for detailed explanation of these parameters.).
#' @return The function returns a list with the following elements: \code{estimates}, a 5-column data frame containing estimates for each date in the window of predictions (up to "now") with corresponding date of case onset, lower and upper bounds of the prediction interval, and the number of cases for that onset date reported up to `now`. If quantiles is not NULL added columns will report the estimates for the requested quantiles;  \code{estimates.inflated}, a Tx4 data frame containing estimates inflated by the proportion_reported for each date in the time series (up to "now") with corresponding date of case onset, lower and upper bounds of the prediction interval, and the number of cases for that onset date reported up to `now`. If quantiles is not NULL added columns will report the inflated estimates for the requested quantiles; \code{nowcast.post.samples}, vector of 10,000 samples from the posterior predictive distribution of the nowcast, and \code{params.post}, a 10,000xN dataframe containing 10,000 posterior samples for the "N" parameters specified in specs[["param_names"]]. See McGough et al. 2019 (https://www.biorxiv.org/content/10.1101/663823v1) for detailed explanation of parameters.
#' @examples
#' # Load the data
#' data(denguedat)
#' # Perform stratified 'NobBS' assuming Poisson distribution, vague priors, and default
#' # specifications.
#' nowcast <- NobBS.strat(denguedat, as.Date("1990-02-05"),units="1 week",onset_date="onset_week",
#' report_date="report_week",strata="gender")
#' nowcast$estimates
#' @export
#' @import coda
#' @import rjags
#' @importFrom rlang sym
#' @importFrom dplyr last
#' @importFrom dplyr select
#' @importFrom dplyr starts_with
#' @importFrom dplyr group_by
#' @importFrom dplyr n
#' @importFrom dplyr summarise
#' @importFrom dplyr left_join
#' @importFrom magrittr %>%
#' @importFrom stats median quantile update
#'
#' @section Notes:
#' 'NobBS' requires that JAGS (Just Another Gibbs Sampler) is downloaded to the system.
#' JAGS can be downloaded at <http://mcmc-jags.sourceforge.net/>.

NobBS.strat <- function(data, now, units, onset_date, report_date, strata, moving_window=NULL, max_D=NULL, cutoff_D=NULL, add_dow_cov=FALSE, quiet=TRUE,
                        proportion_reported=1,
                        specs=list(
                          dist=c("Poisson","NB"),
                          alpha1.mean.prior=0,
                          alpha1.prec.prior=0.001,
                          alphat.shape.prior=0.001,
                          alphat.rate.prior=0.001,
                          beta.priors=NULL,
                          gamma.mean.prior=rep(0,6),
                          gamma.prec.prior=rep(0.25,6),
                          param_names=NULL,
                          conf=0.95,
                          quantiles=c(0.025,0.25,0.5,0.75,0.975),
                          dispersion.prior=NULL,
                          nAdapt=1000,
                          nChains=1,
                          nBurnin=1000,
                          nThin=1,
                          nSamp=10000)) {

  # Check that "now" is entered as a Date
  if(inherits(now, "Date")==FALSE){
    stop("'Now' argument must be of datatype Date (as.Date)")
  }

  # Check that "now" is possible in the sequence of reporting data
  if(dplyr::last(seq(min(data[,onset_date]),now,by=units))!=now){
    stop("The date `now` is not possible to estimate: the possible nowcast dates are seq(unique(data[,onset_date])[1],now,by=units).")
  }

  # Print date
  message(paste("Computing a nowcast for ",now))
  # Define "T", the length of dates between the first date of data and "now", making sure that "T" is unaffected by skipped-over dates in the time series
  # If the moving window is specified, "T" considers only the dates within the moving window; otherwise considers all historical data
  now.T <- ifelse(is.null(moving_window),length(seq(min(data[,onset_date]),as.Date(now),by=units)),
                  moving_window)

  # Define the strata
  strat <- unique(data[,strata])
  S <- length(strat)

  # Check the default arguments
  if (is.null(moving_window)) {
    moving_window <- now.T
  }
  if (is.null(max_D)) {
    max_D <- now.T-1 # ifelse(is.null(moving_window),now.T-1,moving_window-1)
  }
  if (is.null(cutoff_D)) {
    cutoff_D <- TRUE
  }
  if(quiet==TRUE){
    progress.bar <- "none"
  }
  if(quiet==FALSE){
    progress.bar <- "text"
  }

  # Check that proportion_reported is between 0,1
  if (proportion_reported > 1 | proportion_reported<=0){
    stop("The proportion_reported must be a number between (0,1].")
  }

  # Manipulate the control arguments
  if ("Poisson"%in%(specs[["dist",exact=TRUE]])) { # if no distribution specified, take Poisson as default
    specs$dist <- "Poisson"
  }
  if (is.null(specs[["dist",exact=TRUE]])) {
    specs$dist <- "Poisson"
  }
  if (is.null(specs[["alpha1.mean.prior",exact=TRUE]])) {
    specs$alpha1.mean.prior <- 0
  }
  if (is.null(specs[["alpha1.prec.prior",exact=TRUE]])) {
    specs$alpha1.prec.prior <- 0.001
  }
  if (is.null(specs[["alphat.shape.prior",exact=TRUE]])) {
    specs$alphat.shape.prior <- 0.001
  }
  if (is.null(specs[["alphat.rate.prior",exact=TRUE]])) {
    specs$alphat.rate.prior <- 0.001
  }
  if (is.null(specs[["beta.priors",exact=TRUE]])) {
    specs$beta.priors <- rep(0.1, times=(max_D)+1)
  }
  if (is.null(specs[["gamma.mean.prior",exact=TRUE]])) {
    specs$gamma.mean.prior <- rep(0,6)
  }
  if (is.null(specs[["gamma.prec.prior",exact=TRUE]])) {
    specs$gamma.prec.prior <- rep(0.25,6)
  }
  if (is.null(specs[["param_names",exact=TRUE]])&(specs[["dist"]]=="Poisson")) {
    specs$param_names <- c( "lambda","alpha","beta.logged","tau2.alpha","sum.n")
  }
  if (is.null(specs[["param_names",exact=TRUE]])&(specs[["dist"]]=="NB")) {
    specs$param_names <- c( "lambda","alpha","beta.logged","tau2.alpha","sum.n","r")
  }
  if (is.null(specs[["conf",exact=TRUE]])) {
    specs$conf <- 0.95
  }
  if (is.null(specs[["quantiles",exact=TRUE]])) {
    specs$quantiles <- c(0.025,0.25,0.5,0.75,0.975)
  }
  if (is.null(specs[["dispersion.prior",exact=TRUE]])&(specs[["dist"]]=="NB")) {
    specs$dispersion.prior <- c(0.001,0.001)
  }
  if (is.null(specs[["nAdapt",exact=TRUE]])) {
    specs$nAdapt <- 1000
  }
  if (is.null(specs[["nChains",exact=TRUE]])) {
    specs$nChains <- 1
  }
  if (is.null(specs[["nBurnin",exact=TRUE]])) {
    specs$nBurnin <- 1000
  }
  if (is.null(specs[["nThin",exact=TRUE]])) {
    specs$nThin <- 1
  }
  if (is.null(specs[["nSamp",exact=TRUE]])) {
    specs$nSamp <- 10000
  }

  # Warnings
  if(max_D>(moving_window-1)){
    stop("Maximum delay cannot be greater than the length of the moving window minus 1 time unit")
  }

  if(add_dow_cov==TRUE) {
    if(units!="1 day")
      stop("Day-of-week covariates can be added only with daily time unit")
  }

  # Prep the data: filter only to observable cases reported at or before "now"
  unit.num <- switch(units, "1 day"=1,"1 week"=7)
  w.days <- max((moving_window-1)*unit.num,(now.T-1)*unit.num) # moving window converted to days

  sel_rows <- (data[,onset_date]<=now) & (data[,onset_date]>=now-w.days) & (data[,report_date]<=now) & (data[,report_date]>=now-w.days)
  realtime.data <- subset(data,sel_rows)
  realtime.data$week.t <- (as.numeric(realtime.data[,onset_date]-min(realtime.data[,onset_date]))/unit.num)+1
  realtime.data$delay <- as.numeric(realtime.data[,report_date]-realtime.data[,onset_date])/unit.num

  if(cutoff_D==FALSE){
    realtime.data$delay <- ifelse(realtime.data$delay>=max_D,max_D,realtime.data$delay)
  }

  if(length(unique(realtime.data$week.t))!=now.T){
    warning("Warning! The line list has zero case reports for one or more possible onset dates at one or more delays. Proceeding under the assumption that the true number of cases at the associated delay(s) and week(s) is zero.")
  }

  if(!('batched' %in% names(realtime.data))) {
    realtime.data$batched <- FALSE
  }

  # Build the reporting triangle, fill with NAs where unobservable
  reporting.triangle.strat <- array(NA, dim=c(now.T,length(strat),(max_D+1)))
  batch_cases <- matrix(0, nrow=now.T, ncol=length(strat))
  for (t in 1:now.T) {
    for(s in 1:length(strat)){
      batch_cases[t,s] <- nrow(realtime.data[which(realtime.data$week.t==t & realtime.data[,strata] == strat[s] & realtime.data$batched==TRUE),])
      for (d in 0:max_D) {
        reporting.triangle.strat[t, s, (d + 1)] <- nrow(realtime.data[which(realtime.data$week.t == t &
                                                                              realtime.data$delay == d &
                                                                              realtime.data[,strata] == strat[s] &
                                                                              realtime.data$batched==FALSE), ])
        if (now.T < (t + d)) {
          reporting.triangle.strat[t, s, (d + 1)] <- NA
        }
      }
    }
  }

  if(add_dow_cov==TRUE) {
    X <- matrix(0, nrow = now.T, ncol = 7)
    onset_dates <- seq(now-now.T+1,now,by=1)
    row_index <- 1
    for (t in 1:now.T) {
      # Determine the day of the week as an index (1 = Monday, ..., 7 = Sunday)
      day_index <- as.numeric(format(onset_dates[t], "%u"))  # Returns 1 for Monday to 7 for Sunday
      # Set the corresponding day index in the matrix to 1
      X[row_index, day_index] <- 1
      # Move to the next row
      row_index <- row_index + 1
    }
  }

  # Run the JAGS model

  if(specs[["dist"]]=="Poisson"){
    params=c( "lambda","alpha","beta.logged","tau2.alpha","n","sum.n","sum.lambda")
  }
  if(specs[["dist"]]=="NB"){
    params=c( "lambda","alpha","beta.logged","tau2.alpha","n","sum.n","sum.lambda","r")
  }
  nAdapt = specs[["nAdapt"]] #default = 1000
  nChains = specs[["nChains"]] # default=1
  nBurnin = specs[["nBurnin"]] # default=1000
  nThin = specs[["nThin"]] # default=1
  nKeep = specs[["nSamp"]] # default=10,000
  nIter = nKeep * nThin

  if(specs[["dist"]]=="Poisson"){
    dataList = list(Today = now.T,
                    D = max_D,
                    S = length(strat),
                    n = reporting.triangle.strat,
                    m = batch_cases,
                    alpha1.mean.prior=specs$alpha1.mean.prior,
                    alpha1.prec.prior=specs$alpha1.prec.prior,
                    alphat.rate.prior=specs$alphat.rate.prior,
                    alphat.shape.prior=specs$alphat.shape.prior,
                    beta.priors=specs$beta.priors)
  }

  if(specs[["dist"]]=="NB"){
    dataList = list(Today = now.T,
                    D = max_D,
                    S = length(strat),
                    n = reporting.triangle.strat,
                    m = batch_cases,
                    alpha1.mean.prior=specs$alpha1.mean.prior,
                    alpha1.prec.prior=specs$alpha1.prec.prior,
                    alphat.rate.prior=specs$alphat.rate.prior,
                    alphat.shape.prior=specs$alphat.shape.prior,
                    beta.priors=specs$beta.priors,
                    dispersion.prior.shape=specs$dispersion.prior[1],
                    dispersion.prior.rate=specs$dispersion.prior[2])
  }

  if(add_dow_cov==TRUE){
    params=c(params,'gamma')
    specs$param_names <- c(specs$param_names,'gamma')
    dataList$gamma.mean.prior=specs$gamma.mean.prior
    dataList$gamma.prec.prior=specs$gamma.prec.prior
    dataList$X <- X
  }

  if(add_dow_cov==TRUE){
    JAGSmodPois_strat <- (system.file("JAGS", "nowcastPois_strata_withDowCov.txt", package="NobBS"))
    JAGSmodNB_strat <- (system.file("JAGS", "nowcastNB_strata_withDowCov.txt", package="NobBS"))
  } else {
    JAGSmodPois_strat <- (system.file("JAGS", "nowcastPois_strata.txt", package="NobBS"))
    JAGSmodNB_strat <- (system.file("JAGS", "nowcastNB_strata.txt", package="NobBS"))
  }

  nowcastmodel = jags.model(
    file = ifelse(specs[["dist"]]=="Poisson",JAGSmodPois_strat,JAGSmodNB_strat), # test "nowcastPois_strata.txt"
    data = dataList,
    n.chains = nChains,
    n.adapt = nAdapt,
    inits=list(.RNG.seed=1,.RNG.name="base::Super-Duper"),
    quiet=quiet)

  update( object = nowcastmodel, n.iter = nBurnin , progress.bar = progress.bar)

  lambda.output = coda.samples(
    model = nowcastmodel,
    variable.names =  if("sum.n"%in%specs$param_names) c(specs$param_names) else c(specs$param_names,"sum.n"),
    n.iter = nIter,
    thin = nThin,
    quiet=quiet,
    progress.bar=progress.bar)

  mymod.mcmc <- as.mcmc(lambda.output)
  mymod.dat <- as.data.frame(as.matrix(mymod.mcmc))

  # Extract all hindcasts and 95% credible intervals
  t.extract <- (now.T-(now.T-1)):(now.T) # nowcast all weeks up through the present

  estimates <- array(NA, dim = c(now.T,4,S), dimnames=list(NULL,c("estimate","lower","upper","stratum"),strat))
  for(v in t.extract){
    for(s in 1:S){
      estimates[v,1,s] <- median(mymod.dat[, grep(paste("sum.n[",v,",",s,"]",sep=""), colnames(mymod.dat), fixed=TRUE)])
      estimates[v,2,s] <- quantile((mymod.dat[, grep(paste("sum.n[",v,",",s,"]",sep=""), colnames(mymod.dat), fixed=TRUE)]),probs = c((1-specs$conf)/2,1-((1-specs$conf)/2)))[1]
      estimates[v,3,s] <- quantile((mymod.dat[, grep(paste("sum.n[",v,",",s,"]",sep=""), colnames(mymod.dat), fixed=TRUE)]),probs = c((1-specs$conf)/2,1-((1-specs$conf)/2)))[2]
      estimates[v,4,s] <- strat[s]
    }
  }

  # Estimates inflated by proportion reported
  estimates.inflated <- array(NA, dim = c(now.T,4,S), dimnames=list(NULL,c("estimate","lower","upper","stratum"),strat))
  for(v in t.extract){
    for(s in 1:S){
      estimates.inflated[v,1,s] <- median(mymod.dat[, grep(paste("sum.n[",v,",",s,"]",sep=""), colnames(mymod.dat), fixed=TRUE)])/proportion_reported
      estimates.inflated[v,2,s] <- quantile((mymod.dat[, grep(paste("sum.n[",v,",",s,"]",sep=""), colnames(mymod.dat), fixed=TRUE)]),probs = c((1-specs$conf)/2,1-((1-specs$conf)/2)))[1]/proportion_reported
      estimates.inflated[v,3,s] <- quantile((mymod.dat[, grep(paste("sum.n[",v,",",s,"]",sep=""), colnames(mymod.dat), fixed=TRUE)]),probs = c((1-specs$conf)/2,1-((1-specs$conf)/2)))[2]/proportion_reported
      estimates.inflated[v,4,s] <- strat[s]
    }
  }

  estimates <- data.frame(apply(estimates,2,rbind))
  estimates.inflated <- data.frame(apply(estimates.inflated,2,rbind))

  estimates[,c(1:3)] <-  sapply(estimates[,c(1:3)], function(x) as.numeric(as.character(x)))
  estimates.inflated[,c(1:3)] <-  sapply(estimates.inflated[,c(1:3)], function(x) as.numeric(as.character(x)))

  if(!is.null(specs$quantiles)) {
    estimates.quantiles <- array(NA, dim = c(now.T,length(specs$quantiles)+1,S), dimnames=list(NULL,c(paste0('q_',specs$quantiles),"stratum"),strat))
    estimates.quantiles.inflated <- array(NA, dim = c(now.T,length(specs$quantiles)+1,S), dimnames=list(NULL,c(paste0('q_',specs$quantiles),"stratum"),strat))
    for(v in t.extract){
      for(s in 1:S){
        q.vals = quantile((mymod.dat[, grep(paste("sum.n[",v,",",s,"]",sep=""), colnames(mymod.dat), fixed=TRUE)]),probs = specs$quantiles)
        estimates.quantiles[v,,s] <- c(q.vals,strat[s])
        estimates.quantiles.inflated[v,,s] <- c(q.vals/proportion_reported,strat[s])
      }
    }

    estimates.quantiles <- data.frame(apply(estimates.quantiles,2,rbind))
    estimates.quantiles.inflated <- data.frame(apply(estimates.quantiles.inflated,2,rbind))

    estimates.quantiles[,c(1:length(specs$quantiles))] <- sapply(estimates.quantiles[,c(1:length(specs$quantiles))], function(x) as.numeric(as.character(x)))
    estimates.quantiles.inflated[,c(1:length(specs$quantiles))] <- sapply(estimates.quantiles.inflated[,c(1:length(specs$quantiles))], function(x) as.numeric(as.character(x)))

    estimates <- cbind(estimates,estimates.quantiles)
    estimates.inflated <- cbind(estimates.inflated,estimates.quantiles.inflated)
  }

  # Combine nowcast estimates with: dates, number of cases reported at each date and strata
  reported <- data.frame(
    realtime.data %>%
      dplyr::group_by(!!sym(onset_date),!!sym(strata)) %>%
      dplyr::summarise(n.reported=dplyr::n())
  )
  names(reported)[1] <- "onset_date"
  names(reported)[2] <- "stratum"

  estimates <- data.frame(estimates, onset_date=rep(seq(as.Date(now)-w.days,as.Date(now),by=units),times=S)) %>%
    dplyr::left_join(reported,by=c("onset_date","stratum"))

  estimates.inflated <- data.frame(estimates.inflated, onset_date=rep(seq(as.Date(now)-w.days,as.Date(now),by=units),times=S)) %>%
    dplyr::left_join(reported,by=c("onset_date","stratum"))

  t <- now.T

  parameter_extract <- matrix(NA, nrow=10000)

  if("lambda"%in%specs$param_names){
    parameter_extract <- cbind(parameter_extract,mymod.dat %>% dplyr::select(starts_with(paste("lambda[",t,",",sep=""))))
  }
  if("beta.logged"%in%specs$param_names){
    betas.logged<- matrix(NA,nrow=10000,ncol=(max_D+1))
    dimnames(betas.logged) = list(NULL,c(paste("Beta",c(0:max_D))))
    for(d in 0:max_D){
      betas.logged[,(d+1)] <- (mymod.dat %>% dplyr::select(starts_with(paste("beta.logged[",(d+1),"]",sep=""))))[,1]
    }
    parameter_extract <- cbind(parameter_extract,betas.logged)
  }
  if("alpha"%in%specs$param_names){
    parameter_extract <- cbind(parameter_extract,mymod.dat %>% dplyr::select(starts_with(paste("alpha[",t,sep=""))))
  }
  if("tau2.alpha"%in%specs$param_names){
    parameter_extract <- cbind(parameter_extract,mymod.dat %>% dplyr::select(starts_with("tau2.alpha")))
  }

  if("gamma"%in%specs$param_names){
    parameter_extract <- cbind(parameter_extract,mymod.dat %>% dplyr::select(starts_with("gamma")))
  }

  nowcast.post.samps <- (mymod.dat %>% dplyr::select(starts_with(paste("sum.n[",t,sep=""))))[,1]

  list(estimates=estimates,estimates.inflated=estimates.inflated,
       nowcast.post.samps=nowcast.post.samps,params.post=parameter_extract[,2:ncol(parameter_extract)])

}
