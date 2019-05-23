#' Produce smooth Bayesian nowcasts of incomplete, time-stamped reporting data.
#' 
#' @param data A time series of reporting data in line list format (one row per case), with a column \code{onset_date} indicating date of case onset, and a column \code{report_date} indicating date of case report.
#' @param now An object of class \code{Date} indicating the date at which to perform the nowcast.
#' @param units Time scale of reporting. Options: "1 day", "1 week".
#' @param onset_date In quotations, the name of the column of class \code{Date} designating the date of case onset. e.g. "onset_week"
#' @param report_date In quotations, he name of the column of class \code{Date} designating the date of case report. e.g. "report_week"
#' @param moving_window Size of moving window for estimation of cases. The moving window size should be specified in the same date units as the reporting data (i.e. specify 7 to indicate 7 days, 7 weeks, etc). Default: NULL, i.e. takes all historical dates into consideration.
#' @param max_D Maximum possible delay observed or considered for estimation of the delay distribution. Default: (length of unique dates in time series)-1
#' @param cutoff_D Consider only delays d<=\code{max_D}? Default: TRUE. If \code{cutoff_D=TRUE}, delays beyond \code{max_D} are ignored. If \code{cutoff_D=FALSE}, \code{max_D} is interpreted as delays>=\code{max_D} but within the moving window given by \code{moving_window}.
#' @param specs A list with arguments specifying the Bayesian model used: \code{dist} (Default: "Poisson"), \code{beta.priors} (Default: 0.1 for each delay d), \code{nSamp} (Default: 10000), \code{nBurnin} (Default: 1000), \code{nAdapt} (Default: 1000), \code{nChains} (Default: 1), \code{nChains} (Default: 1), \code{nThin} (Default: 1), \code{alphat.shape.prior} (Default: 0.001), \code{alphat.rate.prior} (Default: 0.001), \code{alpha1.mean.prior} (Default: 0), \code{alpha1.prec.prior} (Default: 0.001), \code{dispersion.prior} (Default: NULL, i.e. no dispersion. Otherwise, enter c(shape,rate) for a Gamma distribution.), \code{conf} (Default: 0.95), \code{param_names} (Default: NULL, i.e. output for all parameters is provided)
#' @return \code{NobBS} returns the list \code{nowcast_results} with the following elements: \code{estimates}, a Tx3 data frame containing estimates for each date in the time series (up to "now") with corresponding prediction intervals; \code{nowcast.post.samples}, vector of 10,000 samples from the posterior predictive distribution of the nowcast, and \code{params.post}, a 10,000x5 dataframe containing 10,000 posterior samples for each of the following parameters: \code{beta0}, \code{beta1}, \code{beta2}, \code{alpha[T]}, \code{tau2-alpha}.
#' @examples
#' Load the data
#' data(denguedat)
#' Perform default NobBS assuming Poisson distribution, vague priors, and default specifications.
#' NobBS(denguedat, as.Date("1990-10-01"),units="1 week",onset_date="onset_week",report_date="report_week")
#' @export
#' @import coda
#' @import rjags
#' @importFrom dplyr select

NobBS <- function(data, now, units, onset_date, report_date, moving_window=NULL, max_D=NULL, cutoff_D=NULL, 
                  specs=list(
                    dist=c("Poisson","NB"),
                    alpha1.mean.prior=0,
                    alpha1.prec.prior=0.001,
                    alphat.shape.prior=0.001,
                    alphat.rate.prior=0.001,
                    beta.priors=NULL,
                    param_names=NULL,
                    conf=0.95,
                    dispersion.prior=NULL,
                    nAdapt=1000,
                    nChains=1,
                    nBurnin=1000,
                    nThin=1,
                    nSamp=10000)) {
  # Print date
  print(paste("Computing a nowcast for ",now))
  # Define "T", the length of dates between the first date of data and "now", making sure that "T" is unaffected by skipped-over dates in the time series
  now.T <- length(seq(min(data[,onset_date]),as.Date(now),by=units))
  
  # Check the default arguments
  if (is.null(max_D)) {
    max_D <- now.T-1
  }
  if (is.null(cutoff_D)) {
    cutoff_D <- TRUE
  }
  
  # Manipulate the control arguments
  if ("Poisson"%in%(specs[["dist",exact=TRUE]])) { # if no distribution specified, take Poisson as default
    specs$dist <- "Poisson"
  } 
  if (is.null(specs[["beta.priors",exact=TRUE]])) {
    specs$beta.priors <- rep(0.1, times=(max_D)+1)
  }
  if (is.null(specs[["param_names",exact=TRUE]])&(specs[["dist"]]=="Poisson")) {
    specs$param_names <- c( "lambda","alpha","pi.logged","tau2.alpha","n","pi","sum.n","sum.lambda")
  }
  if (is.null(specs[["param_names",exact=TRUE]])&(specs[["dist"]]=="NB")) {
    specs$param_names <- c( "lambda","alpha","pi.logged","tau2.alpha","n","pi","sum.n","sum.lambda","r")
  }
  if (is.null(specs[["dispersion.prior",exact=TRUE]])&(specs[["dist"]]=="NB")) {
    specs$dispersion.prior <- c(0.001,0.001)
  }
  
  # Prep the data: filter only to observable cases reported at or before "now"
  unit.num <- switch(units, "1 day"=1,"1 week"=7)
  w.days <- max((moving_window-1)*unit.num,(now.T-1)*unit.num) # moving window converted to days
  
  realtime.data <- subset(data,(data[,onset_date]<=now) & (data[,onset_date]>=now-w.days) & (data[,report_date]<=now) & (data[,report_date]>=now-w.days))
  realtime.data$week.t <- (as.numeric(realtime.data$onset_week-min(realtime.data$onset_week))/unit.num)+1
  realtime.data$delay <- as.numeric(realtime.data[,report_date]-realtime.data[,onset_date])/unit.num
  
  if(cutoff_D==FALSE){
    realtime.data$delay <- ifelse(realtime.data$delay>=max_D,max_D,realtime.data$delay)
  }
  
  if(range(realtime.data$week.t)[2]!=now.T){
    print("stop!")
    break
  }
  
  # Build the reporting triangle, fill with NAs where unobservable
  reporting.triangle <- matrix(NA, nrow=now.T,ncol=(max_D+1))
  
  for(t in 1:now.T){
    for(d in 0:max_D){
      reporting.triangle[t,(d+1)] <- nrow(realtime.data[which(realtime.data$week.t==t & realtime.data$delay==d),])
      if(now.T < (t+d)){
        reporting.triangle[t,(d+1)] <- NA 
      }
    }
  }
  
  # Run the JAGS model
  
  if(specs[["dist"]]=="Poisson"){
    params=c( "lambda","alpha","pi.logged","tau2.alpha","n","pi","sum.n","sum.lambda")
  }
  if(specs[["dist"]]=="NB"){
    params=c( "lambda","alpha","pi.logged","tau2.alpha","n","pi","sum.n","sum.lambda","r")
  }
  nAdapt = specs[["nAdapt"]] #default = 1000
  nChains = specs[["nChains"]] # default=1
  nBurnin = specs[["nBurnin"]] # default=1000
  nThin = specs[["nThin"]] # default=1
  nKeep = specs[["nSamp"]] # default=10,000
  nIter = nKeep * nThin
  
  if(specs[["dist"]]=="Poisson"){
    dataList = list(T = now.T,
                    D = max_D,
                    n = reporting.triangle,
                    alpha1.mean.prior=specs$alpha1.mean.prior,
                    alpha1.prec.prior=specs$alpha1.prec.prior,
                    alphat.rate.prior=specs$alphat.rate.prior,
                    alphat.shape.prior=specs$alphat.shape.prior,
                    beta.priors=specs$beta.priors)
  }
  
  if(specs[["dist"]]=="NB"){
    dataList = list(T = now.T,
                    D = max_D,
                    n = reporting.triangle,
                    alpha1.mean.prior=specs$alpha1.mean.prior,
                    alpha1.prec.prior=specs$alpha1.prec.prior,
                    alphat.rate.prior=specs$alphat.rate.prior,
                    alphat.shape.prior=specs$alphat.shape.prior,
                    beta.priors=specs$beta.priors,
                    dispersion.prior.shape=specs$dispersion.prior[1],
                    dispersion.prior.rate=specs$dispersion.prior[2])
  }
  
  JAGSmodPois <- (system.file("JAGS", "nowcastPois.txt", package="NobBS")) # file.path(path.package('NobBS'),"nowcastPois.txt")
  JAGSmodNB <- (system.file("JAGS", "nowcastNB.txt", package="NobBS")) #file.path(path.package('NobBS'),"nowcastNB.txt")
  
  nowcastmodel = jags.model(
    file = ifelse(specs[["dist"]]=="Poisson",JAGSmodPois,JAGSmodNB),
    data = dataList,
    n.chains = nChains,
    n.adapt = nAdapt,
    inits=list(.RNG.seed=1,.RNG.name="base::Super-Duper"))
  
  update( object = nowcastmodel, n.iter = nBurnin ) 
  
  lambda.output = coda.samples(
    model = nowcastmodel,
    variable.names =  specs$param_names, #c("lambda","alpha","pi.logged","tau2.alpha","n","pi","sum.nt","sum.lambda"), # extract only lambda[t,d] of week now w/ 0 delay   ## #params,
    n.iter = nIter,
    thin = nThin)
  
  mymod.mcmc <- as.mcmc(lambda.output)
  mymod.dat <- as.data.frame(as.matrix(mymod.mcmc))
  
  # Extract all hindcasts and 95% credible intervals
  t.extract <- (now.T-(now.T-1)):(now.T) # nowcast all weeks up through the present
  
  estimates <- matrix(NA, ncol=3, nrow=now.T,dimnames=list(NULL,c("estimate","lower","upper")))
  for(v in t.extract){
    estimates[v,1] <- median(mymod.dat[, grep(paste("sum.n[",v,"]",sep=""), colnames(mymod.dat), fixed=T)])
    estimates[v,2] <- quantile((mymod.dat[, grep(paste("sum.n[",v,"]",sep=""), colnames(mymod.dat), fixed=T)]),probs = c((1-specs$conf)/2,1-((1-specs$conf)/2)))[1]
    estimates[v,3] <- quantile((mymod.dat[, grep(paste("sum.n[",v,"]",sep=""), colnames(mymod.dat), fixed=T)]),probs = c((1-specs$conf)/2,1-((1-specs$conf)/2)))[2]
  }
  
  estimates <- as.data.frame(estimates)
  
  t <- now.T 
  
  pi.logged.td1 <- mymod.dat %>% dplyr::select(select_vars(names(mymod.dat),starts_with(paste("pi.logged[1]",sep=""))))
  pi.logged.td2 <- mymod.dat %>% dplyr::select(select_vars(names(mymod.dat),starts_with(paste("pi.logged[2]",sep=""))))
  pi.logged.td3 <- mymod.dat %>% dplyr::select(select_vars(names(mymod.dat),starts_with(paste("pi.logged[3]",sep=""))))
  alpha.last <- mymod.dat %>% dplyr::select(select_vars(names(mymod.dat),starts_with(paste("alpha[",t,sep=""))))
  tau2.alpha <- mymod.dat %>% dplyr::select(select_vars(names(mymod.dat),starts_with("tau2.alpha")))
  
  parameter_extract <- cbind(pi.logged.td1,pi.logged.td2,pi.logged.td3,
                             alpha.last,tau2.alpha)
  
  nowcast.post.samps <- (mymod.dat %>% dplyr::select(select_vars(names(mymod.dat),starts_with(paste("sum.n[",t,sep="")))))[,1]
  
  nowcast_results <<- list(estimates=estimates, nowcast.post.samps=nowcast.post.samps,params.post=parameter_extract)
  
}