source('~/R/lib_2020/rpact_functions.r')

sim_surv1 <- function(n_total, readoutTime,  m_b=10, hr=trt_hr, group_text=''){
  # Create a data frame with the subject IDs and treatment covariate
  cov <- data.frame(id = 1:n_total,
                    treatmentGroup = 0)
  cov$treatmentGroup[sample(n_total, size = n_total/2)] <- 1   
  
  # Simulate the event times
    dat <- simsurv(dist= 'exponential', lambdas = log(2)/m_b, 
                    
                   betas = c(treatmentGroup = log( hr)), 
                   x = cov, 
                   maxt = readoutTime) 
  
  # Merge the simulated event times onto covariate data frame
  dat <- merge(cov, dat)
  dat$arm <- factor(dat$treatmentGroup, levels=c(0, 1), labels = c('B','A'))
  dat$grp <- group_text 
  colnames(dat) <- gsub('eventtime','tte', colnames(dat))
  colnames(dat) <- gsub('status','event', colnames(dat))
  dat
}

sim_rpact2 <- function(n_total, readoutTime,  m_b=10, hr=trt_hr, group_text='',  ai= seq(1, 10,by=1), epr_sim=0.8 ){
  # n_total is the number of patients in one subgroup , 
  
  if(n_total == 0) return (NULL);
  
  d1 <- NULL
  j <- 1
  while( (is.null(d1) || nrow(d1) < n_total) && j < 10){
    j <- j+1
    if( j > 2) print(paste('iteration',j,'current ai is',paste(ai,collapse = '.')))
  trial1 <- getSimulationSurvival(lambda2 = log(2) / m_b, lambda1 = log(2) / ( m_b / hr ), dropoutRate1 = .05, dropoutRate2 = .05, 
            dropoutTime = 24, 
            accrualTime = (1:length(ai)) -1 , accrualIntensity = ai  , plannedEvents = round( n_total * epr_sim, 0),  maxNumberOfSubjects = n_total,    maxNumberOfRawDatasetsPerStage = 1, longTimeSimulationAllowed=T) # simulate a trial with EPR 90%
  
  d1 <- as.data.frame( getRawData(trial1) )
  ai <- ai/j
  }
  
  if(nrow(d1) < n_total){
    stop( paste( 'required number of events reached despite reducing enrollment speed to', paste(ai, collapse = '.')))
    }
 # print(head(d1))
  d1$arm <- factor(d1$treatmentGroup, levels = c(2,1), labels = c('B','A'))
  d1$grp <- group_text
  colnames(d1) <- gsub('timeUnderObservation','tte', colnames(d1))
  if( d1[1,'observationTime'] < readoutTime) {
    print( paste('read out time at', readoutTime,'which is later than the study follow up time:', d1[1,'observationTime']))
    return (NULL)
  }  
  
  cut_early(d1, readoutTime)
	
}

qad_hr <- function(x) { # quick and dirty estimates
 ( x[1,'Expected'] / x[-1,'Expected'] ) / 
    ( x[1,'Observed'] / x[-1,'Observed'] )
  
}

# modify the original print.survdiff to output statistics
summary.survdiff <- function(x, digits = max(options()$digits - 4, 3), ...) {

    saveopt <-options(digits=digits)
    on.exit(options(saveopt))

    if (!inherits(x, 'survdiff'))
	stop("Object is not the result of survdiff")
    if (!is.null(cl<- x$call)) {
	# cat("Call:\n")
	# dput(cl)
	# cat("\n")
	}

    omit <- x$na.action
    if (length(omit)) cat("n=", sum(x$n), ", ", naprint(omit),
					  ".\n\n", sep='')

    if (length(x$n)==1)  {
	z <- sign(x$exp - x$obs) * sqrt(x$chisq)
	temp <- c(x$obs, x$exp, z, signif(pchisq(x$chisq, 1, lower.tail=FALSE),
                                          digits))
	names(temp) <- c("Observed", "Expected", "Z", "p")
	print(temp)
	}
    else {
	if (is.matrix(x$obs)){
	    otmp <- apply(x$obs,1,sum)
	    etmp <- apply(x$exp,1,sum)
	    }
	else {
	    otmp <- x$obs
	    etmp <- x$exp
	    }
	df <- (sum(1*(etmp>0))) -1
	temp <- cbind(x$n, otmp, etmp, ((otmp-etmp)^2)/ etmp,
					 ((otmp-etmp)^2)/ diag(x$var), diag(x$var), 
					 pchisq(x$chisq, df, lower.tail=FALSE))
	dimnames(temp) <- list(names(x$n), c("N", "Observed", "Expected",
				  "(O-E)^2/E", "(O-E)^2/V",'V','P'))
	temp <- cbind(temp, qad_hr=qad_hr(temp))
	return(temp)
# 	cat("\n Chisq=", format(round(x$chisq,1)),
# 		 " on", df, "degrees of freedom, p=",
# 		 format.pval(pchisq(x$chisq, df, lower.tail=FALSE)),
#             "\n")
       }
#     invisible(x)
    }


ana1 <- function(d, form  ){
  if(nrow(d)== 0) return(matrix(nrow = 1, ncol=4))
  stopifnot(all(c('tte','event','arm') %in% colnames(d)))
  stopifnot(sum( d[,'arm'] %in% 'B') > 0)
  
  s1 <- summary( survdiff( as.formula(form)  , data =   d) )
	 
#   if(grepl('+', form)) {
#      c(#  cox.table2(ds1[sel, ], tte='tte', cens='event', form='arm'),
# 	   n = nrow(d), nevent = sum(d[, 'event']), p=s1[1,'P'])
#   }
#   else{  
	   c(#  cox.table2(ds1[sel, ], tte='tte', cens='event', form='arm'),
	   n = nrow(d), nevent = sum(d[, 'event']), 
	   `O-E` = s1['arm=B','Observed'] - s1['arm=B','Expected'],
	   V = s1['arm=B','V'], p=s1[1,'P'], qad_hr=s1[1,'qad_hr']) 
#  }
}


check_para <- function(p) {
  print(unlist(p)) 
  if (is.null(p$fractionActiveArm)) p$fractionActiveArm <- 0.5 # handle unequal randomization
  
  stopifnot(all(sapply(p$n, function(x) {
    x1 <- x * p$fractionActiveArm 
    abs(as.integer(x) - x) < 1e-6  && x >= 2  &&  abs(as.integer(x1) - x1) < 1e-6 
    })  )  )
  stopifnot(length(p$n) == length(p$hr))
  stopifnot(length(p$n) == length(p$med))
  if(!is.null(p$nEvent)) stopifnot(sum(p$n) >= p$nEvent)
}


collapse_grp <- function(dsIn, minEvent = 10, rule = c('collapseGroup','removeStratum') ) {
  stopifnot(all(c('grp','event') %in% colnames(dsIn)))
  stopifnot(length(rule) == 1)
    e1 <- with(dsIn, tapply(event,grp, sum))
    if(any(e1) < minEvent){
    if(rule == 'collapseGroup'){
      group2Collapse <- names(e1)[e1 < minEvent]
      if(sum(e1[group2Collapse]) >= minEvent) {
        dsIn[dsIn$grp %in% group2Collapse,'grp'] <- 'MCollapse'
      }else{
        e2 <- e1[setdiff(names(e1), group2Collapse)]
        nextSmallest <- names(e2)[which.min(e2)]
        dsIn[dsIn$grp %in% group2Collapse,'grp'] <- nextSmallest
      }
    }
      
    if(rule == 'removeStratum') { # a very specific rule to remove the 2nd group; not general at all
      dsIn[dsIn$grp %in% c('M1','M3'),'grp'] <- c('M13') 
      dsIn[dsIn$grp %in% c('M2','M4'),'grp'] <- c('M24')
      dsIn[dsIn$grp %in% c('M5','M7'),'grp'] <- c('M57') 
      dsIn[dsIn$grp %in% c('M6','M8'),'grp'] <- c('M68')
    }
  }
    
    dsIn
  }


# leave collapse_grp() alone for imp150 simulation. the new function collapse_grp2() will be a general function
collapse_grp2 <- function(dsIn, minN = 10, rule = c('collapseGroup','removeStratum'), sumFun=sum,  
                          strat.order ) 
  {
  # if sumFun is sum, the decision to collapse strata is based on number of events
  # if sumfun is length, the decision to collapse is based on number of patients
  
  
  # strat.order gives a vector of strata (M1 etc) with an order to removed, eg(c('M1','M2','M3'))
  # markerList defines each stratum. e.g. c('M1-M2-M3-M4-','M1+M2-M3-M4-', ...), it must be in the same order
  
  stopifnot(all(c('grp','event') %in% colnames(dsIn)))
  stopifnot(length(rule) == 1)
  dsIn$grp <- as.character(dsIn$grp)
    e1 <- with(dsIn, tapply(event,grp, sumFun))
    if(rule == 'collapseGroup'){
          
        if(any(e1) < minN){
    
          group2Collapse <- names(e1)[e1 < minN]
          if(sum(e1[group2Collapse]) >= minN) {
            dsIn[dsIn$grp %in% group2Collapse,'grp'] <- 'MCollapse'
          }else{
            e2 <- e1[setdiff(names(e1), group2Collapse)]
            nextSmallest <- names(e2)[which.min(e2)]
            dsIn[dsIn$grp %in% group2Collapse,'grp'] <- nextSmallest
          }
        }
    }
      
    if(rule == 'removeStratum') {  
      
      grp <- dsIn$grp
      e2 <- e1
      j <- 1
      while(any(e2  < minN)) {
           grp <- gsub(paste(strat.order[j],'(\\+|\\-)', sep =''),'', grp)
           e2 <-   tapply(dsIn$event, grp, sumFun)
           j <- j+1
      }
      dsIn$grp <- grp
    }
 
    
    dsIn
}




sim4 <- function( p  , time , maxPms = 10, rampUpTime, estimatePrognostic=F ) {
# sim4 is similar to sim3, but it adds rules to collapse stratification groups: if the number of event is fewer than minE, the group will be collapsed
  
  #ai0 <- with(p, rep(1:maxPms, each=enrollmentSpeadFactor) * dup )
  if(missing(rampUpTime)) rampUpTime <- round( time /5 )
  ai0 <- with(p, seq(1, maxPms, length.out = rampUpTime) * dup )
  ds1 <- do.call(rbind, lapply(1:length(p$n), function(i) {
    	ds1a <-  sim_rpact2(n_total = p$n[i] * p$dup ,    ai=ai0 * p$n[i]/sum(p$n),   readoutTime=time, m_b = p$med[i], hr = p$hr[i], group_text=paste('M',i, sep='') )
    	
    		if( is.null(ds1a)  ) return(rep(NA, 11))
    	  ds1a$dup <- NA
    	  for (a in unique(ds1a$treatmentGroup)) {
		        ds1a[ds1a$treatmentGroup == a,'dup'] <- sample( rep(1:p$dup, p$n[i]/2) )
		
	      } # enforce the stratified randomization
    	  
    	  ds1a
    }) )
   
  
	stopifnot(nrow(ds1) == with(p, sum(n) *dup))
	c1 <- with(ds1, table(treatmentGroup,  dup, grp, useNA = 'ifany'))
	for (i in 1:length(p$n)) {
   stopifnot(all(as.vector(c1[,, i]) == p$n[i] /2  ))
	}
	 
	res0 <- t ( sapply(split(ds1, paste(ds1$dup)), ana1, form= 'Surv(tte, event) ~ arm') )
	res2 <-  t ( sapply(split(ds1, paste(ds1$dup)), function(d) {
	  d2 <- collapse_grp(d,minEvent=p$minE, rule=p$collapseRule)
 
	  c( ana1( d = d2 , form= 'Surv(tte, event) ~ arm + strata(grp)'), n_strata = length(unique(d2$grp)))
	  }))
	res_s <- t ( sapply(split(ds1, paste(ds1$dup)), ana1, form='Surv(tte, event) ~ arm + strata(grp)') )
	
	colnames(res2) <- paste(colnames(res2),'c', sep='_')
	colnames(res_s) <- paste(colnames(res_s), 's', sep = '_')
	#with(ds1, tapply(tte, list(arm, grp), median))
	
 stopifnot(all(row.names(res0) == row.names(res_s)))
 
	cbind(res0[, c('n','nevent','p','qad_hr')], res_s[, c('p_s','qad_hr_s')], res2[, c('p_c','qad_hr_c','n_strata_c')], delta_p = res0[,'p'] - res_s[,'p_s'], delta_p2 = res2[,'p_c'] - res_s[,'p_s'] )

}



sim3 <- function( p  , time , maxPms = 10, rampUpTime,epr_sim=.8, estimatePrognostic=F) {
    # a larger enrollmentSpeadFactor means slower wrap up
  #enrollmentSpeadFactor = 2
  #ai0 <- with(p, rep(1:maxPms, each=enrollmentSpeadFactor) * dup )
  if(missing(rampUpTime)) rampUpTime <- round( time /5 )
  ai0 <- with(p, seq(1, maxPms, length.out = rampUpTime) * dup )
  ds1 <- do.call(rbind, lapply(1:length(p$n), function(i) {
    	ds1a <-  sim_rpact2(n_total = p$n[i] * p$dup ,    ai=ai0 * p$n[i]/sum(p$n),   readoutTime=time, m_b = p$med[i], hr = p$hr[i], group_text=paste('M',i, sep=''), epr_sim=epr_sim )
    	
    		if( is.null(ds1a)  ) return(rep(NA, 10))
    	  ds1a$dup <- NA
    	  for (a in unique(ds1a$treatmentGroup)) {
		        ds1a[ds1a$treatmentGroup == a,'dup'] <- sample( rep(1:p$dup, p$n[i]/2) )
		
	      } # enforce the stratified randomization
    	  
    	  ds1a
    }) )
   
  
	stopifnot(nrow(ds1) == with(p, sum(n) *dup))
	c1 <- with(ds1, table(treatmentGroup,  dup, grp, useNA = 'ifany'))
	for (i in 1:length(p$n)) {
   stopifnot(all(as.vector(c1[,, i]) == p$n[i] /2  ))
	}
	 
	res0 <- t ( sapply(split(ds1, paste(ds1$dup)), ana1, form= 'Surv(tte, event) ~ arm') )
	res1 <- t ( sapply(split(ds1, paste(ds1$dup, ds1$grp)), ana1, form= 'Surv(tte, event) ~ arm') )
	res1 <- as.data.frame(res1)
	nv <- strsplit(row.names(res1), split =  ' ')
	res1$dup <- sapply(nv,'[',1)
	res1$grp <- sapply(nv,'[',2)
	colnames(res1) <- make.names(colnames(res1))
	res1b <-  dcast(res1, dup ~ grp, value.var='O.E')
	colnames(res1b) <- paste('O-E', colnames(res1b),sep ='_')
	res_s <- t ( sapply(split(ds1, paste(ds1$dup)), ana1, form='Surv(tte, event) ~ arm + strata(grp)') )
	
	colnames(res_s) <- paste(colnames(res_s), 's', sep = '_')
	#with(ds1, tapply(tte, list(arm, grp), median))
	
 stopifnot(all(row.names(res0) == row.names(res_s)))
 stopifnot(all(row.names(res0) == res1b[,'O-E_dup']))
	cbind(res0[, c('n','nevent','p','qad_hr')], res_s[, c('p_s','qad_hr_s')], delta_p = res0[,'p'] - res_s[,'p_s'], res1b)

}


sim3b <- function( p     ) {
  # this function no longer depends on an external custom function. 
  # use simsurv (instead of rpact) to generate survival data
  # try to control the number of event in each simulated  data
 stopifnot(p$dup == 1)
   # Simulate the event times
  j <- 0;
  ds1 <- NULL
  
  if (is.null(p$fractionActiveArm)) p$fractionActiveArm <- 0.5 # handle unequal randomization
  
  readoutTime <- round( -log(1 - p$nEvent / sum(p$n)) / (log(2)/ min(p$med)), 1) # exponential distribution SDF

  while(j <= 10 && (is.null(ds1) || ( (!is.null(ds1)) && sum(ds1[,'event']) < pL$nEvent))) {
 
    if(j > 0) {
       print(paste('iteration', j  , 'with read out time at', readoutTime ))
    }
    
    ds1 <- do.call(rbind, lapply(1:length(p$n), function(i) {
      
    # Create a data frame with the subject IDs and treatment covariate
        cov <- data.frame(id = 1:p$n[i],
                          treatmentGroup = 0)
        cov$treatmentGroup[sample(p$n[i], size = p$n[i] * p$fractionActiveArm)] <- 1   

      dat  <- simsurv(dist= 'exponential', lambdas = log(2)/p$med[i], 
                   betas = c(treatmentGroup = log( p$hr[i])), 
                   x = cov, 
                   maxt = readoutTime) 
      
      # rename variables
      dat$accrualTime <- 0  # all patients were enrolled instantaneously
      dat$observationTime <- readoutTime
      dat$survivalTime <- dat$eventtime
      dat$dropoutTime <-  readoutTime * 2 # just make the dropout far away
        colnames(dat) <- gsub('status','event', colnames(dat))
  
        dat$grp <- paste('M',i, sep='')
    # Merge the simulated event times onto covariate data frame
  dat <- merge(cov, dat)
    dat$treatmentGroup <- dat$treatmentGroup + 1 # 1 for control and 2 for active; which is opposite to rpact
      
      
   
    dat
  }))
    
    j  <- j+1
    readoutTime  <- readoutTime + j * min(p$med)
     
    
    
  }
  
   
  # cut event
  if(sum(ds1[,'event']) > p$nEvent){
       ds1 <- cut2event(ds1, p$nEvent) 
  }else{
    ds1$'timeUnderObservation' <- ds1$eventtime
  }
   
  ds1$arm <- factor(ds1$treatmentGroup, levels=c(1, 2), labels = c('B','A'))
  ds1$grp <- factor(ds1$grp, levels=paste('M', 1:length(p$n), sep=''))
  colnames(ds1) <- gsub('timeUnderObservation','tte', colnames(ds1))

   
 
    ds1$dup <- 1 
  
	stopifnot(nrow(ds1) == with(p, sum(n) *dup))
	c1 <- with(ds1, table(treatmentGroup,  dup, grp, useNA = 'ifany'))
	for (i in 1:length(p$n)) {
   #stopifnot(all(as.vector(c1[,, i]) == p$n[i]  * p$fractionActiveArm ))
	  stopifnot( abs( c1['2',, i]  - p$n[i]  * p$fractionActiveArm ) < 1e-6 )
	  stopifnot( abs( c1['1',, i]  - p$n[i]  * ( 1 - p$fractionActiveArm ) ) < 1e-6)
	}
	 
	

	
	res0 <- t ( sapply(split(ds1, paste(ds1$dup)), ana1, form= 'Surv(tte, event) ~ arm') )
	res1 <- t ( sapply(split(ds1, paste(ds1$dup, ds1$grp)), ana1, form= 'Surv(tte, event) ~ arm') )
	res1 <- as.data.frame(res1)
	nv <- strsplit(row.names(res1), split =  ' ')
	res1$dup <- sapply(nv,'[',1)
	res1$grp <- sapply(nv,'[',2)
	colnames(res1) <- make.names(colnames(res1))
	res1b <-  dcast(res1, dup ~ grp, value.var='O.E')
	colnames(res1b) <- paste('O-E', colnames(res1b),sep ='_')
	
	if(length(unique(ds1$grp)) ==2 ){
	
	res_p <- t ( sapply(split(ds1, paste(ds1$dup, ds1$arm)), function(d, form= 'Surv(tte, event) ~ grp') {
	  s1 <- summary( survdiff( as.formula(form)  , data =   d) )
	  hr <- s1[1, 'qad_hr']
	}))
	
	colnames(res_p) <- paste('prog_hr', colnames(res_p),sep='_')
	}else{
	  res_p <- NA
	}
	
	res_s <- t ( sapply(split(ds1, paste(ds1$dup)), ana1, form='Surv(tte, event) ~ arm + strata(grp)') )
	
	colnames(res_s) <- paste(colnames(res_s), 's', sep = '_')
	#with(ds1, tapply(tte, list(arm, grp), median))
	
 stopifnot(all(row.names(res0) == row.names(res_s)))
 stopifnot(all(row.names(res0) == res1b[,'O-E_dup']))
	cbind(res0[, c('n','nevent','p','qad_hr'), drop=F],
	      res_s[, c('p_s','qad_hr_s'), drop=F], delta_p = res0[,'p'] - res_s[,'p_s'], res1b, res_p)

}




sim4b <- function( p  , cm ,collapseRule ,stratumOrder  ) {
  # modified by sim3b(), it mimics sim4() to allow different rules to collapse stratum
  # cm, a data frame that defines strata, comes from the function constructMedians()
  stopifnot(!is.null(p$fractionActiveArm ))
  
 stopifnot(p$dup == 1)
   # Simulate the event times
  j <- 0;
  ds1 <- NULL
  readoutTime <- round( -log(1 - p$nEvent / sum(cm$n)) / (log(2)/ min(cm$med)), 1) # exponential distribution SDF

  while(j <= 10 && (is.null(ds1) || ( (!is.null(ds1)) && sum(ds1[,'event']) < pL$nEvent))) {
 
    if(j > 0) {
       print(paste('iteration', j  , 'with read out time at', readoutTime ))
    }
    
    ds1 <- do.call(rbind, lapply(1:nrow(cm), function(i) {
      
    # Create a data frame with the subject IDs and treatment covariate
        cov <- data.frame(id = 1:cm[i,'n'],
                          treatmentGroup = 0)
        cov$treatmentGroup[sample(cm[i,'n'], size = cm[i,'n'] * p$fractionActiveArm )] <- 1   

      dat  <- simsurv(dist= 'exponential', lambdas = log(2)/cm[i,'med'], 
                   betas = c(treatmentGroup = log( cm[i,'hr'])), 
                   x = cov, 
                   maxt = readoutTime) 
      
      # rename variables
      dat$accrualTime <- 0  # all patients were enrolled instantaneously
      dat$observationTime <- readoutTime
      dat$survivalTime <- dat$eventtime
      dat$dropoutTime <-  readoutTime * 2 # just make the dropout far away
        colnames(dat) <- gsub('status','event', colnames(dat))
  
        dat$grp <- cm[i,'stratum']
        
    # Merge the simulated event times onto covariate data frame
  dat <- merge(cov, dat)
    dat$treatmentGroup <- dat$treatmentGroup + 1 # 1 for control and 2 for active; which is opposite to rpact
      
      
   
    dat
  }))
    
    j  <- j+1
    readoutTime  <- readoutTime + j * min(cm[,'med'])
     
    
    
  }
  
   
  # cut event
  if(sum(ds1[,'event']) > p$nEvent){
       ds1 <- cut2event(ds1, p$nEvent) 
  }else{
    ds1$'timeUnderObservation' <- ds1$eventtime
  }
   
  ds1$arm <- factor(ds1$treatmentGroup, levels=c(1, 2), labels = c('B','A'))
  ds1$grp <- factor(ds1$grp, levels= cm[,'stratum'])
  colnames(ds1) <- gsub('timeUnderObservation','tte', colnames(ds1))

  ds1$dup <- 1 
  
	stopifnot(nrow(ds1) == sum(cm[,'n']) *p$dup)
	c1 <- with(ds1, table(treatmentGroup,  dup, grp, useNA = 'ifany'))
	for (i in 1:nrow(cm) ) {
   # stopifnot(all(as.vector(c1[,, i]) == cm[i,'n'] /2  ))
	  stopifnot( abs( c1['2',, i]  - cm[i,'n']  * p$fractionActiveArm ) < 1e-6 )
	  stopifnot( abs( c1['1',, i]  - cm[i,'n']  * ( 1 - p$fractionActiveArm ) ) < 1e-6)
 
	
	}
	 
	res0 <- t ( sapply(split(ds1, paste(ds1$dup)), ana1, form= 'Surv(tte, event) ~ arm') )
	 	res2 <-  t ( sapply(split(ds1, paste(ds1$dup)), function(d) {
  	  d2 <- collapse_grp2(d, minN=p$minE, rule= collapseRule, sumFun = sum, strat.order =stratumOrder )
   
  	  c( ana1( d = d2 , form= 'Surv(tte, event) ~ arm + strata(grp)'), n_strata = length(unique(d2$grp)), minStrataSize= min(table(d2$grp)) )
	  }))
	 	
	 	res3 <-  t ( sapply(split(ds1, paste(ds1$dup)), function(d) {
  	  d2 <- collapse_grp2(d, minN=p$minE, rule= collapseRule, sumFun = length, strat.order =stratumOrder )
   
  	  c( ana1( d = d2 , form= 'Surv(tte, event) ~ arm + strata(grp)'), n_strata = length(unique(d2$grp)) , minStrataSize= min(table(d2$grp)) )
	  }))
	 	
	 	colnames(res2) <- paste(colnames(res2),'cE', sep='_')
	  colnames(res3) <- paste(colnames(res3),'cN', sep='_')
	  
	res_s <- t ( sapply(split(ds1, paste(ds1$dup)), ana1, form='Surv(tte, event) ~ arm + strata(grp)') )
	
	colnames(res_s) <- paste(colnames(res_s), 's', sep = '_')
	#with(ds1, tapply(tte, list(arm, grp), median))
	
 stopifnot(all(row.names(res0) == row.names(res_s)))

	cbind(res0[, c('n','nevent','p','qad_hr'), drop=F],
	      res_s[, c('p_s','qad_hr_s'), drop=F], 
	      res2[, c('p_cE','qad_hr_cE','n_strata_cE', 'minStrataSize_cE'), drop=F], 
	      res3[, c('p_cN','qad_hr_cN','n_strata_cN', 'minStrataSize_cN'), drop=F],
	      delta_p = res0[,'p'] - res_s[,'p_s'],
	      delta_pE = res2[,'p_cE'] - res_s[,'p_s'],
	      delta_pN = res3[,'p_cN'] - res_s[,'p_s'])

}



constructMedians <- function( bigN = 1e6, overallMed = 5, freq=c(M1=0.45, M2=0.2), hr = c(0.7, 0.4), check_medians =F) {
  # this function tries to construct subgroup hazard so that the hr 

f1 <- as.matrix( expand.grid(lapply(freq, function(x) c(x, 1-x))) ) # subgroup frequencies
f1.name <- NULL
if(!is.null(names(freq))) {
  f1.name <-   expand.grid(lapply(names(freq), function(x) c(paste(x,'+', sep=''), paste(x,'-', sep=''))), stringsAsFactors = F  )   # subgroup frequencies
  colnames(f1.name) <- names(freq)
   
}  
f2 <- exp( rowSums(log(f1)) )
stopifnot(abs( sum(f2)  - 1) < 1e-5)
c1 <- as.matrix( expand.grid(lapply(freq, function(x) c(x, x-1))) ) # coefficients
r1 <- log(2)/ overallMed * exp( c1 %*% matrix(log(hr), ncol =1) ) # rate of exponential distribution

stopifnot(length(f2) == length(r1))

s1 <- lapply(1:length(f2), function(i)  rexp(round(bigN * f2[i], 0), rate=r1[i]))

#names(s1) <- apply(f1, 1, paste, collapse=',')


r2 <- as.vector( r1 * median(unlist(s1)) / overallMed)



s1 <- lapply(1:length(f2), function(i)  rexp(round(bigN * f2[i], 0), rate=r2[i]))

if(!is.null(f1.name)) {
  
  f1.name$rate <- r2
  f1.name$sim_med <- sapply(s1, median)
  names(r2) <- names(f2) <- apply(f1.name[,names(freq)], 1, paste, collapse='')
  
  
  if(check_medians) {
    
    plot(density(unlist(s1)), col='black', main= paste( 'overall distribution from simulation, with median=', round(median(unlist(s1)),2) ))
lines(density(rexp(10000, rate=log(2)/overallMed)), col='red')
legend('topright', lty=1, legend = c('mixed sim','standard exp'), col=c('black','red'))
    
    for( marker in names(freq)){
      print(marker)
      # individual groups: not done
       for ( j in grep(paste(marker,'+', sep=''), as.character(f1.name[[marker]]), fixed = T )) {
         
         cond1 <- paste(f1.name[j, setdiff(names(freq), marker  ),drop=T], collapse = ',')
         print(paste('conditional on', cond1))
         k <- which( apply(f1.name[ , setdiff(names(freq), marker  ),drop=F], 1, paste, collapse = ',') == cond1)
         temp <- sapply(k, function(k1) {if(k1 != j) print(f1.name[j,'rate'] / f1.name[k1,'rate'])})
       }
      
      # collapse other groups
      print(paste( 'marginal hr',
      
      median(do.call(c, s1[which(f1.name[[marker]] %in% paste(marker,'-',sep=''))])) / 
      median(do.call(c, s1[which(f1.name[[marker]] %in% paste(marker,'+',sep=''))]))  
      ) )
    }
  }
}  

 


invisible(  list(overallM = median(unlist(s1)), rate=r2, freq= f2) )
}

