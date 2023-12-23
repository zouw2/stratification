# para <- c('27519', '60,60',  '20', '0.8,0.8', '27.7,17.2', '50', 'oneFactor')
 #para <- c('27528', '60,60', '20', '0.8,0.8' ,'27.7,17.2', '90', '1F_scen3b')
#para <- c('27528', '400', '20', '0.8' ,'10', '320', '1F_scen3b')
# para <-c('27519', '12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,14,14,14,14,14,14,14,14,14,14,10,10', '40', '0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737', '10,8.79,7.73,6.79,5.97,5.25,4.61,4.05,3.56,3.13,2.75,2.42,2.13,1.87,1.64,1.44,1.27,1.12,0.98,0.86,0.76,0.67,0.59,0.51,0.45,0.4,0.35,0.31,0.27,0.24,0.21,0.18', '320', 'rep_aka_3')
 para <- commandArgs(trailingOnly = TRUE)

# strat1b_aka_collapse_neighbor: add an additional function to collapse 2 neighbor strata (in terms of prognostsis ) in analytical stage
 
# strat1b: this program changes the simulation engine to simsurv
 
print(para)

home  = '~/projects/stratification'
seed.base = as.integer(para[1])

parseP <- function(x) {
  x <- gsub(' ','',x)
  strsplit(x, split= ',')[[1]]
}

pL <- list(n = as.integer(parseP(para[2])), 
           dup = 1,
           hr = as.numeric(parseP(para[4])),
           med = as.numeric(parseP(para[5])),
           nEvent = as.integer(para[6])
           )
           

run_name  = para[7]

source('~/R/lib_2020/cox.table.r')
source('~/R/lib_2020/output_functions.r')
source('~/R/lib_2020/yyboxplot.r')
 

source(paste(home,'strat_functions.r', sep = '/'))

library(reshape2)
 library(simsurv)
library(survival)

collapse_neighor_strata <- function (  grp ) {
  id <- gsub('^M','', grp)
  stopifnot(as.character(as.integer(id)) == id)
  id <- as.integer(id)
  stopifnot(!any(is.na(id)))
  id2 <- round( id / 2 + 1e-3) *2
  factor( paste('M', id2, sep=''), levels=paste('M', sort(unique(id2)), sep =''))
}

# modified from function sim3b(), cn stands for collapse neighbor
sim3b_cn <- function( p     ) { 
  # this function no longer depends on an external custom function. 
  # use simsurv (instead of rpact) to generate survival data
  # try to control the number of event in each simulated  data
 stopifnot(p$dup == 1)
   # Simulate the event times
  j <- 0;
  ds1 <- NULL
  readoutTime <- round( -log(1 - p$nEvent / sum(p$n)) / (log(2)/ min(p$med)), 1) # exponential distribution SDF

  while(j <= 10 && (is.null(ds1) || ( (!is.null(ds1)) && sum(ds1[,'event']) < pL$nEvent))) {
 
    if(j > 0) {
       print(paste('iteration', j  , 'with read out time at', readoutTime ))
    }
    
    ds1 <- do.call(rbind, lapply(1:length(p$n), function(i) {
      
    # Create a data frame with the subject IDs and treatment covariate
        cov <- data.frame(id = 1:p$n[i],
                          treatmentGroup = 0)
        cov$treatmentGroup[sample(p$n[i], size = p$n[i]/2)] <- 1   

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
   stopifnot(all(as.vector(c1[,, i]) == p$n[i] /2  ))
	}
 
	ds1$grp  <- collapse_neighor_strata( as.character(ds1$grp) )

	
	res0 <- t ( sapply(split(ds1, paste(ds1$dup)), ana1, form= 'Surv(tte, event) ~ arm') )
	res1 <- t ( sapply(split(ds1, paste(ds1$dup, ds1$grp)), ana1, form= 'Surv(tte, event) ~ arm') )
	res1 <- as.data.frame(res1)
	nv <- strsplit(row.names(res1), split =  ' ')
	res1$dup <- sapply(nv,'[',1)
	res1$grp <- sapply(nv,'[',2)
	colnames(res1) <- make.names(colnames(res1))
	res1b <-  dcast(res1, dup ~ grp, value.var='O.E')
	colnames(res1b) <- paste('O-E', colnames(res1b),sep ='_')
	
	if(length(unique(ds1$grp)) > 1){
	
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

check_para(pL)

set.seed(seed.base)
r1 <- do.call(rbind, lapply(1: as.integer(para[3]), function(x) { 
 # print(x)
#  if(x == 11) debug(sim3b)
   sim3b_cn(p=pL  ) }))


 

outF <- paste(home,'reports', run_name, sep = '/')
if(!file.exists(outF))   dir.create(outF, recursive =T)

save(r1, file = paste(outF, '/', seed.base,'.rdata', sep = ''))
