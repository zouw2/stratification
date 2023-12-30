# para <- c('27519', '60,60',  '20', '0.8,0.8', '27.7,17.2', '50', 'oneFactor')
 #para <- c('27528', '60,60', '20', '0.8,0.8' ,'27.7,17.2', '90', '1F_scen3b')
#para <- c('27528', '400', '20', '0.8' ,'10', '320', '1F_scen3b')
# para <-c('27519', '12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,14,14,14,14,14,14,14,14,14,14,10,10', '40', '0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737,0.737', '10,8.79,7.73,6.79,5.97,5.25,4.61,4.05,3.56,3.13,2.75,2.42,2.13,1.87,1.64,1.44,1.27,1.12,0.98,0.86,0.76,0.67,0.59,0.51,0.45,0.4,0.35,0.31,0.27,0.24,0.21,0.18', '320', 'rep_aka_3')

# para <- c('27519', '57,15,30,3,24,6,12,3,36,18,18,3,3,3,9,6', '20', '0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65,0.65', '2.83,4.02,3.46,4.91,6.31,8.95,7.7,10.93,3.82,5.43,4.67,6.63,8.51,12.08,10.4,14.75', '192', 'ivg011_sim1','1/3')
 para <- commandArgs(trailingOnly = TRUE)

# this program changes the simulation engine to simsurv
 
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

if(!is.na(para[8]) && length(para[8])) pL$fractionActiveArm <- eval(parse(text= para[8]))

source('~/R/lib_2020/cox.table.r')
source('~/R/lib_2020/output_functions.r')
source('~/R/lib_2020/yyboxplot.r')
 

source(paste(home,'strat_functions.r', sep = '/'))

library(reshape2)
 library(simsurv)
library(survival)


check_para(pL)

set.seed(seed.base)
r1 <- do.call(rbind, lapply(1: as.integer(para[3]), function(x) { 
 # print(x)
#  if(x == 11) debug(sim3b)
   sim3b(p=pL  ) }))



outF <- paste(home,'reports', run_name, sep = '/')
if(!file.exists(outF))   dir.create(outF, recursive =T)

save(r1, file = paste(outF, '/', seed.base,'.rdata', sep = ''))
