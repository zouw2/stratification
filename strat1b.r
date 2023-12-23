# para <- c('27519', '60,60',  '20', '0.8,0.8', '27.7,17.2', '50', 'oneFactor')
 #para <- c('27528', '60,60', '20', '0.8,0.8' ,'27.7,17.2', '90', '1F_scen3b')
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
