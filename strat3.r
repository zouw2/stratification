# para <- c('27519', '60,60',  '20', '0.8,0.8', '27.7,17.2', '50', 'oneFactor')
 #para <- c('27519', '56,92,10,18,166,276,32,52','20','0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8' ,'23.2,19,16.2,13.3,16.2,13.3,11.4,9.3' ,'53.5' ,'imp150_sim2','10', 'collapseGroup')

# this program introduce rules to collapse small subgroups

 para <- commandArgs(trailingOnly = TRUE)

print(para)

home  = '~/projects/stratification'
seed.base = as.integer(para[1])

parseP <- function(x) {
  x <- gsub(' ','',x)
  strsplit(x, split= ',')[[1]]
}

pL <- list(n = as.integer(parseP(para[2])), 
           dup = as.integer(para[3]),
           hr = as.numeric(parseP(para[4])),
           med = as.numeric(parseP(para[5])), 
           minE = as.integer(para[8]), 
           collapseRule = para[9])
           
readoutTime <- as.numeric(para[6])
run_name  = para[7]

source('~/R/lib_2020/cox.table.r')
source('~/R/lib_2020/output_functions.r')
source('~/R/lib_2020/yyboxplot.r')
source('~/R/lib_2020/rpact_functions.r')

source(paste(home,'strat_functions.r', sep = '/'))

library(reshape2)
 library(rpact)
library(survival)


check_para(pL)

set.seed(seed.base)
r1 <- sim4(p=pL, time = readoutTime, estimatePrognostic = F, maxPms = 20, rampUpTime = 30  ) 
outF <- paste(home,'reports', run_name, sep = '/')
if(!file.exists(outF))   dir.create(outF, recursive =T)

save(r1, file = paste(outF, '/', seed.base,'.rdata', sep = ''))
