# para <- c('27519', '60,60',  '20', '0.8,0.8', '27.7,17.2', '50', 'oneFactor')
 #para <- c('27519',  '5', '192' , 'ivg011_sim2','6,20,2','2/3' )

# 3c: try various minE values
# 3b: introduces rules to collapse small subgroups

 para <- commandArgs(trailingOnly = TRUE)

print(para)

home  = '~/projects/stratification'
seed.base = as.integer(para[1])

parseP <- function(x) {
  x <- gsub(' ','',x)
  strsplit(x, split= ',')[[1]]
}

pL <- list(  
           dup = 1,
           nEvent = as.integer(para[3])  
            )

minE   = as.integer(parseP(para[5]))

if(length( minE) == 3)  minE <- seq(from = minE[1], to =minE[2], by = minE[3])
minE <- rep(minE, each = as.integer(para[2]))
 
run_name  = para[4]
if(!is.na(para[6]) && length(para[6])) pL$fractionActiveArm <- eval(parse(text= para[6]))

print(unlist(pL))

source('~/R/lib_2020/cox.table.r')
source('~/R/lib_2020/output_functions.r')
source('~/R/lib_2020/yyboxplot.r')
source('~/R/lib_2020/rpact_functions.r')

source(paste(home,'strat_functions.r', sep = '/'))

library(reshape2)
 library(simsurv)
library(survival)

load(paste( home,'imvigor011_sim_setup1.Rdata', sep ='/'))
 

pL_chk <- pL
for (v in c('med','n','hr'))  pL_chk[[v]] <- cm2[[v]]
check_para(pL_chk)


feature.table$var <- gsub('STRAT','M', feature.table$var)
feature.table <- feature.table[order( pmin(feature.table$freq, 1- feature.table$freq) ), ]


set.seed(seed.base)
r1 <-  do.call(rbind, lapply(minE, function(m)  {
  pL$minE <- m
  r2 <- sim4b(p=pL, cm=cm2,  'collapseGroup' , stratumOrder= feature.table$var )
  cbind(r2, minE=m)
  }))

# set.seed(seed.base)
# r2 <- do.call(rbind,  lapply(1: as.integer(para[2]), function(x) sim4b(p=pL, cm=cm2,  'removeStratum' , stratumOrder= feature.table$var )))


outF <- paste(home,'reports', run_name, sep = '/')
if(!file.exists(outF))   dir.create(outF, recursive =T)

save(r1,  file = paste(outF, '/', seed.base,'.rdata', sep = ''))
