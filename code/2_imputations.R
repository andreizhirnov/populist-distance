library(mice)
library(ranger)

rm(list=ls())
gc()

set.seed(123)

## define constants
# number of imputations
ni <- 20L
methodi <- 'rf'
aux <- c("party","wt","wf.sex","wf.edu","wf.age")

## read data
sams0 <- readRDS("rdata/full_sample.rds") 
mods <- read.csv("misc/variables.csv") 
reverse <- with(mods, setNames(reverse=='Y', st.var))
 
sams <- lapply(sams0, function(d) {
  vs <- grep("^prec",colnames(d),value=TRUE)
  ## trim columns 
  cmi <- apply(is.na(d[,vs]), 2L, mean)
  vs <- vs[which(cmi < .2)]
  ## trim rows 
  rmi <- apply(is.na(d[,vs]), 1L, mean) 
  d <- d[which(rmi < 0.2), c("respid","cntry", vs, aux)]
  ## designate ordered factors
  for (v in vs) {
    if (reverse[v]) {
      d[[v]] <-ordered(d[[v]],levels=5:1)
    } else {
      d[[v]]<-ordered(d[[v]],levels=1:5) 
    }
  }
  d 
})

# Imputation of missing data
df.imp <- lapply(sams, function(d) {  
  d$respid <- d$cntry <- NULL
  d |> mice(seed=123, m=ni, method=methodi)
})

# Select preferred and save
saveRDS(sams,'rdata/full_sample_grid.RDS')
saveRDS(df.imp,'rdata/full_sample_imputed.RDS')