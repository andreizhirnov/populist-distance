library(survey)
library(dplyr)
library(tidyr)
library(openxlsx)

set.seed(123)
options(stringsAsFactors = FALSE)
rm(list=ls())
gc() 

### a few constants
names(countries) <- countries <- c("DE","IT","NL","SE")
wei.vars <- c("wf.sex","wf.edu","wf.age")

### read data
mods <- read.csv("misc/variables.csv")
dims <- unique(mods$dim)
labs <- openxlsx::read.xlsx("rdata/labels_vaa2019.xlsx", sheet=1L)
labs <- setNames(labs[[1]],labs[[2]])
# VAA user data with responses to statements
sample0 <- readRDS("rdata/user_q_file.rds")
# VAA user data with demographic characteristic
dem0 <- readRDS("rdata/user_file.rds")
# Political preferences (vote intent)
pdata0 <- readRDS("rdata/user_ptv_file.rds")
# Party data
parties0 <- read.csv("rdata/party_positions_2019_2020.csv")
# Census benchmarks
fl <- load(file="misc/Census benchmarks.RData")
 
### fix country codes and individual IDs 
countrycodes <- c(fra="FR",ger="DE",ita="IT",nl="NL",swe="SE")
samples <- sample0 |> mutate(cntry=countrycodes[country],
                             respid = paste0("U", respid)) |>
  filter(cntry %in% countries) |> 
  select(cntry, respid, starts_with("precede")) |>
  pivot_longer(cols=starts_with("precede"), names_to="st.var", values_to="val", values_drop_na = TRUE)
samples <- split(samples, samples$cntry)

dem <- dem0 |> mutate(cntry=countrycodes[country],
                       respid = paste0("U", respid)) |>
  filter(cntry %in% countries)

### clean up
rm(sample0, dem0)

### combine parties data -- add Swedish parties separately
parties <- parties0 |> filter(cntry %in% countries & cntry !="SE") |>
  select(cntry, party, starts_with("precede")) |>
  pivot_longer(cols=starts_with("precede"), names_to="st.var", values_to="val", values_drop_na = TRUE)

parties.SE <- openxlsx::read.xlsx("rdata/party_coding_ep2019_Sweden.xlsx",sheet=1L)
parties.SE$cntry<-"SE"
parties.SE <- within(parties.SE, {
  party <- trimws(sapply(strsplit(Party,"\\."),"[",2L))
  st.code <- trimws(sapply(strsplit(Statement,"\\."),"[",1))
  st.var <- c('07'='precede0106',
              '03'='precede0503196',
              '32'='precede0503251',
              '04'='precede0246',
              '08'='precede0503223',
              '30'='precede0503225',
              '17'='precede0238',
              '19'='precede0503234',
              '22'='precede0503257',
              '35'='precede0503224',
              '06'='precede0503197',
              '20'='precede0503231',
              '21'='precede0503233',
              '41'='precede0503199',
              '11'='precede0503192',
              '05'='precede0216',
              '33'='precede0102',
              '10'='precede0503203',
              '16'='precede0103',
              '15'='precede0105',
              '14'='precede0503194')[st.code] 
  val <- match(Coding, c("Do not agree at all","Tend not to agree","Neither or","Tend to agree", "Agree completely"))
})

if(any(is.na(parties.SE$val))) cat("Missing values in the SE data")
parties.SE.prep <- parties.SE |> 
  filter(!is.na(st.var)) |>
  select(cntry, party, st.var, val)

parties <- rbind(parties, parties.SE.prep) |> semi_join(mods, by=join_by(cntry, st.var))

### remove observations with all missing values 
bulk <- lapply(countries, function(cn) {
  t <- samples[[cn]] |> inner_join(mods, by=join_by(cntry, st.var))
# make sure there is at least 1 item on each dimension 
  nmi.cat <- t |> group_by(respid) |> 
    summarize(n=n_distinct(dim, na.rm=TRUE)) |> 
    filter(n==length(dims))
  semi_join(t, nmi.cat, by=join_by(respid))
})

# add political preferences
ivc <- pdata0 |> filter(cntry %in% countries & type=="int_EP") |>
  select(cntry, respid, party) 
 
if (any(duplicated(ivc$respid))) cat("There are duplicates in the vote choice data!")
 
# extract demographic information for weights
respids <- unique(do.call('c', lapply(bulk, `[[`, 'respid')))
respids <- intersect(respids, ivc$respid)

dem <- dem |> filter(respid %in% respids) |> 
  mutate(
    wf.sex = case_when(
      female==1 ~ 'F',
      female==0 ~ 'M',
      male==1 ~ 'M',
      TRUE ~ 'UNK' 
    ),
    age = coalesce(age, 2019-birthyear),
    wf.age = case_when(
      age>=65 ~ "Y_GE65",
      age>=50 ~ "Y50-64",
      age>=30 ~ "Y30-49",
      age>=15 ~ "Y15-29",
      TRUE ~ 'UNK'
    ),
    wf.edu = case_when(
      university==1 ~ "ED5:6",
      university==0 ~ "ED1:4" ,
      TRUE ~ 'UNK' 
    )
  )

allmi2 <- apply(dem[,wei.vars], 1L, function(x) all(x=="UNK"))
wei.data <- dem[!allmi2, c("respid","cntry",wei.vars)]
# this information will also be needed for computing constituency weights
saveRDS(wei.data,'rdata/wei_data.RDS')

## Compute sampling weights
wt.template0 <- emp.select.bmks
weights <-  lapply(countries, function(x) { 
  wd <- subset(wei.data, cntry==x) 
  wt.template <- lapply(wei.vars, function(j) {
    ini <- wt.template0[[x]][[j]]
    ini <- ini[ini[[j]]!="UNK",]
    unk <- mean(wd[[j]]=="UNK")
    unk.df <- data.frame(k="UNK",value=unk)
    colnames(unk.df)[1] <- j
    ini$value <- (1-unk)*ini$value/sum(ini$value)
    if (unk>0) ini <- rbind(ini, unk.df)
    ini
  })
  sur <- svydesign(ids=~respid, probs=1, data= wd)
  sur.r <- rake(sur, sample.margins = lapply(paste0("~",wei.vars), function(j) as.formula(j)), population.margins = wt.template, control = list(maxit = 200))
  wt <- weights(sur.r)
  wt <- wt/mean(wt)
  wt[wt>5] <- 5
  data.frame(respid=sur.r$cluster$respid, wt=wt/mean(wt))
})
weights <- do.call("rbind", weights)
rownames(weights) <- NULL

### combine information
combined <- lapply(countries, function(cn) {
  u <- bulk[[cn]] |> select(cntry, respid, st.var, val) |> 
    inner_join(ivc, by=join_by(cntry, respid)) |> inner_join(weights, by=join_by(respid))
  p <- parties |> filter(cntry==cn) |>
    mutate(respid = party, wt = 1)
  x <- rbind(u, p) |> pivot_wider(names_from = st.var, values_from=val) |>
    left_join(wei.data, by=join_by(cntry, respid)) |>
     mutate(party = trimws(party))
  if (cn=='SE') {
    x |> mutate(party = case_match(party,
                                   'Kristdemokraterna' ~ 'KD', 
                                   'Miljopartiet' ~ 'MP',
                                   'Moderaterna' ~ 'M',
                                   c('Socialdemokraterna','S') ~ 'SAP',
                                   'Sverigedemokraterna' ~ 'SD',
                                   'V?nsterpartiet' ~ 'V',
                                   .default=party))
  } else if (cn=='IT') {
    x |> mutate(party = case_match(party,
                                  'Sinistra' ~ 'SI',
                                  .default=party))
  } else if (cn=='NL') {
    x |> mutate(party = case_match(party,
                                   'Denk' ~ 'DENK',
                                   .default=party))
  } else if (cn=='DE') {
    x |> mutate(party = case_match(party,
                                   'Animal Rights Party' ~ 'TSP',
                                   'Pirates' ~ 'Piraten',
                                   .default=party))
  }
  })

### same the data
saveRDS(combined,"rdata/full_sample.rds")

lapply(combined, nrow)