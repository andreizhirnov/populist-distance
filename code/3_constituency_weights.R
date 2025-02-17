library(survey)
library(dplyr)
library(tidyr)


rm(list=ls())
gc()

set.seed(123)

## load data
li0 <- readRDS('rdata/full_sample_grid.RDS')
ess0 <- readstata13::read.dta13("rdata/ESS9e03_2.dta", convert.factors = FALSE )

## constants
vote.vars <- c("prtvede2", "prtvtcit","prtvtcse","prtvtgnl")
names(wei.vars) <- wei.vars <- c("wf.sex","wf.age","wf.edu")
names(countries) <- countries <- names(li0)

### clean up the list
temp.bycp <- lapply(countries, function(cn) {
  x <- li0[[cn]] |> 
    mutate(cntry_party = paste0(cntry, ":", party)) |>
    mutate(cntry_party = case_match(cntry_party,
                             'SE:S' ~ 'SE:SAP',
                             'IT:Sinistra' ~ 'IT:SI',
                             'NL:Denk' ~ 'NL:DENK',
                             .default = cntry_party)) |>
    filter(grepl("^U", respid))
  split(x, x$cntry_party) 
  }) 
names(temp.bycp) <- NULL
temp.bycp <- do.call('c', temp.bycp)


### prepare ESS-based benchmarks
ptran <- read.table(header=TRUE, sep="|", text="
var|val|party|cntry_party
prtvede2|1|CDU/CSU|DE:CDU/CSU
prtvede2|2|SPD|DE:SPD
prtvede2|3|Linke|DE:Linke
prtvede2|4|Greens|DE:Greens
prtvede2|5|FDP|DE:FDP
prtvede2|6|AfD|DE:AfD
prtvede2|7|OTH|
prtvede2|8|OTH|
prtvede2|9|OTH|
prtvtcit|8|FI|IT:FI
prtvtcit|1|PD|IT:PD
prtvtcit|7|M5S|IT:M5S
prtvtcit|9|Lega|IT:Lega
prtvtcit|10|FdI|IT:FdI
prtvtcit|14|OTH|
prtvtcit|2|PlusEuropa|IT:PlusEuropa
prtvtcit|12|OTH|
prtvtcit|6|OTH|
prtvtcit|13|OTH|
prtvtcit|5|OTH|
prtvtcit|3|Verde|IT:Verde
prtvtcit|4|OTH|
prtvtcit|11|OTH|
prtvtgnl|6|D66|NL:D66
prtvtgnl|1|VVD|NL:VVD
prtvtgnl|4|SP|NL:SP
prtvtgnl|5|CDA|NL:CDA
prtvtgnl|10|OTH|
prtvtgnl|2|PvdA|NL:PvdA
prtvtgnl|8|GL|NL:GL
prtvtgnl|16|OTH|
prtvtgnl|14|OTH|
prtvtgnl|3|PVV|NL:PVV
prtvtgnl|7|CU/SGP|NL:CU/SGP
prtvtgnl|9|CU/SGP|NL:CU/SGP
prtvtgnl|11|OTH|
prtvtgnl|13|FvD|NL:FvD
prtvtgnl|12|OTH|
prtvtgnl|17|OTH|
prtvtcse|6|SAP|SE:SAP
prtvtcse|5|M|SE:M
prtvtcse|7|V|SE:V
prtvtcse|4|MP|SE:MP
prtvtcse|1|C|SE:C
prtvtcse|3|KD|SE:KD
prtvtcse|9|SD|SE:SD
prtvtcse|2|L|SE:L
prtvtcse|10|OTH|
prtvtcse|8|OTH|
")

ess.p.l <- ess0 |> 
  filter(cntry %in% countries) |>
  select(idno, cntry, any_of(vote.vars)) |>
  pivot_longer(cols=any_of(vote.vars), names_to='var', values_to='val', values_drop_na = TRUE) |> 
  inner_join(ptran, by=join_by(var, val)) |>
  filter(nchar(cntry_party)>2) |> 
  select(cntry, idno, cntry_party)

if(any(duplicated(ess.p.l$idno))) cat("There are duplicates in the vote records!")

ess <- ess0 |> 
  subset(cntry %in% countries, select=c("idno","cntry", "gndr", "agea", "eisced" )) |>
  mutate(
    wf.sex = case_match(gndr, 1 ~ 'M', 2 ~ 'F'),
    wf.age = case_when(
      agea>=65 ~ "Y_GE65",
      agea>=50 ~ "Y50-64",
      agea>=30 ~ "Y30-49",    
      agea>=15 ~ "Y15-29",
      TRUE ~ NA
    ),
    wf.edu = case_match(eisced, 0:5 ~ "ED1:4", 6:7 ~"ED5:6")
    ) |>
  inner_join(ess.p.l, by = join_by(idno, cntry)) |>
  group_by(cntry_party) |> mutate(nv = n()) |> filter(nv>=10) |> ungroup()

names(grid.cp) <- grid.cp <- sort(unique(ess$cntry_party)) 
wt.template0 <- lapply(grid.cp, function(x) {
  lapply(wei.vars, function(w) {  
      ess |> 
      subset(cntry_party==x & !is.na(ess[[w]]))|>
      group_by(pick({{w}})) |>  
      summarise(value=n()) |> mutate(value=value/sum(value, na.rm=TRUE)) |> as.data.frame()
  })
})
 
## add weights
temp.wtd <- lapply(names(temp.bycp), function(x) {  
  wd <- temp.bycp[[x]]
  wd$wt <- 1
  if (!x %in% grid.cp) return(wd)
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
  sur <- svydesign(ids=~respid,data= wd, weights=~wt)
  sur.r <- rake(sur, sample.margins = lapply(paste0("~",wei.vars), function(j) as.formula(j)), population.margins = wt.template, control = list(maxit = 200))
  wt <- weights(sur.r)
  wt <- wt/mean(wt)
  wt[wt>3] <- 3
  qq <- setNames(wt/mean(wt), sur.r$cluster$respid)
  return(within(wd, wt <- qq[respid]))
})

temp.wtd <- lapply(temp.wtd, subset, select=c('respid','cntry','cntry_party', 'wt')) 
temp.wtd$make.row.names <- FALSE
temp.wtd <- do.call("rbind", temp.wtd)

saveRDS(temp.wtd, "data/weights_by_cp.rds")