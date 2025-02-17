library(dplyr)
library(tidyr)
library(lavaan) 
library(mice)

rm(list=ls())
gc()

set.seed(123)

## extract the names of observed variables from a lavaan model
extract_ovs <- function(model) {
  t <- lavaan::lavaanify(model)
  unique(subset(t, op=='=~', select='rhs', drop=TRUE))
} 

## read data 
li <- readRDS('rdata/full_sample_grid.RDS')
df.imp <- readRDS('rdata/full_sample_imputed.RDS')
mods <- read.csv("misc/variables.csv") 

## set constants
names(countries) <- countries <- names(li)

## CFA models
models <- mods |> filter(use_agg=='Y') |> group_by(cntry, dim) |>
  summarize(text = paste(st.var, collapse=" + ")) |>
  mutate(text = paste0(dim, ' =~ ', text, ';')) |>
  summarize(text = paste(text, collapse=" "))
models <- with(models, setNames(text, cntry))

## extract complete datasets
completes <- lapply(countries, function(cn) {
  lapply(seq_len(df.imp[[cn]]$m), function(i) {
    u <- mice::complete(df.imp[[cn]], action=i)
    u$wt <- li[[cn]]$wt
    u
  })
})

ests <- lapply(countries, function(cn) {
    lapply(completes[[cn]], function(u) {
      cfa(models[[cn]], ordered=TRUE, sampling.weights="wt", data=u)
    })
})
saveRDS(ests,'data/CFA_estimates.RDS')

## export estimates
s_alias <- mods |> with(setNames(alias, st.var))

qq <- lapply(countries, function(cn) {
  g <- lapply(ests[[cn]], function(e) { 
    u <- subset(parTable(e), op %in% c('~~','|','=~'), select=c("lhs","op","rhs","est","se"))
    s <- subset(standardizedSolution(e), op %in% c('~~','|','=~'), select=c("lhs","op","rhs","est.std","se"))
    colnames(s)[4L] <- 'est'
    u$type <- 'raw'
    s$type <- 'std'
    v <- rbind(u,s)
    within(v, {
      var <- se^2
      op <- ifelse(op=='~~', 'var/cov', ifelse(op=='|','threshold','loading')) 
      })
  }) 
  do.call('rbind', g) |> 
    group_by(lhs, op, rhs, type) |>
    summarize(est=mean(est, na.rm=TRUE),
              var=mean(var, na.rm=TRUE)) |>
    ungroup() |>
    mutate(
      se = sqrt(var), 
      value = paste0(formatC(est, 3L, format='f', flag='+'), ' (',
                      formatC(se, 3L, format='f'), ')'),
      variables = ifelse(op=='var/cov',
                          ifelse(lhs==rhs,
                                 ifelse(lhs %in% names(s_alias), s_alias[lhs], lhs),
                                 paste0(lhs,'-',rhs)),
                          ifelse(op=='threshold',
                                 paste0(s_alias[lhs],": ", rhs),
                                 paste0(s_alias[rhs]," on ", lhs)
                          )),
      cntry = cn
      ) |> 
    select(cntry, op, type, value, variables) |>  
    pivot_wider(names_from=type, values_from=value) 
})
est_out <- do.call('rbind', qq) |> arrange(cntry, op, variables)
write.csv(est_out, "output/CFA_estimates.csv", row.names=FALSE, na="")

## inspect the first set of estimates
first <- lapply(ests,'[[',1L)
lapply(first, fitMeasures, fit.measures=c("cfi", "rmsea",'srmr'))
lapply(first, lavInspect, what="cor.lv")

## generate predicted values for latent variables
preds <- lapply(countries, function(cn) {
  lapply(ests[[cn]], lavPredict, type ="lv", method='ML')
})

issues <- mods |> filter(nchar(issue_desc)>0) |>
  select(cntry, st.var, issue) |> rename('cat'='issue')

positions <- lapply(countries, function(cn) {
  # issue positions
  s1 <- lapply(completes[[cn]], function(d) { 
      d |> 
        select(starts_with('precede')) |>
        apply(MARGIN=2L, FUN=as.numeric) |>
        bind_cols(li[[cn]][,c('respid','cntry','party')]) |>
        pivot_longer(cols=starts_with('precede'), 
                   names_to='st.var',
                   values_to='pos', 
                   values_drop_na = TRUE) |>
        inner_join(issues, by=join_by(cntry, st.var)) |>
        mutate( 
          pos = (pos-1)/4,
          cntry_party = paste0(cntry, ':', party)) |>
        group_by(cntry, cntry_party, respid, cat) |>
        summarize(pos=mean(pos, na.rm=TRUE)) |>
        ungroup()
  })
  # aggregate positions
  s2 <- lapply(preds[[cn]], function(d) { 
    vs <- colnames(d)
    d |> bind_cols(li[[cn]][,c('respid','cntry','party')]) |>
       pivot_longer(cols=any_of(vs), 
                    names_to='cat', 
                    values_to="pos", 
                    values_drop_na = TRUE) |>
      mutate(cntry_party = paste0(cntry, ':', party)) |>
      select(cntry, cntry_party, respid, cat, pos)
  })
  # average over imputations and standardize
  do.call('rbind',c(s1,s2)) |>
    group_by(cntry, cntry_party, cat, respid) |>
    summarize(pos=mean(pos, na.rm=TRUE)) |>
    group_by(cat) |>
    mutate(pos = (pos-mean(pos, na.rm=TRUE))/sd(pos, na.rm=TRUE)) |>
    ungroup()
})

positions_c <- do.call('rbind', positions)

positions_p <- subset(positions_c, subset=!is.na(pos) & !grepl('^U', respid),
                      select=c("cntry_party","cat","pos")) 
positions_u <- subset(positions_c, subset=!is.na(pos) & grepl('^U', respid))
 
saveRDS(positions_p,"data/party_positions.rds")
saveRDS(positions_u,"data/user_positions.rds")
