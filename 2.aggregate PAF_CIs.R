################################################################
# Harriet Rumgay - International Agency for Research on Cancer
# 24 July 2024
################################################################
#
# Script to aggregate country PAF for smokeless tobacco and areca 
# nut consumption to world regions


#### set up script ####

# load following libraries
library(tidyverse)
library(data.table)
library(xlsx)
library(Rcan)


# set working directory
setwd("~/STAN_PAF")


##### load data #####

# load Globocan 2022 estimates
DATA_GLOB	<- as.data.table(read.csv("_data/Globocan2022.csv")) %>% as.data.table()

# create total all ages combined
DATA_GLOB %>% 
  group_by(country_code,cancer_code,sex,type) %>% 
  mutate(total=sum(cases)) -> DATA_GLOB

# create specific subregion groupings
DATA_GLOB %>% 
  group_by(type,sex, age, cancer_code) %>% 
  filter(type==0) %>% # filter to cancer incidence
  mutate(cases=case_when(country_label=="Micronesia/Polynesia"~sum(cases[country_label%in%c("Micronesia/Polynesia","Melanesia")]),
                         country_label=="Caribbean"~sum(cases[country_label%in%c("Caribbean","Central America")]),
                         TRUE~cases),
         py=case_when(country_label=="Micronesia/Polynesia"~sum(py[country_label%in%c("Micronesia/Polynesia","Melanesia")]),
                      country_label=="Caribbean"~sum(py[country_label%in%c("Caribbean","Central America")]),
                            TRUE~py),
         total=case_when(country_label=="Micronesia/Polynesia"~sum(total[country_label%in%c("Micronesia/Polynesia","Melanesia")]),
                         country_label=="Caribbean"~sum(total[country_label%in%c("Caribbean","Central America")]),
                          TRUE~total),
         country_label=case_when(country_label=="Micronesia/Polynesia"~"Melanesia, Micronesia and Polynesia",
                                 country_label=="Caribbean"~"Caribbean and Central America",
                                 TRUE~country_label)) -> DATA_GLOB

# create set of oral cancer incidence data per world region
DATA_GLOB %>% 
  filter(country_code>899, cancer_code==1) %>% 
  dplyr::select(-country_label) -> datreg

# create set of population numbers per country, age, sex, and aggregated region
DATA_GLOB %>% 
  filter(cancer_code %in% c(1)) %>%
  ungroup() %>% 
  dplyr::select(country_code,country_label,sex,age,py) -> pops

# create regions code
DATA_GLOB %>% 
  filter(cancer_code %in% c(1),country_code>=900) %>% 
  ungroup() %>% 
  dplyr::select(country_code,country_label) %>% unique() -> regscode

# import region groupings per country
regs <-read.csv("_data/globocan_country_region_groupings.csv") %>% 
  as.data.table()

rm(DATA_GLOB)

# load Melanesia Micronesia Polynesia population proportions
micro <-read.csv("_data/MMP pops.csv") %>% 
  as.data.table()

# import HDI groupings
hdi <-read.csv("_data/HDI/HDI_2021.csv") %>% 
  as.data.table() %>% 
  mutate(hdigrp = case_when(HDI2021 >= 0.8 ~ "Very high HDI",
                            HDI2021 < 0.8 & HDI2021 >= 0.7 ~ "High HDI",
                            HDI2021 < 0.7 & HDI2021 >= 0.55 ~ "Medium HDI",
                            HDI2021 < 0.55 & HDI2021 >= 0~ "Low HDI",
                            TRUE ~ "Missing HDI")) %>% 
  dplyr::select(country_code, hdigrp)


#### merge region groupings with paf ####

# run function to merge country AFs by age per simulation with regions and HDI groups
CIs <- function(sim){
  simname <- paste0( "_results/",(sim), "_list.RData") # should load as t.sim2
  load(simname)
  
  t.sim <- lapply(1:10000, function(zzz){
    print(zzz)
    
    t <- t.sim2[[zzz]]
    t<- as.data.table(t)
    
  t %>%
    mutate(cases=ifelse(is.na(cases),0,cases)) %>% 
    ungroup() %>% 
    dplyr::select(country_code,cancer_code,sex,age,cases,total,af,sim.n,simaf,attr,simattr,simtotattr,simtotpaf) %>% 
    left_join(regs) %>% 
    left_join(hdi) -> paf
   
  paf
  
  })
  
  # split simulations into 5 groups to break up vector size
  t<-t.sim
  t.sim <- t[c(1:2000)]
  save(t.sim, file = "_results/aggr_simpaf_list.1.RData")
  t.sim <- t[c(2001:4000)]
  save(t.sim, file = "_results/aggr_simpaf_list.2.RData")
  t.sim <- t[c(4001:6000)]
  save(t.sim, file = "_results/aggr_simpaf_list.3.RData")
  t.sim <- t[c(6001:8000)]
  save(t.sim, file = "_results/aggr_simpaf_list.4.RData")
  t.sim <- t[c(8001:10000)]
  save(t.sim, file = "_results/aggr_simpaf_list.5.RData")
}
CIs(sim = "pafCIs2")

  

#### create aggregated totals ####

# run function to create aggregates by region group
CIsreg <- function(sim){
  
  for (i in 1:5){
    
  simname <- paste0( "_results/",(sim), "_list.",i,".RData") # should load as t.sim
  load(simname)

  t.sim3 <- lapply(1:2000, function(zzz){
    print(paste0(zzz, ".", i))
    
    t <- t.sim[[zzz]]
    paf<- as.data.table(t)
    
  # make Melanesia, Micronesia and Polynesia (MMP) total from MMP country totals >900 country_code - apply to remaining cases in region
  paf %>% 
    filter(sex!=0) %>% 
    dplyr::select(-prop_code,-prop_label) %>% 
    filter(subregion_label_3%in%c("Melanesia, Micronesia and Polynesia") | country_code%in%c(258,316,16,896,583,580,584)) %>% 
    left_join(micro %>% dplyr::select(country_code,prop_code,prop_label,prop_py2)) %>%
    mutate(across(af:simaf, ~replace_na(.x, 0))  #if NA or INF replace with 0
    ) %>%
    group_by(prop_code, cancer_code, age, sex, sim.n) %>% 
    mutate(cases.c=sum(cases,na.rm=T),
           total.c=sum(total,na.rm=T),
           af=mean(af,na.rm=T),
           simaf=mean(simaf,na.rm=T)) %>% 
    filter(country_code%in%c(258,540,316)) %>% 
    mutate(across(web_country_label:unterm_country_label,~prop_label),
           country_code=prop_code,
           across(c(cases, total,attr:simtotpaf),~0)) %>% 
    dplyr::select(-prop_py2) %>% 
    dplyr::select(-cases,-total) %>% 
    left_join(datreg %>% ungroup() %>% 
                dplyr::select(-type,-py,-cancer_label)) %>%
    mutate(cases=(cases-cases.c),
           total=(total-total.c),
           attr=af*cases,
           simattr=simaf*cases) -> pafmmp
  
  # add MMP total to paf of other countries
  paf %>% 
    filter(sex!=0) %>% 
    rbind(pafmmp %>% dplyr::select(-cases.c,-total.c)) -> pafs_subregion
  
  # create subregion totals
  pafs_subregion %>% 
    group_by(subregion_label_3, cancer_code, sex, age, sim.n) %>% # aggregate country-level results to region
    mutate(casesgrp = sum(cases, na.rm=T),
           attr = sum(attr, na.rm=T),
           simattr= sum(simattr, na.rm=T),
           af = attr/casesgrp,
           simaf = simattr/casesgrp,
           country_label=subregion_label_3,
           country_code=subreg_code_GCO) %>% 
    ungroup() %>% 
    dplyr::select(continent_label,country_label,country_code,cancer_code,sex, age, sim.n, casesgrp, attr, af, simattr, simaf) %>% 
    unique() -> pafs_subregion
  
  # create totals for both sexes combined
  pafs_subregion %>%
    bind_rows(pafs_subregion %>%
                ungroup() %>%
                group_by(country_label,country_code,cancer_code,age,sim.n) %>%
                mutate(sex=0,
                       sex_label="Both sexes",
                       casesgrp = sum(casesgrp, na.rm=T),
                       attr = sum(attr, na.rm=T),
                       simattr= sum(simattr, na.rm=T),
                       af = attr/casesgrp,
                       simaf = simattr/casesgrp) %>% unique()) -> pafs_subregion

  # apply reg PAF to globocan region estimates to sum to world - by age
  pafs_subregion %>% 
    left_join(datreg) %>% 
    mutate(attr=cases*af,
           simattr=cases*simaf) %>% 
    dplyr::select(-casesgrp) -> pafs_subregion
  
  pafs_subregion
  
  })
  
  savename <- paste0("_results/aggr_simpaf_subreg_list.",i,".RData")
  save(t.sim3, file = savename)
  
  }
}
CIsreg(sim = "aggr_simpaf")

# run function to create aggregates by WHO region group
CIswho <- function(sim){
  
  for (i in 1:5){
    simname <- paste0( "_results/",(sim), "_list.",i,".RData") # should load as t.sim
  
  load(simname)
  
  t.sim3 <- lapply(1:2000, function(zzz){
    print(paste0(zzz,".",i))
    
    t <- t.sim[[zzz]]
    paf<- as.data.table(t)
  
  # create WHO region totals
  paf %>% 
    group_by(who_region_label, cancer_code, sex, age, sim.n) %>% # aggregate country-level results to region
    mutate(casesgrp = sum(cases, na.rm=T),
           attr = sum(attr, na.rm=T),
           simattr= sum(simattr, na.rm=T),
           af = attr/casesgrp,
           simaf = simattr/casesgrp,
           country_label=who_region_label,
           country_code=case_when(who_region_label=="AFRO"~991,
                                  who_region_label=="AMRO"~992,
                                  who_region_label=="EMRO"~993,
                                  who_region_label=="EURO"~994,
                                  who_region_label=="SEARO"~995,
                                  who_region_label=="WPRO"~996)) %>% 
    ungroup() %>% 
    dplyr::select(who_region_label,country_code,cancer_code,sex, age, sim.n, casesgrp, attr, af, simattr, simaf) %>% 
    unique() -> pafs_who
  
  # apply WHO reg PAF to globocan region estimates to sum to world - by age
  pafs_who %>% 
    left_join(datreg) %>% 
    mutate(attr=cases*af,
           simattr=cases*simaf) %>% 
    dplyr::select(-casesgrp) -> pafs_who
  
  pafs_who
  
  })
  
  savename <- paste0("_results/aggr_simpaf_who_list.",i,".RData")
  
  save(t.sim3, file = savename)
  }
}
CIswho(sim = "aggr_simpaf")


# run function to create aggregated world totals
CIsreg2 <- function(sim){
  
  for (i in 1:5){
    
  simname <- paste0( "_results/",(sim), "_list.",i,".RData") # should load as t.sim3
  load(simname)
  
  t.sim4 <- lapply(1:2000, function(zzz){
    print(paste0(zzz,".",i))
    
    t <- t.sim3[[zzz]]
    pafs_subregion<- as.data.table(t)
    
  # create world totals
  pafs_subregion %>% 
    group_by(cancer_code, sex, age, sim.n) %>% 
    mutate(casesgrp = sum(cases, na.rm=T),
           attr = sum(attr, na.rm=T),
           simattr= sum(simattr, na.rm=T),
           af = attr/casesgrp,
           simaf = simattr/casesgrp,
           country_label="World",
           country_code=900) %>% 
    dplyr::select(country_label,country_code,cancer_code,sex, age, sim.n, casesgrp, attr, af, simattr, simaf) %>% 
    unique() %>% 
    left_join(datreg %>% 
                filter(country_code==900)) %>% 
    mutate(attr=cases*af,
           simattr=cases*simaf) %>% 
    dplyr::select(-casesgrp) -> pafs_world
  pafs_world
  })
  
  savename <- paste0("_results/aggr_simpaf_world_list.",i,".RData")
  
  save(t.sim4, file = savename)
  
  }
}
CIsreg2(sim = "aggr_simpaf_subreg")

rm(CIs, CIsreg,CIsreg2, CIswho)

#### calculate aggregate age standardised rates (ASRs) ####

# run function to create total for all ages and calculate ASRs per country, HDI, income group
CIsasr <- function(sim){
  
  for (i in 1:5) {
  simname <- paste0( "_results/",(sim), "_list.",i,".RData") # should load as t.sim
  load(simname)
  print(i)
  
  #t.sim <- as.data.table(t.sim)
  paf <- rbindlist(t.sim)
  
  # PAF and ASR by country
  paf %>%
    mutate(attr=ifelse(is.na(attr),0,attr),
           simattr=ifelse(is.na(simattr),0,simattr)) %>% 
    left_join(pops %>%
                dplyr::select(-country_label)) -> paf
  paf %>% 
    filter(sim.n==i*2000) %>% 
    csu_asr(
      var_age = "age",
      var_cases = "attr",
      var_py = "py",
      group_by = c("country_code","country_label", "sex","cancer_code","total"),
      var_age_group = "country_code",
      var_asr = "asr"
    ) %>% 
    left_join(paf %>% 
                csu_asr(
                  var_age = "age",
                  var_cases = "simattr",
                  var_py = "py",
                  group_by = c("country_code","country_label", "sex","cancer_code","sim.n"),
                  var_age_group = "country_code",
                  var_asr = "simasr"
                )) %>% 
    mutate(af=attr/total,
           simaf=simattr/total) %>% 
    filter(country_code!=158)-> pafasr
  
  savename <- paste0("_results/sim_pafasr_list.",i,".RData")
  save(pafasr, file = savename)
  
  # PAF and ASR by HDI
  paf %>%
    filter(country_code!=158) %>% 
    group_by(hdigrp,sex, age, cancer_code,sim.n) %>% 
    mutate(total=sum(total),
           attr=ifelse(is.na(attr),0,attr),
           simattr=ifelse(is.na(simattr),0,simattr)) -> paf_hdi
  paf_hdi %>% 
    filter(sim.n==i*2000) %>% 
    csu_asr(
      var_age = "age",
      var_cases = "attr",
      var_py = "py",
      group_by = c("hdigrp", "sex","cancer_code","total"),
      var_age_group = "hdigrp",
      var_asr = "asr"
    ) %>% 
    left_join(paf_hdi %>% 
                csu_asr(
                  var_age = "age",
                  var_cases = "simattr",
                  var_py = "py",
                  group_by = c("hdigrp", "sex","cancer_code","sim.n"),
                  var_age_group = "hdigrp",
                  var_asr = "simasr"
                )) %>% 
    mutate(af=attr/total,
           simaf=simattr/total,
           country_code=case_when(hdigrp=="Low HDI"~984,
                                  hdigrp=="Medium HDI"~983,
                                  hdigrp=="High HDI"~982,
                                  hdigrp=="Very high HDI"~981,
                                  hdigrp=="Missing HDI"~1001),
           country_label=hdigrp) -> pafasr_hdi
  
  savename <- paste0("_results/sim_pafasr_hdi_list.",i,".RData")
  save(pafasr_hdi, file = savename)
  
  # PAF and ASR by income level
  paf %>%
    filter(country_code!=158) %>% 
    group_by(income_group,sex,age,cancer_code, sim.n) %>%
    mutate(total=sum(total),
           attr=ifelse(is.na(attr),0,attr),
           simattr=ifelse(is.na(simattr),0,simattr))-> paf_inc
  paf_inc %>%
    filter(sim.n==i*2000) %>% 
    csu_asr(
      var_age = "age",
      var_cases = "attr",
      var_py = "py",
      group_by = c("income_group", "sex","cancer_code","total"),
      var_age_group = "income_group",
      var_asr = "asr"
    ) %>% 
    left_join(paf_inc %>% 
                csu_asr(
                  var_age = "age",
                  var_cases = "simattr",
                  var_py = "py",
                  group_by = c("income_group", "sex","cancer_code","sim.n"),
                  var_age_group = "income_group",
                  var_asr = "simasr"
                )) %>% 
    mutate(af=attr/total,
           simaf=simattr/total,
           country_label=income_group,
           country_code=case_when(country_label=="Low income"~1002,
                                  country_label=="Lower middle income"~1003,
                                  country_label=="Upper middle income"~1004,
                                  country_label=="High income"~1005)) -> pafasr_income
  
  savename <- paste0("_results/sim_pafasr_income_list.",i,".RData")
  save(pafasr_income, file = savename)
  
  }
  
}
CIsasr(sim = "aggr_simpaf")

# run function to create total for all ages and calculate ASRs per world region, continent
CIsasr_reg <- function(sim){
  
  for (i in 1:5) {
  simname <- paste0( "_results/",(sim), "_list.",i,".RData") # should load as t.sim3
  load(simname)
  
  #t.sim3 <- as.data.table(t.sim3)
  pafs_subregion <- rbindlist(t.sim3)
  print(i)
  # PAF and ASR by subregion
  pafs_subregion %>%
    mutate(attr=ifelse(is.na(attr),0,attr),
           simattr=ifelse(is.na(simattr),0,simattr))-> pafs_subregion
  pafs_subregion %>% 
    filter(sim.n==i*2000) %>% 
    csu_asr(
      var_age = "age",
      var_cases = "attr",
      var_py = "py",
      group_by = c("continent_label","country_code","country_label", "sex","cancer_code","total"),
      var_age_group = "country_code",
      var_asr = "asr"
    ) %>% 
    left_join(pafs_subregion %>% 
                csu_asr(
                  var_age = "age",
                  var_cases = "simattr",
                  var_py = "py",
                  group_by = c("continent_label","country_code","country_label", "sex","cancer_code","sim.n"),
                  var_age_group = "country_code",
                  var_asr = "simasr"
                )) %>% 
    mutate(af=attr/total,
           simaf=simattr/total) -> pafasr_subreg
  
  savename <- paste0("_results/sim_pafasr_subreg_list.",i,".RData")
  save(pafasr_subreg, file = savename)
  
  # PAF and ASR by continent
  pafs_subregion %>%
    group_by(continent_label,sex, age, cancer_code,sim.n) %>% 
    mutate(total=sum(total),
           attr=ifelse(is.na(attr),0,attr),
           simattr=ifelse(is.na(simattr),0,simattr)) -> pafs_subregion
  pafs_subregion %>% 
    filter(sim.n==i*2000) %>% 
    csu_asr(
      var_age = "age",
      var_cases = "attr",
      var_py = "py",
      group_by = c("continent_label", "sex","cancer_code","total"),
      var_age_group = "continent_label",
      var_asr = "asr"
    ) %>% 
    left_join(pafs_subregion %>% 
                csu_asr(
                  var_age = "age",
                  var_cases = "simattr",
                  var_py = "py",
                  group_by = c("continent_label", "sex","cancer_code","sim.n"),
                  var_age_group = "continent_label",
                  var_asr = "simasr"
                )) %>% 
    mutate(af=attr/total,
           simaf=simattr/total,
           country_code = case_when(continent_label=="Africa"~903,
                                    continent_label=="Latin America and the Caribbean"~904,
                                    continent_label=="North America"~905,
                                    continent_label=="Europe"~908,
                                    continent_label=="Oceania"~909,
                                    continent_label=="Asia"~935),
           country_label=continent_label)-> pafasr_cont
  
  savename <- paste0("_results/sim_pafasr_cont_list.",i,".RData")
  save(pafasr_cont, file = savename)
  }
}
CIsasr_reg(sim = "aggr_simpaf_subreg")

# run function to create total for all ages and calculate ASRs for world total
CIsasr_world <- function(sim){
  for (i in 1:5) {
    
  simname <- paste0( "_results/",(sim), "_list.",i,".RData") # should load as t.sim4
  load(simname)
  
  #t.sim4 <- as.data.table(t.sim4)
  pafs_world <- rbindlist(t.sim4)
  print(i)
  # PAF and ASR world
  pafs_world %>%
    mutate(attr=ifelse(is.na(attr),0,attr),
           simattr=ifelse(is.na(simattr),0,simattr))-> pafs_world
  pafs_world %>% 
    filter(sim.n==i*2000) %>% 
    csu_asr(
      var_age = "age",
      var_cases = "attr",
      var_py = "py",
      group_by = c("sex","cancer_code","total"),
      var_age_group = "sex",
      var_asr = "asr"
    ) %>% 
    left_join(pafs_world %>% 
                csu_asr(
                  var_age = "age",
                  var_cases = "simattr",
                  var_py = "py",
                  group_by = c("sex","cancer_code","sim.n"),
                  var_age_group = "sex",
                  var_asr = "simasr"
                )) %>% 
    mutate(af=attr/total,
           simaf=simattr/total,
           country_code=900,
           country_label="World") -> pafasr_world
  
  
  savename <- paste0("_results/sim_pafasr_world_list.",i,".RData")
  save(pafasr_world, file = savename)
  }
}
CIsasr_world(sim = "aggr_simpaf_world")

# run function to create total for all ages and calculate ASRs per WHO region
CIsasr_who <- function(sim){
  for (i in 1:5){
    
  simname <- paste0( "_results/",(sim), "_list.",i,".RData") # should load as t.sim3
  load(simname)
  
  #t.sim3 <- as.data.table(t.sim3)
  pafs_who <- rbindlist(t.sim3)
  print(i)
  # PAF and ASR world
  pafs_who %>%
    mutate(attr=ifelse(is.na(attr),0,attr),
           simattr=ifelse(is.na(simattr),0,simattr))-> pafs_who
  pafs_who %>% 
    filter(sim.n==i*2000) %>% 
    csu_asr(
      var_age = "age",
      var_cases = "attr",
      var_py = "py",
      group_by = c("who_region_label","country_code", "sex","cancer_code","total"),
      var_age_group = "country_code",
      var_asr = "asr"
    ) %>% 
    left_join(pafs_who %>% 
                csu_asr(
                  var_age = "age",
                  var_cases = "simattr",
                  var_py = "py",
                  group_by = c("who_region_label","country_code", "sex","cancer_code","sim.n"),
                  var_age_group = "country_code",
                  var_asr = "simasr"
                )) %>% 
    mutate(af=attr/total,
           simaf=simattr/total,
           country_label=who_region_label) -> pafasr_who
  
  
  savename <- paste0("_results/sim_pafasr_who_list.",i,".RData")
  save(pafasr_who, file = savename)
  }
}
CIsasr_who(sim = "aggr_simpaf_who")



#### combine country and region results ####

## load all 5 sets of 2000 simulations per vector
load("_results/sim_pafasr_list.1.RData")
pafasr1 <- pafasr
load("_results/sim_pafasr_list.2.RData")
pafasr1 <- pafasr1 %>% rbind(pafasr)
load("_results/sim_pafasr_list.3.RData")
pafasr1 <- pafasr1 %>% rbind(pafasr)
load("_results/sim_pafasr_list.4.RData")
pafasr1 <- pafasr1 %>% rbind(pafasr)
load("_results/sim_pafasr_list.5.RData")
pafasr <- pafasr1 %>% rbind(pafasr)
rm(pafasr1)

load("_results/sim_pafasr_world_list.1.RData")
pafasr_world1 <- pafasr_world
load("_results/sim_pafasr_world_list.2.RData")
pafasr_world1 <- pafasr_world1 %>% rbind(pafasr_world)
load("_results/sim_pafasr_world_list.3.RData")
pafasr_world1 <- pafasr_world1 %>% rbind(pafasr_world)
load("_results/sim_pafasr_world_list.4.RData")
pafasr_world1 <- pafasr_world1 %>% rbind(pafasr_world)
load("_results/sim_pafasr_world_list.5.RData")
pafasr_world <- pafasr_world1 %>% rbind(pafasr_world)
rm(pafasr_world1)

load("_results/sim_pafasr_cont_list.1.RData")
pafasr_cont1 <- pafasr_cont
load("_results/sim_pafasr_cont_list.2.RData")
pafasr_cont1 <- pafasr_cont1 %>% rbind(pafasr_cont)
load("_results/sim_pafasr_cont_list.3.RData")
pafasr_cont1 <- pafasr_cont1 %>% rbind(pafasr_cont)
load("_results/sim_pafasr_cont_list.4.RData")
pafasr_cont1 <- pafasr_cont1 %>% rbind(pafasr_cont)
load("_results/sim_pafasr_cont_list.5.RData")
pafasr_cont <- pafasr_cont1 %>% rbind(pafasr_cont)
rm(pafasr_cont1)

load("_results/sim_pafasr_subreg_list.1.RData")
pafasr_subreg1 <- pafasr_subreg
load("_results/sim_pafasr_subreg_list.2.RData")
pafasr_subreg1 <- pafasr_subreg1 %>% rbind(pafasr_subreg)
load("_results/sim_pafasr_subreg_list.3.RData")
pafasr_subreg1 <- pafasr_subreg1 %>% rbind(pafasr_subreg)
load("_results/sim_pafasr_subreg_list.4.RData")
pafasr_subreg1 <- pafasr_subreg1 %>% rbind(pafasr_subreg)
load("_results/sim_pafasr_subreg_list.5.RData")
pafasr_subreg<- pafasr_subreg1 %>% rbind(pafasr_subreg)
rm(pafasr_subreg1)

load("_results/sim_pafasr_who_list.1.RData")
pafasr_who1 <- pafasr_who
load("_results/sim_pafasr_who_list.2.RData")
pafasr_who1 <- pafasr_who1 %>% rbind(pafasr_who)
load("_results/sim_pafasr_who_list.3.RData")
pafasr_who1 <- pafasr_who1 %>% rbind(pafasr_who)
load("_results/sim_pafasr_who_list.4.RData")
pafasr_who1 <- pafasr_who1 %>% rbind(pafasr_who)
load("_results/sim_pafasr_who_list.5.RData")
pafasr_who <- pafasr_who1 %>% rbind(pafasr_who)
rm(pafasr_who1)

load("_results/sim_pafasr_hdi_list.1.RData")
pafasr_hdi1 <- pafasr_hdi
load("_results/sim_pafasr_hdi_list.2.RData")
pafasr_hdi1 <- pafasr_hdi1 %>% rbind(pafasr_hdi)
load("_results/sim_pafasr_hdi_list.3.RData")
pafasr_hdi1 <- pafasr_hdi1 %>% rbind(pafasr_hdi)
load("_results/sim_pafasr_hdi_list.4.RData")
pafasr_hdi1 <- pafasr_hdi1 %>% rbind(pafasr_hdi)
load("_results/sim_pafasr_hdi_list.5.RData")
pafasr_hdi <- pafasr_hdi1 %>% rbind(pafasr_hdi)
rm(pafasr_hdi1)

load("_results/sim_pafasr_income_list.1.RData")
pafasr_income1 <- pafasr_income
load("_results/sim_pafasr_income_list.2.RData")
pafasr_income1 <- pafasr_income1 %>% rbind(pafasr_income)
load("_results/sim_pafasr_income_list.3.RData")
pafasr_income1 <- pafasr_income1 %>% rbind(pafasr_income)
load("_results/sim_pafasr_income_list.4.RData")
pafasr_income1 <- pafasr_income1 %>% rbind(pafasr_income)
load("_results/sim_pafasr_income_list.5.RData")
pafasr_income <- pafasr_income1 %>% rbind(pafasr_income)
rm(pafasr_income1)

# combine totals to create one data.table for all groups
pafasr %>% 
  rbind(pafasr_world) %>% 
  rbind(pafasr_cont %>% 
          dplyr::select(-continent_label) %>% 
          filter(country_code!=905)) %>% # remove North America
  rbind(pafasr_subreg %>% 
          dplyr::select(-continent_label)) %>% 
  rbind(pafasr_who %>% 
          dplyr::select(-who_region_label)) %>% 
  rbind(pafasr_hdi %>% 
          dplyr::select(-hdigrp)) %>% 
  rbind(pafasr_income %>% 
          dplyr::select(-income_group)) -> pafasr_all


#### get 95% CIs for attributable cases, PAF and ASR per group ####

pafasr_all %>% 
  group_by(country_code,country_label,sex,cancer_code, total, attr, af, asr) %>% 
  mutate(af=ifelse(is.na(af),0,af),
         simaf=ifelse(is.na(simaf),0,simaf)) %>% 
  summarise(attrLCI=quantile(simattr, probs = 0.025),      # calculate 95% CIs of attributable cases
            attrUCI=quantile(simattr, probs = 0.975),
            afLCI=quantile(simaf, probs = 0.025),      # calculate 95% CIs of attributable cases
            afUCI=quantile(simaf, probs = 0.975),
            asrLCI=quantile(simasr, probs = 0.025),      # calculate 95% CIs of attributable cases
            asrUCI=quantile(simasr, probs = 0.975)) -> pafasr_all



#### save results ####
# paf results
write.csv(pafasr_all,"_results/SLT_PAF_asr_all_CIs.csv",row.names=FALSE)


