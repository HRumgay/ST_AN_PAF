################################################################
# Harriet Rumgay - International Agency for Research on Cancer
# 24 July 2024
################################################################
#
# Script to calculate PAF for smokeless tobacco and areca nut 
# consumption per age and country


#### set up script ####

# load following libraries
library(tidyverse)
library(data.table)
library(xlsx)

# set working directory
setwd("~/STAN_PAF")

set.seed(999) # set random sampling

##### load data #####

# SLT prevalence by country, sex, age
dat <-read.csv("_data/Data_Smokeless_focus20240718.csv") %>% 
  dplyr::rename(sex=sex_code) %>% 
  mutate(country_code=as.numeric(country_code),
         country_iso3=wb_country_abv) %>% 
  as.data.table()

# RRs to use for each country and cancer site
rrs <-read.csv("_data/RRs selection_17.05.23.csv") %>% 
  as.data.table() %>% 
  dplyr::select(-source:-reason)

# load Globocan 2022 estimates
DATA_GLOB	<- as.data.table(read.csv("_data/Globocan2022.csv")) %>% as.data.table()

# create total all ages combined
DATA_GLOB %>% 
  group_by(country_code,cancer_code,sex,type) %>% 
  mutate(total=sum(cases,na.rm=T)) -> DATA_GLOB

# select incidence and cancer types
DATA_GLOB %>% 
  filter(type==0,cancer_code %in% c(1,3,5,6,13),sex!=0,country_code<900) %>%
  dplyr::select(-py,-type) -> can

# missing countries selection and countries to impute from
prev.miss <- as.data.table(read.csv("_data/missing.countries.prevalence.csv"))


#### data preparation ####

# add area name to prevalence data for world regions and filter to ages 15-75
dat %>% 
  mutate(area_name = case_when(area_name==""~sub_continent_name,
                               TRUE~area_name)) %>%
  dplyr::select(-wb_country_abv:-pucem,-who_member_state:-country_label,-country_iso3_2:-select_chk_d) %>% 
  filter(age>=4, age<=16) %>% 
  group_by(un_name, country_iso3,sex_name) %>% 
  mutate(maxage = (max(age)),
         minage = min(age)) -> dat2

# fill sample size if missing, calculate standard error and standard deviation
dat2 %>% 
  mutate(sample_size=ifelse(is.na(sample_size),2000,sample_size),
         prev=prev/100,
         prev_lci=prev_lci/100,
         prev_uci=prev_uci/100,
         se=case_when(!is.na(mean_se)~mean_se/100,
                      !is.na(prev_lci)&!is.na(prev_uci)~((prev_uci-prev_lci)/(2*1.96)),
                      is.na(mean_se)&is.na(prev_uci)~sqrt((prev*(1-prev)/sample_size))),
         sd=se*sqrt(sample_size)) -> dat2 


# expand age groups to all ages 15-75
dat2 %>% 
  bind_rows(dat2 %>% 
              filter(minage>5) %>% 
              filter(age==minage)%>% 
              mutate(age=5)) %>% 
  bind_rows(dat2 %>% 
              filter(minage>4) %>% 
              filter(age==minage)%>% 
              mutate(age=4)) %>% 
  bind_rows(dat2 %>% 
              filter(maxage<11) %>% 
              filter(age==maxage)%>% 
              mutate(age=11)) %>%
  bind_rows(dat2 %>% 
              filter(maxage<12) %>% 
              filter(age==maxage)%>% 
              mutate(age=12)) %>%
  bind_rows(dat2 %>% 
              filter(maxage<13) %>% 
              filter(age==maxage)%>% 
              mutate(age=13)) %>%
  bind_rows(dat2 %>% 
              filter(maxage<14) %>% 
              filter(age==maxage)%>% 
              mutate(age=14)) %>%
  bind_rows(dat2 %>% 
              filter(maxage<15) %>% 
              filter(age==maxage)%>% 
              mutate(age=15)) %>%
  bind_rows(dat2 %>% 
              filter(maxage<16) %>% 
              filter(age==maxage)%>% 
              mutate(age=16)) %>% 
  mutate(p.age = age,    # prevalence age group vs globocan age group
         age = age+2) -> dat2


# merge with fill countries file
dat2 %>% 
  ungroup() %>% 
  dplyr::select(-country_code) %>% 
  full_join(prev.miss %>% 
              dplyr::select(-region_label:-who_region_label) %>% 
              rename(country_iso3=ISO3)) %>% 
  mutate(country_label=ifelse(is.na(country_label),un_name,country_label),
         missing=ifelse(is.na(missing),0,missing)) %>% 
  as.data.table() -> d2

dtemp <- d2[country_code==4,.(age,sex,missing)][,missing:=1]

# impute gaps for missing countries based on average prevalence of fill countries
dlist <- lapply(unique(prev.miss[missing==1]$country_code), function(REG){ 
  dtemp %>% 
    full_join(d2 %>% 
                filter(country_code==REG) %>% 
                dplyr::select(-age,-sex)) %>% 
    rbind(d2 %>% 
            filter(missing==0)) %>% 
    group_by(sex,age) %>% 
    mutate(across(c(fill.1.code,fill.2.code,fill.3.code,fill.4.code,fill.5.code), ~replace_na(.x, 0)),
           prev = case_when(missing==1~mean(prev[country_code%in%c(fill.1.code,fill.2.code,fill.3.code,fill.4.code,fill.5.code)],na.rm=T),
                            TRUE~prev),
           se = case_when(missing==1~mean(se[country_code%in%c(fill.1.code,fill.2.code,fill.3.code,fill.4.code,fill.5.code)],na.rm=T),
                            TRUE~se)) %>% 
    filter(country_code==REG) -> d3
  d3
})
d3 <- do.call(rbind.data.frame, dlist)

d3 %>% 
  rbind(d2 %>% filter(missing==0)) -> dat2
rm(dat,d2,d3,dtemp,dlist)

# merge with relative risk and cancer incidence 
dat2 %>% 
  dplyr::select(country_code,country_label,sex,age,prev,se) %>% 
  left_join(rrs %>% #select(-country_label) %>% 
              mutate(RR=case_when(LCI<1~1,       # replace nonsig RRs with 1
                                  TRUE~RR))) %>% # match RRs
  full_join(can %>% dplyr::select(-country_label)) -> dat3      # match incidence data


#### calculate PAFs ####

# calculate attributable fractions and attributable cases
dat3 %>% 
  mutate(af = (prev*(RR-1))/(1+(prev*(RR-1))),
         attr = af*cases) %>% 
  group_by(country_code,sex, cancer_label,cancer_code) %>% 
  mutate(totattr = sum(attr, na.rm=T),
         totpaf = (totattr/total)) -> paf
# filter to oral cancer only
paf %>% 
  filter(cancer_code==1)-> paf

# add country codes for countries without codes
paf %>%
  mutate(country_code=case_when(country_label=="American Samoa"~16,
                                country_label=="Pohnpei"~583,
                                country_label=="Kosrae"~896,
                                country_label=="Northern Mariana Islands"~580,
                                country_label=="Marshall Islands"~584,
                                TRUE~country_code)) -> paf

#### create simulations of PAF ####

# set number of simulations
nnn <- 10000
# set up simulations template
sim <- data.table(c(1:nnn,1:nnn),c(rep(1,nnn),rep(2,nnn)))
names(sim) <- c("sim.n","sex")

pafCIs <- as.data.table(paf)
# remove unnecessary variables
pafCIs$country_label  <- NULL;pafCIs$region_label  <- NULL; pafCIs$subregion_label  <- NULL; pafCIs$WHO_subregion  <- NULL; pafCIs$cancer_label <- NULL
# merge paf file with simulations template
pafCIs <- merge(pafCIs,sim, allow.cartesian = T)

# get random samples of prevalence and RRs
pafCIs %>%
  group_by(country_code, sex,age,cancer_code) %>%
  mutate(simprev = rnorm(nnn,prev,se),   # using se in rnorm to get random samples of mean
         simrr = exp(rnorm(nnn,log(RR),((log(UCI)-log(LCI))/(2*1.96))))) -> pafCIs

# calculate attributable fraction and attributable cases per simulation per age and all ages combined
pafCIs %>%
  mutate(simaf = (simprev*(simrr-1))/(1+simprev*(simrr-1)),
         simattr = ifelse(is.na(simaf*cases),0,simaf*cases)) %>%
  group_by(country_code,sex, cancer_code,sim.n) %>%
  # totals for all ages combined
  mutate(simtotattr = sum(simattr, na.rm=T), 
         simtotpaf = ifelse(is.na(simtotattr/total),0,simtotattr/total)) %>% 
  data.table() -> pafCIs

save(pafCIs, file = "_results/pafCIs.RData")
rm(dat2, dat3, DATA_GLOB,nnn,sim,rrs,can)

# create list for simulations functions
pafCIs.l <- lapply(1:10000, function(sim){
  print(sim)
  t<-pafCIs[sim.n==sim,]
})
save(pafCIs.l, file = "_results/pafCIs_list.RData")
rm(pafCIs,pafCIs.l,clab,clib.list,prev.miss)
  
#### continue simulations functions ####

# create both sexes totals
CIs2 <- function(sim){
  simname <- paste0( "_results/",(sim), "_list.RData") # should load as pafCIs
  load(simname)
  
  t.sim2 <- lapply(1:10000, function(zzz){
    print(zzz)
    
    t <- t.sim[[zzz]]
    t<- as.data.table(t)
  
  # create totals for both sexes combined
  t %>% 
    group_by(country_code,cancer_code,age,sim.n) %>%
    mutate(sex=0,
           sex_label="Both sexes",
           cases=sum(cases, na.rm=T),
           total=sum(total, na.rm=T),
           attr=sum(attr, na.rm=T),
           af=attr/cases,
           totattr=sum(totattr, na.rm=T),
           totpaf=totattr/total,
           simattr=sum(simattr, na.rm=T),
           simtotattr=sum(simtotattr, na.rm=T),
           simaf=simattr/cases,
           simtotpaf=ifelse(is.na(simtotattr/total),0,simtotattr/total)) %>%
    ungroup() %>%
    dplyr::select(-prev,-se,-simprev,-simrr) -> pafCIs2
  
  #print("filter to both sexes only")
  pafCIs2 %>% 
    unique() -> pafCIs2
  
  #print("add to other m f estimates")
  t %>% 
    bind_rows(pafCIs2) -> pafCIs2
  })
  
  save(t.sim2, file = "_results/pafCIs2_list.RData")
}
CIs2(sim = "pafCIs")


load("_results/pafCIs2_list.RData")

pafCIs2 <- rbindlist(t.sim2)
rm(t.sim2)

# calculate CIs of attributable cases by age
pafCIs2 %>% 
  group_by(country_code,sex,age,cancer_code, cases,attr) %>% 
  summarise(attrLCI=quantile(simattr, probs = 0.025),      # calculate 95% CIs of attributable cases
            attrUCI=quantile(simattr, probs = 0.975)) -> pafCIs.a
# calculate CIs of attributable cases for all ages combined
pafCIs2 %>% 
  group_by(country_code,sex,cancer_code, total, totattr) %>% 
  summarise(totattrLCI=quantile(simtotattr, probs = 0.025),      # calculate 95% CIs of attributable cases
            totattrUCI=quantile(simtotattr, probs = 0.975)) -> pafCIsattr
# calculate CIs of PAFs for all ages combined
pafCIs2 %>% 
  group_by(country_code,sex,cancer_code,total, totpaf) %>% 
  summarise(totpafLCI=quantile(simtotpaf, probs = 0.025),      # calculate 95% CIs of attributable cases
            totpafUCI=quantile(simtotpaf, probs = 0.975)) -> pafCIsaf
rm(pafCIs2)


# merge 3 new dataframes
pafCIs.a %>% 
  left_join(pafCIsattr) %>% 
  left_join(pafCIsaf) -> pafCIs

rm(pafCIs.a,pafCIsaf,pafCIsattr)


#### create final dataset of PAF and attributable cases per country ####

# filter to all ages combined
pafCIs %>% 
  left_join(paf %>% 
              filter(!is.na(ISO3)) %>% 
              ungroup() %>% 
              dplyr::select(region_label,WHO_subregion,subregion_label,ISO3,country_label, country_code) %>% 
              unique()) %>%
  filter(age>5) %>% 
  ungroup() %>% 
  dplyr::select(region_label,WHO_subregion,subregion_label,ISO3,country_label, country_code,sex,cancer_code,totpaf,totpafLCI,totpafUCI,totattr,totattrLCI,totattrUCI,total) %>% 
  unique() -> allpaf
allpaf <- data.table(allpaf)


#### save results ####
# paf results
write.csv(allpaf,"_results/SLT_PAF_CIs.csv",row.names=FALSE)
# age-specific paf results
write.csv(pafCIs,"_results/SLT_PAF_CIs_agespec.csv",row.names=FALSE)



