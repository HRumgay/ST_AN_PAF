################################################################
# Harriet Rumgay - International Agency for Research on Cancer
# 24 July 2024
################################################################
#
# Script to calculate PAF for smokeless tobacco by product type in
# Bangladesh, India, Pakistan & Papua New Guinea


#### set up script ####

# load following libraries
library(tidyverse)
library(data.table)
library(Rcan)


set.seed(999) # set random sampling

# set working directory
setwd("~/STAN_PAF")


##### load data #####

# ST/AN prevalence by country, sex, age
dat <-read.csv("_data/smokeless tobacco product prevalence.csv") %>% 
   rename(sex=sex_code) %>% 
   as.data.table()

# RRs to use for each country and cancer site
rrs <-read.csv("_data/SLT_product type RRs_19.02.24.csv") %>% 
  as.data.table() %>% 
  select(product.type.name.in.survey,product_code,RR,LCI,UCI,category,hex_cat,name_fig) %>% 
  filter(!is.na(LCI)) %>% 
  mutate(RR=as.numeric(RR))

# load Globocan 2022 estimates
DATA_GLOB	<- as.data.table(read.csv("_data/Globocan2022.csv")) %>% as.data.table()
DATA_GLOB %>% 
  group_by(country_code,cancer_code,sex,type) %>% 
  mutate(total=sum(cases,na.rm=T)) -> DATA_GLOB

# select incidence, oral cancer, and Bangladesh, India, Pakistan & Papua New Guinea
DATA_GLOB %>% 
  filter(type==0,cancer_code==1,country_code%in%c(50,356,586,598)) %>%
  ungroup() %>% 
  select(-type) -> can


#### data preparation ####

# select variables needed
dat %>% 
  rename(country_label=country_name,
         year=survey_year_start) %>% 
  select(year,country_label,sex,start_age,end_age,product_code,type,type_short:end_ci,sample_size,standard.error) -> dat2

# tidy CI and age group labels
dat2 %>% 
  filter((!(start_age==15&end_age==100)&!(start_age==15&end_age==64)&sex!=0)|(product_code%in%c("p13","p14","p15")&(start_age==15&end_age==100))) %>%
  mutate(prevalence=ifelse(prevalence==".",0,as.numeric(prevalence)),
         start_ci=ifelse(prevalence==".",0,as.numeric(start_ci)),
         end_ci=ifelse(prevalence==".",0,as.numeric(end_ci)),
         sample_size=as.numeric(gsub(",", "", sample_size)),
         prev.age=case_when(start_age==15~4,
                            start_age==25~6,
                            start_age==35~8,
                            start_age==45~10,
                            start_age==55~12,
                            start_age==65~14)) -> dat2

# fill sample size if missing, calculate standard error and standard deviation
dat2 %>% 
  mutate(sample_size=ifelse(is.na(sample_size),2000,sample_size),
         prev = prevalence/100,
         prev_lci = start_ci/100,
         prev_uci = end_ci/100,
         se=case_when(!is.na(prev_lci)&!is.na(prev_uci)~((prev_uci-prev_lci)/(2*1.96)),
                      !is.na(standard.error)~standard.error,
                      TRUE~sqrt((prev*(1-prev)/sample_size))),
         sd=se*sqrt(sample_size)) -> dat2 


# expand age groups to all ages 15-75
dat2 %>% 
  bind_rows(dat2 %>% 
              filter(prev.age==4) %>% 
              mutate(prev.age=5)) %>% 
  bind_rows(dat2 %>% 
              filter(prev.age==6) %>% 
              mutate(prev.age=7)) %>% 
  bind_rows(dat2 %>% 
              filter(prev.age==6&!str_detect(country_label,"Papua")) %>% 
              mutate(prev.age=8)) %>%
  bind_rows(dat2 %>% 
              filter((prev.age==6&!str_detect(country_label,"Papua"))|prev.age==8) %>% 
              mutate(prev.age=9)) %>% 
  bind_rows(dat2 %>% 
              filter(prev.age==10) %>% 
              mutate(prev.age=11)) %>% 
  bind_rows(dat2 %>% 
              filter(prev.age==10&!str_detect(country_label,"Papua")) %>% 
              mutate(prev.age=12)) %>% 
  bind_rows(dat2 %>% 
              filter((prev.age==10&!str_detect(country_label,"Papua"))|prev.age==12)%>% 
              mutate(prev.age=13)) %>% 
  bind_rows(dat2 %>% 
              filter(prev.age==12)%>% 
              mutate(prev.age=14)) %>% 
  bind_rows(dat2 %>% 
              filter(prev.age%in%c(12,14)) %>% 
              mutate(prev.age=15)) %>% 
  bind_rows(dat2 %>% 
              filter(prev.age%in%c(12,14)) %>% 
              mutate(prev.age=16)) %>% 
  mutate(age = prev.age+2) -> dat2
 
# add extra age groups for p13,p14,p15 where data not provided by age group
dat2 %>% 
  filter(product_code%in%c("p13","p14","p15"),prev.age!=5) %>% 
  group_by(year,sex,product_code,type_short) %>% 
  slice(rep(1:n(),each=13)) %>%
  mutate(age=row_number(),
         age=age+5,
         prev.age=age-2)%>% 
  bind_rows(dat2 %>% 
              filter(!product_code%in%c("p13","p14","p15"))) -> dat2

# add age groups 1 to 5 with prev=0 to include in ASR calculation
dat2 %>% 
  bind_rows(dat2 %>% 
              filter(age==6) %>% 
              group_by(year,sex,product_code,type,type_short,country_label) %>% 
              slice(rep(1:n(),each=5)) %>% 
              mutate(age=row_number()) %>% 
              mutate(across(prevalence:sd,~0))) -> dat2

# merge with relative risk and cancer incidence 
dat2 %>% 
  left_join(rrs %>% 
              select(product_code:UCI) %>% unique()) %>% 
  left_join(can %>% select(-py) %>% filter(sex!=0)) -> dat3

#### calculate PAFs ####

# calculate attributable fractions and attributable cases
dat3 %>% 
  mutate(af = (prev*(RR-1))/(1+prev*(RR-1)),
         attr = af*cases) %>%
  group_by(year,country_code,sex, product_code, type_short) %>% 
  mutate(totattr = sum(attr, na.rm=T),
         totpaf = (totattr/total)) %>% 
  select(year,country_label,country_code,sex,age,product_code,type,type_short,prev,se,RR,LCI,UCI,cancer_code,cases,total,af,attr,totattr,totpaf)-> paf


#### create simulations of PAF ####

# set number of simulations
nnn<-10000
# set up simulations template
sim <- data.table(c(1:nnn,1:nnn),c(rep(1,nnn),rep(2,nnn)))
names(sim) <- c("sim.n","sex")

pafCIs <- as.data.table(paf)
# merge paf file with simulations template
pafCIs <- merge(pafCIs,sim, allow.cartesian = T)

# get random samples of prevalence and RRs
pafCIs %>%
  group_by(year,country_code,product_code, type_short,sex,age) %>% 
  mutate(simprev = rnorm(nnn,prev,se),
         simrr = exp(rnorm(nnn,log(RR),((log(UCI)-log(LCI))/(2*1.96))))) -> pafCIs

# calculate attributable fraction and attributable cases per simulation per age and all ages combined
pafCIs %>% 
  mutate(simaf = (simprev*(simrr-1))/(1+simprev*(simrr-1)),
         simattr = ifelse(is.na(simaf*cases),0,simaf*cases),
         simaf = ifelse(is.na(simaf),0,simaf),
         simattr = ifelse(is.na(simattr),0,simattr)) %>% 
  group_by(year,country_code,product_code, type_short,sex,sim.n) %>% 
  mutate(simtotattr = sum(simattr, na.rm=T),
         simtotpaf = ifelse(is.na(simtotattr/total),0,simtotattr/total)) -> pafCIs

# create totals for both sexes combined
pafCIs %>% 
  bind_rows(pafCIs %>% 
              group_by(year,country_code,product_code, type_short,age,sim.n) %>% 
              mutate(sex=0,
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
              select(-prev,-se,-simprev,-simrr) %>% unique()) -> pafCIs

save(pafCIs, file = "_results/pafCIs_products.RData")

#### calculate age standardised rates (ASRs) ####

# calculate ASRs for simulations
pafCIs %>% 
  full_join(can %>% ungroup() %>% select(country_code,country_label,sex,age,py) %>% unique()) -> pafCIs
pafCIs%>% 
  filter(sim.n==1) %>% 
  csu_asr(
    var_age = "age",
    var_cases = "attr",
    var_py = "py",
    group_by = c("country_code","country_label","year","type","product_code","type_short", "sex","total"),
    var_age_group = "country_code",
    var_asr = "asr",
  ) %>% 
  left_join(pafCIs %>% 
              csu_asr(
                var_age = "age",
                var_cases = "simattr",
                var_py = "py",    
                #first_age = 6,
                #last_age = 18,
                group_by = c("country_code","country_label","year","type","product_code","type_short", "sex","sim.n"),
                var_age_group = "country_code",
                var_asr = "simasr",
              )) %>% 
  mutate(af=attr/total,
         simaf=simattr/total) -> pafasrCI


# calculate CIs of attributable cases for all ages combined
pafCIs %>% 
  group_by(year,country_label,country_code,product_code, type,type_short, sex,total, totattr) %>% 
  summarise(totattrLCI=quantile(simtotattr, probs = 0.025),      # calculate 95% CIs of attributable cases
            totattrUCI=quantile(simtotattr, probs = 0.975)) -> pafCIsattr
# calculate CIs of PAFs for all ages combined
pafCIs %>% 
  group_by(year,country_label,country_code,product_code,type, type_short,sex,total, totpaf) %>% 
  summarise(totpafLCI=quantile(simtotpaf, probs = 0.025),      # calculate 95% CIs of attributable cases
            totpafUCI=quantile(simtotpaf, probs = 0.975)) -> pafCIsaf
# calculate CIs of ASRs for all ages combined
pafasrCI %>% 
  group_by(year,country_label,country_code,product_code,type, type_short,sex,total, asr) %>% 
  summarise(asrLCI=quantile(simasr, probs = 0.025),      # calculate 95% CIs of attributable cases
            asrUCI=quantile(simasr, probs = 0.975)) -> pafasrCI

# merge 3 new dataframes
pafCIsattr %>% 
  left_join(pafCIsaf) %>% 
  left_join(pafasrCI)-> pafCIs
rm(pafCIsaf,pafCIsattr,pafasrCI)


#### save results ####
# paf results
write.csv(pafCIs,"_results/SLT_PAF_products_CIs.csv",row.names=FALSE)


