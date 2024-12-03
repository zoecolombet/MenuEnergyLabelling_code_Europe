

rm(list=ls())

##doing mc monte carlo
library(dplyr)
library(data.table)
library(readxl)
library(fst)
library(gamlss)
library(demography)
library(xlsx)
library(foreach)
library(qs) # save R objects to disk
library(haven)
library(doParallel)
library(gamlss.dist)
library(nlme)
library(parallel)
library(bigstatsr)


setDTthreads(1L, restore_after_fork = FALSE)
threads_fst(1L, reset_after_fork = FALSE)

setwd("~/Simulation_Germany/Germany")

#cREATION OF THE POPULATION
# we only need year 2022 (init year) and ages 30-90
# Set simulation parameters
montecarlo <- 200      # Number of iterations
cluster_number <- 4  # Number of CPU cores to use
init_year <- 22L     # initial year
sim_hor <- 20L       # simulation horizon
min_age <- 30L       # min age for the population / need to be the age of interest (especially regarding the RR) minus lagtime. For us here, RR start a 40, we have 6y lag time, we might want to start at 34. For comparibility with UK model, lets keep 30
max_age <- 89L       # max age for the population
#should be 94  // issue with 90 as for IMD its 90+
SES <- 1:3           #SES education level
scale_factor <- 500  # factor to downscale the country population #was 1000 before
SSB_direct_effect <- FALSE # If false only consider the indirect effect through BMI

#DECISION regarding the coverage of the policy
#assuming that the policy is not at all implemented in Germany for now: 0% of the businesses have the policy
label.policy.coverage.large <- 0.09  #implementing the policy only in large business, that are 9% of the number of food outlets, with none having label already = 9% coverage
label.policy.coverage.full <- 1      #policy implemented in every businesses: 100% of the food sector, with none having label already = 100% coverage
label.policy.coverage.large.turnover <- 0.21  #implementing the policy only in large business, that are 21% of the food sector turnover,with none having label already = 21% coverage
label.policy.coverage.full.turnover <- 1      #policy implemented in every businesses: 100% of the food sector, with 0 already having labels,with none having label already = 100% coverage
SSB.policy.coverage <- 1                      #the SSB will impact every SSB so every consummer will be impacted


#database needed
#population 
onspop_import <- read.csv("./Population_Mortality/pop_estimates_education.csv")
setDT(onspop_import)
#WITH NO MISSING DATA
#summary(onspop_import)
# table(onspop_import$sex, useNA="a")
onspop_import <- onspop_import[sex!="",]
# table(onspop_import$edu_cat, useNA="a")
onspop_import <- onspop_import[year>=2000+init_year & year<=2000+init_year+sim_hor & age>=min_age & age<=max_age,]
onspop_import <- onspop_import[, .(year,
                                   sex=as.factor(ifelse(sex=="male","1","2")),
                                   age,
                                   SES=ifelse(edu_cat=="low","1",
                                              ifelse(edu_cat=="middle","2","3")),
                                   pops)]
onspop_import <- dcast(onspop_import, sex + age + SES ~ year, value.var = "pops")


#extracting pop for further calculation on the results
onspop_import_number <- onspop_import

year_columns <- as.character(2022:2041)
pop_SES <- onspop_import_number[ 
  , lapply(.SD, sum, na.rm = TRUE),  
  by = SES,                         
  .SDcols = year_columns            
]

pop_SES <- pop_SES[,
                   .(total_pop = rowSums(.SD, na.rm = TRUE)),
                   by = SES, 
                   .SDcols = year_columns]

print(pop_SES)

#mortality
onsmort_import2 <- read.csv("./Population_Mortality/mort_forecast.csv")
setDT(onsmort_import2)
#WITH NO MISSING DATA
#summary(onsmort_import2)
# table(onsmort_import2$sex, useNA="a")
table(onsmort_import2$type)
onsmort_import2 <- onsmort_import2[year>2000,]
onsmort_import2 <- onsmort_import2[, .(sex=as.factor(ifelse(sex=="male","1","2")),
                                       reg, age,year,
                                       type=ifelse(type=="haem_stroke","haem",
                                                   ifelse(type=="isch_stroke","isch",type)),
                                       mx_total,
                                       mx_total_se = ((log(mx_total_99) - log(mx_total))/qnorm(0.995) + (log(mx_total_1) - log(mx_total))/qnorm(0.005))/2)]
onsmort_import2 <- dcast(onsmort_import2, sex + age + year ~ type, value.var = c("mx_total","mx_total_se"))
onsmort_import2 <- onsmort_import2[,.(sex,age,year,mx_chd=mx_total_chd,mx_chd_se=mx_total_se_chd,mx_isch=mx_total_isch,mx_isch_se=mx_total_se_isch,mx_haem=mx_total_haem,mx_haem_se=mx_total_se_haem)]

##AS WE DONT HAVE MORTALITY BY EDUCATION LEVEL, LETS ASSUME THAT THE RATE IS EQUALLY DIVIDED, AND LATER IN THE CODE, WE WILL APPLY A RISK OF DEATH BY EDUCATION
# CKNote: This is a rate so it is not divided. We assume initially  that it is the same for all SES levels
test <- CJ(sex = c("1","2"), age = 20:90,year = 2001:2052,SES=1:3)
onsmort_import2 <- merge(test, onsmort_import2, by=c("sex","age","year"), all.x=T)
momo <- c("mx_chd","mx_chd_se","mx_isch","mx_isch_se","mx_haem","mx_haem_se")
onsmort_import <- copy(onsmort_import2)
onsmort_import <- onsmort_import[, (momo) := lapply(.SD, function(x) x/1), .SDcols = momo]

onsmort_import[, year:=year-2000]
rm(onsmort_import2)

####AS WE DONT HAVE THE MORTALITY RATES / NUMBER by education level in Germany, we are using the Age-adjusted mortality rate ratios for education category, ages 30–59 years, 2013
#from the Grigoriev paper (https://bmjopen.bmj.com/content/9/10/e028001)
# and we assume that these difference by education are the same for every type of death, including CHD, Stroke (isch and haem)
rr.education.table <- CJ(sex = c("1","2"), SES=as.factor(1:3))
rr.education.table[, `:=` (rr.education_e =1, rr.education_l=1, rr.education_u=1)]
rr.education.table[SES==1 & sex==1, `:=` (rr.education_e =1.69, rr.education_l=1.64, rr.education_u=1.75)]
rr.education.table[SES==1 & sex==2, `:=` (rr.education_e =1.65, rr.education_l=1.57, rr.education_u=1.74)]
rr.education.table[SES==2 & sex==1, `:=` (rr.education_e =1.35, rr.education_l=1.31, rr.education_u=1.40)]
rr.education.table[SES==2 & sex==2, `:=` (rr.education_e =1.23, rr.education_l=1.16, rr.education_u=1.29)]

rr.education.table[, `:=` (rr.education_se=((log(rr.education_u) - log(rr.education_e))/1.96 + (log(rr.education_l) - log(rr.education_e))/-1.96)/2)]

#head(rr.education.table)


#GAMLSS bmi
tbl_bmi_withoutNRJtable <- read_fst("./Exposures/bmi_table.fst", as.data.table = TRUE)
tbl_bmi_withoutNRJtable <- tbl_bmi_withoutNRJtable[year>2000 & year<2045,]
tbl_bmi_withoutNRJtable <- tbl_bmi_withoutNRJtable[, `:=` (year=year-2000,
                                                           sex=as.factor(fifelse(sex=="male","1","2")),
                                                           SES=fifelse(edu_cat=="low","1",
                                                                      fifelse(edu_cat=="middle","2","3")))]
qread("./Exposures/bmi_model.qs")$family
#TO DO
#need to change the distribution name within the code cerca line 135 (qBCPEo) in this case


#GAMLSS nrj
##FOR GERMANY WE DONT HAVE THE ENERGY FROM OOH, so we are importing energy and we will have a ratio
tbl_Energy_out_kcal_table <- read_fst("./Exposures/energy_table.fst", as.data.table = TRUE)
tbl_Energy_out_kcal_table <- tbl_Energy_out_kcal_table[year>2000 & year<2045,]
tbl_Energy_out_kcal_table <- tbl_Energy_out_kcal_table[, `:=` (
  year = year - 2000,
  sex = as.factor(fifelse(sex ==
                            "male", "1", "2")),
  SES = fifelse(edu_cat ==
                  "low", "1", fifelse(edu_cat == "middle", "2", "3"))
)]
qread("./Exposures/energy_model.qs")$family
#TO DO
#need to change the distribution name within the code cerca line 166 (qexGAUS in this case)

############################################################################
############ NEED TO CHANGE ASAP
#for the country already having the fraction included, fraction_OOH=1
#based on Orfanos 2009 - Eating out of home: energy, macro- and micronutrient intakes in 10 European countries. The European Prospective Investigation into Cancer and Nutrition
#the EPIC study give the proportion of energy from out of home in 2 center in germany in 2000s
#              Heindelberg   Potsdam      Heindelberg      Potsdam
#                  Men         Men           Women          Women
# 35-44y          33.9        30.7           24.5           30.0
# 45-54y          33.0        26.1           23.8           24.3
# 55-64y          24.6        17.6           17.8           15.0
# 65-74y           NA         13.4            NA            9.8
#for now, waiting for better, we are taking the mean between both center
fraction_OOH <- CJ(sex = c("1","2"), age=min_age:max_age)
fraction_OOH[, `:=` (fraction_OOH = fifelse(
  sex == 1 & age < 45,
  0.323,
  fifelse(
    sex == 1 & age >= 45 & age < 55,
    0.2955,
    fifelse(
      sex == 1 & age >= 55 & age < 65,
      0.211,
      fifelse(
        sex == 1 & age >= 65,
        0.134,
        fifelse(
          sex == 2 & age < 45,
          0.2725,
          fifelse(
            sex == 2 & age >= 45 & age < 55,
            0.2405,
            fifelse(
              sex == 2 & age >= 55 & age < 65,
              0.164,
              fifelse(sex ==
                        2 & age >= 65, 0.098, NA_real_)
            )
          )
        )
      )
    )
  )
))] 
# summary(fraction_OOH$fraction_OOH)




#GAMLSS SSB
tbl_SSB_round_table <- read_fst("./Exposures/ssb_table.fst", as.data.table = TRUE)
tbl_SSB_round_table <- tbl_SSB_round_table[year>2000 & year<2045,]
tbl_SSB_round_table <- tbl_SSB_round_table[, `:=` (
  year = year - 2000,
  sex = as.factor(fifelse(sex ==
                            "male", "1", "2")),
  SES = fifelse(edu_cat ==
                  "low", "1", fifelse(edu_cat == "middle", "2", "3"))
)]
qread("./Exposures/ssb_model.qs")$family
#TO DO
#need to change the distribution name within the code cerca line 200 (qGG)

###It's SSB intakes, including the diet SSB
#So we are predicting the % of SSB being diet to be able to estimate the SSB non-diet
tbl_SSB_diet_prop_table <- read_fst("./Exposures/ssb_diet_prop.fst", as.data.table = TRUE)
tbl_SSB_diet_prop_table <- tbl_SSB_diet_prop_table[, `:=` (sex=as.factor(fifelse(sex=="men","1","2")))]
#tbl_SSB_diet_prop_table[,table(sex)]


####MIGHT CHANGE IS YOU HAVE CHANGING THE NUMBER OF SCENARIOS
###In this model, we have for now 35 scenarios (Menu: 22; SSB: 13)
# 35 scenarios, 3 SES, 3 diseases (chd,haem,isch)
# 2 baseline, 3 SES, 3 diseases
# obesity prevalence: 35 scenarios, 3 SES, 1 disease 
# 1 baseline, 3 SES, 1 disease
poss_tot <- (35*3*3) + (35*3*1) + (2*3*3) + (1*3*1)





##START OF THE MODEL
#######################################MONTECARLO
#Number of iterations (for testing our model, put 1:2)

mat3 <- FBM(montecarlo, (poss_tot*sim_hor))
cl <- parallel::makeCluster(cluster_number, outfile="") ##HERE you might want to change the number of cluster // change the number of cores
#detectCores(all.tests = FALSE, logical = TRUE)
doParallel::registerDoParallel(cl)
time.begin <- Sys.time()


system.time(
  tmp <- foreach(sim = 1:montecarlo, .combine = 'rbind', .packages=c("data.table","fst","gamlss")) %dopar% {
    set.seed(sim)
    setDTthreads(1L, restore_after_fork = FALSE)
    threads_fst(1L, reset_after_fork = FALSE)
    
    time.old <- Sys.time()
    
    # Import ONS population estimates 
    onspop <- onspop_import
    
    sp <- onspop[between(as.integer(age), min_age, max_age), 
                 .(year = init_year, age = as.integer(age), sex = factor(sex), SES = factor(SES), popsize = `2022`)]
    head(sp)
    # We scale down the population by the scale factor. We will scale up by factor the results (for now ignoring the truncation error)
    sp[, popsize := round(popsize/scale_factor)]
    nrow(sp)
    head(sp)
    
    # Let's create the simulants (person-year table)
    sp <- sp[rep(1:.N, times = popsize), .(pid = .I, year, age, sex, SES)] # rep(1:3, times = 1:3)
    nrow(sp)
    head(sp)
    
    # Now lets project their life course for the next 20 years (simulation horizon)
    sp <- rbindlist(rep(list(sp), sim_hor), idcol = "time")
    sp[, `:=` (year = year + time - 1L, 
               age = age + time - 1L,
               time = NULL)]
    sp[, table(age, year)]   # this is the close cohort
    sp <- sp[age <= max_age] # Some pruning
    setkey(sp, pid, year)  
    
    # The structure above is very useful in dynamic microsimulation, each row is a person year
    
    # We will create an open cohort 
    # With our previous approach, we have no one aged 20 in 2021 or 20-21 in 22
    # We need a new cohort of 20yo simulants enter the model every year
    # Let's get back to the ONS pop size estimates file
    onspop <- melt(onspop, id.vars = c("age", "sex","SES"),
                   variable.name = "year", value.name = "popsize")
    setDT(onspop)
    onspop[, lapply(.SD, class)]
    onspop[, year := as.integer(as.character(year)) - 2000L]
    onspop[, popsize := as.integer(as.character(popsize))]
    
    # Let's consider the cohort of 30yo entering the model in years post 2022
    sp2 <- onspop[as.integer(age) == min_age & 
                    between(year, init_year + 1L, init_year + sim_hor - 1L),
                  .(year, age = as.integer(age), sex = sex, SES=factor(SES), popsize = round(popsize / scale_factor))] #was 1000 before
    
    # Form the person-year table
    sp2 <- sp2[rep(1:.N, times = popsize), .(pid = .I + max(sp$pid), year, age, sex, SES)]
    
    # Create their lifecourse as with the original cohort
    sp2 <- rbindlist(rep(list(sp2), sim_hor), idcol = "time")
    sp2[, `:=` (year = year + time - 1L, 
                age = age + time - 1L,
                time = NULL)]
    setkey(sp2, pid, year) 
    
    # some pruning is necessary
    sp2[, table(age, year)]
    sp2 <- sp2[year < init_year + sim_hor ]
    sp2[, table(age, year)]
    
    # And we bind the 2 cohorts
    sp <- rbind(sp, sp2)
    setkey(sp, pid, year) 
    sp[, table(age, year)]
    
    # Now I can start adding attributes for my simulants. The order we add the parameters
    # reflects the causality network assumptions. 
    # Here I just estimate the bmi I will make the rank stability assumption
    tt <- sp[, .("pid" = unique(pid))]
    #set.seed(42L)
    tt[, rn_bmi := runif(.N) * 0.995] #test here to find good number to avoid long tail
    sp[tt, on = "pid", rn_bmi := rn_bmi]
    head(sp)
    tbl <- tbl_bmi_withoutNRJtable
    head(tbl)
    sp[tbl, on = c("year","age","sex","SES"), `:=` (mu = i.mu, sigma = i.sigma, nu = i.nu, tau = i.tau)]
    head(sp)
    sp[, bmi := qBCPEo(rn_bmi, mu, sigma, nu, tau)] #need to use the quantile function here! qXXX. Also CHANGE the distribution accroding to the GAMLSS fitting
    # plot(density(sp$bmi))
    summary(sp$bmi) #14 to 50, mean at 28, seems high CHRIS
    sp[, c("rn_bmi", "mu", "sigma", "nu", "tau") := NULL]
    head(sp)
    #BMI is the best way of doing that, instead of height and weight as modelling weight with height is super complicated
    #Here I need to do a hard boundary as I have a specific range for Energy
    sp[bmi > 55, bmi := 55]
    sp[bmi < 13, bmi := 13]
    
    sp[, `:=` (bmi_grp="<20")]
    sp[(bmi>=20 & bmi<25), `:=` (bmi_grp="20-25")]
    sp[(bmi>=25 & bmi<30), `:=` (bmi_grp="25-30")]
    sp[(bmi>=30), `:=` (bmi_grp=">=30")]
    head(sp)
    sp[year==22, prop.table(table(bmi_grp))] 
    # bmi_grp / good, in average, 19% of German pop is obese (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9520353/)
    # <20             >=30      20-25      25-30 
    #0.03928694 0.18613503 0.35180847 0.42276955  #edited by Edi
    
    #having energy intakes, also predicted by the bmi
    tt <- sp[, .("pid" = unique(pid))]
    tt[, rn_nrj := runif(.N) * 0.997] #test here to find good number to avoid long tail
    sp[tt, on = "pid", rn_nrj := rn_nrj]
    #head(sp)
    tbl <- tbl_Energy_out_kcal_table
    tbl[,`:=` (bmi2=round(bmi,1))]
    #table(tbl$bmi2)
    #head(tbl)
    sp[,`:=` (bmi2=round(bmi,1))] 
    sp[tbl, on = c("year","age","sex","SES","bmi2"), `:=` (mu = i.mu, sigma = i.sigma, nu = i.nu)]
    #head(sp)
    sp[, Energy_out_kcal := qexGAUS(rn_nrj, mu, sigma, nu)] #need to use the quantile function here! qXXX. Also CHANGE the distribution according to the GAMLSS fitting
    #plot(density(sp$Energy_out_kcal))
    #summary(sp$Energy_out_kcal) 
    sp[, c("rn_nrj", "mu", "sigma", "nu") := NULL]
    #head(sp)
    
    
    ###now, we are estimating the energy from out of home
    #for the country already having the fraction included, fraction_OOH=1
    sp[fraction_OOH, on=c("sex","age"), fraction_OOH := i.fraction_OOH]
    sp[, Energy_out_kcal:=Energy_out_kcal*fraction_OOH] 
    
    
    #having SSB intakes, also predicted by the bmi
    tt <- sp[, .("pid" = unique(pid))]
    tt[, rn_ssb := runif(.N) * 0.997] #test here to find good number to avoid long tail
    sp[tt, on = "pid", rn_ssb := rn_ssb]
    #head(sp)
    tbl <- tbl_SSB_round_table
    tbl[,`:=` (bmi2=round(bmi,1))]
    #head(tbl)
    #sp[,`:=` (sex=factor(sex))]
    ###ASK CHRIS IF PUTTING BMI AS AN INTEGER IS THE BEST THING TO DO
    #no, we should have integer just for the matching, so create another variable with integer that we delete after the merging
    sp[tbl, on = c("year","age","sex","SES","bmi2"), `:=` (mu = i.mu, sigma = i.sigma, nu = i.nu)]
    #head(sp)
    sp[, SSB_all := qGG(rn_ssb, mu, sigma, nu)] #need to use the quantile function here! qXXX. Also CHANGE the distribution accroding to the GAMLSS fitting
    #plot(density(sp$SSB))
    #summary(sp$SSB) 
    sp[, c("rn_ssb", "mu", "sigma", "nu", "bmi2") := NULL]
    #head(sp)
    
    #having % of SSB being diet
    tbl <- tbl_SSB_diet_prop_table
    sp[tbl, on = c("age","sex"), `:=` (diet_prop=i.diet_prop)]
    sp[, SSB_diet:=SSB_all*diet_prop]
    sp[, SSB:=SSB_all-SSB_diet]
    
    
    
    # Now I will model mortality 
    #Artilce ERFC 2011 - Separate and combined associations of body-mass index and abdominal adiposity with cardiovascular disease: collaborative analysis of 58 prospective studies
    # https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(11)60105-0/fulltext#secd19223949e1020
    # Dans le doc https://els-jbs-prod-cdn.jbs.elsevierhealth.com/cms/attachment/6d6bb978-0d17-4017-9845-2903d3b68778/gr2.gif
    # In people with BMI of 20 kg/m² or higher: HRs per 1 SD higher baseline values : 4·56 kg/m² higher BMI
    # CHD
    # 40-59 years : 1.41 (1.30-1.53)
    # 60-69 years : 1.23 (1.15-1.31)
    # 70+ years   : 1.12 (1.05-1.19)
    # 
    # Ischaemic stroke
    # 40-59 years : 1.34 (1.21-1.48)
    # 60-69 years : 1.22 (1.13-1.31)
    # 70+ years   : 1.08 (0.99-1.18)
    
    sp[, `:=` (rr.obesity.chd_e =1, rr.obesity.chd_l =1, rr.obesity.chd_u= 1,
               rr.obesity.isch_e=1, rr.obesity.isch_l=1, rr.obesity.isch_u=1,
               rr.obesity.haem_e=1, rr.obesity.haem_l=1, rr.obesity.haem_u=1)]
    sp[(age>=40 & age<60), `:=` (rr.obesity.chd_e= 1.41, rr.obesity.chd_l= 1.30, rr.obesity.chd_u= 1.53,
                                 rr.obesity.isch_e=1.34, rr.obesity.isch_l=1.21, rr.obesity.isch_u=1.48,
                                 rr.obesity.haem_e=1,    rr.obesity.haem_l=1,    rr.obesity.haem_u=1)]
    sp[(age>=60 & age<70), `:=` (rr.obesity.chd_e= 1.23, rr.obesity.chd_l= 1.15, rr.obesity.chd_u= 1.31,
                                 rr.obesity.isch_e=1.22, rr.obesity.isch_l=1.13, rr.obesity.isch_u=1.31,
                                 rr.obesity.haem_e=1,    rr.obesity.haem_l=1,    rr.obesity.haem_u=1)]
    sp[(age>=70), `:=` (rr.obesity.chd_e= 1.12, rr.obesity.chd_l= 1.05, rr.obesity.chd_u= 1.19,
                        rr.obesity.isch_e=1.08, rr.obesity.isch_l=0.99, rr.obesity.isch_u=1.18,
                        rr.obesity.haem_e=1,    rr.obesity.haem_l=1,    rr.obesity.haem_u=1)]
    sp[, table(rr.obesity.chd_e)]
    #sp[, table(rr.obesity.chd_e,age)]
    sp[, table(rr.obesity.isch_e)]
    sp[, `:=` (rr.obesity.chd_se=((log(rr.obesity.chd_u) - log(rr.obesity.chd_e))/1.96 + (log(rr.obesity.chd_l) - log(rr.obesity.chd_e))/-1.96)/2,
               rr.obesity.isch_se=((log(rr.obesity.isch_u) - log(rr.obesity.isch_e))/1.96 + (log(rr.obesity.isch_l) - log(rr.obesity.isch_e))/-1.96)/2,
               rr.obesity.haem_se=((log(rr.obesity.haem_u) - log(rr.obesity.haem_e))/1.96 + (log(rr.obesity.haem_l) - log(rr.obesity.haem_e))/-1.96)/2)]
    head(sp)
    
    
    ########################################################################################################################################################################
    ########################################################################################################################################################################
    ########################################################################################################################################################################
    ##NEED TO BE CHANGE ACCORDING TO DECISION (and make sure its mortality and not risk on incidence / prevalence)
    # ALSO NEED TO SEE IF CHANGE FOR SEX, AGE, BMI AND SES
    #FOR NOW SAME RR FOR AGE SEX BMI AND SES
    ########################################################################################################################################################################
    ########################################################################################################################################################################
    ########################################################################################################################################################################
    # CVD Mortality RR for SSB intake (Santos et al (2022) https://doi.org/10.1016/j.clnesp.2022.08.021
    ### SSB intake increased the risk of coronary heart disease (RR = 1.15; 9% C.I. 1.06–1.25), and stroke (RR = 1.10; 9% C.I. 1.01–1.19) 
    #in adults after adjustment for all potential confounders
    # sp[, `:=` (rr.SSB.sugar.chd_e =1.15, rr.SSB.sugar.chd_l =1.06, rr.SSB.sugar.chd_u= 1.25,
    #            rr.SSB.sugar.isch_e=1.10, rr.SSB.sugar.isch_l=1.01, rr.SSB.sugar.isch_u=1.19,   #stroke=isch+haem? so risk for isch is the same that for stroke? same for haem?
    #            rr.SSB.sugar.haem_e=1.10, rr.SSB.sugar.haem_l=1.01, rr.SSB.sugar.haem_u=1.19)]  #stroke=isch+haem? so risk for isch is the same that for stroke? same for haem?
    #sp[, table(rr.SSB.sugar.chd_e)]
    #sp[, `:=` (rr.SSB.sugar.chd_se =((log(rr.SSB.sugar.chd_u) - log(rr.SSB.sugar.chd_e))/1.96 + (log(rr.SSB.sugar.chd_l) - log(rr.SSB.sugar.chd_e))/-1.96)/2,
    #           rr.SSB.sugar.isch_se=((log(rr.SSB.sugar.isch_u) - log(rr.SSB.sugar.isch_e))/1.96 + (log(rr.SSB.sugar.isch_l) - log(rr.SSB.sugar.isch_e))/-1.96)/2,
    #           rr.SSB.sugar.haem_se=((log(rr.SSB.sugar.haem_u) - log(rr.SSB.sugar.haem_e))/1.96 + (log(rr.SSB.sugar.haem_l) - log(rr.SSB.sugar.haem_e))/-1.96)/2)]
    #head(sp)
    
    #Micha et al, etable 5 - 1serving of SSB (8oz or 227mL) per day, increase the chd risk, BMI adjusted differently by age
    # 25-34 years	: 1.33 (1.19; 1.47)	
    # 35-44 years : 1.31 (1.18; 1.45)	
    # 45-54 years : 1.26 (1.15; 1.37)	
    # 55-64 years : 1.21 (1.13; 1.30)
    # 65-74 years	: 1.17 (1.10; 1.24)
    # 75+ years   : 1.09 (1.06; 1.13)
    sp[, `:=` (rr.SSB.intake.BMIadj.chd_e =1.33, rr.SSB.intake.BMIadj.chd_l=1.19, rr.SSB.intake.BMIadj.chd_u=1.47,
               rr.SSB.intake.BMIadj.isch_e=1, rr.SSB.intake.BMIadj.isch_l=1, rr.SSB.intake.BMIadj.isch_u=1,   
               rr.SSB.intake.BMIadj.haem_e=1, rr.SSB.intake.BMIadj.haem_l=1, rr.SSB.intake.BMIadj.haem_u=1)]
    sp[(age>=35 & age<45), `:=` (rr.SSB.intake.BMIadj.chd_e =1.31, rr.SSB.intake.BMIadj.chd_l =1.18, rr.SSB.intake.BMIadj.chd_u= 1.45,
                                 rr.SSB.intake.BMIadj.isch_e=1,    rr.SSB.intake.BMIadj.isch_l=1,    rr.SSB.intake.BMIadj.isch_u=1,   
                                 rr.SSB.intake.BMIadj.haem_e=1,    rr.SSB.intake.BMIadj.haem_l=1,    rr.SSB.intake.BMIadj.haem_u=1)]
    sp[(age>=45 & age<55), `:=` (rr.SSB.intake.BMIadj.chd_e =1.26, rr.SSB.intake.BMIadj.chd_l =1.15, rr.SSB.intake.BMIadj.chd_u= 1.37,
                                 rr.SSB.intake.BMIadj.isch_e=1,    rr.SSB.intake.BMIadj.isch_l=1,    rr.SSB.intake.BMIadj.isch_u=1,   
                                 rr.SSB.intake.BMIadj.haem_e=1,    rr.SSB.intake.BMIadj.haem_l=1,    rr.SSB.intake.BMIadj.haem_u=1)]
    sp[(age>=55 & age<65), `:=` (rr.SSB.intake.BMIadj.chd_e =1.21, rr.SSB.intake.BMIadj.chd_l =1.13, rr.SSB.intake.BMIadj.chd_u= 1.30,
                                 rr.SSB.intake.BMIadj.isch_e=1,    rr.SSB.intake.BMIadj.isch_l=1,    rr.SSB.intake.BMIadj.isch_u=1,   
                                 rr.SSB.intake.BMIadj.haem_e=1,    rr.SSB.intake.BMIadj.haem_l=1,    rr.SSB.intake.BMIadj.haem_u=1)]
    sp[(age>=65 & age<75), `:=` (rr.SSB.intake.BMIadj.chd_e =1.17, rr.SSB.intake.BMIadj.chd_l =1.10, rr.SSB.intake.BMIadj.chd_u= 1.24,
                                 rr.SSB.intake.BMIadj.isch_e=1,    rr.SSB.intake.BMIadj.isch_l=1,    rr.SSB.intake.BMIadj.isch_u=1,   
                                 rr.SSB.intake.BMIadj.haem_e=1,    rr.SSB.intake.BMIadj.haem_l=1,    rr.SSB.intake.BMIadj.haem_u=1)]
    sp[(age>=75), `:=` (rr.SSB.intake.BMIadj.chd_e =1.09, rr.SSB.intake.BMIadj.chd_l =1.06, rr.SSB.intake.BMIadj.chd_u= 1.13,
                        rr.SSB.intake.BMIadj.isch_e=1, rr.SSB.intake.BMIadj.isch_l=1, rr.SSB.intake.BMIadj.isch_u=1,   
                        rr.SSB.intake.BMIadj.haem_e=1, rr.SSB.intake.BMIadj.haem_l=1, rr.SSB.intake.BMIadj.haem_u=1)]
    sp[, table(rr.SSB.intake.BMIadj.chd_e)]
    sp[, `:=` (rr.SSB.intake.BMIadj.chd_se=((log(rr.SSB.intake.BMIadj.chd_u) - log(rr.SSB.intake.BMIadj.chd_e))/1.96 + (log(rr.SSB.intake.BMIadj.chd_l) - log(rr.SSB.intake.BMIadj.chd_e))/-1.96)/2,
               rr.SSB.intake.BMIadj.isch_se=((log(rr.SSB.intake.BMIadj.isch_u) - log(rr.SSB.intake.BMIadj.isch_e))/1.96 + (log(rr.SSB.intake.BMIadj.isch_l) - log(rr.SSB.intake.BMIadj.isch_e))/-1.96)/2,
               rr.SSB.intake.BMIadj.haem_se=((log(rr.SSB.intake.BMIadj.haem_u) - log(rr.SSB.intake.BMIadj.haem_e))/1.96 + (log(rr.SSB.intake.BMIadj.haem_l) - log(rr.SSB.intake.BMIadj.haem_e))/-1.96)/2)]
    head(sp)
    
    #as sometime we don't have mortality by education, we are adding a ratio (define at the top of the code)
    sp <- sp[rr.education.table, on=c("SES","sex")]
    
    
    # Simulation: RR 
    k <- runif(1)
    sp[, `:=` (
      rr.obesity.chd = qlnorm(
        k,
        meanlog = log(rr.obesity.chd_e),
        sdlog = rr.obesity.chd_se
      ),
      rr.obesity.isch = qlnorm(
        k,
        meanlog = log(rr.obesity.isch_e),
        sdlog = rr.obesity.isch_se
      ),
      rr.obesity.haem = qlnorm(
        k,
        meanlog = log(rr.obesity.haem_e),
        sdlog = rr.obesity.haem_se
      )
    )]   
    sp[rr.obesity.chd < 1, rr.obesity.chd := 1]
    sp[rr.obesity.isch < 1, rr.obesity.isch := 1]
    sp[rr.obesity.haem < 1, rr.obesity.haem := 1]
    
   # k <- runif(1)
   #  sp[, `:=` (
   #    rr.SSB.sugar.chd = qlnorm(
   #      k,
   #      meanlog = log(rr.SSB.sugar.chd_e),
   #      sdlog = rr.SSB.sugar.chd_se
   #   ),
   #    rr.SSB.sugar.isch = qlnorm(
   #      k,
   #      meanlog = log(rr.SSB.sugar.isch_e),
   #      sdlog = rr.SSB.sugar.isch_se
   #    ),
   #    rr.SSB.sugar.haem = qlnorm(
   #      k,
   #      meanlog = log(rr.SSB.sugar.haem_e),
   #      sdlog = rr.SSB.sugar.haem_se
   #    )
   #  )]
   #  sp[rr.SSB.sugar.chd < 1,  rr.SSB.sugar.chd := 1]
   #  sp[rr.SSB.sugar.isch < 1, rr.SSB.sugar.isch := 1]
   #  sp[rr.SSB.sugar.haem < 1, rr.SSB.sugar.haem := 1]
    
    k <- runif(1)
    sp[, `:=` (
      rr.SSB.intake.BMIadj.chd = qlnorm(
        k,
        meanlog = log(rr.SSB.intake.BMIadj.chd_e),
        sdlog = rr.SSB.intake.BMIadj.chd_se
      ),
      rr.SSB.intake.BMIadj.isch = qlnorm(
        k,
        meanlog = log(rr.SSB.intake.BMIadj.isch_e),
        sdlog = rr.SSB.intake.BMIadj.isch_se
      ),
      rr.SSB.intake.BMIadj.haem = qlnorm(
        k,
        meanlog = log(rr.SSB.intake.BMIadj.haem_e),
        sdlog = rr.SSB.intake.BMIadj.haem_se
      )
    )]
    sp[rr.SSB.intake.BMIadj.chd < 1,  rr.SSB.intake.BMIadj.chd := 1]
    sp[rr.SSB.intake.BMIadj.isch < 1, rr.SSB.intake.BMIadj.isch := 1]
    sp[rr.SSB.intake.BMIadj.haem < 1, rr.SSB.intake.BMIadj.haem := 1]
    
    k <- runif(1)
    sp[, `:=` (rr.education = qlnorm(k,meanlog = log(rr.education_e), sdlog = rr.education_se))]
    sp[rr.education < 1, rr.education := 1]          
               
    
    sp[, c("rr.obesity.chd_e", "rr.obesity.chd_l", "rr.obesity.chd_u", "rr.obesity.chd_se",
           "rr.obesity.isch_e", "rr.obesity.isch_l", "rr.obesity.isch_u", "rr.obesity.isch_se",
           "rr.obesity.haem_e", "rr.obesity.haem_l","rr.obesity.haem_u", "rr.obesity.haem_se",
           "rr.SSB.intake.BMIadj.chd_e", "rr.SSB.intake.BMIadj.chd_l", "rr.SSB.intake.BMIadj.chd_u", "rr.SSB.intake.BMIadj.chd_se",
           "rr.SSB.intake.BMIadj.isch_e", "rr.SSB.intake.BMIadj.isch_l", "rr.SSB.intake.BMIadj.isch_u", "rr.SSB.intake.BMIadj.isch_se",
           "rr.SSB.intake.BMIadj.haem_e", "rr.SSB.intake.BMIadj.haem_l","rr.SSB.intake.BMIadj.haem_u", "rr.SSB.intake.BMIadj.haem_se",
           "rr.education_e", "rr.education_l","rr.education_u", "rr.education_se") := NULL]
    #View(sp)
    
    
    #We are assuming that we have a 5y lag time between the bmi and the RR see https://doi.org/10.1136/bmj.i2793 
    #to be easier here, we are creating a bmi_lagged to obtain the rr 
    #Calculating the RR according to BMI
    setkey(sp, pid, year) # NECESSARY!!!
    sp[, bmi_lagged := shift(bmi, 5), by = pid]
    #check <- sp[,.(pid,year,bmi,bmi_lagged)]
    
    sp[, `:=` (rr.obesity.chd.indiv=rr.obesity.chd^((bmi_lagged-20)/4.56),  #BMIss bring back to 4.56 to be comparative with rr.Obesity
               rr.obesity.isch.indiv=rr.obesity.isch^((bmi_lagged-20)/4.56),  
               rr.obesity.haem.indiv=rr.obesity.haem^((bmi_lagged-20)/4.56))]
    sp[bmi_lagged<20, `:=` (rr.obesity.chd.indiv=1,  
                            rr.obesity.isch.indiv=1,  
                            rr.obesity.haem.indiv=1)]
    
    
    ##################DO YOU WANT A LAG TIME FOR RR SSB BMIADJ? NEED TO MAKE THE DECISION AND PUT THE LAG TIME HERE IF NEEDED (NEED TO SEARCH IN LITERATURE AND DISCUSS WITH CHRIS AND MARTIN)
    #do we need a lag time for SSB consumption on risk????? IF YES, NEED TO HAVE A SSB_lagged
    setkey(sp, pid, year) #(added by Edi)
    sp[, SSB_lagged:= shift(SSB, 5), by = pid] 
    sp[, `:=` (rr.SSB.intake.BMIadj.chd.indiv=rr.SSB.intake.BMIadj.chd^((SSB_lagged)/227.3045),  #227.3045 as the risk is for 8oz or 20 grams of sugar (risk increase for each 8oz or 20 grams of sugar)
               rr.SSB.intake.BMIadj.isch.indiv=rr.SSB.intake.BMIadj.isch^((SSB_lagged)/227.3045),  #it is the same when this SSB in ml transformed into sugar --> SSB*20/227.3045, and risk calculated based on sugar --> (SSB*20/227.3045)/20 --> SSB/227.3045
               rr.SSB.intake.BMIadj.haem.indiv=rr.SSB.intake.BMIadj.haem^((SSB_lagged)/227.3045))]
    sp[SSB_lagged<227.3045, `:=` (rr.SSB.intake.BMIadj.chd.indiv=1,  
                                  rr.SSB.intake.BMIadj.isch.indiv=1,  
                                  rr.SSB.intake.BMIadj.haem.indiv=1)]
    
    ###################CHIRS/KARL NEED TO CHECK. FOR NOW, AS WE DONT HAVE MORTALITY BY EDUCATION IN GERMANY, 
    #from table 2 here, we know the age-adjusted mortality rate ratios for education in Germany  (https://bmjopen.bmj.com/content/9/10/e028001)
    # so, for now, I am calculating the risk for chd for each individual based on their bmi (with the lag) as above 
    #and to add the mortality differences by education I multiply this risk by the mortality rate ratios for education form the paper, and then I do the parf.
    ##DOES THAT IS OK????
    # CKNote: this is likely wrong
    
    if (SSB_direct_effect) {
    sp[, `:=` (
      rr.obesity.chd.indiv.educ = rr.obesity.chd.indiv * rr.education * rr.SSB.intake.BMIadj.chd.indiv,
      rr.obesity.isch.indiv.educ = rr.obesity.isch.indiv * rr.education * rr.SSB.intake.BMIadj.isch.indiv,
      rr.obesity.haem.indiv.educ = rr.obesity.haem.indiv * rr.education * rr.SSB.intake.BMIadj.haem.indiv
    )]
    } else {
      sp[, `:=` (
        rr.obesity.chd.indiv.educ = rr.obesity.chd.indiv * rr.education,
        rr.obesity.isch.indiv.educ = rr.obesity.isch.indiv * rr.education,
        rr.obesity.haem.indiv.educ = rr.obesity.haem.indiv * rr.education
      )]
    }
    
    # adding RR =1 for those with NA for BMI/SSB lag as they are within the lag time (have no risk) #added by Edi
    sp[is.na(bmi_lagged), rr.obesity.chd.indiv.educ := 1]
    sp[is.na(bmi_lagged), rr.obesity.isch.indiv.educ := 1]
    sp[is.na(bmi_lagged), rr.obesity.haem.indiv.educ := 1]
               

    mrtlparf1 <-
      sp[, .(parf_obesity_chd = 1 - 1 / (sum(rr.obesity.chd.indiv.educ) / .N )), keyby = .(year, age, sex, SES)] #dont need bmi grp here/ if strong trend no year to avoid double counting, for us Chris say keep year
    head(mrtlparf1)
    mrtlparf2 <-
      sp[, .(parf_obesity_isch = 1 - 1 / (sum(rr.obesity.isch.indiv.educ) / .N )), keyby = .(year, age, sex, SES)]
    mrtlparf3 <-
      sp[, .(parf_obesity_haem = 1 - 1 / (sum(rr.obesity.haem.indiv.educ) / .N )), keyby = .(year, age, sex, SES)]
    
    mrtlparf <- mrtlparf1[mrtlparf2, on = .NATURAL,]
    mrtlparf <- mrtlparf[mrtlparf3, on = .NATURAL,]
    head(mrtlparf) #parf=0 until year 2027 in line with lag time
    rm(mrtlparf1,mrtlparf2,mrtlparf3) 
    
    
    # Get mortality rates from ons
    onsmort <- onsmort_import
    
    # mx = central rate of mortality for 1
    k=runif(1)
    onsmort <- setDT(onsmort)[, `:=` (chd.rate=qlnorm(k,meanlog = log(mx_chd), sdlog = mx_chd_se),
                                      isch.rate=qlnorm(k,meanlog = log(mx_isch), sdlog = mx_isch_se),
                                      haem.rate=qlnorm(k,meanlog = log(mx_haem), sdlog = mx_haem_se))]
    
    
    onsmort <- onsmort[,c("year","age","SES","sex","chd.rate","isch.rate","haem.rate")]
    
    yeye <- onsmort[between(as.integer(age), min_age, max_age) &
                      between(year, init_year, init_year + sim_hor - 1L), .(
                        year,
                        age = as.integer(age),
                        sex = factor(sex),
                        SES = factor(SES),
                        #here you should divide by 100,000 if the rate is for 100,000. mr = mr/1e5
                        chd.rate,
                        isch.rate,
                        haem.rate
                      )]
    
    
    mrtlparf[yeye, on = c("year", "age", "sex", "SES"), `:=` (chd.rate = i.chd.rate,
                                                              isch.rate = i.isch.rate,
                                                              haem.rate = i.haem.rate)]
    
    #mortality rate for unexposed bmi  
    mrtlparf[, `:=` (p0_bmi_mrtl_chd  = chd.rate  * (1 - parf_obesity_chd),# p0 is the mortality rate of the unexposed (bmi<20)
                     p0_bmi_mrtl_isch = isch.rate * (1 - parf_obesity_isch), 
                     p0_bmi_mrtl_haem = haem.rate * (1 - parf_obesity_haem))] 
    head(mrtlparf)
    sp[mrtlparf, on = c("year","age","sex","SES"), `:=` (p0_bmi_mrtl_chd = i.p0_bmi_mrtl_chd, p0_bmi_mrtl_isch = i.p0_bmi_mrtl_isch, p0_bmi_mrtl_haem = i.p0_bmi_mrtl_haem)]
    

    #we are estimating the mortality
    sp[, `:=` (mr_chd  = p0_bmi_mrtl_chd  * rr.obesity.chd.indiv.educ,
               mr_isch = p0_bmi_mrtl_isch * rr.obesity.isch.indiv.educ,
               mr_haem = p0_bmi_mrtl_haem * rr.obesity.haem.indiv.educ)]
    #no need to do mr_chd_ssb  = p0_SSB_intake_mrtl_chd  * rr.SSB.intake.BMIadj.chd.indiv 
    #as they should give the same results (see test below)
    
    # Let's test it, should be equal
    # sp[mrtlparf, on = c("year", "age",  "sex", "SES"), `:=` (mr_chd_ons = i.chd.rate, mr_isch_ons = i.isch.rate, mr_haem_ons = i.haem.rate)]
    # sp[, .("Our deaths" = sum(mr_chd, na.rm=T), "ONS deaths" = sum(mr_chd_ons, na.rm=T)), keyby = year] # our expected number of deaths by year vs the ONS one
    # sp[, .("Our deaths" = sum(mr_isch, na.rm=T), "ONS deaths" = sum(mr_isch_ons, na.rm=T)), keyby = year] # our expected number of deaths by year vs the ONS one
    # sp[, .("Our deaths" = sum(mr_haem, na.rm=T), "ONS deaths" = sum(mr_haem_ons, na.rm=T)), keyby = year] # our expected number of deaths by year vs the ONS one
    # sp[, .("Our deaths" = sum(mr_chd+mr_isch+mr_haem, na.rm=T), "ONS deaths" = sum(mr_chd_ons+mr_isch_ons+mr_haem_ons, na.rm=T)), keyby = year] # our expected number of deaths by year vs the ONS one
    # sp[, .("Our deaths" = sum(mr_chd+mr_isch+mr_haem, na.rm=T), "ONS deaths" = sum(mr_chd_ons+mr_isch_ons+mr_haem_ons, na.rm=T))] # our expected number of deaths by year vs the ONS one
    # mrtlparf[, .("Our deaths" = sum(chd.rate+isch.rate+haem.rate, na.rm=T)), keyby = year]
    # mrtlparf[, .("Our deaths" = sum(chd.rate+isch.rate+haem.rate, na.rm=T))] 
    
    # ###TEST to be sure that mr_chd obtain with bmi is the same as the one obtain with ssb
    # sp[, `:=` (mr_chd_ssb  = p0_SSB_intake_mrtl_chd  * rr.SSB.intake.BMIadj.chd.indiv.educ,
    #            mr_isch_ssb = p0_SSB_intake_mrtl_isch * rr.SSB.intake.BMIadj.isch.indiv.educ,
    #            mr_haem_ssb = p0_SSB_intake_mrtl_haem * rr.SSB.intake.BMIadj.haem.indiv.educ)]
    # #everything should be equal
    # sp[, .("Our deaths" = sum(mr_chd_ssb), "ONS deaths" = sum(mr_chd_ons)), keyby = year] # our expected number of deaths by year vs the ONS one
    # sp[, .("Our deaths" = sum(mr_isch_ssb), "ONS deaths" = sum(mr_isch_ons)), keyby = year] # our expected number of deaths by year vs the ONS one
    # sp[, .("Our deaths" = sum(mr_haem_ssb), "ONS deaths" = sum(mr_haem_ons)), keyby = year] # our expected number of deaths by year vs the ONS one
    # sp[, .("Our deaths" = sum(mr_chd_ssb+mr_isch_ssb+mr_haem_ssb), "ONS deaths" = sum(mr_chd_ons+mr_isch_ons+mr_haem_ons)), keyby = year] # our expected number of deaths by year vs the ONS one
    # sp[, .("Our deaths" = sum(mr_chd_ssb+mr_isch_ssb+mr_haem_ssb), "ONS deaths" = sum(mr_chd_ons+mr_isch_ons+mr_haem_ons))] # our expected number of deaths by year vs the ONS one
    # 
    # sp[, .("Our deaths SSB" = sum(mr_chd_ssb), "Our deaths BMI" = sum(mr_chd)), keyby = year] 
    # sp[, .("Our deaths SSB" = sum(mr_isch_ssb), "Our deaths BMI" = sum(mr_isch)), keyby = year] 
    # sp[, .("Our deaths SSB" = sum(mr_haem_ssb), "Our deaths BMI" = sum(mr_haem)), keyby = year] 
    # sp[, .("Our deaths SSB" = sum(mr_chd_ssb+mr_isch_ssb+mr_haem_ssb), "Our deaths BMI" = sum(mr_chd+mr_isch+mr_haem)), keyby = year] 
    # sp[, .("Our deaths SSB" = sum(mr_chd_ssb+mr_isch_ssb+mr_haem_ssb), "Our deaths BMI" = sum(mr_chd+mr_isch+mr_haem))] 
    # ##end of test
    
    
    #####################CHRIS NEED TO CHECK
    #for now, before CHRIS check, I consider that for the baseline, I can have the same expected/baseline number of death from both pathways (obesity and SSBBMIadj)
    #my concern CHRIS, is the competitive risk maybe
    
    # Let's see who dies
    sp[, rn_mrtl := runif(.N)] # Keep this as well
    pid_toremove <- integer(0) # a place holder
    for (i in sort(unique(sp$year))) {
      #print(paste0("Year: ", i))
      
      # rebalance the mortality rates to account for simulants removed from the
      # population when dead
      if (length(pid_toremove) > 0L) {
        correction_factor <- sp[year == i, sum(mr_chd)]/sp[year == i & !pid %in% pid_toremove, sum(mr_chd)] 
        sp[year == i, mr_chd := mr_chd * correction_factor] 
        #print(paste0("correction factor: ", correction_factor))
      }
      
      sp[year == i & !pid %in% pid_toremove, dead_chd := rn_mrtl < mr_chd] 
      pid_toremove <- sp[year <= i & (dead_chd), pid]
      #print(paste0("PID to remove: ", paste(pid_toremove, collapse = ",")))
    }
    sp[, longdead_chd := cumsum(dead_chd), by = pid]
    sp[dead_chd == FALSE & longdead_chd == 1, dead_chd := NA] # Best for later calculations.
    # So dead == FALSE mean alive, dead == TRUE means died that year,
    # and dead == NA means died in a previous year
    # dcast(sp, year~dead_chd) # Note NA some years is decreasing because simulants get removed from the dataset because their age is > 89
    # table(sp$(dead_chd)
    rm(correction_factor)
    
    for (i in sort(unique(sp$year))) {
      #print(paste0("Year: ", i))
      
      # rebalance the mortality rates to account for simulants removed from the
      # population when dead
      if (length(pid_toremove) > 0L) {
        correction_factor <- sp[year == i, sum(mr_isch)]/sp[year == i & !pid %in% pid_toremove, sum(mr_isch)] 
        sp[year == i, mr_isch := mr_isch * correction_factor]
        #print(paste0("correction factor: ", correction_factor))
      }
      
      sp[year == i & !pid %in% pid_toremove, dead_isch := rn_mrtl < mr_isch] 
      pid_toremove <- sp[year <= i & (dead_isch), pid]
      #print(paste0("PID to remove: ", paste(pid_toremove, collapse = ",")))
    }
    sp[, longdead_isch := cumsum(dead_isch), by = pid]
    sp[dead_isch == FALSE & longdead_isch == 1, dead_isch := NA] # Best for later calculations.
    # So dead == FALSE mean alive, dead == TRUE means died that year,
    # and dead == NA means died in a previous year
    #dcast(sp, year~dead_isch) # Note NA some years is decreasing because simulants get removed from the dataset because their age is > 89
    rm(correction_factor)
    
    for (i in sort(unique(sp$year))) {
      #print(paste0("Year: ", i))
      
      # rebalance the mortality rates to account for simulants removed from the
      # population when dead
      if (length(pid_toremove) > 0L) {
        correction_factor <- sp[year == i, sum(mr_haem)]/sp[year == i & !pid %in% pid_toremove, sum(mr_haem)] 
        sp[year == i, mr_haem := mr_haem * correction_factor] 
        #print(paste0("correction factor: ", correction_factor))
      }
      
      sp[year == i & !pid %in% pid_toremove, dead_haem := rn_mrtl < mr_haem] 
      pid_toremove <- sp[year <= i & (dead_haem), pid]
      #print(paste0("PID to remove: ", paste(pid_toremove, collapse = ",")))
    }
    sp[, longdead_haem := cumsum(dead_haem), by = pid]
    sp[dead_haem == FALSE & longdead_haem == 1, dead_haem := NA] # Best for later calculations.
    # So dead == FALSE mean alive, dead == TRUE means died that year,
    # and dead == NA means died in a previous year
    #dcast(sp, year~dead_haem) # Note NA some years is decreasing because simulants get removed from the dataset because their age is > 24
    rm(correction_factor)
    
    
    ##############################################################################
    #                      DEFINING ALL THE MODEL SCENARIOS                      #																				
    ##############################################################################
    # S0. no policy 
    
    ##############################################################################
    ##############################################################################
    ## MENU CALORIE LABELLING
    ##############################################################################
    ##############################################################################
    #Let's see the effect of kcal menu labelling
    #Energy kcal labelling decrease Energy (Crockett 2018): reduction of 47 kcal in energy purchased (MD:-46.72 kcal, 95% CI: -78.35, -15.10, N = 1877).
    #same effect in all SES, age, sex, BMI
    sp[, policy.energy.effect_sensi := rnorm(1,mean = -47, sd = (-15-(-78))/3.92)]
    sp[policy.energy.effect_sensi > 0, policy.energy.effect_sensi := 0]
    
    #Food labeling decreased consumer intakes of energy by 7.3% (95% CI:-10.1%, -4.4%) ; Shangguan 2019
    sp[, policy.energy.effect := rnorm(1,mean = -0.073, sd = (-0.044-(-0.101))/3.92)] #Updated on 2024.10.31 as the main effect
    sp[policy.energy.effect > 0, policy.energy.effect := 0]
    
    #Reformulation: Zlatevskaa
    #Email from Zlatevska saying that Estimate-15.34 [-23.1780;-7.5150]
    sp[, policy.energy.refeffect := rnorm(1,mean = -15, sd = (-7.5150-(-23.1780))/3.92)]
    sp[policy.energy.refeffect > 0, policy.energy.refeffect := 0]
    
    #Reformulation from the US studies https://doi.org/10.1136/bmjopen-2022-063614 , https://doi.org/10.1161/CIRCOUTCOMES.119.006313 #Updated on 2024.10.31
    sp[, policy.energy.refeffectUS := rnorm(1,mean = -0.05, sd = (0-0)/3.92)] #No CIs
    sp[policy.energy.refeffectUS > 0, policy.energy.refeffectUS := 0]
    
    sp[, policy.energy.refeffectUS.sensi := rnorm(1,mean = -0.05, sd = (-0.0425-(-0.0575))/3.92)] #No CIs
    sp[policy.energy.refeffectUS.sensi > 0, policy.energy.refeffectUS.sensi := 0]
    
    # S1.LABEL Energy current coverage 
    # Let's see what would happen if we put the policy in large food outlet (coverage of 0.07)x
    #with 26.5% compensation (we don't have info on SD)
    #assumption same compensation in all SES, age, sex, BMI
    sp[, sc1_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.large)) + ((Energy_out_kcal + ((Energy_out_kcal*policy.energy.effect) * (1-0.265))) * label.policy.coverage.large)]
    # S1.mincomp.LABEL Energy full coverage : people affected by the MenuLabelling - 11% compensation
    sp[, sc1mincomp_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.large)) + ((Energy_out_kcal + ((Energy_out_kcal*policy.energy.effect) * (1-0.11))) * label.policy.coverage.large)]
    # S1.maxcomp.LABEL Energy full coverage : people affected by the MenuLabelling - 42% compensation
    sp[, sc1maxcomp_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.large)) + ((Energy_out_kcal + ((Energy_out_kcal*policy.energy.effect) * (1-0.42))) * label.policy.coverage.large)]
    
    # S3.LABEL Energy full coverage with 26.5% compensation
    sp[, sc3_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.full)) + ((Energy_out_kcal + ((Energy_out_kcal*policy.energy.effect) * (1-0.265))) * label.policy.coverage.full)]
    # S3.mincomp.LABEL Energy full coverage : people affected by the MenuLabelling - 11% compensation
    sp[, sc3mincomp_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.full)) + ((Energy_out_kcal + ((Energy_out_kcal*policy.energy.effect) * (1-0.11))) * label.policy.coverage.full)]
    # S3.maxcomp.LABEL Energy full coverage : people affected by the MenuLabelling - 42% compensation
    sp[, sc3maxcomp_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.full)) + ((Energy_out_kcal + ((Energy_out_kcal*policy.energy.effect) * (1-0.42))) * label.policy.coverage.full)]
    
    # S2.LABEL Reformulation of the offer with no compensation = the retailers will decrease 5% each meal, current implementation
    sp[, sc2_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.large)) + ((Energy_out_kcal + ((Energy_out_kcal*policy.energy.refeffectUS) * (1-0))) * label.policy.coverage.large)]
    # S2.LABEL Reformulation of the offer with no compensation = the retailers will decrease 5% each meal with +- 15% margin errors, current implementation
    sp[, sc2sensi_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.large)) + ((Energy_out_kcal + ((Energy_out_kcal*policy.energy.refeffectUS.sensi) * (1-0))) * label.policy.coverage.large)]
    # S4.LABEL Reformulation of the offer with no compensation = the retailers will decrease 5% each meal, full implementation
    sp[, sc4_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.full)) + ((Energy_out_kcal + ((Energy_out_kcal*policy.energy.refeffectUS) * (1-0))) * label.policy.coverage.full)]
    # S4.LABEL Reformulation of the offer with no compensation = the retailers will decrease 5% each meal with +- 15% margin errors, full implementation
    sp[, sc4sensi_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.full)) + ((Energy_out_kcal + ((Energy_out_kcal*policy.energy.refeffectUS.sensi) * (1-0))) * label.policy.coverage.full)]
    
    
    #COMBINED
    # S5.LABEL Conso - Compensation + Reformulation, current coverage
    sp[, sc5_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.large)) + ((Energy_out_kcal + (((Energy_out_kcal*policy.energy.effect) * (1-0.265)) + (Energy_out_kcal*policy.energy.refeffectUS))) * label.policy.coverage.large)]
    # S5.mincomp.LABEL 11% compensation
    sp[, sc5mincomp_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.large)) + ((Energy_out_kcal + (((Energy_out_kcal*policy.energy.effect) * (1-0.11)) + (Energy_out_kcal*policy.energy.refeffectUS))) * label.policy.coverage.large)]
    # S5.maxcomp.LABEL 42% compensation
    sp[, sc5maxcomp_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.large)) + ((Energy_out_kcal + (((Energy_out_kcal*policy.energy.effect) * (1-0.42)) + (Energy_out_kcal*policy.energy.refeffectUS))) * label.policy.coverage.large)]
    
    # S6.LABEL Conso - Compensation + Reformulation, full coverage
    sp[, sc6_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.full)) + ((Energy_out_kcal + (((Energy_out_kcal*policy.energy.effect) * (1-0.265)) + (Energy_out_kcal*policy.energy.refeffectUS))) * label.policy.coverage.full)]
    # S6.mincomp.LABEL 11% compensation
    sp[, sc6mincomp_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.full)) + ((Energy_out_kcal + (((Energy_out_kcal*policy.energy.effect) * (1-0.11)) + (Energy_out_kcal*policy.energy.refeffectUS))) * label.policy.coverage.full)]
    # S6.maxcomp.LABEL 42% compensation
    sp[, sc6maxcomp_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.full)) + ((Energy_out_kcal + (((Energy_out_kcal*policy.energy.effect) * (1-0.42)) + (Energy_out_kcal*policy.energy.refeffectUS))) * label.policy.coverage.full)]
    
    #sensitivity analysis, reduction of -47 kcal in energy intake for the consumer
    sp[, sc1sensi_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.large)) + ((Energy_out_kcal + ((policy.energy.effect_sensi) * (1-0.265))) * label.policy.coverage.large)]
    sp[, sc3sensi_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.full)) + ((Energy_out_kcal + ((policy.energy.effect_sensi) * (1-0.265))) * label.policy.coverage.full)]
    sp[, sc5sensi_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.large)) + ((Energy_out_kcal + (((policy.energy.effect_sensi) * (1-0.265)) + (Energy_out_kcal*policy.energy.refeffectUS))) * label.policy.coverage.large)]
    sp[, sc6sensi_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.full)) + ((Energy_out_kcal + (((policy.energy.effect_sensi) * (1-0.265)) + (Energy_out_kcal*policy.energy.refeffectUS))) * label.policy.coverage.full)]
    
    #sensitivity analysis, using turnover instead of number of outlets
    sp[, sc1sensiT_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.large.turnover)) + ((Energy_out_kcal + ((Energy_out_kcal*policy.energy.effect) * (1-0.265))) * label.policy.coverage.large.turnover)]
    sp[, sc3sensiT_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.full.turnover)) + ((Energy_out_kcal + ((Energy_out_kcal*policy.energy.effect) * (1-0.265))) * label.policy.coverage.full.turnover)]
    sp[, sc5sensiT_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.large.turnover)) + ((Energy_out_kcal + (((Energy_out_kcal*policy.energy.effect) * (1-0.265)) + (Energy_out_kcal*policy.energy.refeffectUS))) * label.policy.coverage.large.turnover)]
    sp[, sc6sensiT_label_nrj := (Energy_out_kcal*(1-label.policy.coverage.full.turnover)) + ((Energy_out_kcal + (((Energy_out_kcal*policy.energy.effect) * (1-0.265)) + (Energy_out_kcal*policy.energy.refeffectUS))) * label.policy.coverage.full.turnover)]
    
    
    #List all your scenarios here!
    lili_label <- c("sc1_label","sc1mincomp_label","sc1maxcomp_label",
                    "sc3_label","sc3mincomp_label","sc3maxcomp_label",
                    "sc2_label","sc2sensi_label","sc4_label", "sc4sensi_label",
                    "sc5_label","sc5mincomp_label","sc5maxcomp_label",
                    "sc6_label","sc6mincomp_label","sc6maxcomp_label",
                    "sc1sensi_label","sc3sensi_label","sc5sensi_label","sc6sensi_label",
                    "sc1sensiT_label","sc3sensiT_label","sc5sensiT_label","sc6sensiT_label")
    
    #change in nrj
    #sc1_label_nrjdif        := ifelse(sc1_label_nrj > Energy_out_kcal, 0, sc1_label_nrj - Energy_out_kcal)]
    fnrj <- function(x,y) ifelse(x > y, 0, x - y)
    
    lili_nrj <- paste0(lili_label,"_nrj")
    lili_nrjdif <- paste0(lili_label,"_nrjdif")
    sp[, (lili_nrjdif)  := lapply(.SD,fnrj,Energy_out_kcal), .SDcols = lili_nrj]
    rm(fnrj,lili_nrj)
    
    
    #test <- sp[, sc1_label_nrj:sc6sensi_label_nrjdif]
    #test <- sp[,.(Energy_out_kcal,sc1_label_nrj,sc1_label_nrjdif)]
    #test <- sp[,.(Energy_out_kcal,sc6sensi_label_nrj,sc6sensi_label_nrjdif)]
    
    
    #linking the change in nrj and the change in weight (kg)
    #for men   = 17.7 * ((sc1_label_nrjdif*(4.2/1000))/PAL) #1 Kilocalorie = (4.2/1000) MJ
    #for women = 20.7 * ((sc1_label_nrjdif*(4.2/1000))/PAL)
    #Assuming a PAL of 1.5
    
    fbwss <- function(x,y) ifelse(y=="1", 17.7 * ((x*(4.2/1000))/1.5), 20.7 * ((x*(4.2/1000))/1.5)) #1 Kilocalorie = (4.2/1000) MJ
    lili_BWss <- paste0(lili_label,"_BWss")
    
    sp[, (lili_BWss)  := lapply(.SD,fbwss,sex), .SDcols = lili_nrjdif]
    
    #we have an absolute change in weight and we cannot convert in a relative one, so we Assume that everyone has the same height. 
    #We then calculate their weight given BMI and height (we assumed a height of 1.6m)
    sp[, wt_assumption := (bmi)*(1.6*1.6)]
    
    #we then apply the absolute difference to their calculated weight, and from their, calculate the relative difference
    fbwtrd <- function(x,y) (x+y)/(y)
    lili_wtrd <- paste0(lili_label,"_wtrd")
    
    sp[, (lili_wtrd)  := lapply(.SD,fbwtrd,wt_assumption), .SDcols = lili_BWss]
    
    
    
    #Then we apply that relative change to the baseline BMI
    #All the calculations are independent of height so id doesn't matter what value will you use.
    #sc1_label_bmi = bmi*sc1_label_wtrd
    fbmi <- function(x,y) y*x
    lili_bmi <- paste0(lili_label,"_bmi")
    
    sp[, (lili_bmi)  := lapply(.SD,fbmi,bmi), .SDcols = lili_wtrd]  #we are working on bmi and not bmi_lagged!! this is because we assumed that change in nrj change your bmi very quickly while your change in bmi take around 6y to change your RR
    #tt <- sp[,.(bmi,sc1_label_wtrd,sc1_label_bmi)]
    #summary(sp$bmi)
    #summary(sp$sc1_label_bmi)
    ###TEST FOR THE MODEl
    #10% DECREASE IN BMI
    #sp[, sc1_label_bmi_test := bmi*0.90)]
    #test <- sp[,.(sc1_label_nrjdif, sc1_label_BWss, bmi, wt_assumption, sc1_label_wtrd, sc1_label_bmi )]
    sp[, c(lili_BWss,lili_wtrd) := NULL] 
    
    
    
    # hist(sp$Energy_out_kcal, 10, freq = FALSE, col = rgb(0,0,1,1/4), xlim = c(0,2000), ylim = c(0, 0.002))  
    # hist(sp$sc1_label_nrj, 10, freq = FALSE, col = rgb(1,0,0,1/4), xlim = c(0.1,0.3), add = TRUE) 
    # legend('topright',c('Base-case','Traffic light mandatory'),
    #        fill = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), bty = 'n',
    #        border = NA)
    # 
    # plot(density(sp$Energy_out_kcal))
    # lines(density(sp$sc1_label_nrj), col=("purple"))
    # 
    # plot(density(sp$bmi))
    # lines(density(sp$sc1_label_bmi), col=("green"))
    
    
    ##############################################################################
    ##############################################################################
    ## SSB tax scenarios
    ##############################################################################
    ############################################################################## 
    
    #SSB levy effect Afhsin et al. (2017)     http://dx.doi.org/10.1371/journal.pone.0172277
    #each 10% increase in SSB tax was associated with a decline in SSB intake of 6.74% (95% CI = 3.08%, 10.40%). 
    sp[, policy.SSB.intake.effect.afshin := rnorm(1,mean = -0.0674, sd = (-0.0308-(-0.1040))/3.92)]
    sp[policy.SSB.intake.effect.afshin > 0, policy.SSB.intake.effect.afshin := 0]
    
    #SSB levy effect Teng et al. (2019)       https://doi.org/10.1111/obr.12868 
    #10% increase in SSB tax was associated with a decline in purchases and dietary intake of 10.0% (95% CI= 5.0% to 14.7%)
    sp[, policy.SSB.intake.effect.teng := rnorm(1,mean = -0.10, sd = (-0.05-(-0.147))/3.92)]
    sp[policy.SSB.intake.effect.teng > 0, policy.SSB.intake.effect.teng := 0]
    
    #SSB tax effect Andreyeva et al. (2022)   https://pubmed.ncbi.nlm.nih.gov/35648398/ Edited on 2024.10.28
    #Calculating the effect based on pass through rate 82% (95% CI = 66% to 98%) and price elasticity −1·59 (95% CI: [−2·11, −1·08]) 
    
    sp[, policy.SSB.pass.through := rnorm(1,mean = 0.82, sd = (0.98-0.66)/3.92)]
    sp[policy.SSB.pass.through < 0,  policy.SSB.pass.through := 0]
    
    sp[, policy.SSB.price.elasticity := rnorm(1,mean = -1.59, sd = (-1.08-(-2.11))/3.92)]
    sp[policy.SSB.price.elasticity > 0,  policy.SSB.price.elasticity := 0]
    
    #SSB tax effect Andreyeva et al. (2022)   https://pubmed.ncbi.nlm.nih.gov/35648398/ Edited on 2024.10.28
    #10% increase in SSB tax was associated with a decline in purchases and dietary intake
    sp[, policy.SSB.intake.effect.andreyava := 0.1 * policy.SSB.pass.through * policy.SSB.price.elasticity]
    sp[policy.SSB.intake.effect.andreyava > 0,  policy.SSB.intake.effect.andreyava := 0]
    
    #20% increase in SSB tax was associated with a decline in purchases and dietary intake 
    sp[, policy.SSB.intake.effect.andreyava2 := 0.2 * policy.SSB.pass.through * policy.SSB.price.elasticity]
    sp[policy.SSB.intake.effect.andreyava2 > 0,  policy.SSB.intake.effect.andreyava2 := 0]
    
    #30% increase in SSB tax was associated with a decline in purchases and dietary intake
    sp[, policy.SSB.intake.effect.andreyava3 := 0.3 * policy.SSB.pass.through * policy.SSB.price.elasticity]
    sp[policy.SSB.intake.effect.andreyava3 > 0,  policy.SSB.intake.effect.andreyava3 := 0]
    
    #Reformulation effect #edited by Edi 30% lower sugar content due to SSB taxes without changing consumption (as in Karl's paper)
    #see below 
    
    # ###############################################
    # ############### NEED TO DELETE IF YOU ARE NOT DOING SUGAR PATHWAY. 
    # ##############################################
    ##for now, only bmi pathway with this scenario, not the sugar or ssb intake!!!!!!!!!!!!!!!!!
    #As for now we are not having SSB_sugar_intake, we will assume linearity and that 20g sugar per 8 fluid ounces of SSB and 8 fluid ounces (≈ 227.3045ml) of SSB
    sp[, SSB_sugar := SSB*20/227.3045] #this will be used for later calculation.
    
    #####NOT WORKING FOR NOW!!!!!!!!!!!! FOR NOW, THE SUGAR PATHWAY DOES NOT EXIST
    #JUST FOR THE TEST, TO CHANGE
    #SDIL effect : the sugar intake from SSB will be reduce by X %, so my new sugar intake is intake*(1-%red)
    #assuming that every SSB is concerned by the policy and that every people in the sample will be too
    #assuming no compensation
    #sp[, sc1_SSB_sugar := (SSB_sugar*(1-SSB.policy.coverage)) + ((SSB_sugar + ((SSB_sugar*policy.SSB.sugar.effect) * (1-0))) * SSB.policy.coverage)] 
    
    #effect of reformulation (sugar content is reduced by 30% per 20 grams in 227.3045 ml --> 14 grams in 227.3045 ml) #edited by Edi
    sp[, sc1_SSB_sugar := ((SSB*(1-SSB.policy.coverage)) + ((SSB + ((SSB*0) * (1-0))) * SSB.policy.coverage)*14/227.3045)] 
    
    # change in sugar consumption due to reformulation
    sp[, sc1_SSB_sugardif := sc1_SSB_sugar - SSB_sugar]
    
    
    
    ### new SSB intake after the policy and CONVERTED INTO SUGAR WITHOUT REFORMULATION
    # S2.SSB - 10% tax - intake Afshin Full implementation
    #a 10% ssb tax will reduce the consumption of SSB
    #assuming no compensation 
    sp[, sc2_SSB := ((SSB*(1-SSB.policy.coverage)) + ((SSB + ((SSB*policy.SSB.intake.effect.afshin) * (1-0))) * SSB.policy.coverage)*20/227.3045)]
    
    # S3.SSB - 10% tax - intake TENG Full implementation
    #assuming no compensation    
    sp[, sc3_SSB := ((SSB*(1-SSB.policy.coverage)) + ((SSB + ((SSB*policy.SSB.intake.effect.teng) * (1-0))) * SSB.policy.coverage)*20/227.3045)]
    
    # S4.SSB - 10% tax - intake Andreyeva Full implementation (added by Edi)
    #assuming no compensation    
    sp[, sc4_SSB := ((SSB*(1-SSB.policy.coverage)) + ((SSB + ((SSB*policy.SSB.intake.effect.andreyava) * (1-0))) * SSB.policy.coverage)*20/227.3045)]
    
    # S5.SSB - 20% tax - intake Andreyeva Full implementation (added by Edi)
    #assuming no compensation    
    sp[, sc5_SSB := ((SSB*(1-SSB.policy.coverage)) + ((SSB + ((SSB*policy.SSB.intake.effect.andreyava2) * (1-0))) * SSB.policy.coverage)*20/227.3045)]
    
    # S6.SSB - 30% tax - intake Andreyeva Full implementation (added by Edi)
    #assuming no compensation    
    sp[, sc6_SSB := ((SSB*(1-SSB.policy.coverage)) + ((SSB + ((SSB*policy.SSB.intake.effect.andreyava3) * (1-0))) * SSB.policy.coverage)*20/227.3045)]
    
    
    #List all your scenarios here!
    #####BEWARE, FOR NOW SC1 IS ABOUT SSB SUGAR NOT SSB SO NOT LISTED HERE
    lili_ssbI <- c("sc2_SSB","sc3_SSB","sc4_SSB","sc5_SSB","sc6_SSB")
    
    # change in SSB consumption MIGHT NEED TO BE CHANGE AND DISCUSS!!!!!!!!!!!!!
    #sc2_SSBdif := sc2_SSB-SSB]
    fssb <- function(x,y) ifelse(x > y, 0, x - y)
    
    lili_ssbI_intakedif <- paste0(lili_ssbI,"dif")
    sp[, (lili_ssbI_intakedif)  := lapply(.SD,fssb,SSB_sugar), .SDcols = lili_ssbI] #here I used SSB_sugar instead of SSB (by Edi)
    rm(fssb,lili_ssbI_intakedif)
    
    
    
    
    #added by Edi
    #Consumption and Reformulation (sugar content is reduced by 30% per 20 grams in 227.3045 ml --> 14 grams in 227.3045 ml)
    sp[, sc2sensi_SSB := ((SSB*(1-SSB.policy.coverage)) + ((SSB + ((SSB*policy.SSB.intake.effect.afshin) * (1-0))) * SSB.policy.coverage)*14/227.3045)]
    
    # S3.SSB - 10% tax - intake TENG Full implementation
    #assuming no compensation    
    sp[, sc3sensi_SSB := ((SSB*(1-SSB.policy.coverage)) + ((SSB + ((SSB*policy.SSB.intake.effect.teng) * (1-0))) * SSB.policy.coverage)*14/227.3045)]
    
    # S4.SSB - 10% tax - intake Andreyeva Full implementation (added by Edi)
    #assuming no compensation    
    sp[, sc4sensi_SSB := ((SSB*(1-SSB.policy.coverage)) + ((SSB + ((SSB*policy.SSB.intake.effect.andreyava) * (1-0))) * SSB.policy.coverage)*14/227.3045)]
    
    # S5.SSB - 20% tax - intake Andreyeva Full implementation (added by Edi)
    #assuming no compensation    
    sp[, sc5sensi_SSB := ((SSB*(1-SSB.policy.coverage)) + ((SSB + ((SSB*policy.SSB.intake.effect.andreyava2) * (1-0))) * SSB.policy.coverage)*14/227.3045)]
    
    # S6.SSB - 30% tax - intake Andreyeva Full implementation (added by Edi)
    #assuming no compensation    
    sp[, sc6sensi_SSB := ((SSB*(1-SSB.policy.coverage)) + ((SSB + ((SSB*policy.SSB.intake.effect.andreyava3) * (1-0))) * SSB.policy.coverage)*14/227.3045)]
    
    #List all your scenarios here!
    #####BEWARE, FOR NOW SC1 IS ABOUT SSB SUGAR NOT SSB SO NOT LISTED HERE
    lili_ssbsensi <- c("sc2sensi_SSB","sc3sensi_SSB","sc4sensi_SSB","sc5sensi_SSB","sc6sensi_SSB")
    
    # change in SSB consumption MIGHT NEED TO BE CHANGE AND DISCUSS!!!!!!!!!!!!!
    #sc2_SSBdif := sc2_SSB-SSB]
    fssb <- function(x,y) ifelse(x > y, 0, x - y)
    
    lili_ssbsensi_intakedif <- paste0(lili_ssbsensi,"dif")
    sp[, (lili_ssbsensi_intakedif)  := lapply(.SD,fssb,SSB_sugar), .SDcols = lili_ssbsensi] #edited by Edi
    rm(fssb,lili_ssbsensi_intakedif)
    
    
    
    
    #####MODELLING THE BMI PATHWAY
    # SSB TAX POLICY > CHANAGE IN SSB > CHANGE IN SUGAR > CHANGE IN BMI > CHANGE IN RISK FROM BMI
    
    #need to calculate the new BMI 
    #if SSB scenario, change in sugar lead to change in BMI. 
    #Micha paper. assuming linearity and 20g sugar per 8 fluid ounces of SSB and 8 fluid ounces (≈ 227.3045ml) of SSB
    #Meta-analysis of 3 cohort studies: 0.10kg/m2 (0.05–0.15) increased BMI (baseline BMI <25 per 8 oz/d) , 0.23 (0.14–0.32) increased BMI (baseline BMI >=25 Per 8 oz/d) 
    #effect diff if BMI<25 or >=25
    #for scenario 2.SSB.sensi, I WANT TO use Nguyen et al (2023) https://doi.org/10.1016/j.ajcnut.2022.11.008  0.42-kg (95% CI: 0.26 kg, 0.58 kg; P < 0.01) // from prospective cohort study
    #BUT I DONT HAVE A DOSE RESPONSE
    
    
    #scenario looking at change from sugar
    r_new_bmi <- runif(1, min = 0, max = 1) ###HERE CHRIS, DOES THAT IS THE CORRECT WAY TO DO IT?
    eff <- qnorm(r_new_bmi,mean = 0.1/20, sd = ((0.15/20)-(0.05/20))/3.92)
    if (eff < 0 ) eff <- 0
    sp[, `:=` (sc1_SSB_sugar_bmi = bmi + (sc1_SSB_sugardif * eff),
               sc2_SSB_bmi = bmi + (sc2_SSBdif * eff),
               sc3_SSB_bmi = bmi + (sc3_SSBdif * eff),
               sc4_SSB_bmi = bmi + (sc4_SSBdif * eff),
               sc5_SSB_bmi = bmi + (sc5_SSBdif * eff),
               sc6_SSB_bmi = bmi + (sc6_SSBdif * eff),
               
               sc2sensi_SSB_bmi = bmi + (sc2sensi_SSBdif * eff),
               sc3sensi_SSB_bmi = bmi + (sc3sensi_SSBdif * eff),
               sc4sensi_SSB_bmi = bmi + (sc4sensi_SSBdif * eff),
               sc5sensi_SSB_bmi = bmi + (sc5sensi_SSBdif * eff),
               sc6sensi_SSB_bmi = bmi + (sc6sensi_SSBdif * eff))]
    
    eff <- qnorm(r_new_bmi,mean = 0.23/20, sd = ((0.32/20)-(0.14/20))/3.92)
    if (eff < 0 ) eff <- 0
    sp[bmi >=25, `:=` (sc1_SSB_sugar_bmi = bmi + (sc1_SSB_sugardif * eff),
                       sc2_SSB_bmi = bmi + (sc2_SSBdif * eff),
                       sc3_SSB_bmi = bmi + (sc3_SSBdif * eff),
                       sc4_SSB_bmi = bmi + (sc4_SSBdif * eff),
                       sc5_SSB_bmi = bmi + (sc5_SSBdif * eff),
                       sc6_SSB_bmi = bmi + (sc6_SSBdif * eff),
                       
                       sc2sensi_SSB_bmi = bmi + (sc2sensi_SSBdif * eff),
                       sc3sensi_SSB_bmi = bmi + (sc3sensi_SSBdif * eff),
                       sc4sensi_SSB_bmi = bmi + (sc4sensi_SSBdif * eff),
                       sc5sensi_SSB_bmi = bmi + (sc5sensi_SSBdif * eff),
                       sc6sensi_SSB_bmi = bmi + (sc6sensi_SSBdif * eff))]
    
    
    # test <- sp[bmi < sc6sensi_SSB_bmi,] #should never happen
    
    
    # hist(sp$SSB, 10, freq = FALSE, col = rgb(0,0,1,1/4), xlim = c(0,2000), ylim = c(0, 1e-3))
    # hist(sp$sc2_SSB, 10, freq = FALSE, col = rgb(1,0,0,1/4), xlim = c(0.1,0.3), add = TRUE)
    # legend('topright',c('Base-case','SSB 10% tax'),
    #        fill = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), bty = 'n',
    #        border = NA)
    # 
    # plot(density(sp$SSB))
    # lines(density(sp$sc2_SSB), col=("purple"))
    # 
    # plot(density(sp$bmi))
    # lines(density(sp$sc2_SSB_bmi), col=("green"))
    
    
    
    #List all your scenarios here!
    lili_1 <- c("sc1_label","sc1mincomp_label","sc1maxcomp_label",
                "sc3_label","sc3mincomp_label","sc3maxcomp_label",
                "sc2_label","sc2sensi_label","sc4_label", "sc4sensi_label",
                "sc5_label","sc5mincomp_label","sc5maxcomp_label",
                "sc6_label","sc6mincomp_label","sc6maxcomp_label",
                "sc1sensi_label","sc3sensi_label","sc5sensi_label", "sc6sensi_label",
                "sc1sensiT_label","sc3sensiT_label","sc5sensiT_label", "sc6sensiT_label",
                "sc1_SSB_sugar","sc2_SSB","sc3_SSB", "sc4_SSB", "sc5_SSB", "sc6_SSB",
                                "sc2sensi_SSB","sc3sensi_SSB","sc4sensi_SSB","sc5sensi_SSB","sc6sensi_SSB") #edited by Edi
    
    
    #obesity prevalence
    lili_bmi <- paste0(lili_1,"_bmi")
    lili_bmi_grp <- paste0(lili_1,"_bmi_grp")
    
    fobprev <- function(x) fifelse(x<20,"<20",
                                  fifelse(x>=20 & x<25,"20-25",
                                         fifelse(x>=25 & x<30,"25-30",
                                                fifelse(x>=30,">=30", NA_character_))))
    sp[, (lili_bmi_grp)  := lapply(.SD,fobprev), .SDcols = lili_bmi]
    #check
    # tt <- sp[,.(sc1_label_bmi,sc1_label_bmi_grp,sc2_SSB_bmi,sc2_SSB_bmi_grp)]
    # tt[,table(sc1_label_bmi_grp)]
    # tt[,table(sc2_SSB_bmi_grp)]
    
    rm(lili_bmi_grp,lili_bmi,lili_nrjdif,lili_BWss,lili_wtrd,fbmi,fbwtrd,fobprev,fbwss)
    
    #head(sp[,.(bmi, sc1_label_bmi, sc1_SSB_sugar_bmi, sc2_SSB_bmi, sc2sensi_SSB_bmi)]) 
    ##reduction in BMI by SSB < kcals policies
    
    #We are assuming that we have a 5y lag time between the bmi and the RR see https://doi.org/10.1136/bmj.i2793 
    #to be simpler, we will lag the bmi 
    setkey(sp, pid, year) # NECESSARY!!!
    
    ##LOOP NOT WORKING ON SHIFT FOR NOW, NEED TO DO IT BY HAND
    ##REMEMBER TO ADD NEW SCENARIOS IF NEEDED!!!!!!!!!!!!!!!!!!!!!!!!!!!
    sp[, `:=` (sc1_label_bmi_lagged         = shift(sc1_label_bmi, 5),    #BMIss bring back to 4.56 to be comparative with rr.Obesity
               sc1mincomp_label_bmi_lagged  = shift(sc1mincomp_label_bmi, 5),
               sc1maxcomp_label_bmi_lagged  = shift(sc1maxcomp_label_bmi, 5),
               sc3_label_bmi_lagged         = shift(sc3_label_bmi, 5),
               sc3mincomp_label_bmi_lagged  = shift(sc3mincomp_label_bmi, 5),
               sc3maxcomp_label_bmi_lagged  = shift(sc3maxcomp_label_bmi, 5),
               sc2_label_bmi_lagged         = shift(sc2_label_bmi, 5),
               sc2sensi_label_bmi_lagged    = shift(sc2sensi_label_bmi, 5),
               sc4_label_bmi_lagged         = shift(sc4_label_bmi, 5),
               sc4sensi_label_bmi_lagged    = shift(sc4sensi_label_bmi, 5),
               sc5_label_bmi_lagged         = shift(sc5_label_bmi, 5),
               sc5mincomp_label_bmi_lagged  = shift(sc5mincomp_label_bmi, 5),
               sc5maxcomp_label_bmi_lagged  = shift(sc5maxcomp_label_bmi, 5),
               sc6_label_bmi_lagged         = shift(sc6_label_bmi, 5),
               sc6mincomp_label_bmi_lagged  = shift(sc6mincomp_label_bmi, 5),
               sc6maxcomp_label_bmi_lagged  = shift(sc6maxcomp_label_bmi, 5),
               sc1sensi_label_bmi_lagged    = shift(sc1sensi_label_bmi, 5),
               sc3sensi_label_bmi_lagged    = shift(sc3sensi_label_bmi, 5),
               sc5sensi_label_bmi_lagged    = shift(sc5sensi_label_bmi, 5),
               sc6sensi_label_bmi_lagged    = shift(sc6sensi_label_bmi, 5),          
               sc1sensiT_label_bmi_lagged    = shift(sc1sensiT_label_bmi, 5),
               sc3sensiT_label_bmi_lagged    = shift(sc3sensiT_label_bmi, 5),
               sc5sensiT_label_bmi_lagged    = shift(sc5sensiT_label_bmi, 5),
               sc6sensiT_label_bmi_lagged    = shift(sc6sensiT_label_bmi, 5),
               #SSB (indirect)
               sc1_SSB_sugar_bmi_lagged    = shift(sc1_SSB_sugar_bmi, 5),
               sc2_SSB_bmi_lagged          = shift(sc2_SSB_bmi, 5),
               sc3_SSB_bmi_lagged          = shift(sc3_SSB_bmi, 5),
               sc4_SSB_bmi_lagged          = shift(sc4_SSB_bmi, 5),
               sc5_SSB_bmi_lagged          = shift(sc5_SSB_bmi, 5),
               sc6_SSB_bmi_lagged          = shift(sc6_SSB_bmi, 5),
              
               sc2sensi_SSB_bmi_lagged          = shift(sc2sensi_SSB_bmi, 5), #edited by Edi #identified a mistake #17.10.2024 and have changed it.
               sc3sensi_SSB_bmi_lagged          = shift(sc3sensi_SSB_bmi, 5),
               sc4sensi_SSB_bmi_lagged          = shift(sc4sensi_SSB_bmi, 5),
               sc5sensi_SSB_bmi_lagged          = shift(sc5sensi_SSB_bmi, 5),
               sc6sensi_SSB_bmi_lagged          = shift(sc6sensi_SSB_bmi, 5)), by = pid] 
    #check <- sp[,.(pid,year,bmi,bmi_lagged,sc1mincomp_label_bmi,sc1mincomp_label_bmi_lagged)]
    
    #SSB (direct)
    sp[, `:=` (sc1_SSB_sugar_lagged   = shift(sc1_SSB_sugar,5), #edited by Edi
               sc2_SSB_lagged         = shift(sc2_SSB, 5),
               sc3_SSB_lagged         = shift(sc3_SSB, 5),
               sc4_SSB_lagged         = shift(sc4_SSB, 5), 
               sc5_SSB_lagged         = shift(sc5_SSB, 5),
               sc6_SSB_lagged         = shift(sc6_SSB, 5),
               
               sc2sensi_SSB_lagged         = shift(sc2sensi_SSB, 5),
               sc3sensi_SSB_lagged         = shift(sc3sensi_SSB, 5),
               sc4sensi_SSB_lagged         = shift(sc4sensi_SSB, 5), 
               sc5sensi_SSB_lagged         = shift(sc5sensi_SSB, 5),
               sc6sensi_SSB_lagged         = shift(sc6sensi_SSB, 5)), by = pid ]
  
    
    
    
    #####MODELLING THE BMI PATHWAY
    # KCAL LABEL POLICY > CHANGE IN KCAL > CHANGE IN BMI > CHANGE IN RISK FROM BMI
    # SSB TAX POLICY > CHANAGE IN SSB > CHANGE IN SUGAR > CHANGE IN BMI > CHANGE IN RISK FROM BMI
    lili_bmilag <- paste0(lili_1,"_bmi_lagged")
    sc_name_all <- lili_1
    sc_name_SSB <- grep("_SSB", sc_name_all, value = TRUE)
    sc_name_label <- grep("_label$", sc_name_all, value = TRUE)

    lili_lagbmichd <- paste0(lili_1,"_lag.bmi.chd")
    lili_lagbmiisch <- paste0(lili_1,"_lag.bmi.isch")
    lili_lagbmihaem <- paste0(lili_1,"_lag.bmi.haem")
    
    #then we calculate the risk with the lagged bmi to take into account the 6y lag time
    frr.indiv.bmi <- function(x,y) y^((x-20)/4.56) #BMIss bring back to 4.56 to be comparative with rr.Obesity
    sp[, (lili_lagbmichd)  := lapply(.SD,frr.indiv.bmi,rr.obesity.chd), .SDcols = lili_bmilag]
    sp[, (lili_lagbmiisch) := lapply(.SD,frr.indiv.bmi,rr.obesity.isch), .SDcols = lili_bmilag]
    sp[, (lili_lagbmihaem) := lapply(.SD,frr.indiv.bmi,rr.obesity.haem), .SDcols = lili_bmilag]
    #for info, the formula for sc1_label for example then is 
    #sc1_label_lag.bmi.chd = rr.obesity.chd^((sc1_label_bmi_lagged-20)/4.56)
    rm(frr.indiv.bmi)
    
    # direct effect of ssb
    lili_ssbp <- c("sc1_SSB_sugar","sc2_SSB","sc3_SSB", "sc4_SSB", "sc5_SSB", "sc6_SSB",
                   "sc2sensi_SSB","sc3sensi_SSB","sc4sensi_SSB","sc5sensi_SSB","sc6sensi_SSB")
    
    # CKNote: Why above doesn't include the sc2sensi_SSB_lag.chd etc. scenarios???
    # The following chunk is to be replaced if the direct effect is implemented
    sp[, `:=` (
      sc1_SSB_sugar_lag.chd = 1,
      sc1_SSB_sugar_lag.isch = 1,
      sc1_SSB_sugar_lag.haem = 1,
      sc2_SSB_lag.chd = 1,
      sc2_SSB_lag.isch = 1,
      sc2_SSB_lag.haem = 1,
      sc3_SSB_lag.chd = 1,
      sc3_SSB_lag.isch = 1,
      sc3_SSB_lag.haem = 1,
      sc4_SSB_lag.chd = 1,
      sc4_SSB_lag.isch = 1,
      sc4_SSB_lag.haem = 1,
      sc5_SSB_lag.chd = 1,
      sc5_SSB_lag.isch = 1,
      sc5_SSB_lag.haem = 1,
      sc6_SSB_lag.chd = 1,
      sc6_SSB_lag.isch = 1,
      sc6_SSB_lag.haem = 1,
      sc2sensi_SSB_lag.chd = 1,
      sc2sensi_SSB_lag.isch = 1,
      sc2sensi_SSB_lag.haem = 1,
      sc3sensi_SSB_lag.chd = 1,
      sc3sensi_SSB_lag.isch = 1,
      sc3sensi_SSB_lag.haem = 1,
      sc4sensi_SSB_lag.chd = 1,
      sc4sensi_SSB_lag.isch = 1,
      sc4sensi_SSB_lag.haem = 1,
      sc5sensi_SSB_lag.chd = 1,
      sc5sensi_SSB_lag.isch = 1,
      sc5sensi_SSB_lag.haem = 1,
      sc6sensi_SSB_lag.chd = 1,
      sc6sensi_SSB_lag.isch = 1,
      sc6sensi_SSB_lag.haem = 1
    )]
    
    ##see if you want to add a lag time
    lili_ssbplag <- paste0(lili_ssbp,"_lagged")
    lili_lagssbpchd <- paste0(lili_ssbp,"_lag.chd")
    lili_lagssbpisch <- paste0(lili_ssbp,"_lag.isch")
    lili_lagssbphaem <- paste0(lili_ssbp,"_lag.haem")
    #then we calculate the risk with the new ssb intakes using lag time 
    frr.indiv.ssb <- function(x,y) y^((x)/20) #227.3045 as the risk is for 8oz (risk increase for each 8oz) #edited by Edi as we converted into sugar (20 grams)
    sp[, (lili_lagssbpchd)  := lapply(.SD,frr.indiv.ssb,rr.SSB.intake.BMIadj.chd), .SDcols = lili_ssbplag]
    sp[, (lili_lagssbpisch) := lapply(.SD,frr.indiv.ssb,rr.SSB.intake.BMIadj.isch), .SDcols = lili_ssbplag]
    sp[, (lili_lagssbphaem) := lapply(.SD,frr.indiv.ssb,rr.SSB.intake.BMIadj.haem), .SDcols = lili_ssbplag]
    #for info, the formula for sc1_label for example then is 
    #sc2_SSB.chd = rr.SSB.intake.BMIadj.chd^((SSB)/227.3045)
  
    rm(frr.indiv.ssb)
    
    # TODO delete unnecessary columns with scenario parameters thar are not
    # needed as soon as we get the individual level risks
    
    
    
    #sp_check <- sp[, . (sc1_label_bmi_lagged, sc1_label_lag.bmi.chd, rr.obesity.chd) ]
    #View (sp_check)
    
    #check <- (sp[,.(age, sc1_label_bmi_lagged, sc1_label_lag.bmi.chd, sc1_SSB_sugar_bmi_lagged, sc1_SSB_sugar_lag.bmi.chd)])
    
    #the people bmi<20 have no risk, not a protective risk
    sp[sc1_label_bmi_lagged < 20, `:=` (
      sc1_label_lag.bmi.chd = 1,
      sc1_label_lag.bmi.isch = 1,
      sc1_label_lag.bmi.haem = 1
    )]
    sp[sc1mincomp_label_bmi_lagged < 20, `:=` (
      sc1mincomp_label_lag.bmi.chd = 1,
      sc1mincomp_label_lag.bmi.isch = 1,
      sc1mincomp_label_lag.bmi.haem = 1
    )]
    sp[sc1maxcomp_label_bmi_lagged < 20, `:=` (
      sc1maxcomp_label_lag.bmi.chd = 1,
      sc1maxcomp_label_lag.bmi.isch = 1,
      sc1maxcomp_label_lag.bmi.haem = 1
    )]
    sp[sc3_label_bmi_lagged < 20, `:=` (
      sc3_label_lag.bmi.chd = 1,
      sc3_label_lag.bmi.isch = 1,
      sc3_label_lag.bmi.haem = 1
    )]
    sp[sc3mincomp_label_bmi_lagged < 20, `:=` (
      sc3mincomp_label_lag.bmi.chd = 1,
      sc3mincomp_label_lag.bmi.isch = 1,
      sc3mincomp_label_lag.bmi.haem = 1
    )]
    sp[sc3maxcomp_label_bmi_lagged < 20, `:=` (
      sc3maxcomp_label_lag.bmi.chd = 1,
      sc3maxcomp_label_lag.bmi.isch = 1,
      sc3maxcomp_label_lag.bmi.haem = 1
    )]
    sp[sc2_label_bmi_lagged < 20, `:=` (
      sc2_label_lag.bmi.chd = 1,
      sc2_label_lag.bmi.isch = 1,
      sc2_label_lag.bmi.haem = 1
    )]
    sp[sc2sensi_label_bmi_lagged < 20, `:=` (
      sc2sensi_label_lag.bmi.chd = 1,
      sc2sensi_label_lag.bmi.isch = 1,
      sc2sensi_label_lag.bmi.haem = 1
    )]
    sp[sc4_label_bmi_lagged < 20, `:=` (
      sc4_label_lag.bmi.chd = 1,
      sc4_label_lag.bmi.isch = 1,
      sc4_label_lag.bmi.haem = 1
    )]
    sp[sc4sensi_label_bmi_lagged < 20, `:=` (
      sc4sensi_label_lag.bmi.chd = 1,
      sc4sensi_label_lag.bmi.isch = 1,
      sc4sensi_label_lag.bmi.haem = 1
    )]
    sp[sc5_label_bmi_lagged < 20, `:=` (
      sc5_label_lag.bmi.chd = 1,
      sc5_label_lag.bmi.isch = 1,
      sc5_label_lag.bmi.haem = 1
    )]
    sp[sc5mincomp_label_bmi_lagged < 20, `:=` (
      sc5mincomp_label_lag.bmi.chd = 1,
      sc5mincomp_label_lag.bmi.isch = 1,
      sc5mincomp_label_lag.bmi.haem = 1
    )]
    sp[sc5maxcomp_label_bmi_lagged < 20, `:=` (
      sc5maxcomp_label_lag.bmi.chd = 1,
      sc5maxcomp_label_lag.bmi.isch = 1,
      sc5maxcomp_label_lag.bmi.haem = 1
    )]
    sp[sc6_label_bmi_lagged < 20, `:=` (
      sc6_label_lag.bmi.chd = 1,
      sc6_label_lag.bmi.isch = 1,
      sc6_label_lag.bmi.haem = 1
    )]
    sp[sc6mincomp_label_bmi_lagged < 20, `:=` (
      sc6mincomp_label_lag.bmi.chd = 1,
      sc6mincomp_label_lag.bmi.isch = 1,
      sc6mincomp_label_lag.bmi.haem = 1
    )]
    sp[sc6maxcomp_label_bmi_lagged < 20, `:=` (
      sc6maxcomp_label_lag.bmi.chd = 1,
      sc6maxcomp_label_lag.bmi.isch = 1,
      sc6maxcomp_label_lag.bmi.haem = 1
    )]
    sp[sc1sensi_label_bmi_lagged < 20, `:=` (
      sc1sensi_label_lag.bmi.chd = 1,
      sc1sensi_label_lag.bmi.isch = 1,
      sc1sensi_label_lag.bmi.haem = 1
    )]
    sp[sc3sensi_label_bmi_lagged < 20, `:=` (
      sc3sensi_label_lag.bmi.chd = 1,
      sc3sensi_label_lag.bmi.isch = 1,
      sc3sensi_label_lag.bmi.haem = 1
    )]
    sp[sc5sensi_label_bmi_lagged < 20, `:=` (
      sc5sensi_label_lag.bmi.chd = 1,
      sc5sensi_label_lag.bmi.isch = 1,
      sc5sensi_label_lag.bmi.haem = 1
    )]
    sp[sc6sensi_label_bmi_lagged < 20, `:=` (
      sc6sensi_label_lag.bmi.chd = 1,
      sc6sensi_label_lag.bmi.isch = 1,
      sc6sensi_label_lag.bmi.haem = 1
    )]
    sp[sc1sensiT_label_bmi_lagged < 20, `:=` (
      sc1sensiT_label_lag.bmi.chd = 1,
      sc1sensiT_label_lag.bmi.isch = 1,
      sc1sensiT_label_lag.bmi.haem = 1
    )]
    sp[sc3sensiT_label_bmi_lagged < 20, `:=` (
      sc3sensiT_label_lag.bmi.chd = 1,
      sc3sensiT_label_lag.bmi.isch = 1,
      sc3sensiT_label_lag.bmi.haem = 1
    )]
    sp[sc5sensiT_label_bmi_lagged < 20, `:=` (
      sc5sensiT_label_lag.bmi.chd = 1,
      sc5sensiT_label_lag.bmi.isch = 1,
      sc5sensiT_label_lag.bmi.haem = 1
    )]
    sp[sc6sensiT_label_bmi_lagged < 20, `:=` (
      sc6sensiT_label_lag.bmi.chd = 1,
      sc6sensiT_label_lag.bmi.isch = 1,
      sc6sensiT_label_lag.bmi.haem = 1
    )]
    #SSB
    sp[sc1_SSB_sugar_bmi_lagged < 20, `:=` (
      sc1_SSB_sugar_lag.bmi.chd = 1,
      sc1_SSB_sugar_lag.bmi.isch = 1,
      sc1_SSB_sugar_lag.bmi.haem = 1
    )]
    sp[sc2_SSB_bmi_lagged < 20, `:=` (
      sc2_SSB_lag.bmi.chd = 1,
      sc2_SSB_lag.bmi.isch = 1,
      sc2_SSB_lag.bmi.haem = 1
    )]
    sp[sc3_SSB_bmi_lagged < 20, `:=` (
      sc3_SSB_lag.bmi.chd = 1,
      sc3_SSB_lag.bmi.isch = 1,
      sc3_SSB_lag.bmi.haem = 1
    )]
    sp[sc4_SSB_bmi_lagged < 20, `:=` (
      sc4_SSB_lag.bmi.chd = 1,
      sc4_SSB_lag.bmi.isch = 1,
      sc4_SSB_lag.bmi.haem = 1
    )]
    sp[sc5_SSB_bmi_lagged < 20, `:=` (
      sc5_SSB_lag.bmi.chd = 1,
      sc5_SSB_lag.bmi.isch = 1,
      sc5_SSB_lag.bmi.haem = 1
    )]
    sp[sc6_SSB_bmi_lagged < 20, `:=` (
      sc6_SSB_lag.bmi.chd = 1,
      sc6_SSB_lag.bmi.isch = 1,
      sc6_SSB_lag.bmi.haem = 1
    )]
    
    sp[sc2sensi_SSB_bmi_lagged < 20, `:=` (
      sc2sensi_SSB_lag.bmi.chd = 1,
      sc2sensi_SSB_lag.bmi.isch = 1,
      sc2sensi_SSB_lag.bmi.haem = 1
    )]
    sp[sc3sensi_SSB_bmi_lagged < 20, `:=` (
      sc3sensi_SSB_lag.bmi.chd = 1,
      sc3sensi_SSB_lag.bmi.isch = 1,
      sc3sensi_SSB_lag.bmi.haem = 1
    )]
    sp[sc4sensi_SSB_bmi_lagged < 20, `:=` (
      sc4sensi_SSB_lag.bmi.chd = 1,
      sc4sensi_SSB_lag.bmi.isch = 1,
      sc4sensi_SSB_lag.bmi.haem = 1
    )]
    sp[sc5sensi_SSB_bmi_lagged < 20, `:=` (
      sc5sensi_SSB_lag.bmi.chd = 1,
      sc5sensi_SSB_lag.bmi.isch = 1,
      sc5sensi_SSB_lag.bmi.haem = 1
    )]
    sp[sc6sensi_SSB_bmi_lagged < 20, `:=` (
      sc6sensi_SSB_lag.bmi.chd = 1,
      sc6sensi_SSB_lag.bmi.isch = 1,
      sc6sensi_SSB_lag.bmi.haem = 1
    )] #added by Edi
    
    #for SSB direct effect RR = 1 for people with < 20 grams of sugar
    sp[sc1_SSB_sugar_lagged < 20, `:=` (
      sc1_SSB_sugar_lag.chd = 1,
      sc1_SSB_sugar_lag.isch = 1,
      sc1_SSB_sugar_lag.haem = 1
    )]
    sp[sc2_SSB_lagged < 20, `:=` (
      sc2_SSB_lag.chd = 1,
      sc2_SSB_lag.isch = 1,
      sc2_SSB_lag.haem = 1
    )]
    sp[sc3_SSB_lagged < 20, `:=` (
      sc3_SSB_lag.chd = 1,
      sc3_SSB_lag.isch = 1,
      sc3_SSB_lag.haem = 1
    )]
    sp[sc4_SSB_lagged < 20, `:=` (
      sc4_SSB_lag.chd = 1,
      sc4_SSB_lag.isch = 1,
      sc4_SSB_lag.haem = 1
    )]
    sp[sc5_SSB_lagged < 20, `:=` (
      sc5_SSB_lag.chd = 1,
      sc5_SSB_lag.isch = 1,
      sc5_SSB_lag.haem = 1
    )]
    sp[sc6_SSB_lagged < 20, `:=` (
      sc6_SSB_lag.chd = 1,
      sc6_SSB_lag.isch = 1,
      sc6_SSB_lag.haem = 1
    )]
    
    sp[sc2sensi_SSB_lagged < 20, `:=` (
      sc2sensi_SSB_lag.chd = 1,
      sc2sensi_SSB_lag.isch = 1,
      sc2sensi_SSB_lag.haem = 1
    )]
    sp[sc3sensi_SSB_lagged < 20, `:=` (
      sc3sensi_SSB_lag.chd = 1,
      sc3sensi_SSB_lag.isch = 1,
      sc3sensi_SSB_lag.haem = 1
    )]
    sp[sc4sensi_SSB_lagged < 20, `:=` (
      sc4sensi_SSB_lag.chd = 1,
      sc4sensi_SSB_lag.isch = 1,
      sc4sensi_SSB_lag.haem = 1
    )]
    sp[sc5sensi_SSB_lagged < 20, `:=` (
      sc5sensi_SSB_lag.chd = 1,
      sc5sensi_SSB_lag.isch = 1,
      sc5sensi_SSB_lag.haem = 1
    )]
    sp[sc6sensi_SSB_lagged < 20, `:=` (
      sc6sensi_SSB_lag.chd = 1,
      sc6sensi_SSB_lag.isch = 1,
      sc6sensi_SSB_lag.haem = 1
    )]

    
    # RR product
    if (SSB_direct_effect) {
      for (j in sc_name_label) { # No SSB change hence rr is as sc0
        for (ds in c("chd", "isch", "haem")) {
          all_rr <- paste0("all_rr_", ds, "_", j) # product of RR
          bmi_rr <- paste0(j, "_lag.bmi.", ds)
          ssb_rr <- paste0("rr.SSB.intake.BMIadj.", ds, ".indiv") # sc0 risk
          sp[, (all_rr) := rr.education * get(ssb_rr) * get(bmi_rr)]
          sp[, (bmi_rr) := NULL]
        }
      }
      
      for (j in sc_name_SSB) { # No SSB change hence rr is as sc0
        for (ds in c("chd", "isch", "haem")) {
          all_rr <- paste0("all_rr_", ds, "_", j) # product of RR
          bmi_rr <- paste0(j, "_lag.bmi.", ds)
          ssb_rr <- paste0(j, "_lag.", ds)
          sp[, (all_rr) := rr.education * get(ssb_rr) * get(bmi_rr)]
          sp[, c(bmi_rr, ssb_rr) := NULL]
        }
      }
      sp[, c(
        "rr.education",
        "rr.SSB.intake.BMIadj.chd.indiv",
        "rr.SSB.intake.BMIadj.isch.indiv",
        "rr.SSB.intake.BMIadj.haem.indiv"
      ) := NULL]
    } else { # If no direct SSB effect
      for (j in sc_name_label) { # No SSB change hence rr is as sc0
        for (ds in c("chd", "isch", "haem")) {
          all_rr <- paste0("all_rr_", ds, "_", j) # product of RR
          bmi_rr <- paste0(j, "_lag.bmi.", ds)
          sp[, (all_rr) := rr.education * get(bmi_rr)]
          sp[, (bmi_rr) := NULL]
        }
      }
      
      for (j in sc_name_SSB) { # No SSB change hence rr is as sc0
        for (ds in c("chd", "isch", "haem")) {
          all_rr <- paste0("all_rr_", ds, "_", j) # product of RR
          bmi_rr <- paste0(j, "_lag.bmi.", ds)
          ssb_rr <- paste0(j, "_lag.", ds) # Used only to delete the col
          sp[, (all_rr) := rr.education * get(bmi_rr)]
          sp[, c(bmi_rr, ssb_rr) := NULL]
        }
      }
      
      sp[, c(
        "rr.education",
        "rr.SSB.intake.BMIadj.chd.indiv",
        "rr.SSB.intake.BMIadj.isch.indiv",
        "rr.SSB.intake.BMIadj.haem.indiv"
      ) := NULL]
    }

    
    #cleaning var
    #sp[, c(lili_bmilag) := NULL] #deactivated by Edi
    
    ##TEST FOR A LOOP, NOT WORKING YET
    # test <- copy(sp)
    # test <- test[,.(sc1_label_bmi_lagged,sc1_label_lag.bmi.chd,sc1mincomp_label_bmi_lagged,sc1mincomp_label_lag.bmi.chd)]
    
    # lili_bmilagTEST <- c("sc1_label_bmi_lagged","sc1mincomp_label_bmi_lagged")
    # lili_lagchdTEST <- c("sc1_label_lag.chdTEST22","sc1mincomp_label_lag.chdTEST22")
    # lili_lagchdTEST <- c("sc1_label_lag.chd","sc1mincomp_label_lag.chd")
    #working but I can not put a list in y (instead of sc1_label_bmi_lagged in the current example...)
    #test[, (lili_lagchdTEST)  := lapply(.SD, function(x,y) ifelse(y<20 & !is.na(y),1,x), sc1_label_bmi_lagged), .SDcols = lili_lagchdTEST]
    #tt <- test[sc1_label_lag.chdTEST!=sc1_label_lag.chd,] #should be nrow=12995
    
    #not working at all
    #test[, mapply(function(X,Y) {(lili_lagchdTEST) := lapply(X, function(x,y) ifelse(y<20 & !is.na(y),1,x), Y)}, X=lili_lagchdTEST, Y=lili_bmilagTEST)]
    
    
    
    ###################CHIRS/KARL NEED TO CHECK. FOR NOW, AS WE DONT HAVE MORTALITY BY EDUCATION IN GERMANY, 
    #from table 2 here, we know the age-adjusted mortality rate ratios for education in Germany  (https://bmjopen.bmj.com/content/9/10/e028001)
    # so, for now, I am calculating the risk for chd for each individual based on their bmi (with the lag) as above 
    #and to add the mortality differences by education I multiply this risk by the mortality rate ratios for education form the paper, and then I do the parf.
    ##DOES THAT IS OK???? # CK: Nope, it wasn't. 

    
    # adding RR =1 for those with NA for BMI/SSB lag as they within the lag time (have no risk) #added by Edi

    if (SSB_direct_effect) {
      all_rr_chd_label <- paste0("all_rr_chd_",sc_name_label)
      all_rr_haem_label <- paste0("all_rr_haem_",sc_name_label)
      all_rr_isch_label <- paste0("all_rr_isch_",sc_name_label)
    
      for (i in seq_along(lili_bmilag)) {
      sp[is.na(get(lili_bmilag[i])), (all_rr_chd_label[i]) := 1]
       }
      for (i in seq_along(lili_bmilag)) {
      sp[is.na(get(lili_bmilag[i])), (all_rr_haem_label[i]) := 1]
       }
      for (i in seq_along(lili_bmilag)) {
      sp[is.na(get(lili_bmilag[i])), (all_rr_isch_label[i]) := 1] 
       }  
    
      all_rr_chd_SSB <- paste0("all_rr_chd_",sc_name_SSB)
      all_rr_haem_SSB <- paste0("all_rr_haem_",sc_name_SSB)
      all_rr_isch_SSB <- paste0("all_rr_isch_",sc_name_SSB)
    
      for (i in seq_along(lili_ssbplag)) {
      sp[is.na(get(lili_ssbplag[i])), (all_rr_chd_SSB[i]) := 1]
       }
      for (i in seq_along(lili_ssbplag)) {
      sp[is.na(get(lili_ssbplag[i])), (all_rr_haem_SSB[i]) := 1]
       }
      for (i in seq_along(lili_ssbplag)) {
      sp[is.na(get(lili_ssbplag[i])), (all_rr_isch_SSB[i]) := 1]   
       }  
      
    } else {
      all_rr_chd_all <- paste0("all_rr_chd_",sc_name_all)
      all_rr_haem_all <- paste0("all_rr_haem_",sc_name_all)
      all_rr_isch_all <- paste0("all_rr_isch_",sc_name_all)
    
      for (i in seq_along(lili_bmilag)) {
      sp[is.na(get(lili_bmilag[i])), (all_rr_chd_all[i]) := 1]
        }
      for (i in seq_along(lili_bmilag)) {
      sp[is.na(get(lili_bmilag[i])), (all_rr_haem_all[i]) := 1]
       }
      for (i in seq_along(lili_bmilag)) {
      sp[is.na(get(lili_bmilag[i])), (all_rr_isch_all[i]) := 1] 
       }  
    }
    
    
    # fill all NAs in numeric cols with 1. This excludes the "dead_chd"  "dead_isch" "dead_haem"
    tt <- sp[, sapply(.SD, anyNA)]
    ttt <- sp[, sapply(.SD, is.numeric)]
    cn <- names(sp)[as.logical(tt*ttt)]
    cn <- grep("^dead_", cn, value = TRUE, invert = TRUE) # extra safety to exclude dead cols
    setnafill(sp, "const", 1, cols = cn)
    
    

    # Estimate mrtl in scenario through bmi
    lili_mr_bmi_chd <- paste0(sc_name_all,"_mr_bmi_chd")
    lili_mr_bmi_isch <- paste0(sc_name_all,"_mr_bmi_isch")
    lili_mr_bmi_haem <- paste0(sc_name_all,"_mr_bmi_haem")
    
    
    sp[, (lili_mr_bmi_chd) := lapply(.SD, `*`, p0_bmi_mrtl_chd),
       .SDcols = paste0("all_rr_chd_", sc_name_all)]
    sp[, paste0("all_rr_chd_", sc_name_all) := NULL]
    
    sp[, (lili_mr_bmi_isch) := lapply(.SD, `*`, p0_bmi_mrtl_isch),
       .SDcols = paste0("all_rr_isch_", sc_name_all)]
    sp[, paste0("all_rr_isch_", sc_name_all) := NULL]
    
    sp[, (lili_mr_bmi_haem) := lapply(.SD, `*`, p0_bmi_mrtl_haem),
       .SDcols = paste0("all_rr_haem_", sc_name_all)]
    sp[, paste0("all_rr_haem_", sc_name_all) := NULL]
    
    #for info, the formula for sc1_label for example then is 
    #sc1_label_mr_bmi_chd = p0_bmi_mrtl_chd * sc1_label_lag.bmi.chd.educ
    
    #as <20  have rr of 1 already no need to do >=20 here, just p0*rr
    
    #tt<-sp[,.(pid,year,rr.obesity.chd,bmi,bmi_lagged,rr.obesity.chd.indiv,mr_chd,sc3_label_bmi,sc3_label_bmi_lagged,sc3_label_lag.bmi.chd,p0_bmi_mrtl_chd,sc3_label_mr_bmi_chd)]
    
    #check <- (sp[,.(age, sc1_label_bmi_lagged, sc1_label_lag.bmi.chd, sc1_label_lag.bmi.chd.educ, sc1_SSB_sugar_bmi_lagged, sc1_SSB_sugar_lag.bmi.chd, sc1_SSB_sugar_lag.bmi.chd.educ, p0_bmi_mrtl_chd, sc1_label_mr_bmi_chd, sc1_SSB_sugar_mr_bmi_chd)])
    
    
    
    #####MODELLING THE SSB > RRAdjBMI PATHWAY
    # ##FOR SSB ON CVD - BMI adjusted
    # SSB TAX POLICY > CHANGE IN SSB > CHANGE IN RISK FROM BMI
    
    #List all your scenarios here! Only the one using the SSB direct pathway (no labels)
    ##NOTE I DO NOT LIST sc1_SSB_sugar AS ONLY CHANGE SUGAR, NOT INTAKE!! NEED TO CHANGE IF YOU ARE ADDING PATHWAY MAYBE??? NEED TO CHECK
    
    ##################DO YOU WANT A LAG TIME FOR RR SSB BMIADJ? NEED TO MAKE THE DECISION AND PUT THE LAG TIME HERE IF NEEDED (NEED TO SEARCH IN LITERATURE AND DISCUSS WITH CHRIS AND MARTIN)
    #do we need a lag time for SSB consumption on risk????? IF YES, NEED TO HAVE A SSB_lagged
    
    ## adding lag time (edited by Edi)

    
    

    
    ###################CHIRS/KARL NEED TO CHECK. FOR NOW, AS WE DONT HAVE MORTALITY BY EDUCATION IN GERMANY, 
    #from table 2 here, we know the age-adjusted mortality rate ratios for education in Germany  (https://bmjopen.bmj.com/content/9/10/e028001)
    # so, for now, I am calculating the risk for chd for each individual based on their bmi (with the lag) as above 
    #and to add the mortality differences by education I multiply this risk by the mortality rate ratios for education form the paper, and then I do the parf.
    ##DOES THAT IS OK???? CKNote: It wasn't, deleted
    

    

    
    #check <- (sp[,.(age, rr.obesity.chd, sc1_label_bmi_lagged, sc1_label_lag.bmi.chd, sc1_label_lag.bmi.chd.educ, 
    # sc2_SSB_bmi_lagged, sc2_SSB_lag.bmi.chd, sc2_SSB_lag.bmi.chd.educ, p0_bmi_mrtl_chd, 
    # sc1_label_mr_bmi_chd, sc2_SSB_mr_bmi_chd, sc2_SSB_lagged, sc2_SSB_lag.chd, sc2_SSB_lag.chd.educ, p0_SSB_intake_mrtl_chd, sc2_SSB_mr_I_chd)])
    

    
    
    #################################################################################################################################################
    #################################################################################################################################################
# Simulate mortality

    for (h in sc_name_all) {
      pid_toremove <- integer(0) # a place holder
      
      hchd <- paste0(h,"_mr_bmi_chd")
      hh <- sub("_mr_", "_dead_", hchd)
      hhh <- sub("_mr_", "_longdead_", hchd)
      # hh <- paste(sub("_mr.*", "", h),"_dead_bmi_chd",sep="")
      # hhh <- paste(sub("_mr.*", "", h),"_longdead_bmi_chd",sep="")
      for (i in sort(unique(sp$year))) {
        #print(paste0("Year: ", i))
        
        # rebalance the mortality rates to account for simulants removed from the
        # population when dead
        if (length(pid_toremove) > 0L) {
          correction_factor <- as.numeric(sp[year == i, lapply(.SD,sum), .SDcols=hchd]/sp[year == i & !pid %in% pid_toremove, lapply(.SD,sum), .SDcols=hchd])
          sp[year == i, (hchd) := lapply(.SD,"*",correction_factor), .SDcols=hchd]
          #print(paste0("correction factor: ", correction_factor))
        }
        
        sp[year == i & !pid %in% pid_toremove, paste(hh) := lapply(.SD,">",rn_mrtl), .SDcols=hchd]
        pid_toremove <- sp[year <= i & (sp[[hh]]), pid]
        #print(paste0("PID to remove: ", paste(pid_toremove, collapse = ",")))
      }
      sp[, (hhh) := lapply(.SD,cumsum), by = pid, .SDcols=hh]
      sp[sp[[hh]] == FALSE & sp[[hhh]] == 1, paste(hh) := NA] # Best for later calculations.
      sp[, (hhh) := NULL]
      # So dead == FALSE mean alive, dead == TRUE means died that year,
      # and dead == NA means died in a previous year
      #dcast(sp, year~sp[[hh]]) # Note NA some years is decreasing because simulants get removed from the dataset because their age is > 89
      #table(sp$sc1_label_dead_bmi_chd)
      
      hisch <- paste0(h,"_mr_bmi_isch")
      hh <- sub("_mr_", "_dead_", hisch)
      hhh <- sub("_mr_", "_longdead_", hisch)
      # hh <- paste(sub("_mr.*", "", h),"_dead_bmi_isch",sep="")
      # hhh <- paste(sub("_mr.*", "", h),"_longdead_bmi_isch",sep="")
      for (i in sort(unique(sp$year))) {
        #print(paste0("Year: ", i))
        
        # rebalance the mortality rates to account for simulants removed from the
        # population when dead
        if (length(pid_toremove) > 0L) {
          correction_factor <- as.numeric(sp[year == i, lapply(.SD,sum), .SDcols=hisch]/sp[year == i & !pid %in% pid_toremove, lapply(.SD,sum), .SDcols=hisch])
          sp[year == i, (hisch) := lapply(.SD,"*",correction_factor), .SDcols=hisch]
          #print(paste0("correction factor: ", correction_factor))
        }
        
        sp[year == i & !pid %in% pid_toremove, paste(hh) := lapply(.SD,">",rn_mrtl), .SDcols=hisch]
        pid_toremove <- sp[year <= i & (sp[[hh]]), pid]
        #print(paste0("PID to remove: ", paste(pid_toremove, collapse = ",")))
      }
      sp[, c(hhh) := lapply(.SD,cumsum), by = pid, .SDcols=hh]
      sp[sp[[hh]] == FALSE & sp[[hhh]] == 1, paste(hh) := NA] # Best for later calculations.
      sp[, (hhh) := NULL]
      
      hhaem <- paste0(h,"_mr_bmi_haem")
      hh <- sub("_mr_", "_dead_", hhaem)
      hhh <- sub("_mr_", "_longdead_", hhaem)
      # hh <- paste(sub("_mr.*", "", h),"_dead_bmi_haem",sep="")
      # hhh <- paste(sub("_mr.*", "", h),"_longdead_bmi_haem",sep="")
      for (i in sort(unique(sp$year))) {
        #print(paste0("Year: ", i))
        
        # rebalance the mortality rates to account for simulants removed from the
        # population when dead
        if (length(pid_toremove) > 0L) {
          correction_factor <- as.numeric(sp[year == i, lapply(.SD,sum), .SDcols=hhaem]/sp[year == i & !pid %in% pid_toremove, lapply(.SD,sum), .SDcols=hhaem])
          sp[year == i, (hhaem) := lapply(.SD,"*",correction_factor), .SDcols=hhaem]
          #print(paste0("correction factor: ", correction_factor))
        }
        
        sp[year == i & !pid %in% pid_toremove, paste(hh) := lapply(.SD,">",rn_mrtl), .SDcols=hhaem]
        pid_toremove <- sp[year <= i & (sp[[hh]]), pid]
        #print(paste0("PID to remove: ", paste(pid_toremove, collapse = ",")))
      }
      sp[, c(hhh) := lapply(.SD,cumsum), by = pid, .SDcols=hh]
      sp[sp[[hh]] == FALSE & sp[[hhh]] == 1, paste(hh) := NA] # Best for later calculations.
      sp[, (hhh) := NULL]
      
      sp[, c(hchd, hisch, hhaem) := NULL]
    } # end of mortality simulation
    
    sp[, `:=` (dead_bmi_chd=dead_chd,dead_bmi_isch=dead_isch,dead_bmi_haem=dead_haem,
               dead_I_chd=dead_chd,dead_I_isch=dead_isch,dead_I_haem=dead_haem)]
    
    liliok <- paste0(sub(
      "_mr_",
      "_dead_",
      c(lili_mr_bmi_chd, lili_mr_bmi_isch, lili_mr_bmi_haem)
    ))
    
    ## Results (distributional impact)
    lili_totY <- c()
    testdata2 <- data.table()
    top_year <- (init_year+sim_hor)-1
    for (y in init_year:top_year){
      byyear <- sp[year==y,]
      testdata1 <- data.table()
      for (d in c("dead_bmi_chd","dead_bmi_isch","dead_bmi_haem","dead_I_chd","dead_I_isch","dead_I_haem")){
        for (ses in SES){
          byses <- byyear[SES==ses,]
          lili_totY <- c(lili_totY, paste("expected",d, "SES",ses, "Year", y, sep="_"))
          tt <- data.table((sum(byses[[d]]==TRUE,na.rm=T)) * scale_factor)
          names(tt) <- c(paste("expected",d, "SES",ses, "Year", y, sep="_"))
          testdata1 <- cbind(testdata1, tt)
          rm(tt)
          for (e in liliok){
            if (grepl(d, e)==TRUE) {
              #f <- sub('_dead.*', '', e)
              #lili_totY <- c(lili_totY, paste("DPPs", d, f, "SES",ses, "Year", y, sep="_"))
              lili_totY <- c(lili_totY, paste("DPPs", e, "SES",ses, "Year", y, sep="_"))
              tt <-  data.table((sum(byses[[d]]==TRUE,na.rm=T) - sum(byses[[e]]==TRUE,na.rm=T)) * scale_factor)
              #names(tt) <- c(paste("DPPs",d, f, "SES",ses, "Year", y, sep="_"))
              names(tt) <- c(paste("DPPs",e, "SES",ses, "Year", y, sep="_"))
              testdata1 <- cbind(testdata1, tt)
              rm(tt)
              
              # tt <-  data.table((sum(byses[[e]]==TRUE,na.rm=T)) * scale_factor)
              # names(tt) <- c(paste("verif",e, "SES",ses, "Year", y, sep="_"))
              # testdata1 <- cbind(testdata1, tt)
              # rm(tt,f)
              
            }
          }
          rm(byses)
        }
      }
      testdata2 <- cbind(testdata2, testdata1)
      rm(byyear,testdata1)
    }
    rm(y,d,e,ses)
    #saveRDS(lili_totY, file="lili_totY.RData")
    
    ##obesity prevalence
    lili_totOY <- c()
    obesity_prevalence <- data.table(V1=c("<20", ">=30", "20-25", "25-30"))
    top_year <- (init_year+sim_hor)-1
    for (y in init_year:top_year){
      byyear <- sp[year==y,]
      testdata1 <- data.table(V1=c("<20", ">=30", "20-25", "25-30"))
      for (d in c("bmi_grp")){
        for (ses in SES){
          byses <- byyear[SES==ses,]
          lili_totOY <- c(lili_totOY, paste("expected",d, "SES",ses, "Year", y, sep="_"))
          tt <- data.table(prop.table(table(byses[[d]])))
          names(tt) <- c("V1", paste("expected",d, "SES",ses, "Year", y, sep="_"))
          testdata1 <- testdata1[tt, on=c("V1")]
          rm(tt)
          for (e in lili_1){  ##its normal to have lili_1 for obesity!!
            sc <- paste(e, d,sep="_")
            lili_totOY <- c(lili_totOY, paste("expected", d, e, "SES",ses, "Year", y, sep="_"))
            tt <- data.table(prop.table(table(byses[[sc]])))
            names(tt) <- c("V1", paste("expected", d, e, "SES",ses, "Year", y, sep="_"))
            testdata1 <- testdata1[tt, on=c("V1")]
            rm(tt)
          }
          rm(byses)
        }
      }
      obesity_prevalence <- obesity_prevalence[testdata1, on=c("V1")]
      rm(byyear,testdata1)
    }
    
    rm(y,d,e,ses,top_year)
    #saveRDS(lili_totOY, file="lili_totOY.RData")
    
    obesity_prevalence2 <- obesity_prevalence[V1==">=30",]
    obesity_prevalence2<-obesity_prevalence2[,V1:=NULL]
    #saveRDS(obesity_prevalence2, file="obesity_prevalence2.RData")
    
    lili_totYtest <- c(lili_totY, lili_totOY)
    saveRDS(lili_totYtest, file="lili_totYtest.RData")
    
    time.new <- difftime(Sys.time(),time.old, units = "secs")
    time.end <- difftime(Sys.time(),time.begin, units= "mins")
    print(paste0("Simulation: ",sim,"   Time spent: ",round(time.new, digits=2), " seconds    Total time: ", round(time.end, digits=2)," minutes"))
    
    #putting the results in the mat3 for each iterations
    for (nbl in 1:ncol(mat3)){
      if (nbl <= length(lili_totY)){
        mat3[sim,nbl] <- testdata2[[nbl]]
      } else {
        tt <- nbl-length(lili_totY)
        mat3[sim,nbl] <- obesity_prevalence2[[tt]]
      }
    }
    rm(nbl)
  }
)
##############END OF THE monte carlo
parallel::stopCluster(cl)
#mat3[]
#View(mat3)

#####PREPARING THE RESULTS
lili_totYtest <- readRDS("lili_totYtest.RData")
resultsYY2 <- data.table()
for (nbl in 1:(length(lili_totYtest))){
  resultsYY2 <- cbind(resultsYY2,  mat3[,nbl])
}
rm(nbl)
colnames(resultsYY2) <- lili_totYtest
#View(resultsYY2)

resultsYY2$sim <- seq.int(nrow(resultsYY2))

melt_data <- melt(resultsYY2, id = c("sim")) 
melt_data[,`:=` (variable=as.character(variable))]
melt_data[,`:=` (year=substring(variable,nchar(variable)-1,nchar(variable)),
                 SES=substring(variable,nchar(variable)-8,nchar(variable)-8),
                 var2=substring(variable,1,nchar(variable)-14))]
melt_data[,`:=` (year=as.numeric(as.character(year)))]
melt_data[,variable:=NULL]
dcast_data <- dcast(melt_data, sim + year + SES ~var2, mean) 

resultsYY <- dcast_data

####SAVING THE RESULTS
saveRDS(resultsYY, file=paste("results.test", toString(Sys.Date()), "rds", sep = "."))








##############################################################
##############################################################
##RESULTS
##############################################################
##############################################################
# Load results
resultsYY <- readRDS("new_results.test.2024-11-04_direct_200i_500s.rds")

#results for year 2022-2041 AS WE ARE STARTING THE POLICIES IN 2022!!!!!!!!!!!!!!!!!!!!!
resultsYYall <- resultsYY[year >=22] #for 10-year implementation [year >=22 & year <=31]
tail(resultsYYall)
#resultsYYall <- cbind(resultsYYall[,sim:DPPs_sc6sensi_label_dead_bmi_isch],resultsYYall[,expected_dead_I_chd:expected_dead_bmi_isch])

#we don't want results by year or SES
setDT(resultsYYall)
resultsYYall <- resultsYYall[, c("SES","year"):=NULL]
resultsYYall <- aggregate(resultsYYall[,2:ncol(resultsYYall)],by=list(resultsYYall$sim), sum)
setnames(resultsYYall, c("Group.1"), c("sim"))
setDT(resultsYYall)

##NEED TO DO A LOOP LATER
#DOING CVD, which is the sum of CHD+ISCH+HAEM
resultsYYall[, `:=` (expected_dead_bmi_cvd=expected_dead_bmi_chd+expected_dead_bmi_isch+expected_dead_bmi_haem,
                     expected_dead_I_cvd=expected_dead_I_chd+expected_dead_I_isch+expected_dead_I_haem,
                     DPPs_sc1_label_dead_bmi_cvd=(DPPs_sc1_label_dead_bmi_chd+DPPs_sc1_label_dead_bmi_isch+DPPs_sc1_label_dead_bmi_haem),
                     DPPs_sc1mincomp_label_dead_bmi_cvd=(DPPs_sc1mincomp_label_dead_bmi_chd+DPPs_sc1mincomp_label_dead_bmi_isch+DPPs_sc1mincomp_label_dead_bmi_haem),
                     DPPs_sc1maxcomp_label_dead_bmi_cvd=(DPPs_sc1maxcomp_label_dead_bmi_chd+DPPs_sc1maxcomp_label_dead_bmi_isch+DPPs_sc1maxcomp_label_dead_bmi_haem),
                     DPPs_sc3_label_dead_bmi_cvd=(DPPs_sc3_label_dead_bmi_chd+DPPs_sc3_label_dead_bmi_isch+DPPs_sc3_label_dead_bmi_haem),
                     DPPs_sc3mincomp_label_dead_bmi_cvd=(DPPs_sc3mincomp_label_dead_bmi_chd+DPPs_sc3mincomp_label_dead_bmi_isch+DPPs_sc3mincomp_label_dead_bmi_haem),
                     DPPs_sc3maxcomp_label_dead_bmi_cvd=(DPPs_sc3maxcomp_label_dead_bmi_chd+DPPs_sc3maxcomp_label_dead_bmi_isch+DPPs_sc3maxcomp_label_dead_bmi_haem),
                     DPPs_sc2_label_dead_bmi_cvd=(DPPs_sc2_label_dead_bmi_chd+DPPs_sc2_label_dead_bmi_isch+DPPs_sc2_label_dead_bmi_haem),
                     DPPs_sc2sensi_label_dead_bmi_cvd=(DPPs_sc2sensi_label_dead_bmi_chd+DPPs_sc2sensi_label_dead_bmi_isch+DPPs_sc2sensi_label_dead_bmi_haem),
                     DPPs_sc4_label_dead_bmi_cvd=(DPPs_sc4_label_dead_bmi_chd+DPPs_sc4_label_dead_bmi_isch+DPPs_sc4_label_dead_bmi_haem),
                     DPPs_sc4sensi_label_dead_bmi_cvd=(DPPs_sc4sensi_label_dead_bmi_chd+DPPs_sc4sensi_label_dead_bmi_isch+DPPs_sc4sensi_label_dead_bmi_haem),
                     DPPs_sc5_label_dead_bmi_cvd=(DPPs_sc5_label_dead_bmi_chd+DPPs_sc5_label_dead_bmi_isch+DPPs_sc5_label_dead_bmi_haem),
                     DPPs_sc5mincomp_label_dead_bmi_cvd=(DPPs_sc5mincomp_label_dead_bmi_chd+DPPs_sc5mincomp_label_dead_bmi_isch+DPPs_sc5mincomp_label_dead_bmi_haem),
                     DPPs_sc5maxcomp_label_dead_bmi_cvd=(DPPs_sc5maxcomp_label_dead_bmi_chd+DPPs_sc5maxcomp_label_dead_bmi_isch+DPPs_sc5maxcomp_label_dead_bmi_haem),
                     DPPs_sc6_label_dead_bmi_cvd=(DPPs_sc6_label_dead_bmi_chd+DPPs_sc6_label_dead_bmi_isch+DPPs_sc6_label_dead_bmi_haem),
                     DPPs_sc6mincomp_label_dead_bmi_cvd=(DPPs_sc6mincomp_label_dead_bmi_chd+DPPs_sc6mincomp_label_dead_bmi_isch+DPPs_sc6mincomp_label_dead_bmi_haem),
                     DPPs_sc6maxcomp_label_dead_bmi_cvd=(DPPs_sc6maxcomp_label_dead_bmi_chd+DPPs_sc6maxcomp_label_dead_bmi_isch+DPPs_sc6maxcomp_label_dead_bmi_haem),
                     DPPs_sc1sensi_label_dead_bmi_cvd=(DPPs_sc1sensi_label_dead_bmi_chd+DPPs_sc1sensi_label_dead_bmi_isch+DPPs_sc1sensi_label_dead_bmi_haem),
                     DPPs_sc3sensi_label_dead_bmi_cvd=(DPPs_sc3sensi_label_dead_bmi_chd+DPPs_sc3sensi_label_dead_bmi_isch+DPPs_sc3sensi_label_dead_bmi_haem),
                     DPPs_sc5sensi_label_dead_bmi_cvd=(DPPs_sc5sensi_label_dead_bmi_chd+DPPs_sc5sensi_label_dead_bmi_isch+DPPs_sc5sensi_label_dead_bmi_haem),
                     DPPs_sc6sensi_label_dead_bmi_cvd=(DPPs_sc6sensi_label_dead_bmi_chd+DPPs_sc6sensi_label_dead_bmi_isch+DPPs_sc6sensi_label_dead_bmi_haem),
                     DPPs_sc1sensiT_label_dead_bmi_cvd=(DPPs_sc1sensiT_label_dead_bmi_chd+DPPs_sc1sensiT_label_dead_bmi_isch+DPPs_sc1sensiT_label_dead_bmi_haem),
                     DPPs_sc3sensiT_label_dead_bmi_cvd=(DPPs_sc3sensiT_label_dead_bmi_chd+DPPs_sc3sensiT_label_dead_bmi_isch+DPPs_sc3sensiT_label_dead_bmi_haem),
                     DPPs_sc5sensiT_label_dead_bmi_cvd=(DPPs_sc5sensiT_label_dead_bmi_chd+DPPs_sc5sensiT_label_dead_bmi_isch+DPPs_sc5sensiT_label_dead_bmi_haem),
                     DPPs_sc6sensiT_label_dead_bmi_cvd=(DPPs_sc6sensiT_label_dead_bmi_chd+DPPs_sc6sensiT_label_dead_bmi_isch+DPPs_sc6sensiT_label_dead_bmi_haem),
                     
                     #BMI 
                     DPPs_sc1_SSB_sugar_dead_bmi_cvd=(DPPs_sc1_SSB_sugar_dead_bmi_chd+DPPs_sc1_SSB_sugar_dead_bmi_isch+DPPs_sc1_SSB_sugar_dead_bmi_haem),
                     DPPs_sc2_SSB_dead_bmi_cvd=(DPPs_sc2_SSB_dead_bmi_chd+DPPs_sc2_SSB_dead_bmi_isch+DPPs_sc2_SSB_dead_bmi_haem),
                     DPPs_sc3_SSB_dead_bmi_cvd=(DPPs_sc3_SSB_dead_bmi_chd+DPPs_sc3_SSB_dead_bmi_isch+DPPs_sc3_SSB_dead_bmi_haem),
                     DPPs_sc4_SSB_dead_bmi_cvd=(DPPs_sc4_SSB_dead_bmi_chd+DPPs_sc4_SSB_dead_bmi_isch+DPPs_sc4_SSB_dead_bmi_haem),
                     DPPs_sc5_SSB_dead_bmi_cvd=(DPPs_sc5_SSB_dead_bmi_chd+DPPs_sc5_SSB_dead_bmi_isch+DPPs_sc5_SSB_dead_bmi_haem),
                     DPPs_sc6_SSB_dead_bmi_cvd=(DPPs_sc6_SSB_dead_bmi_chd+DPPs_sc6_SSB_dead_bmi_isch+DPPs_sc6_SSB_dead_bmi_haem),
                     
                     DPPs_sc2sensi_SSB_dead_bmi_cvd=(DPPs_sc2sensi_SSB_dead_bmi_chd+DPPs_sc2sensi_SSB_dead_bmi_isch+DPPs_sc2sensi_SSB_dead_bmi_haem), #added by Edi
                     DPPs_sc3sensi_SSB_dead_bmi_cvd=(DPPs_sc3sensi_SSB_dead_bmi_chd+DPPs_sc3sensi_SSB_dead_bmi_isch+DPPs_sc3sensi_SSB_dead_bmi_haem),
                     DPPs_sc4sensi_SSB_dead_bmi_cvd=(DPPs_sc4sensi_SSB_dead_bmi_chd+DPPs_sc4sensi_SSB_dead_bmi_isch+DPPs_sc4sensi_SSB_dead_bmi_haem),
                     DPPs_sc5sensi_SSB_dead_bmi_cvd=(DPPs_sc5sensi_SSB_dead_bmi_chd+DPPs_sc5sensi_SSB_dead_bmi_isch+DPPs_sc5sensi_SSB_dead_bmi_haem),
                     DPPs_sc6sensi_SSB_dead_bmi_cvd=(DPPs_sc6sensi_SSB_dead_bmi_chd+DPPs_sc6sensi_SSB_dead_bmi_isch+DPPs_sc6sensi_SSB_dead_bmi_haem))]

#number of expected death
#should be the same for bmi and I pathway (I just duplicated them for model purpose)
resultsYYall[, as.list(quantile(expected_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))]
resultsYYall[, as.list(quantile(expected_dead_I_cvd, probs = c(0.5, 0.025, 0.975)))]
# resultsYYall[, as.list(quantile(expected_dead_bmi_chd, probs = c(0.5, 0.025, 0.975)))]
# resultsYYall[, as.list(quantile(expected_dead_bmi_isch, probs = c(0.5, 0.025, 0.975)))]
# resultsYYall[, as.list(quantile(expected_dead_bmi_haem, probs = c(0.5, 0.025, 0.975)))]


###############################################################
##      DPPs CVD by scenario
###############################################################
############################
## MENU LABELLING
############################
#main results
signif(resultsYYall[, as.list(quantile(DPPs_sc1_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))],2)
signif(resultsYYall[, as.list(mean(DPPs_sc1_label_dead_bmi_cvd))],2) #calculating the average instead of median for a small effect

signif(resultsYYall[, as.list(quantile(DPPs_sc3_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc2_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(mean(DPPs_sc2_label_dead_bmi_cvd))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc4_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc5_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(mean(DPPs_sc5_label_dead_bmi_cvd))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc6_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)

# % out of predicted deaths
round(resultsYYall[, as.list(quantile((DPPs_sc6_label_dead_bmi_cvd/expected_dead_bmi_cvd*100), probs = c(0.5, 0.025, 0.975)))], 2)

# sensi compensation
signif(resultsYYall[, as.list(quantile(DPPs_sc1mincomp_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(mean(DPPs_sc1mincomp_label_dead_bmi_cvd))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc1maxcomp_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(mean(DPPs_sc1maxcomp_label_dead_bmi_cvd))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc3mincomp_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(quantile(DPPs_sc3maxcomp_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc5mincomp_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(mean(DPPs_sc5mincomp_label_dead_bmi_cvd))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc5maxcomp_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(mean(DPPs_sc5maxcomp_label_dead_bmi_cvd))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc6mincomp_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(quantile(DPPs_sc6maxcomp_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)

# sensi reformulation
signif(resultsYYall[, as.list(quantile(DPPs_sc2sensi_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(mean(DPPs_sc2sensi_label_dead_bmi_cvd))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc4sensi_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)

# sensi -47 kcal
signif(resultsYYall[, as.list(quantile(DPPs_sc1sensi_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(mean(DPPs_sc1sensi_label_dead_bmi_cvd))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc3sensi_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc5sensi_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(mean(DPPs_sc5sensi_label_dead_bmi_cvd))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc6sensi_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)

# turnover
signif(resultsYYall[, as.list(quantile(DPPs_sc1sensiT_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(mean(DPPs_sc1sensiT_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc3sensiT_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(quantile(DPPs_sc5sensiT_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(quantile(DPPs_sc6sensiT_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)


############################
## SSB
############################
signif(resultsYYall[, as.list(quantile(DPPs_sc1_SSB_sugar_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(mean(DPPs_sc1_SSB_sugar_dead_bmi_cvd))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc2_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(mean(DPPs_sc2_SSB_dead_bmi_cvd))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc3_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(mean(DPPs_sc3_SSB_dead_bmi_cvd))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc4_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(mean(DPPs_sc4_SSB_dead_bmi_cvd))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc5_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(mean(DPPs_sc5_SSB_dead_bmi_cvd))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc6_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(mean(DPPs_sc6_SSB_dead_bmi_cvd))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc2sensi_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(quantile(DPPs_sc3sensi_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc4sensi_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(mean(DPPs_sc4sensi_SSB_dead_bmi_cvd))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc5sensi_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(mean(DPPs_sc5sensi_SSB_dead_bmi_cvd))], 2)

signif(resultsYYall[, as.list(quantile(DPPs_sc6sensi_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)))], 2)
signif(resultsYYall[, as.list(mean(DPPs_sc6sensi_SSB_dead_bmi_cvd))], 2)

# % out of predicted deaths
round(resultsYYall[, as.list(quantile((DPPs_sc6sensi_SSB_dead_bmi_cvd/expected_dead_bmi_cvd*100), probs = c(0.5, 0.025, 0.975)))], 2)



####################################################################################
#####################RESULTS ON OBESITY
####################################################################################
#we don't want results by ses
#resultsYYO <- cbind(resultsYY[,sim:SES],resultsYY[,expected_bmi_grp:expected_bmi_grp_sc6sensi_label])
#setDT(resultsYYO)
resultsYYO <- resultsYY[, c("SES"):=NULL]
resultsYYO <- aggregate(resultsYYO[,3:ncol(resultsYYO)],by=list(resultsYYO$sim,resultsYYO$year), mean)
setnames(resultsYYO, c("Group.1","Group.2"), c("sim","year"))
setDT(resultsYYO)

#expected obesity prevalence 
resultsYYO[, as.list(quantile(expected_bmi_grp*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
###############################################################
#difference in prevalence by scenario (-2 means a decrease in 2 percentage point in the obesity prevalence)
###############################################################
############################
## MENU LABELLING
############################
#main results
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc1_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc3_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc2_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc4_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc5_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc6_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]

# sensi compensation
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc1mincomp_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc1maxcomp_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc3mincomp_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc3maxcomp_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc5mincomp_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc5maxcomp_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc6mincomp_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc6maxcomp_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]

# sensi -47 kcal
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc1sensi_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc3sensi_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc5sensi_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc6sensi_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]

# turnover
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc1sensiT_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc3sensiT_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc5sensiT_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc6sensiT_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]


############################
## SSB
############################
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc1_SSB_sugar - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc2_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc3_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc4_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc5_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc6_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]

resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc2sensi_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc3sensi_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc4sensi_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc5sensi_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]
resultsYYO[, as.list(round(quantile((expected_bmi_grp_sc6sensi_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year")]




####################################################################################
####################################################################################
###BY SES
####################################################################################
####################################################################################
# Load results
resultsYY <- readRDS("new_results.test.2024-11-04_direct_200i_500s.rds")

#results for year 2022-2041 AS WE ARE STARTING THE POLICIES IN 2022!!!!!!!!!!!!!!!!!!!!!
resultsYYallSES <- resultsYY[year >=22] #for 10-year implementation [year >=22 & year <=31]
tail(resultsYYallSES)
#resultsYYall <- cbind(resultsYYall[,sim:DPPs_sc6sensi_label_dead_bmi_isch],resultsYYall[,expected_dead_I_chd:expected_dead_bmi_isch])

#we don't want results by year or SES
setDT(resultsYYallSES)
resultsYYallSES <- resultsYYallSES[, c("year"):=NULL]
resultsYYallSES <- aggregate(resultsYYallSES[,3:ncol(resultsYYallSES)],by=list(resultsYYallSES$sim, resultsYYallSES$SES), sum)
setnames(resultsYYallSES, c("Group.1", "Group.2"), c("sim", "SES"))
setDT(resultsYYallSES)

resultsYYallSES[, `:=` (expected_dead_bmi_cvd=expected_dead_bmi_chd+expected_dead_bmi_isch+expected_dead_bmi_haem,
                        expected_dead_I_cvd=expected_dead_I_chd+expected_dead_I_isch+expected_dead_I_haem,
                        DPPs_sc1_label_dead_bmi_cvd=(DPPs_sc1_label_dead_bmi_chd+DPPs_sc1_label_dead_bmi_isch+DPPs_sc1_label_dead_bmi_haem),
                        DPPs_sc1mincomp_label_dead_bmi_cvd=(DPPs_sc1mincomp_label_dead_bmi_chd+DPPs_sc1mincomp_label_dead_bmi_isch+DPPs_sc1mincomp_label_dead_bmi_haem),
                        DPPs_sc1maxcomp_label_dead_bmi_cvd=(DPPs_sc1maxcomp_label_dead_bmi_chd+DPPs_sc1maxcomp_label_dead_bmi_isch+DPPs_sc1maxcomp_label_dead_bmi_haem),
                        DPPs_sc3_label_dead_bmi_cvd=(DPPs_sc3_label_dead_bmi_chd+DPPs_sc3_label_dead_bmi_isch+DPPs_sc3_label_dead_bmi_haem),
                        DPPs_sc3mincomp_label_dead_bmi_cvd=(DPPs_sc3mincomp_label_dead_bmi_chd+DPPs_sc3mincomp_label_dead_bmi_isch+DPPs_sc3mincomp_label_dead_bmi_haem),
                        DPPs_sc3maxcomp_label_dead_bmi_cvd=(DPPs_sc3maxcomp_label_dead_bmi_chd+DPPs_sc3maxcomp_label_dead_bmi_isch+DPPs_sc3maxcomp_label_dead_bmi_haem),
                        DPPs_sc2_label_dead_bmi_cvd=(DPPs_sc2_label_dead_bmi_chd+DPPs_sc2_label_dead_bmi_isch+DPPs_sc2_label_dead_bmi_haem),
                        DPPs_sc4_label_dead_bmi_cvd=(DPPs_sc4_label_dead_bmi_chd+DPPs_sc4_label_dead_bmi_isch+DPPs_sc4_label_dead_bmi_haem),
                        DPPs_sc5_label_dead_bmi_cvd=(DPPs_sc5_label_dead_bmi_chd+DPPs_sc5_label_dead_bmi_isch+DPPs_sc5_label_dead_bmi_haem),
                        DPPs_sc5mincomp_label_dead_bmi_cvd=(DPPs_sc5mincomp_label_dead_bmi_chd+DPPs_sc5mincomp_label_dead_bmi_isch+DPPs_sc5mincomp_label_dead_bmi_haem),
                        DPPs_sc5maxcomp_label_dead_bmi_cvd=(DPPs_sc5maxcomp_label_dead_bmi_chd+DPPs_sc5maxcomp_label_dead_bmi_isch+DPPs_sc5maxcomp_label_dead_bmi_haem),
                        DPPs_sc6_label_dead_bmi_cvd=(DPPs_sc6_label_dead_bmi_chd+DPPs_sc6_label_dead_bmi_isch+DPPs_sc6_label_dead_bmi_haem),
                        DPPs_sc6mincomp_label_dead_bmi_cvd=(DPPs_sc6mincomp_label_dead_bmi_chd+DPPs_sc6mincomp_label_dead_bmi_isch+DPPs_sc6mincomp_label_dead_bmi_haem),
                        DPPs_sc6maxcomp_label_dead_bmi_cvd=(DPPs_sc6maxcomp_label_dead_bmi_chd+DPPs_sc6maxcomp_label_dead_bmi_isch+DPPs_sc6maxcomp_label_dead_bmi_haem),
                        DPPs_sc1sensi_label_dead_bmi_cvd=(DPPs_sc1sensi_label_dead_bmi_chd+DPPs_sc1sensi_label_dead_bmi_isch+DPPs_sc1sensi_label_dead_bmi_haem),
                        DPPs_sc3sensi_label_dead_bmi_cvd=(DPPs_sc3sensi_label_dead_bmi_chd+DPPs_sc3sensi_label_dead_bmi_isch+DPPs_sc3sensi_label_dead_bmi_haem),
                        DPPs_sc5sensi_label_dead_bmi_cvd=(DPPs_sc5sensi_label_dead_bmi_chd+DPPs_sc5sensi_label_dead_bmi_isch+DPPs_sc5sensi_label_dead_bmi_haem),
                        DPPs_sc6sensi_label_dead_bmi_cvd=(DPPs_sc6sensi_label_dead_bmi_chd+DPPs_sc6sensi_label_dead_bmi_isch+DPPs_sc6sensi_label_dead_bmi_haem),
                        DPPs_sc1sensiT_label_dead_bmi_cvd=(DPPs_sc1sensiT_label_dead_bmi_chd+DPPs_sc1sensiT_label_dead_bmi_isch+DPPs_sc1sensiT_label_dead_bmi_haem),
                        DPPs_sc3sensiT_label_dead_bmi_cvd=(DPPs_sc3sensiT_label_dead_bmi_chd+DPPs_sc3sensiT_label_dead_bmi_isch+DPPs_sc3sensiT_label_dead_bmi_haem),
                        DPPs_sc5sensiT_label_dead_bmi_cvd=(DPPs_sc5sensiT_label_dead_bmi_chd+DPPs_sc5sensiT_label_dead_bmi_isch+DPPs_sc5sensiT_label_dead_bmi_haem),
                        DPPs_sc6sensiT_label_dead_bmi_cvd=(DPPs_sc6sensiT_label_dead_bmi_chd+DPPs_sc6sensiT_label_dead_bmi_isch+DPPs_sc6sensiT_label_dead_bmi_haem),
                        
                        #BMI 
                        DPPs_sc1_SSB_sugar_dead_bmi_cvd=(DPPs_sc1_SSB_sugar_dead_bmi_chd+DPPs_sc1_SSB_sugar_dead_bmi_isch+DPPs_sc1_SSB_sugar_dead_bmi_haem),
                        DPPs_sc2_SSB_dead_bmi_cvd=(DPPs_sc2_SSB_dead_bmi_chd+DPPs_sc2_SSB_dead_bmi_isch+DPPs_sc2_SSB_dead_bmi_haem),
                        DPPs_sc3_SSB_dead_bmi_cvd=(DPPs_sc3_SSB_dead_bmi_chd+DPPs_sc3_SSB_dead_bmi_isch+DPPs_sc3_SSB_dead_bmi_haem),
                        DPPs_sc4_SSB_dead_bmi_cvd=(DPPs_sc4_SSB_dead_bmi_chd+DPPs_sc4_SSB_dead_bmi_isch+DPPs_sc4_SSB_dead_bmi_haem),
                        DPPs_sc5_SSB_dead_bmi_cvd=(DPPs_sc5_SSB_dead_bmi_chd+DPPs_sc5_SSB_dead_bmi_isch+DPPs_sc5_SSB_dead_bmi_haem),
                        DPPs_sc6_SSB_dead_bmi_cvd=(DPPs_sc6_SSB_dead_bmi_chd+DPPs_sc6_SSB_dead_bmi_isch+DPPs_sc6_SSB_dead_bmi_haem),
                        
                        DPPs_sc2sensi_SSB_dead_bmi_cvd=(DPPs_sc2sensi_SSB_dead_bmi_chd+DPPs_sc2sensi_SSB_dead_bmi_isch+DPPs_sc2sensi_SSB_dead_bmi_haem), #added by Edi
                        DPPs_sc3sensi_SSB_dead_bmi_cvd=(DPPs_sc3sensi_SSB_dead_bmi_chd+DPPs_sc3sensi_SSB_dead_bmi_isch+DPPs_sc3sensi_SSB_dead_bmi_haem),
                        DPPs_sc4sensi_SSB_dead_bmi_cvd=(DPPs_sc4sensi_SSB_dead_bmi_chd+DPPs_sc4sensi_SSB_dead_bmi_isch+DPPs_sc4sensi_SSB_dead_bmi_haem),
                        DPPs_sc5sensi_SSB_dead_bmi_cvd=(DPPs_sc5sensi_SSB_dead_bmi_chd+DPPs_sc5sensi_SSB_dead_bmi_isch+DPPs_sc5sensi_SSB_dead_bmi_haem),
                        DPPs_sc6sensi_SSB_dead_bmi_cvd=(DPPs_sc6sensi_SSB_dead_bmi_chd+DPPs_sc6sensi_SSB_dead_bmi_isch+DPPs_sc6sensi_SSB_dead_bmi_haem))]

#number of expected death
#should be the same for bmi and I pathway (I just duplicated them for model purpose)
resultsYYallSES[, as.list(signif(quantile(expected_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)), 2)),by=c("SES")]
resultsYYallSES[, as.list(quantile(expected_dead_I_cvd, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
# resultsYYall[, as.list(quantile(expected_dead_bmi_chd, probs = c(0.5, 0.025, 0.975)))]
# resultsYYall[, as.list(quantile(expected_dead_bmi_isch, probs = c(0.5, 0.025, 0.975)))]
# resultsYYall[, as.list(quantile(expected_dead_bmi_haem, probs = c(0.5, 0.025, 0.975)))]


###############################################################
##      DPPs CVD by scenario
###############################################################
############################
## MENU LABELLING
############################
#main results
resultsYYallSES[, as.list(signif(quantile(DPPs_sc1_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)), 2)), by = "SES"]
resultsYYallSES[, as.list(signif(quantile(DPPs_sc3_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)), 2)), by = "SES"]
resultsYYallSES[, as.list(signif(quantile(DPPs_sc2_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)), 2)), by = "SES"]
resultsYYallSES[, as.list(signif(quantile(DPPs_sc4_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)), 2)), by = "SES"]

resultsYYallSES[, as.list(signif(quantile(DPPs_sc5_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)), 2)), by = "SES"]
resultsYYallSES[, as.list(signif(mean(DPPs_sc5_label_dead_bmi_cvd), 2)), by = "SES"]

resultsYYallSES[, as.list(signif(quantile(DPPs_sc6_label_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)), 2)), by = "SES"]

# % out of predicted deaths
resultsYYallSES[, as.list(round(quantile((DPPs_sc6_label_dead_bmi_cvd/expected_dead_bmi_cvd)*100, probs = c(0.5, 0.025, 0.975)), 2)),by=c("SES")]


# % out of number of population
print(pop_SES)
resultsYYallSES[ pop_SES, on = c("SES"), `:=` (pop=total_pop)]
resultsYYallSES[, as.list(mean(pop, na.rm = TRUE)), by = SES]
resultsYYallSES[, `:=` (exp_rate=expected_dead_bmi_cvd/pop*100000,
                        #Combined consumer + reformulation
                        #Calorie labelling
                        DPPs_sc5_label_rate=DPPs_sc5_label_dead_bmi_cvd/pop*100000,
                        DPPs_sc6_label_rate=DPPs_sc6_label_dead_bmi_cvd/pop*100000,
                        
                        #SSB
                        DPPs_sc4sensi_SSB_rate=DPPs_sc4sensi_SSB_dead_bmi_cvd/pop*100000,
                        DPPs_sc5sensi_SSB_rate=DPPs_sc5sensi_SSB_dead_bmi_cvd/pop*100000,
                        DPPs_sc6sensi_SSB_rate=DPPs_sc6sensi_SSB_dead_bmi_cvd/pop*100000)]
                        #exp_sc6_label_rate=(expected_dead_bmi_cvd-DPPs_sc6_label_dead_bmi_cvd)/pop*100000,
                        #exp_sc6sensi_SSB_rate=(expected_dead_bmi_cvd-DPPs_sc6sensi_SSB_dead_bmi_cvd)/pop*100000)]

resultsYYallSES[, as.list(signif(quantile(exp_rate, probs = c(0.5, 0.025, 0.975)), 2)),by=c("SES")] 
resultsYYallSES[, as.list(signif(quantile(DPPs_sc5_label_rate, probs = c(0.5, 0.025, 0.975)), 2)),by=c("SES")] 
resultsYYallSES[, as.list(signif(quantile(DPPs_sc6_label_rate, probs = c(0.5, 0.025, 0.975)), 2)),by=c("SES")] 

resultsYYallSES[, as.list(signif(quantile(DPPs_sc4sensi_SSB_rate, probs = c(0.5, 0.025, 0.975)), 2)),by=c("SES")] 
resultsYYallSES[, as.list(signif(quantile(DPPs_sc5sensi_SSB_rate, probs = c(0.5, 0.025, 0.975)), 2)),by=c("SES")] 
resultsYYallSES[, as.list(signif(quantile(DPPs_sc6sensi_SSB_rate, probs = c(0.5, 0.025, 0.975)), 2)),by=c("SES")] 
#resultsYYallSES[, as.list(signif(quantile(exp_sc6_label_rate, probs = c(0.5, 0.025, 0.975)), 2)),by=c("SES")] 
#resultsYYallSES[, as.list(signif(quantile(exp_sc6sensi_SSB_rate, probs = c(0.5, 0.025, 0.975)), 2)),by=c("SES")] 

resultsYYallSES[, as.list(round(quantile(exp_rate[SES==1]/exp_rate[SES==3], probs = c(0.5, 0.025, 0.975)), 2))] 
sum((resultsYYallSES[SES == 1, exp_rate] / resultsYYallSES[SES == 3, exp_rate]) > 1)/200 # Probability of reducing absolute ineq

#Calorie labelling
resultsYYallSES[, as.list(round(quantile(DPPs_sc5_label_rate[SES==1]/DPPs_sc5_label_rate[SES==3], probs = c(0.5, 0.025, 0.975),na.rm = TRUE), 2))] 
resultsYYallSES[, as.list(round(quantile(DPPs_sc5_label_rate[SES==1]/DPPs_sc5_label_rate[SES==3 & DPPs_sc5_label_rate > 0], probs = c(0.5, 0.025, 0.975),na.rm = TRUE), 2))] 

ratios <- resultsYYallSES[SES == 1, DPPs_sc5_label_rate] / resultsYYallSES[SES == 3, DPPs_sc5_label_rate]
valid_ratios <- ratios[is.finite(ratios)]

quantiles <- quantile(valid_ratios, probs = c(0.5, 0.025, 0.975), na.rm = TRUE)
quantiles_rounded <- round(quantiles, 2)
as.list(quantiles_rounded) # 95% UI are same

length(valid_ratios)
sum(valid_ratios > 1)/length(valid_ratios)

resultsYYallSES[, as.list(round(quantile(DPPs_sc6_label_rate[SES==1]/DPPs_sc6_label_rate[SES==3], probs = c(0.5, 0.025, 0.975),na.rm = TRUE), 2))] 
resultsYYallSES[, as.list(round(quantile(DPPs_sc6_label_rate[SES==1]/DPPs_sc6_label_rate[SES==3 & DPPs_sc6_label_rate > 0], probs = c(0.5, 0.025, 0.975),na.rm = TRUE), 2))] 
sum(
  (resultsYYallSES[SES == 1, DPPs_sc6_label_rate] / resultsYYallSES[SES == 3, DPPs_sc6_label_rate]) > 1 &
    is.finite(resultsYYallSES[SES == 1, DPPs_sc6_label_rate] / resultsYYallSES[SES == 3, DPPs_sc6_label_rate])
) / 200

ratios <- resultsYYallSES[SES == 1, DPPs_sc6_label_rate] / resultsYYallSES[SES == 3, DPPs_sc6_label_rate]
valid_ratios <- ratios[is.finite(ratios)]

quantiles <- quantile(valid_ratios, probs = c(0.5, 0.025, 0.975), na.rm = TRUE)
quantiles_rounded <- round(quantiles, 2)
as.list(quantiles_rounded)  # 95% UI are same

length(valid_ratios)
sum(valid_ratios > 1)/length(valid_ratios)

#SSB
resultsYYallSES[, as.list(round(quantile(DPPs_sc4sensi_SSB_rate[SES==1]/DPPs_sc4sensi_SSB_rate[SES==3], probs = c(0.5, 0.025, 0.975), na.rm = TRUE), 2))] 
resultsYYallSES[, as.list(round(quantile(DPPs_sc4sensi_SSB_rate[SES==1]/DPPs_sc4sensi_SSB_rate[SES==3 & DPPs_sc4sensi_SSB_rate > 0], probs = c(0.5, 0.025, 0.975),na.rm = TRUE), 2))] 
sum(
  (resultsYYallSES[SES == 1, DPPs_sc4sensi_SSB_rate] / resultsYYallSES[SES == 3, DPPs_sc4sensi_SSB_rate]) > 1 &
    is.finite(resultsYYallSES[SES == 1, DPPs_sc4sensi_SSB_rate] / resultsYYallSES[SES == 3, DPPs_sc4sensi_SSB_rate])
) / 200

ratios <- resultsYYallSES[SES == 1, DPPs_sc4sensi_SSB_rate] / resultsYYallSES[SES == 3, DPPs_sc4sensi_SSB_rate]
valid_ratios <- ratios[is.finite(ratios)]

quantiles <- quantile(valid_ratios, probs = c(0.5, 0.025, 0.975), na.rm = TRUE)
quantiles_rounded <- round(quantiles, 2)
as.list(quantiles_rounded) # 95% UI is litte bit different

length(valid_ratios)
sum(valid_ratios > 1)/length(valid_ratios)


resultsYYallSES[, as.list(round(quantile(DPPs_sc5sensi_SSB_rate[SES==1]/DPPs_sc5sensi_SSB_rate[SES==3], probs = c(0.5, 0.025, 0.975), na.rm = TRUE), 2))] 
resultsYYallSES[, as.list(round(quantile(DPPs_sc5sensi_SSB_rate[SES==1]/DPPs_sc5sensi_SSB_rate[SES==3 & DPPs_sc5sensi_SSB_rate > 0], probs = c(0.5, 0.025, 0.975),na.rm = TRUE), 2))] 
sum(
  (resultsYYallSES[SES == 1, DPPs_sc5sensi_SSB_rate] / resultsYYallSES[SES == 3, DPPs_sc5sensi_SSB_rate]) > 1 &
    is.finite(resultsYYallSES[SES == 1, DPPs_sc5sensi_SSB_rate] / resultsYYallSES[SES == 3, DPPs_sc5sensi_SSB_rate])
) / 200

ratios <- resultsYYallSES[SES == 1, DPPs_sc5sensi_SSB_rate] / resultsYYallSES[SES == 3, DPPs_sc5sensi_SSB_rate]
valid_ratios <- ratios[is.finite(ratios)]

quantiles <- quantile(valid_ratios, probs = c(0.5, 0.025, 0.975), na.rm = TRUE)
quantiles_rounded <- round(quantiles, 2)
as.list(quantiles_rounded) # 95% UI is litte bit different

length(valid_ratios)
sum(valid_ratios > 1)/length(valid_ratios)

resultsYYallSES[, as.list(round(quantile(DPPs_sc6sensi_SSB_rate[SES==1]/DPPs_sc6sensi_SSB_rate[SES==3], probs = c(0.5, 0.025, 0.975), na.rm = TRUE), 2))] 
resultsYYallSES[, as.list(round(quantile(DPPs_sc6sensi_SSB_rate[SES==1]/DPPs_sc6sensi_SSB_rate[SES==3 & DPPs_sc6sensi_SSB_rate > 0], probs = c(0.5, 0.025, 0.975),na.rm = TRUE), 2))] 
sum(
  (resultsYYallSES[SES == 1, DPPs_sc6sensi_SSB_rate] / resultsYYallSES[SES == 3, DPPs_sc6sensi_SSB_rate]) > 1 &
    is.finite(resultsYYallSES[SES == 1, DPPs_sc6sensi_SSB_rate] / resultsYYallSES[SES == 3, DPPs_sc6sensi_SSB_rate])
) / 200

ratios <- resultsYYallSES[SES == 1, DPPs_sc6sensi_SSB_rate] / resultsYYallSES[SES == 3, DPPs_sc6sensi_SSB_rate]
valid_ratios <- ratios[is.finite(ratios)]

quantiles <- quantile(valid_ratios, probs = c(0.5, 0.025, 0.975), na.rm = TRUE)
quantiles_rounded <- round(quantiles, 2)
as.list(quantiles_rounded) # 95% UI is litte bit different

length(valid_ratios)
sum(valid_ratios > 1)/length(valid_ratios)

#resultsYYallSES[, as.list(signif(quantile(exp_sc6_label_rate[SES==1]/exp_sc6_label_rate[SES==3], probs = c(0.5, 0.025, 0.975), na.rm = TRUE), 2))] 
#resultsYYallSES[, as.list(signif(quantile(exp_sc6sensi_SSB_rate[SES==1]/exp_sc6sensi_SSB_rate[SES==3], probs = c(0.5, 0.025, 0.975), na.rm = TRUE), 2))] 


all.equal(resultsYYallSES[SES == 1, DPPs_sc6_label_rate],
          resultsYYallSES[, DPPs_sc6_label_rate[SES == 1]])
resultsYYallSES[ (DPPs_sc6_label_rate[SES == 1]) > 0, ]


############################
## SSB
############################
resultsYYallSES[, as.list(signif(quantile(DPPs_sc1_SSB_sugar_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)), 2)), by = "SES"]
resultsYYallSES[, as.list(signif(quantile(DPPs_sc2_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)), 2)), by = "SES"]
resultsYYallSES[, as.list(signif(quantile(DPPs_sc3_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)), 2)), by = "SES"]
resultsYYallSES[, as.list(signif(quantile(DPPs_sc4_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)), 2)), by = "SES"]
resultsYYallSES[, as.list(signif(quantile(DPPs_sc5_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)), 2)), by = "SES"]
resultsYYallSES[, as.list(signif(quantile(DPPs_sc6_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)), 2)), by = "SES"]

resultsYYallSES[, as.list(signif(quantile(DPPs_sc2sensi_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)), 2)), by = "SES"]
resultsYYallSES[, as.list(signif(quantile(DPPs_sc3sensi_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)), 2)), by = "SES"]
resultsYYallSES[, as.list(signif(quantile(DPPs_sc4sensi_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)), 2)), by = "SES"]
resultsYYallSES[, as.list(signif(quantile(DPPs_sc5sensi_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)), 2)), by = "SES"]
resultsYYallSES[, as.list(signif(quantile(DPPs_sc6sensi_SSB_dead_bmi_cvd, probs = c(0.5, 0.025, 0.975)), 2)), by = "SES"]

# % out of predicted deaths
resultsYYallSES[, as.list(round(quantile((DPPs_sc6sensi_SSB_dead_bmi_cvd/expected_dead_bmi_cvd)*100, probs = c(0.5, 0.025, 0.975)), 2)),by=c("SES")]



####################################################################################
#####################RESULTS ON OBESITY
####################################################################################
resultsYYOSES <- resultsYY[year >=22] #for 10-year implementation [year >=22 & year <=31]
tail(resultsYYOSES)
resultsYYOSES <- aggregate(resultsYYOSES[,4:ncol(resultsYYOSES)],by=list(resultsYYOSES$sim,resultsYYOSES$year, resultsYYOSES$SES), mean)
setnames(resultsYYOSES, c("Group.1","Group.2", "Group.3"), c("sim","year", "SES"))
setDT(resultsYYOSES)

#expected obesity prevalence for most deprived (SES==1) and least deprived (SES==5) 
resultsYYOSES[, as.list(quantile(expected_bmi_grp*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
###############################################################
#difference in prevalence by scenario (-2 means a decrease in 2 percentage point in the obesity prevalence)
###############################################################
############################
## MENU LABELLING
############################
#main results
resultsYYOSES[, as.list(round(quantile((expected_bmi_grp_sc1_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year", "SES")]
resultsYYOSES[, as.list(round(quantile((expected_bmi_grp_sc3_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year", "SES")]
resultsYYOSES[, as.list(round(quantile((expected_bmi_grp_sc2_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year", "SES")]
resultsYYOSES[, as.list(round(quantile((expected_bmi_grp_sc4_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year", "SES")]
resultsYYOSES[, as.list(round(quantile((expected_bmi_grp_sc5_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year", "SES")]
resultsYYOSES[, as.list(round(quantile((expected_bmi_grp_sc6_label - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year", "SES")]

############################
## SSB
############################
resultsYYOSES[, as.list(round(quantile((expected_bmi_grp_sc1_SSB_sugar - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year", "SES")]
resultsYYOSES[, as.list(round(quantile((expected_bmi_grp_sc2_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year", "SES")]
resultsYYOSES[, as.list(round(quantile((expected_bmi_grp_sc3_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year", "SES")]
resultsYYOSES[, as.list(round(quantile((expected_bmi_grp_sc4_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year", "SES")]
resultsYYOSES[, as.list(round(quantile((expected_bmi_grp_sc5_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year", "SES")]
resultsYYOSES[, as.list(round(quantile((expected_bmi_grp_sc6_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year", "SES")]

resultsYYOSES[, as.list(round(quantile((expected_bmi_grp_sc2sensi_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year", "SES")]
resultsYYOSES[, as.list(round(quantile((expected_bmi_grp_sc3sensi_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year", "SES")]
resultsYYOSES[, as.list(round(quantile((expected_bmi_grp_sc4sensi_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year", "SES")]
resultsYYOSES[, as.list(round(quantile((expected_bmi_grp_sc5sensi_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year", "SES")]
resultsYYOSES[, as.list(round(quantile((expected_bmi_grp_sc6sensi_SSB - expected_bmi_grp) * 100, probs = c(0.5, 0.025, 0.975)), 2)), by = c("year", "SES")]
######end 






