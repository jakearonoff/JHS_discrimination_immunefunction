library("car")
library("betareg")
library("stargazer")

d <- read.csv("path_to_jhs_data.csv", header = T, check.names = F, sep = ",") # reading JHS data (demographics, 
#                                                                               discrimination, cell counts)

d$nlr <- d$SEGS/d$LYMPHS # generating NLR (SEGS=Neutrophils)
d$lognlr <- log(d$nlr) # highly skewed, so log-transforming 
d$lognlrz <- scale(d$lognlr) # standardizing (mean = 0, SD = 1)
d$midlifetime <- ifelse(d$lifetimeDiscrm == 3 | d$lifetimeDiscrm == 4,1,0) # categorizing into middle lifetime discrim.
d$highlifetime <- ifelse(d$lifetimeDiscrm > 4,1,0) # categorizing into high lifetime discrimination
d$middaily <- ifelse(d$dailyDiscr > 1 & d$dailyDiscr <= 3,1,0) # categorizing into less frequent daily discrimination
d$highdaily <- ifelse(d$dailyDiscr > 3,1,0) # categorizing into more frequent daily discrimination
d$someburden <- ifelse(d$discrmBurden >= 1 & d$discrmBurden <= 2.5,1,0) # categorizing into "some burden"
d$highburden <- ifelse(d$discrmBurden > 2.5,1,0) # categorizing into "high burden"

# stargazer code = generating table for descriptive stats of sample 
stargazer(d[c("age","female", "dailyDiscr","lifetimeDiscrm","discrmBurden","SEGS","LYMPHS","MONOS")], type = "html",
          out = "filepath.html", digits = 2, title = "Select Descriptive Stats",
          covariate.labels = c("Age","Female", "Daily Discrimination","Lifetime Discrimination","Discrimination Burden","Neutrophils","Lymphocytes","Monocytes"),
          summary.stat = c("mean","sd","min","max"))


#####################################################################
## innate vs adaptive cells analysis 
####################################################################

segs <- betareg(SEGS ~ midburden + highburden + middaily + highdaily + midlifetime + highlifetime + age + female + currentSmoker + 
                    alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                    edu3cat + statinMeds + antiArythMeds + 
                    Income + CVDHx, data = d)
summary(segs)
vif(segs) # checking VIF using car package to assess colinearity 


lymphs <- betareg(LYMPHS ~ midburden + highburden + middaily + highdaily + midlifetime + highlifetime + age + female + currentSmoker + 
                      alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                      edu3cat + statinMeds + antiArythMeds + 
                      Income + CVDHx, data = d)
summary(lymphs)
vif(lymphs)


monos <- betareg(MONOS ~ midburden + highburden + middaily + highdaily + midlifetime + highlifetime + age + female + currentSmoker + 
                     alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                     edu3cat + statinMeds + antiArythMeds + 
                     Income + CVDHx, data = d)
summary(monos)
vif(monos)


nlr <- lm(lognlrz ~ midburden + highburden + middaily + highdaily + midlifetime + highlifetime + age + female + currentSmoker + 
                   alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                   edu3cat + statinMeds + antiArythMeds + 
                   Income + CVDHx, data = d)
summary(nlr)
vif(nlr)

# stargazer code generating a table with model outputs 
stargazer(segs, lymphs, monos, nlr,  type = "html", 
          out = "filepath.html", 
          title = "Regression Models Predicting Innate and Adaptive Cells",
          covariate.labels = c("Mid Burden","High Burden","Mid Daily","High Daily","Mid Lifetime",
                               "High Lifetime","Age","Female","Current Smoker?",
                               "Weekly Alc. Cons.","Waist Circ.",
                               "Diastolic BP","LDL","HDL","Trigs","eGFRmdrd", "Education",
                               "Statin Meds?","Anti-Aryth. Meds?","Income",
                               "CVD history?"),
          dep.var.labels = c("Neutrophils", "Lymphocytes", "Monocytes", "NLR"),
          dep.var.caption = "Coefficients with 95% CI",
          digits = 2, ci = T,
          star.char = c("t","*","**"))


########################################################################

## Analysis using estimates of lymphocyte %'s and immune system aging from methylation data

########################################################################

m <- read.csv("path_to_data_file.csv", header = T, check.names = F, sep = ",") # reading methylation-derived data 
#                                                                                (lymphocyte %'s, immune system aging)
m <- merge(d,m, by = "subjid") # merging with other dataset 

m$CD4.naivez <- scale(m$CD4.naive) # standardizing naive CD4 estimates (mean = 0, SD = 1)  
m$CD8.naivez <- scale(m$CD8.naive) # standardizing naive CD8 estimates (mean = 0, SD = 1)  
m$PlasmaBlastz <- scale(m$PlasmaBlast) # standardizing plasma blast estimates (mean = 0, SD = 1)  
m$CD8pCD28nCD45RAnz <- scale(m$CD8pCD28nCD45RAn) # standardizing exhausted CD8T cells estimates (mean = 0, SD = 1) 

# stargazer code generating table with descriptive stats for the subsample 
stargazer(m[,c("age","female", "dailyDiscr","lifetimeDiscrm","discrmBurden","SEGS","LYMPHS","MONOS", "CD4T", "CD8T", 
               "Bcell","CD4.naive","CD8.naive", "PlasmaBlast","CD8pCD28nCD45RAn")], type = "html",
          out = "filepath.html", digits = 3, title = "Select Descriptive Stats",
          covariate.labels = c("Age","Female", "Daily Discrimination","Lifetime Discrimination","Discrimination Burden",
                               "Neutrophils","Lymphocytes","Monocytes","CD4T %","CD8T %","B cell %","CD4 Naive",
                               "CD8 Naive","Plasma Blast","Exhausted T cells"),
          summary.stat = c("mean","sd","min","max"))


# Beta regression models predicting %'s of CD4T, CD8T and B cells
cd4t <- betareg(CD4T ~ midburden + highburden + middaily + highdaily + midlifetime + highlifetime + age + female + 
                  currentSmoker + alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                   edu3cat + statinMeds + antiArythMeds + 
                   Income + CVDHx, data = m)
summary(cd4t)


cd8t <- betareg(CD8T ~ midburden + highburden + middaily + highdaily + midlifetime + highlifetime + age + female + 
                  currentSmoker + alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                   edu3cat + statinMeds + antiArythMeds + 
                   Income + CVDHx, data = m)
summary(cd8t)


bcell <- betareg(Bcell ~ midburden + highburden + middaily + highdaily + midlifetime + highlifetime + age + female + 
                   currentSmoker + alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                    edu3cat + statinMeds + antiArythMeds + 
                    Income + CVDHx, data = m)
summary(bcell)

# stargazer code generating table for beta regression model outputs
stargazer(cd4t,cd8t,bcell,  type = "html", 
          out = "filepath.html", 
          title = "Beta Regression Models Predicting Specific Lymphocytes",
          covariate.labels = c("Mid Burden","High Burden","Mid Daily","High Daily","Mid Lifetime",
                               "High Lifetime","Age","Female","Current Smoker?",
                               "Weekly Alc. Cons.","Waist Circ.",
                               "Diastolic BP","LDL","HDL","Trigs","eGFRmdrd", "Education",
                               "Statin Meds?","Anti-Aryth. Meds?","Income",
                               "CVD history?"),
          dep.var.labels = c("CD4T", "CD8T", "B Cells"),
          dep.var.caption = "Coefficients with 95% CI",
          digits = 2, ci = T,
          star.char = c("t","*","**"))


# Regression models predicting markers of immune system aging 
cd4naive <- lm(CD4.naivez ~ midburden + highburden + middaily + highdaily + midlifetime + highlifetime + age + 
                 female + currentSmoker + alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                  edu3cat + statinMeds + antiArythMeds + 
                  Income + CVDHx, data = m)
summary(cd4naive)


cd8naive <- lm(CD8.naivez ~ midburden + highburden + middaily + highdaily + midlifetime + highlifetime + age + 
                 female + currentSmoker + alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                  edu3cat + statinMeds + antiArythMeds + 
                  Income + CVDHx, data = m)
summary(cd8naive)


exhaustedt <- lm(CD8pCD28nCD45RAnz ~ midburden + highburden + middaily + highdaily + midlifetime + highlifetime + 
                   age + female + currentSmoker + alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                    edu3cat + statinMeds + antiArythMeds + 
                    Income + CVDHx, data = m)
summary(exhaustedt)


plasmab <- lm(PlasmaBlastz ~ midburden + highburden + middaily + highdaily + midlifetime + highlifetime + age + 
                female + currentSmoker + alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                 edu3cat + statinMeds + antiArythMeds + 
                 Income + CVDHx, data = m)
summary(plasmab)

# stargazer code generating table with regression model outputs 
stargazer(cd4naive,cd8naive,exhaustedt, plasmab,  type = "html", 
          out = "filepath.html", 
          title = "OLS Regression Models Predicting Immune Aging",
          covariate.labels = c("Mid Burden","High Burden","Mid Daily","High Daily","Mid Lifetime",
                               "High Lifetime","Age","Female","Current Smoker?",
                               "Weekly Alc. Cons.","Waist Circ.",
                               "Diastolic BP","LDL","HDL","Trigs","eGFRmdrd", "Education",
                               "Statin Meds?","Anti-Aryth. Meds?","Income",
                               "CVD history?"),
          dep.var.labels = c("CD4 Naive", "CD8 Naive", "Exhausted CD8T cells", "Plasma Blasts"),
          dep.var.caption = "Coefficients with 95% CI",
          digits = 2, ci = T,
          star.char = c("t","*","**"))





