library("tidyverse")
library("car") # for checking model VIF
library("betareg") # for implementing beta regression models 
library("stargazer") # for outputting model results into tables 
library("jtools") # for plotting model results


########################################################################
## example of modeling steps for cell counts and NLR, shown here with neutrophils  
########################################################################
neutcountz_unadj <- lm(neutcountz ~ midburden + highburden + middaily + 
                         highdaily + midlifetime + 
                         highlifetime, data = d)
summary(neutcountz_unadj)
vif(neutcountz_unadj)

neutcountz <- lm(neutcountz ~ midburden + highburden + middaily + 
                   highdaily + midlifetime + 
                   highlifetime + age + male + 
                   currentSmoker + alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                   edu3cat + statinMeds + antiArythMeds + 
                   Income + CVDHx, data = d)
summary(neutcountz)
vif(neutcountz)

neutcountz_sex <- lm(neutcountz ~ middaily*male + 
                       highdaily*male + midlifetime*male + 
                       highlifetime*male + midburden*male + highburden*male + age + 
                       currentSmoker + alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                       edu3cat + statinMeds + antiArythMeds + 
                       Income + CVDHx, data = d)
summary(neutcountz_sex)

neutcountz_age <- lm(neutcountz ~ middaily*age + 
                       highdaily*age + midlifetime*age + 
                       highlifetime*age + midburden*age + highburden*age + male + 
                       currentSmoker + alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                       edu3cat + statinMeds + antiArythMeds + 
                       Income + CVDHx, data = d)
summary(neutcountz_age)


########################################################################
## Figure 2  
########################################################################

# specify how the variables should be presented in the plot 
coef_names <- c("Less Frequent Everyday" = "middaily", "More Frequent Everyday" = "highdaily", 
                "Middle Lifetime" = "midlifetime", "Highest Lifetime" = "highlifetime", 
                "Some Burden from Lifetime" = "midburden", "High Burden from Lifetime" = "highburden")

# generate and store the plot
plot <- plot_summs(lymphcount, neutcountz, monocountz, lognlrz, coefs = coef_names, 
                   model.names = c("Lymphocytes", "Neutrophils", "Monocytes", "NLR"),
                   colors = c("CUD"))
# plot settings
apatheme=theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family='Helvetica'),
        legend.title=element_blank(), 
        axis.text=element_text(size=15),
        axis.title=element_text(size = 15),
        legend.text = element_text(size = 15))
# plot
plot  + apatheme + labs(x = "\n Estimates with 95% CI \n ", y = NULL) 


########################################################################
## Example stargazer code for putting model results into tables. I saved the tables in HTML and opened/edited in Word
########################################################################
stargazer(lymphcount, neutcountz,  type = "html", 
          out = "~/Downloads/jhs_lymph_neut.html", 
          title = "OLS Regression Models Predicting Cell Counts",
          digits = 2, ci = T,
          star.char = c("t","*","**"))


########################################################################
## analysis and plot for Figure 3
########################################################################
oldermen <- filter(d, age > 53 & female == 0)
olderwomen <- filter(d, age > 53 & female == 1)
youngerwomen <- filter(d, age < 54 & female == 1)
youngermen <- filter(d, age < 54 & female == 0)

lymphcountom <- lm(lymphcountz ~ midburden + highburden + middaily + 
                   highdaily + midlifetime + 
                   highlifetime + age + 
                   currentSmoker + alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                   edu3cat + statinMeds + antiArythMeds + 
                   Income + CVDHx, data = oldermen)
summary(lymphcountom)

lymphcountow <- lm(lymphcountz ~ midburden + highburden + middaily + 
                     highdaily + midlifetime + 
                     highlifetime + age + 
                     currentSmoker + alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                     edu3cat + statinMeds + antiArythMeds + 
                     Income + CVDHx, data = olderwomen)
summary(lymphcountow)

lymphcountyw <- lm(lymphcountz ~ midburden + highburden + middaily + 
                     highdaily + midlifetime + 
                     highlifetime + age + 
                     currentSmoker + alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                     edu3cat + statinMeds + antiArythMeds + 
                     Income + CVDHx, data = youngerwomen)
summary(lymphcountyw)

lymphcountym <- lm(lymphcountz ~ midburden + highburden + middaily + 
                     highdaily + midlifetime + 
                     highlifetime + age + 
                     currentSmoker + alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                     edu3cat + statinMeds + antiArythMeds + 
                     Income + CVDHx, data = youngermen)
summary(lymphcountym)


plot <- plot_summs(lymphcountom, lymphcountym, lymphcountow, lymphcountyw, coefs = coef_names, 
                   model.names = c("Older Men (n=567)", "Younger Men (n=638)", "Older Women (n=1,094)", "Younger Women (n=1,020)"),
                   colors = c("Qual1"))

plot  + apatheme + labs(x = "\n Estimates with 95% CI \n ", y = NULL)



########################################################################
## Example code for beta regression models predicting cell proportions, here CD4T. 
## The code for modeling is identical as with lm(), just using betareg() function instead. 
########################################################################

cd4_unadj <- betareg(CD4T ~ middaily + 
                 highdaily + midlifetime + 
                 highlifetime + midburden + highburden, data = m)
summary(cd4_unadj)
vif(cd4_unadj)

cd4 <- betareg(CD4T ~ middaily + 
                    highdaily + midlifetime + 
                    highlifetime + midburden + highburden + age + male + 
                    currentSmoker + alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                    edu3cat + statinMeds + antiArythMeds + 
                    Income + CVDHx, data = m)
summary(cd4)
vif(cd4)

cd4_sex <- betareg(CD4T ~ middaily*male + 
                 highdaily*male + midlifetime*male + 
                 highlifetime*male + midburden*male + highburden*male + age + 
                 currentSmoker + alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                 edu3cat + statinMeds + antiArythMeds + 
                 Income + CVDHx, data = m)
summary(cd4_sex)

cd4_age <- betareg(CD4T ~ middaily*age + 
                     highdaily*age + midlifetime*age + 
                     highlifetime*age + midburden*age + highburden*age + male + 
                     currentSmoker + alcw + waist + dbp + ldl + hdl + trigs + eGFRmdrd +
                     edu3cat + statinMeds + antiArythMeds + 
                     Income + CVDHx, data = m)
summary(cd4_age)



