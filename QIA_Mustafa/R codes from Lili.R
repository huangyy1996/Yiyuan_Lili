
#=========================================================================
# ANOVA model 
#  see  https://rcompanion.org/handbook/I_08.html
#  see http://www.sthda.com/english/wiki/two-way-anova-test-in-r;  
#==============================================================================


data=read.csv("C:/Users/zhaolili/Downloads/OJOC May/bone2.csv")

# First: make categorical predictors as factors

data$group=as.factor(data$group)
data$gender=as.factor(data$gender)
model = lm(y ~ group + gender + group:gender,  data = data)

# obtain parameter estiamtes
summary(model)  

# get the type III table as in SAS
anova(model)

# residual plot
plot(model)



## Alternative functions
# run two-way ANOVA model with interaction
res.aov <- aov(y ~ group + gender + group:gender, data = data)
summary(res.aov)

#Tukey multiple pairwise-comparisons  (Multiple comparisons using "multcomp" package or "lsmeans" package)
TukeyHSD(res.aov, which = "group")

library(lsmeans)
lsmeans(res.aov, ~ group:gender, adjust="tukey")       ### Tukey-adjusted comparisons


#======================================================================
# ICC to assess agreeemnt for continuous ata
#======================================================================
library(irr)
library(dplyr)
library(psych)

# row is observation and colomn (humnan vs machine)
icc(ddm, model = "twoway", 
    type = "agreement", unit = "single")


#======================================================================
# kappa to assess agreeemnt for categorical data
#======================================================================
cohen.kappa 