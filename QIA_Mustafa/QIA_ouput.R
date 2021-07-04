#-------------------------QIA project-----------------------
#Dataset: QIA project final-OCT-2020.csv
#Author: Yiyuan Huang
#Date: Oct 5, 2020
#-----------------------------------------------------------
library(dplyr)
#library(qwraps2)
library(arsenal)
library(irr)
library(psych)
library(ggpubr)
library(multcompView)
library(lsmeans)
library(ggplot2)
library(DescTools)

#----------------------------READ DATA--------------------------
QIA <- read.csv("E:/research_lili/QIA/data/QIA project final-OCT-2020.csv", header = T)

#===============================Aim 1====================================================
#Aim 1: Summary statistics for Demographics features:
#Variables: AGE, DIAGNOSIS, ER_HUMAN, PgR_HUMAN, HER2_HUMAN. 
QIA$DIAGNOSIS <- as.factor(QIA$DIAGNOSIS)
df_summary_cont <- QIA %>% select(AGE, ER_HUMAN, PgR_HUMAN, HER2_HUMAN)

#-----------------------------table1-----------------------------------------
my_controls <- tableby.control(
  numeric.stats = c("meansd", "medianq1q3", "range", "Nmiss2"),
  cat.stats = c("countpct", "Nmiss2"),
  stats.labels = list(
    meansd = "Mean (SD)",
    medianq1q3 = "Median (Q1, Q3)",
    range = "Min - Max",
    Nmiss2 = "Missing"
  )
)

table1 <- tableby(~AGE+ ER_HUMAN+ PgR_HUMAN+ HER2_HUMAN+ DIAGNOSIS, data=QIA, control = my_controls)
summary(table1)

ER_HUMAN_CAT = as.factor(QIA$ER_HUMAN>=0.01)
PgR_HUMAN_CAT = as.factor(QIA$PgR_HUMAN>=0.01)
QIA <- QIA %>% mutate(ER_HUMAN_CAT, PgR_HUMAN_CAT) 
levels(QIA$ER_HUMAN_CAT) <- c("Negative", "Positive")
levels(QIA$PgR_HUMAN_CAT) <- c("Negative", "Positive")

my_controls2 <- tableby.control(
  cat.stats = c("countpct", "Nmiss2")
)
table2 <- tableby(~ER_HUMAN_CAT + PgR_HUMAN_CAT, data = QIA, control = my_controls2)
summary(table2)

#===============================AIM 2====================================================
#Aim 2: Correlation between the Human score of (ER_HUMAN, PgR_HUMAN, HER2_HUMAN) 
#       and the Digital image analysis score (ER_COREO, PgR_COREO, HER2_COREO)

#ER_HUMAN vs ER_COREO (both cont --- icc)
corr_ER <- icc(QIA %>% select(ER_HUMAN, ER_COREO), model = "twoway", type = "agreement", unit = "single")

#PgR_HUMAN vs PgR_COREO (both cont --- icc)
corr_PgR <- icc(QIA %>% select(PgR_HUMAN, PgR_COREO), model = "twoway", type = "agreement", unit = "single")

#HER2_HUMAN vs HER2_COREO (both cat --- kappa)
QIA$HER2_COREO <- as.factor(QIA$HER2_COREO)
corr_HER2 <- cohen.kappa(QIA %>% select(HER2_HUMAN, HER2_COREO)) 

StuartTauC(as.numeric(QIA$HER2_HUMAN), as.numeric(QIA$HER2_COREO), conf.level = 0.95)


#===============================AIM 3=================================================
#Aim 3A: Correlation between HER2_COREO (Negative, equivocal and Positive) 
#        and the HER2_FISH (Negative and Positive). 
QIA$HER2_FISH[which(QIA$HER2_FISH=="")] = NA
QIA$HER2_FISH <- as.factor(QIA$HER2_FISH)
QIA$HER2_COREO_CAT <- cut(QIA$HER2_COREO, breaks = c(0,1.5,2.5,3.5),
                          labels = c("Negative", "Equivocal", "Positive"))

table(QIA$HER2_FISH, QIA$HER2_COREO_CAT)

#         Negative Equivocal Positive
#Negative       15       329        3
#Positive        1        23        1

#corr_COREO_FISH <- cohen.kappa(QIA %>% select(HER2_COREO_CAT, HER2_FISH))
StuartTauC(as.numeric(QIA$HER2_FISH), as.numeric(QIA$HER2_COREO_CAT), conf.level = 0.95)


#Aim 3B: Correlation between HER2_COREO (cat) and HER2D17Z1_RATIO (cont).---ANOVA
QIA$LOG_RATIO <- log(QIA$HER2D17Z1_RATIO)
QIA$LOG_COPYNUMBER <- log(QIA$HER2_COPYNUMBER)

#boxplot
ggboxplot(QIA, x = "HER2_COREO_CAT", y = "LOG_RATIO")

#type III anova 
fit_3b <- lm(HER2D17Z1_RATIO ~ HER2_COREO_CAT, data = QIA)
#summary(fit_3b)
p_3b <- anova(fit_3b)[1,5] #p-value for F-test
res_3b <- plot(fit_3b) #model checking plots

#post-hoc analysis
marginal = lsmeans(fit_3b, ~ HER2_COREO_CAT)
pairs(marginal, adjust="tukey")

#Aim 3C: Correlation between HER2_COREO (cat) and HER2_COPYNUMBER (cont). ---ANOVA
ggboxplot(QIA, x = "HER2_COREO_CAT", y = "LOG_COPYNUMBER")

fit_3c <- lm(HER2_COPYNUMBER ~ HER2_COREO_CAT, data = QIA)

p_3b <- anova(fit_3c)[1,5] #p-value for F-test
res_3b <- plot(fit_3c) #model checking plots

marginal = lsmeans(fit_3c, ~ HER2_COREO_CAT)
pairs(marginal, adjust="tukey")


#Aim 3D: Correlation between DIAGNOSIS (CAT) 
#        and HER2_FISH (CAT)/ASCOCAP_Giudline_Groups (CAT).
QIA$ASCOCAP_Giudline_Groups[which(QIA$ASCOCAP_Giudline_Groups == "")] = NA
QIA$ASCOCAP_Giudline_Groups <- as.factor(QIA$ASCOCAP_Giudline_Groups)

corr_DIAG_FISH <- cohen.kappa(QIA %>% select(DIAGNOSIS, HER2_FISH))
corr_DIAG_ASCO <- cohen.kappa(QIA %>% select(DIAGNOSIS, ASCOCAP_Giudline_Groups))

#=======================AIM 4========================================================
#Aim 4A: Discordant between the Human score of (ER_HUMAN, PgR_HUMAN, HER2_HUMAN) 
#and the Digital image analysis score (ER_COREO, PgR_COREO, HER2_COREO). 

#SAME AS BEFORE: CATEGORIZE COREO SCORES
ER_COREO_CAT = as.factor(QIA$ER_COREO>=0.01)
PgR_COREO_CAT = as.factor(QIA$PgR_COREO>=0.01)
QIA <- QIA %>% mutate(ER_COREO_CAT, PgR_COREO_CAT) 
levels(QIA$ER_COREO_CAT) <- c("Negative", "Positive")
levels(QIA$PgR_COREO_CAT) <- c("Negative", "Positive")

corr_ER_CAT <- cohen.kappa(QIA %>% select(ER_HUMAN_CAT, ER_COREO_CAT)) 
corr_PgR_CAT <- cohen.kappa(QIA %>% select(PgR_HUMAN_CAT, PgR_COREO_CAT)) 

#As for HER2_HUMAN VS HER2_COREO, the same as Aim 2

#Aim 4B: Discordant between HER2_COREO and HER2_FISH
#the same as 3A

#======================END OF THE MAIN PART===========================================


##ANOTHER WAY OF MAKING SUMMARY TABLE FOR CAT VAR-----------------------
df_summary_cat <- QIA %>% select(DIAGNOSIS)
summ_cont <- as.data.frame(summary(df_summary_cont))
summ_cat <- table(df_summary_cat)
summ_cat_prop <- as.data.frame(prop.table(summ_cat))
summ_cat <- cbind(as.data.frame(summ_cat), summ_cat_prop)
summ_cat <- summ_cat[,c(1,2,4)]
colnames(summ_cat) <- c("DIAGNOSIS", "Freqency", "Proportion")

#-------------------qwraps2----------------------------------------------
summary_stat <- list(
  "AGE" =
    list(
      "mean (sd)" = ~qwraps2::mean_sd(pop, na_rm = TRUE),
      "median (Q1, Q3)" = ~qwraps2::median_iqr(pop, na_rm = TRUE),
      "min" = ~min(pop, na.rm = TRUE),
      "max" = ~max(pop, na.rm = TRUE),
      "Missing" = ~sum(is.na(pop))
    ),
  "ER_HUMAN" =
    list(
      "mean (sd)" = ~qwraps2::mean_sd(pop, na_rm = TRUE),
      "median (Q1, Q3)" = ~qwraps2::median_iqr(pop, na_rm = TRUE),
      "min" = ~min(pop, na.rm = TRUE),
      "max" = ~max(pop, na.rm = TRUE),
      "Missing" = ~sum(is.na(pop))
    ),
  "PgR_HUMAN" =
    list(
      "mean (sd)" = ~qwraps2::mean_sd(pop, na_rm = TRUE),
      "median (Q1, Q3)" = ~qwraps2::median_iqr(pop, na_rm = TRUE),
      "min" = ~min(pop, na.rm = TRUE),
      "max" = ~max(pop, na.rm = TRUE),
      "Missing" = ~sum(is.na(pop))
    ),
  "HER2_HUMAN" =
    list(
      "mean (sd)" = ~qwraps2::mean_sd(pop, na_rm = TRUE),
      "median (Q1, Q3)" = ~qwraps2::median_iqr(pop, na_rm = TRUE),
      "min" = ~min(pop, na.rm = TRUE),
      "max" = ~max(pop, na.rm = TRUE),
      "Missing" = ~sum(is.na(pop))
    )
)
summary_table(df_summary_cont, summary_stat)
#---------------------------------------------------------------------------------


#---------------------------------------------------------------------------------
#lsmeans boxplot
CLD = cld(marginal,
          alpha   = 0.05,
          Letters = letters,    ### Use lower-case letters for .group
          adjust  = "tukey")    ### Tukey-adjusted comparisons


ggplot(na.omit(QIA %>% select(HER2_COREO_CAT, HER2D17Z1_RATIO)),
       aes(x     = HER2_COREO_CAT,
           y     = fit_3b$fitted.values)) +
  geom_point(shape  = 15,
             size   = 4) +
  geom_errorbar(aes(ymin  =  lower.CL,
                    ymax  =  upper.CL),
                width =  0.2,
                size  =  0.7) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  ylab("Least square mean\n HER2D17Z1_RATIO") +
  geom_text(nudge_x = c(0, 0, 0),
            nudge_y = c(120, 120, 120),
            color   = "black")
