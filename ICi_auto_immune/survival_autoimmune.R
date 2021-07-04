#===================================================================================
# STUDY the effect of Autoimmunes on overall survival
# Method: time-dependent Cox-PH model 

# required data files: base_pat_elig.rds (baseline)
#                      final_diag_elig.rds (time-depnedent covariates: autoinmune)
#                      Autoimmune_ICD10_list_with_my_score.csv (name list of autoimmunes)

# by Yiyuan Huang
# Date:  
#==================================================================================

rm(list=ls())
library(dplyr)
#library(lubricate)
library(survival)
 
#---------------------------------------0. DATA PREPERATION------------------------------------------------------------------------

#read data (in "Ray" folder)
df_base <- get(load("X:/Cancer/Ray/Immuno_6/Output/base_pat_elig.RData"))
df_diag <- get(load("X:/Cancer/Ray/Immuno_6/Output/final_diag_elig.RData"))

#df_base <- readRDS("X:/Cancer/Yiyuan Huang/data/base_pat_elig.rds") # baseline
#df_diag <- readRDS("X:/Cancer/Yiyuan Huang/data/final_diag_elig.rds") # diagnostic recordings
#df_auto <- read.csv("X:/Cancer/Code/Yiyuan Huang/data/Immu6_new_autoimmune_inpat.csv", header = TRUE) #names of autoimmunes

# Use this new file (in "xuan" folder)
df_auto <- read.csv("X:/Cancer/Data/ICD10_code/Autoimmune_ICD10_list_with_my_score.csv", header = TRUE) #names of autoimmunes

# select autoimmune ICD10 based on : doctor_score>0 or My_score==1 (lili)
names_auto  = 
  df_auto %>% 
  #dplyr::filter(doctor_score>0 | My_score==1) %>%
  dplyr::transmute(ICD_CODE = gsub("\\.", "", ICD_CODE)) %>%
  unlist() %>% 
  as.vector()
#char vector of the coding names of autoimmunes

#--------------------------1. SURVIVAL WRT A SPECIFIC AUTOIMMUNE (ONE-BY-ONE)-------------------------------------

#---------------------------1.1 SET THE RANGE---------------------------------------------------------------------

#NEED TO MODIFY THIS BLOCK BASED ON WHETHER IMPAT RESTRICTION IS NEEDED
df_diag <- df_diag %>% 
  select(Patid, Diag,  Fst_Immuno, Fst_Dt, history, inpatient_flag, autoimmune_flag)  %>%
  filter(inpatient_flag==1) #Delete this pipe when inpat restriction is not needed

#df_auto<- df_auto %>% filter(n >= 10) #number of cases >= 10
#names_auto <- df_auto_sub$Diag #vector of the coding names of autoimmunes
df_base <- df_base %>% 
  select(Patid, YMDOD, Lst_followup) %>%
  mutate(death = 1L * !is.na(YMDOD) ) 

#-----------------------------1.2 function for synthesizing data and fitting the model------------------------------

# function designed to implement the task
# input: A character variable denoting the name of the autoimmune
# output: A list contain 4 elements: 
#               name(character): the name of the autoimmune;
#               n(integer): the number of cases used to fit the model; 
#               hr(double): hazard ratio (patients w. autoimmune/patients wo. autoimmune);
#               p(double): p-value of corresponding estimatiion of coefficient.
               
Surv_auto <- function(autoname){
  
  df_diag_temp <- df_diag %>% filter(Diag == autoname)
  
  # filter by the patient WITHOUT the history of the autoimmune
  filter_hist <- df_diag_temp %>%
    group_by(Patid) %>%
    summarize(sum_hist = sum(history)) %>%
    filter(sum_hist == 0)
  
  df_diag_temp <- df_diag_temp %>% filter(Patid %in% filter_hist$Patid)
  #df_base_temp <- df_base %>% filter(Patid %in% filter_hist$Patid)
  
  #merge fst-immuno into base to calculate the time-to-event denoted as "time"
  #generate a dataset with Fst_immuno of all the patients
  fst_im <- df_diag %>% group_by(Patid) %>% summarise(Fst_Immuno=mean(Fst_Immuno)) %>% ungroup()
  #merge Fst_immuno and Lst_followup (in baseline) together
  df_base_temp <- df_base %>%
    left_join(., fst_im, by = "Patid") %>% 
    mutate(time = Lst_followup - Fst_Immuno) %>% 
    select(-Lst_followup, -Fst_Immuno, -YMDOD) %>%
    filter(time != 0) # delete those with same Fst_immuno and Lst_followup               
  ## well-prepared baseline dataset
  
  # merge base and diag together 
  # to calsulate the sequential time of each patient, which is denoted as "day" 
  # Since we just care about the first day when the autoimmune appears so only rows with min(day) within each patient are filtered
  # Generate covariates immu to denote whether the autoimmune appears at the time recorded
  # All the immu=1 in this dataset 
  df_merge <- df_diag_temp %>%
    left_join(., df_base, by ="Patid") %>%
    mutate(day = Fst_Dt - Fst_Immuno ) %>%
    group_by(Patid) %>%
    filter(day == min(day)) %>%
    select(Patid, day) %>%
    mutate(immu=1)   
  
  df_merge$day <- as.numeric(df_merge$day)
  # find those rows whose immu=0
  # These are diveded into 2 cases: (a) patients with different Fst_immuno and Fst_date of autoimmune (in df_merge);
  #                                 (b) patients without autoimmune at all (beyond df_merge).
  
  # find the patients whose autoimmune appears in the day of Fst_immuno 
  #(then all the other subjets in df_merge should have a new row with immu=0 at day=0) (???)
  df_temp1 <- df_merge %>% filter(day == 0)
  # all the other patients in baseline but not in df_temp1 are qualified subjects
  df_temp2 <- df_base_temp %>% 
    filter(!(Patid %in% df_temp1$Patid)) %>%
    mutate(day=0, immu=0) %>%
    select(-death, -time)
  # conbine ummu=0 and ummu=1
  df_seq <- rbind(df_merge, df_temp2) %>% 
    arrange(Patid) %>% 
    arrange(day) %>%
    filter(Patid %in% df_base_temp$Patid)
  ##sequential dataset
  
  ## tmerge part to make tstart and tstop
  df_2 <- tmerge(df_base_temp, df_base_temp, id= Patid, death = event(time, death))
  df_2 <- tmerge(df_2, df_seq, id=Patid, immu=tdc(day, immu))
  
  # fit coxph with time-dependent cov
  fit_coxph <- coxph(Surv(tstart, tstop, death) ~ immu, df_2)
  
  name <- autoname  #the code name of auto immune
  n <- sum(unique(df_base_temp$Patid) %in% filter_hist$Patid) #number of cases in the model
  hr <- exp(fit_coxph$coefficients) #hazard ratio
  p <- summary(fit_coxph)$coefficients[1,5] #p-value
  
  result <- list(name, n, hr, p)
  return(result)
}

#----------------------------------------------1.3 WARNING HANDLING----------------------------------------

catchToList <- function(expr) {
  val <- NULL
  myWarnings <- NA #denote any warnings that may arisedf_diag
  wHandler <- function(w) {
    myWarnings <<- w$message
    invokeRestart("muffleWarning")
  }
  val <- withCallingHandlers(expr, warning = wHandler)
  list(value = val, warnings = myWarnings)
}
#-----------------------------------------------1.4 MAIN PART: LOOP FOR EVERY AUTOIMMUNE-------------------------------------------------

#Main part 1 (w/o inpatient restriction): loop applying function Surv_auto()
#i=1
result_final <- data.frame(name = character(), n = numeric(), hr = numeric(), p = numeric(), warning = character())
for(autoname in names_auto){
#  temp_result <- Surv_auto(autoname)
#  df_result <- data.frame(Diag=temp_result[[1]], n=temp_result[[2]], hr=temp_result[[3]], p=temp_result[[4]])
  temp_result <-catchToList(Surv_auto(autoname))  
  df_result <- data.frame(Diag=temp_result$value[[1]], n=temp_result$value[[2]], hr=temp_result$value[[3]], p=temp_result$value[[4]], warning = temp_result$warnings)
  result_final <- rbind(result_final, df_result)
  print(autoname)
}

#Main part 2 (with npatient restriction): loop applying function Surv_auto()
#i=1
result_final_inpat <- data.frame(name_inpat = character(), n_inpat = numeric(), hr_inpat = numeric(), p_inpat = numeric(), warning_inpat = character())
for(autoname in names_auto){
  #  temp_result <- Surv_auto(autoname)
  #  df_result <- data.frame(Diag=temp_result[[1]], n=temp_result[[2]], hr=temp_result[[3]], p=temp_result[[4]])
  temp_result <- catchToList(Surv_auto(autoname))  
  df_result <- data.frame(Diag_inpat=temp_result$value[[1]], n_inpat=temp_result$value[[2]], hr_inpat=temp_result$value[[3]], p_inpat=temp_result$value[[4]], warning_inpat = temp_result$warnings)
  result_final_inpat <- rbind(result_final_inpat, df_result)
  print(autoname)
}


#----------------------------------1.5 RESULTS PROCESSING--------------------------------------------------------------------

# merge part I and II adding diag names
#result_final=result_final %>% left_join(., df_auto %>% 
#                                 dplyr::select(Diag, rheumo_flag, DESC)) PROBLEM: those columns are in old csv file

df_auto <- df_auto %>% mutate(Diag=names_auto)
# merge result_final (w/o inpat) and result_final_inpat as well as 3 columns (doctor_score, newly-added, My_score) 
result_final_merge <- result_final %>%
  left_join(., result_final_inpat, by <- c("Diag"="Diag_inpat")) %>%
  left_join(., df_auto %>% select(ICD_DESC, doctor_score, newly_added, My_score, Diag))

## filter meaningful rows (autoimmunes)
#merge result
result_filter <- result_final_merge %>% 
  filter(is.na(warning) == 1 | is.na(warning_inpat) == 1) %>% #non-warning
  filter(p <= 0.01 | p_inpat <= 0.01 ) %>% #significant estimate of hazard ratio
  filter(n > 3 | n_inpat > 3) #adequate number of cases

#without inpat
wo_inpat_filter <- result_final %>%
  left_join(., df_auto %>% select(ICD_DESC, doctor_score, newly_added, My_score, Diag)) %>%
  filter(is.na(warning) == 1) %>% #non-warning
  filter(p <= 0.01) %>% #significant estimate of hazard ratio
  filter(n > 3) %>%#adequate number of cases
  select(-warning)

#with inpat
with_inpat_filter <- result_final_inpat %>%
  left_join(., df_auto %>% select(ICD_DESC, doctor_score, newly_added, My_score, Diag), by=c("Diag_inpat"="Diag")) %>%
  filter(is.na(warning_inpat) == 1) %>% #non-warning
  filter(p_inpat <= 0.01 ) %>% #significant hazard ratio
  filter(n_inpat > 3) %>% #adequate number of cases
  select(-warning_inpat)

#save as .csv file
write.csv(result_final_merge, "X:/Cancer/Yiyuan Huang/output/Surv_merge_warning.csv", row.names = F)
write.csv(result_filter, "X:/Cancer/Yiyuan Huang/output/Surv_final.csv", row.names = F)
write.csv(wo_inpat_filter, "X:/Cancer/Yiyuan Huang/output/Surv_final_wo_inpat.csv", row.names = F)
write.csv(with_inpat_filter, "X:/Cancer/Yiyuan Huang/output/Surv_final_with_inpat.csv", row.names = F)

