
rm(list=ls())
library(dplyr)
#library(lubricate)
library(survival)

#---------------------------------------0. DATA PREPERATION------------------------------------------------------------------------

#read data (in "Ray" folder)
df_base <- get(load("X:/Cancer/Ray/Immuno_6/Output/base_pat_elig.RData"))
df_diag <- get(load("X:/Cancer/Ray/Immuno_6/Output/final_diag_elig.RData"))

# Use this new file (in "xuan" folder)
df_auto <- read.csv("X:/Cancer/Data/ICD10_code/Autoimmune_ICD10_list_with_my_score.csv", header = TRUE) #names of autoimmunes

# select autoimmune ICD10 based on : doctor_score>0 or My_score==1 (lili)
names_auto  = 
  df_auto %>% 
  dplyr::filter(My_score>0 | dan_score==1) %>%  # My_score inlucde doctor's score
  dplyr::transmute(ICD_CODE = gsub("\\.", "", ICD_CODE)) %>%
  unlist() %>% 
  as.vector()

#--------------------------2. PRE-INVESTIGATION OF SURVIVAL WRT ANY AUTOIMMUNE----------------------------------------------

# BEFORE YOUR TRY RUN:
# YOU MUST RUN PART 0 FIRST BEFORE YOU RUN PART 2

#--------------------------2.1 NUMBER OF AUTOIMMUNES OF EACH PATIENT AND SUMMARY STATISTICS---------------------------------------------------

#drop all histirucal recordings
df_n_auto_temp <- df_diag %>% 
  group_by(Patid, Diag) %>%
  filter(sum(history) <= 0) 

# get unique autoimmne codes per subject and then count n
df_n_auto <- df_n_auto_temp %>%
  filter(Diag %in% names_auto) %>%
  select(Patid, Diag) %>%
  distinct() %>%
  group_by(Patid) %>%
  summarise(n=n()) %>%
  arrange(desc(n))

#show summary statistics
#summary(df_n_auto)
#hist(df_n_auto$n)   
table(df_n_auto$n)   
#   1   2   3   4   5 
#948  109  30   2   1

#-------------------------2.2 SIMILARITY MATRIX--------------------------------------------------------------------------

# Find the patients with more than 1 autoimmune
df_similar_patid <- df_n_auto %>%
  filter(n != 1)

# Prepare the data we need
df_similar <- df_n_auto_temp %>%
  filter(Patid %in% df_similar_patid$Patid) %>% #keep all the patients with more than 1 autoimmune
  filter(Diag %in% names_auto) %>%
  select(Patid, Diag)

#Initialize the similarity matrix
mat_similar <- matrix(0, length(names_auto), length(names_auto))
colnames(mat_similar) <- names_auto
rownames(mat_similar) <- names_auto

t <- 1
#Loop for each patient
for (pat in df_similar_patid$Patid) {
  df_similar_pat <- df_similar %>%
    filter(Patid == pat) %>% #filter the recordings of this patient 
    distinct()
  d <- df_similar_pat$Diag # the autoimmunes this patient have
  axis <- numeric(0)
  
  #find the position of autoimmunes of this patient in the name list
  for (p in 1: length(d)) {  
    axis <- c(axis, which(names_auto == d[p]))
  }
  
  #find the posision of pairs in the matrix and add 1 to each 
  for(i in 1:length(axis)){
    for (j in i+1:length(axis)) {
      #To make up an upper-triangular matrix, we always access the element with larger column index than row
      a <- min(axis[i], axis[j])
      b <- max(axis[i], axis[j])
      mat_similar[a, b] <- mat_similar[a, b] + 1
    }
  }
  print(t)
  t <- t+1
}

# show the frequencies and positions of all the pairs
# find all the non-zero element and its indexes in the matrix
df_pair <- as.data.frame(which(mat_similar != 0, arr.ind = T)) #get indexes of all the non-zero elements
rownames(df_pair) <- c()
df_pair <- df_pair %>% 
  mutate(freq = mat_similar[which(mat_similar!=0)]) %>%  #frequency of each pair
  mutate(name_1=names_auto[row], name_2=names_auto[col]) %>% #autoimmune code names of each pair
  arrange(desc(freq))

#mat_pat_auto <- matrix(0, length(df_similar_patid$Patid), length(names_auto) + 1)
#colnames(mat_pat_auto) <- c("Patid", names_auto)

#saving results
#write.csv(df_n_auto, "X:/Cancer/Yiyuan Huang/output/freq_auto_by_patient.csv", row.names = F)
write.csv(df_pair, "X:/Cancer/Yiyuan Huang/output/similarity_freq.csv", row.names = F)
saveRDS(mat_similar, "X:/Cancer/Yiyuan Huang/output/similarity_matrix.rds")

#--------------------------3. SURVIVAL WRT ANY AUTOIMMUNE-----------------------------------------------------------------

#prepare baseline dataset
df_base <- df_base %>% 
  select(Patid, YMDOD, Fst_Immuno, Lst_followup) %>%
  mutate(death = 1L * !is.na(YMDOD), time = Lst_followup - Fst_Immuno ) 

#prepare sequential dataset (diag)
#drop all histirucal recordings
df_temp <- df_diag %>% 
  group_by(Patid, Diag) %>%
  filter(sum(history) <= 0) 

df_any_auto <- df_temp %>%
  filter(Diag %in% names_auto) %>%
  select(Patid:Fst_Immuno, inpatient_flag) %>%
  group_by(Patid, Diag) %>%
  filter(Fst_Dt == min(Fst_Dt)) %>% #only keep the patients' first recording of each autoimmune  
  arrange(Fst_Dt) %>%
  arrange(Patid) %>%
  mutate(day = Fst_Dt - Fst_Immuno)

# GENERATE COVARIATE (any_immu)
# initialize the covariate (all 1)
df_any_auto <- df_any_auto %>% mutate(any_immu = 1)
# for patients with more than 1 autoimmune, any_immu can take values from 1 up to number of autoimmunes this patient has 
for (patid in df_any_auto$Patid) {
  df_any_auto$any_immu[which(df_any_auto$Patid==patid)] <- 1:sum(df_any_auto$Patid==patid)
}
df_any_auto <- df_any_auto %>%
  ungroup() %>%
  select(-(Fst_Dt:inpatient_flag))

# find the patients whose autoimmune appears in the day of Fst_immuno 
#(then all the other subjets in df_merge should have a new row with immu=0 at day=0) (???)
df_temp1 <- df_any_auto %>% filter(day == 0)
# all the other patients in baseline but not in df_temp1 are qualified subjects
df_temp2 <- df_base %>% 
  filter(!(Patid %in% df_temp1$Patid)) %>%
  mutate(day = 0, any_immu = 0) %>%
  select(-death, -time) %>%
  select(-(YMDOD:Lst_followup))
# conbine immu=0 and immu=1
df_seq <- rbind(df_any_auto, df_temp2) %>% 
  arrange(day) %>% 
  arrange(Patid)
##sequential dataset 

#delete all duplicate identifiers
df_seq <- df_seq %>%
  group_by(Patid, day) %>%
  filter(any_immu==max(any_immu)) %>%
  ungroup()

## tmerge part to make tstart and tstop (Problem 3 arises here, see details in PROBLEM CHECKING 3)
df_2 <- tmerge(df_base, df_base, id= Patid, death = event(time, death))
df_2 <- tmerge(df_2, df_seq, id=Patid, any_immu=tdc(day, any_immu))



#--------------------------END OF THE MAIN CODES--------------------------------------------------------------------------


#-----------------TODO: Another way to check NUMBER OF CASES--------------------------------------------------------------
df_diag_check <- df_diag %>% 
  group_by(Patid, Diag) %>%
  filter(sum(history) <= 0) 

df_diag_check2 <- df_diag_check %>%
  filter(Diag %in% names_auto) %>%
  select(Patid, Diag) %>%
  distinct() %>%
  group_by(Diag) %>%
  summarise(n=n())


#-----------------------------------------PROBLEM CHECKING----------------------------------
#s1 <- setdiff(names_auto, df_diag_check2$Diag) ## check which elements in names_auto are not qualified autoimmune 

#filter(df_diag, Diag == "H4413") ##for example ("H4413" does not even exist in df_diag)

#names_check <- read.csv("X:/Cancer/Yiyuan Huang/data/Immu6_new_autoimmune_inpat.csv", header = TRUE)$Diag #names of autoimmunes

#sum(s1 %in% names_check)

#-----------------------------------------PROBLEM CHECKING 2----------------------------------------------

filter(result_final_merge, n<n_inpat)
#3 autoimmunes whose n < n_inpat!!! B20, K743, D68311

#TODO, remove p>0.01; warnings; n<=3

#----------------------------------------PROBLEM CHECKING 3----------------------------------------------------

#Run part 0 and part 3: <prepare baseline dataset> before you run the following codes

#df_base: there are 74 patients with more than 1 row 
df_check <- df_base %>% 
  group_by(Patid) %>% 
  summarise(n=n()) %>% 
  filter(n!=1) #problematic patients id

df_prob <- df_base %>% 
  filter(Patid %in% df_check$Patid) #corresponding rows of the problematic patients

