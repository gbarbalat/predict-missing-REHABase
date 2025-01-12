# Clean raw data then run SuperLearner algorithm to predict missing rehab outcome data at t0 based on background factors
# includes multiple imputations of missing data in background factors
# also includes calculation of attrition weights (if removing observations of missing values)

# HDR ----

rm(list=ls())

library(ranger)
library(pROC)
library(mice)
library(dplyr)
library(tidyr)
library(tidyverse)
library(SuperLearner)
library(caret)
library(ck37r)
library(fastshap)
library(ggplot2)
library(patchwork)
library(shapviz)

today <- format(Sys.Date(), "%Y-%m-%d")

plot_var_outcome <- FALSE
bivariate <- FALSE
plot_imp <- FALSE
#Addictions: None:0 ;;; Tox: 1,2,3,7,8,10,11,13,15 ;;; Behav:5,6,9,12,16 ;;;;; Both
# ""0"" = ""Aucun""
# ""1"" = ""Tabac"" ; ""2"" = ""Alcool"" ; ""3"" = ""Cannabis"" ; ""10"" = ""Amphétamines/Cocaïne/Ecstasy"" ; 
# ""7"" = ""Benzodiazépines"" ;""11"" = ""Antalgiques/Antidépresseurs/Anxiolytiques/Hypnotiques/Neuroleptiques"" ; 
# ""13"" = ""Hallucinogènes (LSD\\ Champignons)"" ; ""15"" = ""Protoxyde d'azote"" ; ""8"" = ""Opiacés/Morphiniques"" ;
# ""5"" = ""Jeux d'argent (Casino\\FDJ\\PMU)"" ; ""9"" = ""Internet et Jeux video"" ; ""6"" = ""Achats compulsifs"" ; 
# ""12"" = ""Alimentaire (Boulimia nervosa)"" ; ""16"" = ""Sexuelle/Paraphilique/Pornographique"" ;
# ""99"" = ""Autre"" ; 
# "8888"" = ""Non demandé durant l'entretien""


varOI <- "WEMWBS_TOT_NB"#"EGF" STORI_CONSCIENCE

path="C:/Users/Guillaume/Desktop/All projects/Gal REHABASE Fdam HORIZONS/Explore REHABase/"
#path="/bettik/barbalag/missing_data/"

p_miss_indiv <- 100#pct of missing values to keep
maxit <- 10;#for mice, this has to be increased!!!!
m <- 10;#for mice, this has to be increased!!!!
method <- "method.AUC"
#method <- "method.NNloglik"#Imbalanced dataset so use NLL (neg log like) method.NNloglik
#Stratify on Y when randomly assigning i.i.d. units to validation sets  
V=2#Nfolds for SL - 10 V=20 to 200 see Philipps
nsim=1;#for SHAP
k_outer<- 2 #outer (outside caret) folds 
k_inner <- 2 #inner (caret) folds 
top_n_select <- 20 #display SHAP one way plot
set.seed(123)

# Define the SuperLearner library
source("C:/Users/Guillaume/Desktop/All projects/ML/SL_predict_caret_adaptive.R")
#source("/bettik/barbalag/missing_data/SL_predict_caret_adaptive.R")

SL.library <- list("SL.mean",
                   c("SL.caret.xgboost", "All", "screen.mixed", "screen.glmnet", "screen.randomForest"),
                   c("SL.caret.ranger", "All", "screen.mixed", "screen.glmnet", "screen.randomForest"),
                   c("SL.gam_cts", "All", "screen.mixed", "screen.glmnet", "screen.randomForest"),
                   c("SL.rpart", "All", "screen.mixed", "screen.glmnet", "screen.randomForest"),
                   c("SL.caret.naive_bayes", "All", "screen.mixed", "screen.glmnet", "screen.randomForest"),
                   c("SL.caret.earth", "All", "screen.mixed", "screen.glmnet", "screen.randomForest"),
                   c("SL.glm", "All", "screen.mixed", "screen.glmnet", "screen.randomForest"),
                   c("SL.glmnet_1", "All", "screen.mixed", "screen.glmnet", "screen.randomForest"),
                   c("SL.glmnet_0.5", "All", "screen.mixed", "screen.glmnet", "screen.randomForest"),
                   c("SL.glmnet_0", "All", "screen.mixed", "screen.glmnet", "screen.randomForest")
)

SL.library <- list(c("SL.step.interaction",  "All","screen.mixed", "screen.glmnet", "screen.randomForest"))
SL.library <- list(c("SL.caret.step.interaction",  "screen.mixed", "screen.glmnet", "screen.randomForest"))
SL.library <- list(c("SL.caret.ksvm", "All", "screen.mixed", "screen.glmnet", "screen.randomForest"))
SL.library <- list(c("SL.caret.glmnet", "All", "screen.mixed", "screen.glmnet", "screen.randomForest"))
SL.library <- list(c("SL.caret.earth", "All", "screen.mixed", "screen.glmnet", "screen.randomForest"))
SL.library <- list(c("SL.caret.rpart", "All", "screen.mixed", "screen.glmnet", "screen.randomForest"))
SL.library <- list(c("SL.caret.bartMachine", "All", "screen.mixed", "screen.glmnet", "screen.randomForest"))
SL.library <- list(c("SL.caret.gam_cts",  "All"))
SL.library <- list(c("SL.caret.glm", "All", "screen.mixed", "screen.glmnet", "screen.randomForest"))
SL.library <- list(c("SL.caret.ranger", "screen.randomForest"))
SL.library <- list(c("SL.caret.rpart", "All"))

# First, let's create a function to categorize addiction values
categorize_addiction <- function(x) {
  if (is.na(x) || x == "") return(NA)
  if (x == "0") return("none")
  
  values <- as.numeric(strsplit(x, ",")[[1]])
  tox <- any(values %in% c(1,2,3,7,8,10,11,13,15))
  behav <- any(values %in% c(5,6,9,12,16))
  
  if (tox && behav) return("both")
  if (tox) return("tox")
  if (behav) return("behav")
  return(NA)
}

# Function to categorize medications
categorize_med <- function(data, class_name, med_class) {
  data %>%
    filter(if_any(starts_with("LISTE_CLASSE_TTTPSY_"), ~ . == med_class)) %>%
    dplyr::select(StudySubjectID) %>%
    mutate(!!class_name := 1)
}

#Function to find missing values
try_missing <- function(merged_basicRegp) {
  (miss_var_summary_obs <-naniar::miss_var_summary(merged_basicRegp)); print(miss_var_summary_obs)
  (miss_var_table_obs <- naniar::miss_var_table(merged_basicRegp) %>%
      mutate(pct_miss_in_var=n_miss_in_var*100/nrow(merged_basicRegp), .after=n_miss_in_var))
  print(miss_var_table_obs)
  (miss_case_summary_obs <- naniar::miss_case_summary(merged_basicRegp))
  print(miss_case_summary_obs)
  (miss_case_table_obs <- naniar::miss_case_table(merged_basicRegp) %>%
      mutate(pct_miss_in_case=n_miss_in_case*100/ncol(merged_basicRegp), .after=n_miss_in_case) %>%
      arrange(desc(pct_miss_in_case))
  )
  print(miss_case_table_obs)
  return(miss_case_summary_obs)
}

# Function to calculate FMI
calculate_fmi <- function(imp, formula) {
  m <- imp$m  # number of imputations
  n <- nrow(imp$data)#sample size
  k <-model.matrix(formula,imp$data) %>% colnames %>% length#nb of param to fit
  
  # Fit model to each imputed dataset
  fits <- lapply(1:m, function(i) {
    data <- complete(imp, i)
    fit <- lm(formula, data = data)
    list(coef = fit$coefficients, se = sqrt(diag(vcov(fit))))
  })
  
  # Extract coefficients and standard errors
  coefs <- sapply(fits, function(x) x$coef)
  ses <- sapply(fits, function(x) x$se)
  
  # Calculate components
  theta_bar <- rowMeans(coefs)
  # Calculate summed differences for each coefficient
  B <- sapply(rownames(coefs), function(coef_name) {
    sum((coefs[coef_name,] - theta_bar[coef_name])^2)/(m-1)
  })
  names(B) <- rownames(coefs)  
  W <- rowMeans(ses^2)#averaged within imputation variance
  
  # Calculate FMI
  T <- W + B + B/m # T <- W+B
  lambda <- (B + B/m)/ T; lambda
  riv <- lambda/(1-lambda); riv
  df_old <- (m-1)/lambda^2; df_old
  df_obs <- (((n-k)+1)/((n-k)+3))*(n-k)*(1-lambda); df_obs
  df_adj <- (df_old*df_obs)/(df_old+df_obs);df_adj
  #The difference between lambda and FMI is that FMI is adjusted for the fact that the number of imputed datasets 
  #that are generated is not unlimitedly large. These measures differ for a small value of the df.
  num_fmi <- riv + 2/(df_adj+3)
  den_fmi <- 1+riv
  fmi <- num_fmi/den_fmi; fmi
  
  return(list(lambda,fmi))
}


# PRE-PROCESSING ----
extract02012023=read.csv(paste0(path,"extract02012023.csv"), sep=";",encoding="UTF-8",fileEncoding="latin1")
colnames(extract02012023)
unique(extract02012023$StudySubjectID)
unique(extract02012023$TIME)
unique(extract02012023$SEXE)
unique(extract02012023$LANGMAT)
unique(extract02012023$NIVETUD_cat)
unique(extract02012023$LISTE_CLASSE_DIAG1R)
unique(extract02012023$LISTE_DIAG1R)
unique(extract02012023$LISTE_CLASSE_DIAG2R_1)
unique(extract02012023$NIVETUD_cat)

dx_to_keep <- c("2-Spectre de la SCHIZOPHRENIE","1-Troubles NEURODEVELOPPEMENTAUX",#"16-Troubles ADDICTIFS"  ,                         
             "18-Troubles de la PERSONNALITE","4-Troubles DEPRESSIFS","anxiety disorders" ,                       
             "3-Troubles BIPOLAIRES" 
)
dx_to_keep <- "2-Spectre de la SCHIZOPHRENIE"

#basic set of variables
col_to_keep <- c(
              "CENTRE",
              "AGE_MEDSOC",
              "SEXE",
              "NIVETUD_cat",
              #"LISTE_CLASSE_DIAG1R",#if anal only in one dx group
              "COMOR_PSY", 
              "FGA", "SGA", "noRx",#"ANX", "THY", "ATD", 
              "addict",
              "SOMA_1","TTTSOMA_1",
              #"EGF",	
              #"CGI_SEVERITE",
           
              "SIT_FAM_cat", "ETRE_PARENT",
              "ADRSSR_cat",
              "LGMT_cat",
              "SIT_PRO_cat","RQTH",
              "DUREE_MALADIE",
              #"NBR_HOSPI", "DUREE_HOSPI",
              "TS",#"NBR_TS",#
              "MARGIN", "MARGIN_ACTPASS",
              "ATCD_MEDLEG",
              
              # item based
              # "SQOL18_RELFAM_SATISF","SQOL18_AMIS_SATISF","SQOL18_VIESENTIM_SATISF",#phy, psy, res, sel
              # "SQOL18_AUTONOMIE_SATISF",#"SQOL18_AUTONOMIE_SATISF","SQOL18_RELFAM_SATISF",	"SQOL18_AMIS_SATISF",	"SQOL18_VIESENTIM_SATISF",
              # #"SERS_SC_TOT",	#"SERS_POS",	
              # "IS_BIRCHWOOD_SC_TOT",#"IS_BIRCHWOOD_TTT",
              # "MARS_TOT",
              # "ISMI_SC_TOT",	#"ISMI_ALIEN",	
              # #"STORI_ST_RETABLISSEMENT",
              # "WEMWBS_TOT_NB",
              # #EAS
              
              #"MEMCHIF_MCD",#
              #"MEMCHIF_MCI",	
              
              "missing_out" #outcome
              )

#var to be transformed
old_var <- c(  
               # "SOMA_1",
               # "addict",
               # "ATCD_MEDLEG",
               
               "MARGIN","MARGIN_ACTPASS",
               
               #"DUREE_MALADIE",
               #"NBR_HOSPI", 
               #"DUREE_HOSPI",
               
               "FGA", "SGA", "noRx"#"ANX", "THY", "ATD" 
               
               #"SQOL18_RELFAM_SATISF","SQOL18_AMIS_SATISF","SQOL18_VIESENTIM_SATISF"
               )

#transformed variables ???? whilst passive imputation
new_var <- c(
  
  #"Comorb",#Y/N"SOMA_1","ATCD_MEDLEG","addict",
  
  "MARGIN_hx",#No, Past, Current
  
  #"HDR",#"DUREE_MALADIE", "NBR_HOSPI", "DUREE_HOSPI",
  
  "Rx"#SGA, FGA, Both, Nil
  
  #"SQOL18_REL"#SQOL18_RELFAM_SATISF","SQOL18_AMIS_SATISF","SQOL18_VIESENTIM_SATISF"
  
  
)

#auxiliary variables
aux_var <- c(
  
    # "COMOR_PSY", 
    # "TS","NBR_TS",#"TS_ANNEE1",
    # "SIT_FAM_cat", "ETRE_PARENT",
    # "ADRSSR_cat",
    # "LGMT_cat",
    # "SIT_PRO_cat","RQTH",
    # 
    # "EGF",
    # "CGI_SEVERITE",
    
    # "SERS_POS",	
    # "IS_BIRCHWOOD_TTT",
    # "MARS_TOT",
    
    #"MEMCHIF_MCI"
)

numeric_col <- c("TIME",
              "AGE_MEDSOC",#"nDIAG2R", "nTTTPSY",
              #"nDIAGSOMA",
              "EGF",	"CGI_SEVERITE",
              "NBR_TS",
              "NBR_ENFNTS",	"AGE_DBT_MALADIE",
              "DUREE_MALADIE",
              "NBR_HOSPI", "DUREE_HOSPI",
              "SQOL18_ESTIM2SOI_SATISF",	"SQOL18_RESILIENCE_SATISF",	"SQOL18_AUTONOMIE_SATISF",	"SQOL18_BI1ETREPHYSIK_SATISF",	"SQOL18_RELFAM_SATISF",	"SQOL18_AMIS_SATISF",	"SQOL18_VIESENTIM_SATISF",	"SQOL18_BI1ETREPSYKO_SATISF",	"SQOL18_TOT_SATISF",
              "WEMWBS_TOT_NB",
              "EAS_SOINSPERSO",	"EAS_VIEQUOT",	"EAS_RESSOURCES",	"EAS_RELATIONSEXT",	"EAS_AFFECTIVE",	#"EAS_TOTAL",
              "SERS_POS",	"SERS_NEG",	"SERS_SC_TOT",
              "IS_BIRCHWOOD_SYMPT",	"IS_BIRCHWOOD_MALADIE",	"IS_BIRCHWOOD_TTT",	"IS_BIRCHWOOD_SC_TOT",
              "MARS_TOT",
              "ISMI_ALIEN",	"ISMI_APPROB_STEREO",	"ISMI_EXP_DISCRIM",	"ISMI_RET_SOC",	"ISMI_RESIST_STIGMA",	"ISMI_SC_TOT",
              "STORI_MORATOIRE",	"STORI_CONSCIENCE",	"STORI_PREPARATION",	"STORI_RECONSTRUCTION",	"STORI_CROISSANCE",	#"STORI_ST_RETABLISSEMENT",
              "MEMCHIF_MCD", #"MEMCHIF_MCD_STRD",	"MEMCHIF_EMCD",	
              "MEMCHIF_MCI"#,	"MEMCHIF_MCI_STRD",	"MEMCHIF_EMCI",	
              #"D2R_CCT_NB",	"D2R_CCT_CENT",	"D2R_E_NB",	"D2R_E_CENT",	"D2R_CC_NB",	"D2R_CC_CENT",
              #"ACSO",	"MASC_NB",	"MASC_SIGMA",	"AIHQ_HB_NB",	"AIHQ_HB_SIGMA",	"AIHQ_ATTRIB_RESP_NB",	"AIHQ_ATTRIB_RESP_SIGMA",	"AIHQ_AB_NB",	"AIHQ_AB_SIGMA"
)

all_NA <- c("","8888,00","Non demandé durant l'entretien","Dunno", "8888.00","Non demandé durant entretien",
            "8888,0","8888.0")

## merged_a/b ----
#Outcome of interest at time==0
extract02012023_0 <- extract02012023 %>%
  filter(TIME==0) %>%
  mutate_at(vars(numeric_col), ~ gsub(",",".",.)) %>%
  mutate_at(vars(numeric_col),~ as.numeric(.)) %>% 
  mutate_if(is.character,
            ~if_else(. %in% all_NA, NA, .)
  ) %>%
  mutate_if(is.numeric, ~if_else(. == 8888,NA_real_,.)
  ) %>% 
  mutate(RQTH=case_when(RQTH=="Demande en attente" ~ "Non",
                        is.character(RQTH) ~ as.character(RQTH))
  ) 
colSums(is.na(extract02012023_0 %>% dplyr::select(!!sym(varOI))))
if (length(dx_to_keep==1)) {
  extract02012023_0 <- extract02012023_0 %>%
  filter(LISTE_CLASSE_DIAG1R %in% dx_to_keep)
}
#TTT voir script WM_Stigma.R 
# Categorize NLP (FGA and SGA)

NLP1<-unique(with(extract02012023_0,LISTE_DCI_1[LISTE_CLASSE_TTTPSY_1=="(VI)-NEUROLEPTIQUES (ou Antipsychotiques)"]))
NLP2<-unique(with(extract02012023_0,LISTE_DCI_2[LISTE_CLASSE_TTTPSY_2=="(VI)-NEUROLEPTIQUES (ou Antipsychotiques)"]))
NLP3<-unique(with(extract02012023_0,LISTE_DCI_3[LISTE_CLASSE_TTTPSY_3=="(VI)-NEUROLEPTIQUES (ou Antipsychotiques)"]))
NLP4<-unique(with(extract02012023_0,LISTE_DCI_4[LISTE_CLASSE_TTTPSY_4=="(VI)-NEUROLEPTIQUES (ou Antipsychotiques)"]))
NLP_all<-unique(c(NLP1,NLP2,NLP3, NLP4))

SGA <- c("Aripiprazole", "Palipéridone", "Rispéridone", "Olanzapine", "Clozapine", "Amisulpride",
         "Quétiapine", "VI.13 - Antipsychotiques ATYPIQUES")
FGA <- setdiff(NLP_all, SGA);FGA <- FGA[!is.na(FGA)]

# Create who_SGA and who_FGA
who_SGA <- extract02012023_0 %>%
  filter(if_any(starts_with("LISTE_DCI_"), ~ . %in% SGA)) %>%
  dplyr::select(StudySubjectID) %>%
  mutate(SGA = 1)

who_FGA <- extract02012023_0 %>%
  filter(if_any(starts_with("LISTE_DCI_"), ~ . %in% FGA)) %>%
  dplyr::select(StudySubjectID) %>%
  mutate(FGA = 1)

# Categorize other medications
who_ANX <- categorize_med(extract02012023_0, "ANX", "(I)-ANXIOLYTIQUES") %>%
  bind_rows(categorize_med(extract02012023_0, "ANX", "(II)-HYPNOTIQUES"))

who_ATD <- categorize_med(extract02012023_0, "ATD", "(III)-ANTIDEPRESSEURS")
who_THY <- categorize_med(extract02012023_0, "THY", "(V)-THYMOREGULATEURS (ou Normothymiques)")

# who has no Rx
who_noRx <- extract02012023_0 %>%
    filter(if_any(starts_with("LISTE_CLASSE_TTTPSY_1"), ~ . == "(0)-AUCUN Traitement Psychotrope")) %>%
    dplyr::select(StudySubjectID) %>%
    mutate(noRx = 1)

# Combine all categories
merged_a <- extract02012023_0 %>%
  #left_join(who_ANX, by = "StudySubjectID") %>%
  # left_join(who_THY, by = "StudySubjectID") %>%
  # left_join(who_ATD, by = "StudySubjectID") %>%
  left_join(who_FGA, by = "StudySubjectID") %>%
  left_join(who_SGA, by = "StudySubjectID") %>%
  left_join(who_noRx, by = "StudySubjectID")
  # mutate(across(c(ANX, THY, ATD, FGA, SGA), ~ifelse(is.na(.), 0, .))) 

#ADDICTIONS
# apply categorize_addiction to create the new column
merged_a$addict <- sapply(merged_a$CONSO_ADDICT, categorize_addiction)
tmp <- merged_a %>%
  dplyr::select(CONSO_ADDICT,addict)
unique(merged_a$NBR_TS)

non_na_counts <- colSums(!is.na(merged_a))
# Trier les colonnes par ordre décroissant du nombre de valeurs non-NA
sorted_columns <- sort(non_na_counts, decreasing = TRUE)
# Afficher les résultats
print(sorted_columns)

id <- merged_a$StudySubjectID

# Create a binary outcome variable for missingness
merged_a$missing_out <- as.numeric(is.na(merged_a %>% dplyr::select(!!sym(varOI))))
merged_b <- merged_a %>%
  dplyr::select(col_to_keep) %>%
  mutate_if(is.character,~as.factor(.)) 

#explore NA and gross distributions
merged_b %>% dplyr::select(!where(is.numeric)) %>% map(table, useNA="always")
merged_b %>% dplyr::select(where(is.numeric)) %>% map(summary)

#what is obvious
#regroup centers with small N
#regroup dx
#NA if 0-Diagnostic NON ETABLI/NON CONFIRME; ADRSSR_cat=Autre
#NBR_TS=0 if TS=Non
#remove EAS

merged_basic <- merged_b


if (length(dx_to_keep)!=1) {
#group Centre and diag as there are too many levels 
levels(merged_basic$LISTE_CLASSE_DIAG1R) <- c(levels(merged_basic$LISTE_CLASSE_DIAG1R), "anxiety disorders")
merged_basic$LISTE_CLASSE_DIAG1R[startsWith(as.character(merged_basic$LISTE_CLASSE_DIAG1R), "5") | 
                           startsWith(as.character(merged_basic$LISTE_CLASSE_DIAG1R), "6") | 
                           startsWith(as.character(merged_basic$LISTE_CLASSE_DIAG1R), "7") | 
                           startsWith(as.character(merged_basic$LISTE_CLASSE_DIAG1R), "8") | 
                           startsWith(as.character(merged_basic$LISTE_CLASSE_DIAG1R), "9")] <- "anxiety disorders"
merged_basic$LISTE_CLASSE_DIAG1R[startsWith(as.character(merged_basic$LISTE_CLASSE_DIAG1R), "0")] <- NA
merged_basic$ADRSSR_cat[startsWith(as.character(merged_basic$ADRSSR_cat), "Autre")] <- NA

merged_basic <- merged_basic %>%
  mutate(LISTE_CLASSE_DIAG1R = fct_other(LISTE_CLASSE_DIAG1R, 
                                         keep = dx_to_keep, 
                                         other_level = "Other")) 
}

merged_basic$ADRSSR_cat <- droplevels(merged_basic$ADRSSR_cat)

table(merged_basic$LISTE_CLASSE_DIAG1R, useNA="always")
table(merged_basic$ADRSSR_cat, useNA="always")

merged_basic$NBR_TS <- ifelse(merged_basic$TS == "Non", 0, merged_basic$NBR_TS)
table(merged_basic$NBR_TS, useNA="always")

#create an "other" level for all categories that are less than 50
merged_basic <- merged_basic %>%
  mutate(CENTRE = fct_other(CENTRE,
                            keep = names(table(CENTRE))[table(CENTRE) >= 50]))
table(merged_basic$CENTRE, useNA="always") 

merged_basicRegp <- merged_basic %>%
  # mutate(Comorb = case_when(
  #   SOMA_1 == "Non" & addict == "none" & ATCD_MEDLEG == "Non" ~ "N",
  #   SOMA_1 == "Oui" | (addict != "none" & !is.na(addict)) | ATCD_MEDLEG == "Oui" ~ "Y",
  #   TRUE ~ NA_character_
  # )) %>%
  # mutate(HDR = ifelse(is.na(NBR_HOSPI) | is.na(DUREE_HOSPI) | is.na(DUREE_MALADIE) | DUREE_MALADIE == 0,
  #                     NA,
  #                     NBR_HOSPI * DUREE_HOSPI / DUREE_MALADIE)) %>%
  mutate(Rx = case_when(
     FGA == 1 & is.na(SGA) ~ "FGA",
     SGA == 1 & is.na(FGA) ~ "SGA",
     FGA == 1 & SGA==1 ~ "Both",
     noRx == 1 ~ "Nil",
    TRUE ~ NA_character_
  )) %>%
  mutate(MARGIN_hx= case_when(
    MARGIN=="Non" ~ "Non",
    MARGIN=="Oui" & MARGIN_ACTPASS=="Passée" ~ "Past",
    MARGIN=="Oui" & MARGIN_ACTPASS=="Actuelle" ~ "Current"
  )) %>%
  # mutate(SQOL18_REL = ifelse(is.na(SQOL18_RELFAM_SATISF) | 
  #                       is.na(SQOL18_AMIS_SATISF) | 
  #                       is.na(SQOL18_VIESENTIM_SATISF),
  #                     NA,
  #                     SQOL18_RELFAM_SATISF+SQOL18_AMIS_SATISF+SQOL18_VIESENTIM_SATISF)) %>%
  #  "SQOL18_REL"#SQOL18_RELFAM_SATISF","SQOL18_AMIS_SATISF","SQOL18_VIESENTIM_SATISF"

  #mutate(across(c(Comorb, Rx, MARGIN_hx), ~ factor(.))) %>%
  mutate(across(c(Rx, MARGIN_hx), ~ factor(.))) %>%
  
  dplyr::select(all_of(col_to_keep),-all_of(old_var), all_of(new_var))

#relevel ADRESSR_cat
merged_basicRegp$ADRSSR_cat <- relevel(merged_basicRegp$ADRSSR_cat, ref="Professionnel de santé du secteur public") 

#plot var-outcome (inc.missing) ----
if (plot_var_outcome) {
# Function to create appropriate plot based on variable type
plot_variable <- function(data, x_var, y_var) {
  if (is.factor(data[[x_var]]) || is.character(data[[x_var]])) {
    # For categorical variables, use box plot
    ggplot(data, aes_string(y = x_var, x = y_var)) +
      geom_boxplot() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste(y_var, "vs", x_var), y = x_var, x = y_var)
  } else {
    # For continuous variables, use scatter plot
    ggplot(data, aes_string(y = x_var, x = y_var, group = y_var)) +
      geom_boxplot() +
      labs(title = paste(y_var, "vs", x_var), y = x_var, x = y_var)
  }
}

# Create a list to store all plots
plot_list <- list()

# Loop through all predictor variables
for (var in colnames(merged_basicRegp)) {
  plot_list[[var]] <- plot_variable(merged_basicRegp, var, "missing_out")
}
plot_list[[2]]
for (var in colnames(merged_basicRegp)) {
  plot <- plot_variable(merged_basicRegp, var, "missing_out")
  print(plot)
  #readline(prompt="Press [enter] to continue")
}
}

#select cases ----
#based on n_miss in case and calculate attrition weights 
miss_case_summary_obs <- try_missing(merged_basicRegp)
#eg observations with more than 30% NA
keep_id <- miss_case_summary_obs %>%
  filter(pct_miss<=p_miss_indiv) %>%
  dplyr::select(case) %>% unlist
merged_basicRegp[keep_id,] -> merged_basicRegp

#study NA pattern ----
md_pattern <- md.pattern(merged_basicRegp)
write.csv(md_pattern,"md_pattern.csv")
md.pairs(merged_basicRegp)
(fx <- flux(merged_basicRegp)); print(fx)
plot(fx$influx, fx$outflux, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Influx", ylab = "Outflux", main = "Flux Plot")
text(fx$influx, fx$outflux, row.names(fx), pos = 4, cex = 0.8)

# Perform a dry run ----
dryrun <- mice(merged_basicRegp, maxit = 0, print = FALSE)

# Inspect the results
# 1. Check the imputation methods chosen for each variable
print(dryrun$method)

# 2. Examine the predictor matrix
print(dryrun$predictorMatrix)

# 3. Look at the number of missing values per variable
print(dryrun$nmis)

# 4. Look at the logged events produced by mice(). 
#Remove any constant and collinear variables before imputation.
print(dryrun$loggedEvents)

#create folds ----
custom_folds <- createFolds(merged_basicRegp$missing_out, k = k_outer, returnTrain = TRUE)

#loop over outer folds
for (outer_idx in 1:k_outer) {
  ##run mice with ignore argument- #train will be FALSE, test will be TRUE
  train_data <- merged_basicRegp[custom_folds[[outer_idx]],]
  test_data <- merged_basicRegp[-custom_folds[[outer_idx]],]
  
  logical_vector <- rep(TRUE, nrow(merged_basicRegp))
  logical_vector[custom_folds[[outer_idx]]] <- FALSE
  
  ##Method 1 when ignore=FALSE(training set), imputation model will use the observation
  merged_imputed <- mice(merged_basicRegp, m = m, maxit=maxit, ignore=logical_vector, seed=1, print = FALSE)#
  
  ## Imputation dx ----
  warnings()
  merged_imputed$loggedEvents
  colSums(is.na(merged_imputed %>% complete("long")))
  colSums(is.na(merged_imputed %>% complete(1)))
  
  if (plot_imp) {
    plot(merged_imputed) #convergence
    stripplot(merged_imputed)#values for imputed datasets imputed and non-imputed points
    densityplot(merged_imputed)
  }
  
  
  for (m_idx in 1:m) {
    imp.test <- filter(merged_imputed, logical_vector) %>% complete(m_idx)
    imp.train <- filter(merged_imputed, !logical_vector) %>% complete(m_idx)
    
    ##Method 2 (does not give the same results)
    # mod.train <- mice(train_data, m = 1, print = FALSE, seed = 1); imp.train <- mod.train %>% complete
    # imp.test <- mice.mids(mod.train, newdata = test_data, print=FALSE) %>% complete
    
    ## calculate lambda/fmi ----
    predictors <- paste(imp.train %>% dplyr::select(-missing_out) %>% colnames, collapse = "+")
    formula <- formula(paste0("missing_out ~", predictors))
    calculate_fmi(merged_imputed, formula)
    calculate_fmi(filter(merged_imputed, logical_vector), formula)
    calculate_fmi(filter(merged_imputed, !logical_vector), formula)
    
    
    ## SL ----
    test <- CV.SuperLearner(Y = imp.train$missing_out, 
                            X = imp.train %>% dplyr::select(-missing_out),
                            cvControl =list(V=V, stratifyCV=TRUE),
                            innerCvControl = list(list(V=V, stratifyCV=TRUE)),
                            SL.library = SL.library,
                            family=binomial(),
                            verbose = TRUE, 
                            method = method
    )
    summary(test)
  
    
    ## SHAP ----
    sl_fit <- SuperLearner(Y = imp.train$missing_out, 
                           X = imp.train %>% dplyr::select(-missing_out),
                           cvControl =list(V=V, stratifyCV=TRUE),
                           SL.library = SL.library,
                           family=binomial(),
                           verbose = TRUE, 
                           method = method
    )
    #save.image(file=paste0("SL_anal_", today, ".RData"))
    test_SL_predict <- predict(sl_fit,imp.test)$pred
    
    #calculate AUC on hold out fold
    auc <- roc(imp.test$missing_out, test_SL_predict)$auc
    
    # Calculate accuracy
    predicted_classes <- ifelse(test_SL_predict > 0.5, 1, 0)
    accuracy <- mean(predicted_classes == imp.test$missing_out)
    
    print(paste("AUC:", auc))
    print(paste("Accuracy:", accuracy))
    
    pred_wrapper <- function(object, newdata) {
      predict(object, newdata = newdata)$pred %>% as.vector
    }
    
    Shap_values <- explain(X = imp.train %>% dplyr::select(-missing_out),
                           pred_wrapper = pred_wrapper,
                           object=sl_fit,
                           nsim = nsim)  # Adjust nsim as needed
    #save.image(file=paste0("SHAP_anal_", today, ".RData"))
    shap_df <- as.data.frame(Shap_values)
    
    #Calculate mean absolute SHAP values to determine feature importance:
    mean_shap <- colMeans(abs(shap_df))
    sorted_features <- names(sort(mean_shap, decreasing = TRUE))
    
    #Select the top N most important features (e.g., top 5):
    top_features <- head(sorted_features, top_n_select)
    
    #Create dependence plots for each top feature:
    shap_viz <- shapviz(Shap_values, 
                        X = imp.train %>% dplyr::select(-missing_out)
    )
    
    plots <- lapply(top_features, function(feature) {
      
      dep_plot <- sv_dependence(shap_viz, v = feature, color_var = "auto")
      
      shap_data <- data.frame(
        feature_value = shap_viz$X[[feature]],
        shap_value = as.data.frame(shap_viz$S)[[feature]]
      )
      wrap_plots(dep_plot, ncol = 1)
      
    })
    
    combined_plot <- wrap_plots(plots, ncol = 2)
    print(combined_plot)
    
  }
}
