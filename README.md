# Background
Missing data in observational research is a curse  
Ideally, missing data would be predicted as soon as patient is referred to the service  
Patient comes with background factors: Age Sexe Hx Rx etc ...  
Can these background factors be used to predict the occurence of missing values in rehab outcomes?  
Good Predictions of the occurence of missing values based on background factors may  
1- inform clinicians on whether their next patient may be at risk of missing values, and therefore act accordingly (e.g. by re-calling patients)  
2- Identifying predictors of missingness helps researchers understand potential biases in their data and develop strategies to minimize missing data in future studies  
3- Improving data quality would lead to more accurate treatment effect estimates in research studies  
4- Awareness of missing data patterns and their predictors allows for more nuanced interpretation of study findings, acknowledging potential biases or limitations  
5- Analyzing predictors of missingness can provide valuable information about the profiles of service users and their engagement with rehabilitation programs  

# REHABase 
is a multicentric rehab database in France  
When referred, patients are seen by  
1- the physician who gathers background data (Age, sex, hx, Rx etc ...)  
2- the nurse who collects rehab outcome data at baseline e.g. well-being, quality of life, adherence to treatment, self-esteem, insight, autonomy, self-stigmatisation, recovery  
3- the neuropsychologist who collects neuropsych data e.g. IQ, attention or social cognition tests  
4- Then rehab outcomes (and possibly neuropsych data) are re-assessed after one or many interventions are run and most often at 1 year

# RQ
1- Can SCZ missing rehab outcome data at baseline be predicted based on background factors?  
2- What are the specific influence of these background factors on missingness?  
3- Do the predictive model and variable importance vary per questionnaire/test missing 
We will not be able to run the analysis at t1 as there are too many missing values (>90% - wouldn't make any sense to better understand determinants of missing values when the number of missing values is so high). Thus the analysis will be run at baseline (t0) only.

# Data  
- as collected in REHABase
- background factors
- typical cleaning procedure, e.g. regroup levels that have less than 30 observations in variables
- Variables  
              "CENTRE",# the patient's rehabilitation center  
              "AGE_MEDSOC",#age  
              "SEXE",#sex  
              "NIVETUD_cat",#level of education                
              "COMOR_PSY", #presence or absence of psychiatric comorbity  
              "Rx",# being on FGA, SGA, both or none   
              "addict", #having a behavioural or substance addiction comorbity  
              "SOMA_1", #having a physical comorbity  
              "TTTSOMA_1", #being on treatment for a physical problem  
              "EGF",	#Global Assessment of Functioning score  
              "CGI_SEVERITE", #Clinical Global Impression score  
              "SIT_FAM_cat", #Marital status  
              "ETRE_PARENT",#Being a parent  
              "ADRSSR_cat",#which structure referred this patient to the rehab center  
              "LGMT_cat",#what type of accomodation does the patient lives in  
              "SIT_PRO_cat",#what is their professional status  
              "RQTH",#do they have a current status of disabled worker  
              "DUREE_MALADIE",#how long have they been psychiatrically ill  
              "NBR_HOSPI",#how many times have they been admitted in psychiatry  
              "DUREE_HOSPI",#how long have they been admitted in psychiatry  
              "TS",#have they ever tried to commit suicide?  
              "NBR_TS",#how many times have they tried to commit suicide?  
              "MARGIN_ACTPASS",#have they been marginalized (of no fixed above) currently or in the past  
              "ATCD_MEDLEG",#do they have a forensic history  

# Strategy for missing data
- re- missing data in predictors, we plan on imputing missing data using R package mice
- amount to keep will depend on various parameters such as: outflux-influx plot, fmi/lambda parameters, pct of missing values in observations. Ideally, we'll keep all missing observations.
- For instance, we'll most likely remove variables that are on the right handside corner of the outflux-influx plot; we'll carefully inspect fraction of missing information and lambda parameters so that included  variables do not increase the variability of coefficients over imputed datasets; if missing values occur over a very large proportions of predictors, we'll have to remove the observation itself.
- for instance, if all variables in a participants are missing, we're not going to keep that observation. Probably similarly if only age or sex are non-missing.
- in case where observations are discarded while there are still _some_ predictors (more than simply age or sex), then we'll calculate attrition weights using a SuperLearner prediction of being included in the study given observed variables
- auxiliary variables: we are not planning on using auxiliary variables simply because all variables roughly have the same amount of missing values. None can be useful as auxiliary variable.
- regarding the number of imputations m, we'll use the rule of thumb m=%age of missing data
- regarding the number of iterations maxit, we'll use maxit=20
- regarding the imputation model, we'll use mice standard argument (pmm for continuous variables, logreg for binary variables and polyreg for categorical non-binary variables)
- Importantly, imputation will precede predictive modeling, but to avoid any data leaking, we'll split imputation using training vs. testing sets (using the ignore argument in the mice function)

# Predictive modeling 
- After having imputed missing values, our task will be to predict missingness in the outcome based on the predictors
- We'll use a SuperLearner ensemble model (from the R package SuperLearner), using a broad variety of basis learners: glm, glm with one-way interaction, regression regularization, MARS, random forest, extreme gradient boosting, support vector machine.
- for each basis learner, we'll use the standard model as well as ex-ante screens based on correlation coefficients (p>0.1), regression regularization and random forest importance values (10 best predictors)
- other than that, we'll use the R package caret adaptive scheme to find the set of hyperparameters that minimize a log loss function (tuneLength=10, nfolds=10)
- based on Philipps' (2021) guidelines, we aim to use V=10-20 folds
- we'll use the option stratifyCV=TRUE so that the outcome rate is identical across folds 
- We'll use the CIMENT server of the university of Grenoble, France, hoping that such an aggressive strategy will be handled without any overwhelming of resources
- If resources do become overwhelmed, then we'll decrease V (e.g. to 10 folds), and or decrease the number of algorithms/screens (e.g. we won't use screening algorithms when using regression regularization)
- To better understand whether missing values make the calculation of predictive accuracy and importance values unstable, we'll calculate the sd/variance/range of predictive accuracy and importance values across imputed datasets
- We'll fit our model to a training set (70-80% of the data: 80% if N>2000, 70% if N<2000 - set arbitrarily) and evaluate its performance on the testing set. 
   
# Variable importance
- Using the fastshap R package, SHAP value will be calculated for each fold and each training observation (using nsim=100)
- SHAP plots will be made with the shapviz R package 

# set.seed
We'll specify set.seed=123 each time a seed is requested

# sensitvity analysis  
- With/Without variables with high FMI/lambda (on pre-test: duration and number of psych. admission, duration of illness)
- Variable amount of missing values in the final dataset depending on the parameters mentioned above

# Results
- Calculate fit over 30 imputed datasets
- Calculate SHAP values over 30 imputed datasets
- Detect potential interactions using SHAP values and the shapviz::sv_interaction() function

