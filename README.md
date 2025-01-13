# Background
Missing data in observational research is a curse  
Ideally, missing data would be predicted as soon as patient is referred to the service  
Patient comes with background factors: Age Sexe Hx Rx etc ...  
Can these background factors be used to predict the occurence of missing values?  
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

# RQ
1- Can missing rehab outcome data at baseline be predicted based on background factors?  
2- What are the influence of these factors in missingness?  
3- Do the predictive model and variable importance vary per diagnosis group and questionnaire/test missing  

# Data  
- as collected in REHABase
- background factor
- typical cleaning procedure, e.g. regroup levels that have less than 30 observations in variables
- re- missing data in predictors, we plan on imputing missing data
- amount to keep will depend on various parameters such as: outflux-influx plot, fmi/lambda parameters, pct of missing values in observations. Ideally, we'll keep all missing observations.
- For instance, we'll most likely remove variables that are on the right handside corner of the outflux-influx plot; we'll carefully inspect fraction of missing information and lambda parameters so that included variables do not increase the variability of coefficients over imputed datasets; if missing values occur over a very large proportions of predictors, we'll have to remove the observation itself.
- for instance, if all variables in a participants are missing, we're not going to keep that observation. Same if only age or sex are non-missing.
- in case where observations are discarded while  

# Analysis  

