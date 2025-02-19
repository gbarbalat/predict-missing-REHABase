# Background

*   Missing data is a key issue in observational research.
*   Aim: Predict missing rehab outcomes using background factors at referral.
*   Good predictions can:
    *   **Inform clinicians**: Identify at-risk patients and act proactively.
    *   **Help researchers**: Understand biases and minimize missing data in future.
    *   **Improve data quality**: Lead to more accurate treatment effect estimates.
    *   **Nuanced Interpretation**: Allow for more nuanced interpretation of study findings, acknowledging potential biases or limitations
    *   **Service user profiling**: Analyzing predictors of missingness can provide valuable information about the profiles of service users and their engagement with rehabilitation programs

# REHABase

*   Multicentric French rehab database.
*   Data collected upon referral:
    1.  **Physician**: Background data (age, sex, history, etc.).
    2.  **Nurse**: Rehab outcomes (well-being, QoL, etc.).
    3.  **Neuropsychologist**: Neuropsychological data (IQ, attention, etc.).
*   Rehab outcomes reassessed at follow-up (often 1 year).

# RQ (Research Questions)

1.  Can baseline missing rehab outcome data in SCZ be predicted from background?
2.  What is the specific influence of background factors on missingness?
3.  Does the model/variable importance vary per questionnaire/test?
*   Analysis at baseline (T0) only, due to excessive missing data at T1.

# Data

*   From REHABase.
*   Background factors (see variable list below).
*   Typical cleaning (regrouping rare levels).
*   **Variables**:
    *   `"CENTRE"`, `"AGE_MEDSOC"`, `"SEXE"`, `"NIVETUD_cat"`, `"COMOR_PSY"`, `"Rx"`, `"addict"`, `"SOMA_1"`, `"TTTSOMA_1"`, `"EGF"`, `"CGI_SEVERITE"`, `"SIT_FAM_cat"`, `"ETRE_PARENT"`, `"ADRSSR_cat"`, `"LGMT_cat"`, `"SIT_PRO_cat"`, `"RQTH"`, `"DUREE_MALADIE"`, `"NBR_HOSPI"`, `"DUREE_HOSPI"`, `"TS"`, `"NBR_TS"`, `"MARGIN_ACTPASS"`, `"ATCD_MEDLEG"`

# Initial analysis plan

## Strategy for missing data

*   **Imputation**: Use R `mice` package.
    *   Keep observations (if possible based on outflux-influx plots and the lambda parameter).
    *   Discard observations with near full missingness (if only age and sex are available).
    *   If some predictors remain, calculate attrition weights using SuperLearner.
*   **Auxiliary variables**: Not planned, due to similar missingness across variables.
*   `m` = % missing data. `maxit` = 20. Standard `mice` imputation models.
*   **Split sample**: Imputation will precede predictive modeling and will be done separately in each training vs. testing sets using the ignore argument in the mice function

## Predictive modeling

*   SuperLearner ensemble model (R `SuperLearner` package).
*   Basis learners: glm, glm+interactions, regularized regression, MARS, RF, XGBoost, SVM.
*   Ex-ante screens: correlations, regression regularization, RF importance.
*   `caret` adaptive hyperparameter tuning (tuneLength=10, nfolds=10).
*   V=10-20 folds.
*   `stratifyCV=TRUE`.
*   CIMENT server (Grenoble). Reduce V or algorithms/screens if needed.
*   Assess stability of accuracy/importance by calculating SD across imputed datasets.
*   Training (70-80%) / Testing split.

## Variable importance

*   `fastshap` R package: SHAP values (nsim=100) for each fold and training observation.
*   SHAP plots with `shapviz` R package.

## set.seed

*   `set.seed=123`.

## Results

*   Fit calculated over 30 imputed datasets.
*   SHAP values calculated over 30 imputed datasets.
*   Interaction detection using SHAP values and `shapviz::sv_interaction()`.

## Sensitivity analysis

*   With/Without variables with high FMI/lambda.
*   Selecting predictors based on predictive accuracy and/or SHAP values, e.g. if they vary among imputed databases.

# Changes to the initial analysis plan

*   FMI/lambda not used as criteria to select variables as we did not carry out standard glm analysis.
*   All observations and variables were kept.
*   Questioned 2 variables with large SHAP range (esp. in lower end of the distribution): "NBR_HOSPI", "DUREE_HOSPI". May want to re-do analysis w/o these variables. (Incidentally, those variables also had large FMI/lambda values and other variables had low FMI/lambda values).
*   Due to computational issues, we had to restrict to outer data partition instead of outer cross-validation (CV); training set = 70% of dataset; inner CV based on V=8 folds.
