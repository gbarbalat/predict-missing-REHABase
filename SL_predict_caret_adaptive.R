#happens that screen.glmnet => no variables selected ==> gives an error message! ==> screen.glmnet removed 
SL.library <- list("SL.mean",
                   c("SL.caret.xgboost", "All", "screen.mixed", "screen.randomForest"),
                   c("SL.caret.ranger", "All","screen.mixed",  "screen.randomForest"),
                   #c("SL.rpart", "All","screen.mixed",   "screen.randomForest"),#bratMachine not working with caret, rpart somehow useless 
                   #c("SL.caret.ksvm", "All","screen.mixed",  "screen.randomForest"),#pbs with predictions in caret
                   c("SL.caret.naive_bayes", "All", "screen.mixed", "screen.randomForest"),
                   #c("SL.caret.earth", "All", "screen.mixed", "screen.randomForest"),#only for cont predictors....
                   c("SL.caret.glm", "All", "screen.mixed", "screen.randomForest"),
                   c("SL.caret.step.interaction", "screen.mixed", "screen.randomForest"),#,
                   c("SL.caret.glmnet", "All", "screen.mixed",  "screen.randomForest")
)

# SL fns using caret ----
screen.corP <- function (Y, X, family, obsWeights, id, method = "pearson", minPvalue = 0.2, 
          minscreen = 2, ...) 
{
  listp <- apply(X, 2, function(x, Y, method) {
    ifelse(var(x) <= 0, 1, cor.test(x, y = Y, method = method)$p.value)
  }, Y = Y, method = method)
  whichVariable <- (listp <= minPvalue)
  if (sum(whichVariable) < minscreen) {
    warning("number of variables with p value less than minPvalue is less than minscreen")
    whichVariable[rank(listp) <= minscreen] <- TRUE
  }
  whichVariable
  return(whichVariable)
}

screen.corRank <- function (Y, X, family, method = "pearson", rank = 2, ...) 
{
  listp <- apply(X, 2, function(x, Y, method) {
    ifelse(var(x) <= 0, 1, cor.test(x, y = Y, method = method)$p.value)
  }, Y = Y, method = method)
  whichVariable <- (rank(listp) <= rank)
  return(whichVariable)
}

screen.mixed <- function(Y, X, family, obsWeights, id, method = "pearson", minPvalue = 0.2, minscreen = 2, ...) 
{
  listp <- sapply(X, function(x) {
    if (is.numeric(x)) {
      ifelse(var(x) <= 0, 1, cor.test(x, y = Y, method = method)$p.value)
    } else if (is.factor(x)) {
      chisq.test(table(x, Y))$p.value
    } else {
      #1  # For other types, assume no correlation
      chisq.test(table(x, Y))$p.value
    }
  })
  
  whichVariable <- (listp <= minPvalue)
  if (sum(whichVariable) < minscreen) {
    warning("number of variables with p value less than minPvalue is less than minscreen")
    whichVariable[rank(listp) <= minscreen] <- TRUE
  }
  return(whichVariable)
}

screen.glmnet <- function (Y, X, family, alpha = 1, minscreen = 2, nfolds = 10, 
                           nlambda = 100, ...) 
{
  original_names <- colnames(X)
  
  if (!is.matrix(X)) {
    X_matrix <- model.matrix(~-1 + ., X)
  } else {
    X_matrix <- X
  }
  
  fitCV <- glmnet::cv.glmnet(x = X_matrix, y = Y, lambda = NULL, type.measure = "deviance", 
                             nfolds = nfolds, family = binomial(), alpha = alpha, 
                             nlambda = nlambda)
  
  whichVariable <- (as.numeric(coef(fitCV$glmnet.fit, s = fitCV$lambda.min))[-1] != 0)
  
  if (sum(whichVariable) < minscreen) {
    warning("fewer than minscreen variables passed the glmnet screen, increased lambda to allow minscreen variables")
    sumCoef <- apply(as.matrix(fitCV$glmnet.fit$beta), 2, 
                     function(x) sum((x != 0)))
    newCut <- which.max(sumCoef >= minscreen)
    whichVariable <- (as.matrix(fitCV$glmnet.fit$beta)[, newCut] != 0)
  }
  
  # Map the selected variables back to the original variable names
  if (!is.matrix(X)) {
    selected_vars <- colnames(X_matrix)[whichVariable]
    original_vars <- sapply(strsplit(selected_vars, ":"), function(x) x[1])
    whichVariable_original <- original_names %in% unique(original_vars)
  } else {
    whichVariable_original <- whichVariable
  }
  whichVariable_original
  return(whichVariable_original)
}


screen.randomForest <- function (Y, X, family, nVar = 10, ntree = 1000, mtry = ifelse(family$family == 
                                                                 "gaussian", floor(sqrt(ncol(X))), max(floor(ncol(X)/3), 1)), 
          nodesize = ifelse(family$family == "gaussian", 5, 1), maxnodes = NULL, 
          ...) 
{
  if (family$family == "gaussian") {
    rank.rf.fit <- randomForest::randomForest(Y ~ ., data = X, 
                                              ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                              keep.forest = FALSE, maxnodes = maxnodes)
  }
  if (family$family == "binomial") {
    rank.rf.fit <- randomForest::randomForest(as.factor(Y) ~ 
                                                ., data = X, ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                              keep.forest = FALSE, maxnodes = maxnodes)
  }
  whichVariable <- (rank(-rank.rf.fit$importance) <= nVar)
  return(whichVariable)
}

SL.gam_cts <- function (Y, X, newX, family, obsWeights, deg.gam = 2, cts.num = 11, 
          ...) 
{
  if (!requireNamespace("gam")) {
    stop("SL.gam requires the gam package, but it isn't available")
  }
  if (!"package:gam" %in% search()) 
    attachNamespace("gam")
  if ("mgcv" %in% loadedNamespaces()) 
    warning("mgcv and gam packages are both in use. You might see an error because both packages use the same function names.")
  cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
  if (sum(!cts.x) > 0) {
    gam.model <- as.formula(paste("Y~", paste(paste("s(", 
                                                    colnames(X[, cts.x, drop = FALSE]), ",", deg.gam, 
                                                    ")", sep = ""), collapse = "+"), "+", paste(colnames(X[, 
                                                                                                           !cts.x, drop = FALSE]), collapse = "+")))
  }
  else {
    gam.model <- as.formula(paste("Y~", paste(paste("s(", 
                                                    colnames(X[, cts.x, drop = FALSE]), ",", deg.gam, 
                                                    ")", sep = ""), collapse = "+")))
  }
  if (sum(!cts.x) == length(cts.x)) {
    gam.model <- as.formula(paste("Y~", paste(colnames(X), 
                                              collapse = "+"), sep = ""))
  }
  fit.gam <- gam::gam(gam.model, data = X, family = family, 
                      control = gam::gam.control(maxit = 50, bf.maxit = 50), 
                      weights = obsWeights)
  if (packageVersion("gam") >= "1.15") {
    pred <- gam::predict.Gam(fit.gam, newdata = newX, type = "response")
  }
  else {
    stop("This SL.gam wrapper requires gam version >= 1.15, please update the gam package with 'update.packages('gam')'")
  }
  fit <- list(object = fit.gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.gam_cts")
  return(out)
}

predict.SL.gam_cts <- function (object, newdata, ...) 
{
  .SL.require("gam")
  if (packageVersion("gam") >= "1.15") {
    pred <- gam::predict.Gam(object = object$object, newdata = newdata, 
                             type = "response")
  }
  else {
    stop("This SL.gam wrapper requires gam version >= 1.15, please update the gam package with 'update.packages('gam')'")
  }
  return(pred)
}

SL.caret.glm <- function (Y, X, newX, family, obsWeights, method = "glm", tuneLength = 10, 
                           trControl = caret::trainControl(method = "cv", 
                                                           search = "random",
                                                           summaryFunction = mnLogLoss,
                                                           number = 10, 
                                                           verboseIter = TRUE,
                                                           classProbs = TRUE,
                                                           adaptive = list(min = 2, alpha = 0.05, method = "gls", complete = TRUE)
                           ), 
                          preProcess = NULL, #c("center", "scale"), NULL,
                           metric = ifelse(family$family == "gaussian", "RMSE", "logLoss"), ...) 
{
  
  X <- as.data.frame(X)
  newX <- as.data.frame(newX)
  
  if (family$family == "gaussian") {
    fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                              metric = metric, method = method, tuneLength = tuneLength, 
                              preProcess=preProcess,
                              trControl = trControl, 
                              verbose=FALSE, verbosity=0)
    pred <- predict(fit.train, newdata = newX, type = "raw")
  }
  if (family$family == "binomial") {
    Y.f <- as.factor(Y)
    levels(Y.f) <- c("A0", "A1")
    fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                              metric=metric, method = method, tuneLength = tuneLength, 
                              preProcess=preProcess,
                              trControl = trControl)
    pred <- predict(fit.train, newdata=newX, type = "prob")[, 
                                                            2]
  }
  fit <- list(object = fit.train)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.caret.glm")
  return(out)
}

predict.SL.caret.glm <- function(object, newdata, family, ...) {
  newdata <- as.data.frame(newdata)
  
  if (family$family == "gaussian") {
    pred <- predict(object$object, newdata = newdata, type = "raw")
  } else if (family$family == "binomial") {
    pred <- predict(object$object, newdata = newdata, type = "prob")[, 2]
  } else {
    stop("Unsupported family")
  }
  
  return(pred)
}





SL.caret.gam_cts <- function(Y, X, newX, family, obsWeights, method="gam", tuneLength = 10,
                             deg.gam = "df = 2", cts.num = 11,
                             trControl = caret::trainControl(method = "adaptive_cv", 
                                                             search = "random",
                                                             summaryFunction = mnLogLoss,
                                                             number = 10, 
                                                             verboseIter = TRUE,
                                                             classProbs = TRUE,
                                                             adaptive = list(min = 2, alpha = 0.05, method = "gls", complete = TRUE)
                             ), 
                             metric = ifelse(family$family == "gaussian", "RMSE", "logLoss"), ...) 
{
  # Check for required packages
  if (!requireNamespace("gam")) {
    stop("SL.caret.gam_cts requires the gam package, but it isn't available")
  }
  
  X <- as.data.frame(X)
  newX <- as.data.frame(newX)
  Y.f <- as.factor(Y)
  levels(Y.f) <- c("A0", "A1")
  
  # Identify continuous predictors
  cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
  
  # Create smooth terms for continuous predictors
  for (col in colnames(X)[cts.x]) {
    X[[paste0(col, "_smooth")]] <- gam::s(X[[col]], df = deg.gam)
    newX[[paste0(col, "_smooth")]] <- gam::s(newX[[col]], df = deg.gam)
  }
  
  # Remove original continuous columns from X and newX
  X <- X[, !cts.x]  # Keep only non-continuous predictors
  newX <- newX[, !cts.x]  # Keep only non-continuous predictors
  
  # Combine smooth terms with remaining predictors
  data <- data.frame(Y.f = Y.f, X)
  
  # Create the formula for the GAM model using smooth terms
 # smooth_terms <- paste(paste0(colnames(X), "_smooth"), collapse = "+")
  linear_terms <- paste(colnames(X), collapse = "+")

    gam.model <- as.formula(paste("Y.f ~", linear_terms))
  
  # Train the GAM model using caret
  fit.train <- train(gam.model,
                     data = data,
                     weights = obsWeights, 
                     metric = metric,
                     method = method,
                     tuneLength = tuneLength, 
                     trControl = trControl,
                     ...)
  
  
  pred <- predict(fit.train, newdata = newX, type = "raw")[,2]
  
  fit <- list(object = fit.train)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.caret.gam_cts")
  return(out)
}


predict.SL.caret.gam_cts <- function(object, newdata, family, ...) {
  #newdata <- model.matrix(~. -1, newdata)
  newdata <- as.data.frame(newdata)
  
  if (family$family == "gaussian") {
    pred <- predict(object$object, newdata = newdata, type = "raw")
  } else if (family$family == "binomial") {
    pred <- predict(object$object, newdata = newdata, type = "prob")[, 2]
  } else {
    stop("Unsupported family")
  }
  
  return(pred)
}


SL.caret.step.interaction <- function (Y, X, newX, family, obsWeights, method = "glmStepAIC", tuneLength = 10, 
                          trControl = caret::trainControl(method = "cv", 
                                                          search = "random",
                                                          summaryFunction = mnLogLoss,
                                                          number = 10, 
                                                          verboseIter = TRUE,
                                                          classProbs = TRUE,
                                                          adaptive = list(min = 2, alpha = 0.05, method = "gls", complete = TRUE)
                          ), 
                          metric = ifelse(family$family == "gaussian", "RMSE", "logLoss"), ...) 
{
  
  # X <- model.matrix(~. -1, X)
  # newX <- model.matrix(~. -1, newX)
  X <- as.data.frame(X)
  newX <- as.data.frame(newX)
  
  if (family$family == "gaussian") {
    fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                              metric = metric, method = method, tuneLength = tuneLength, 
                              trControl = trControl)
    pred <- predict(fit.train, newdata = newX, type = "raw")
  }
  if (family$family == "binomial") {
    Y.f <- as.factor(Y)
    levels(Y.f) <- c("A0", "A1")
    fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                              metric=metric, method = method, tuneLength = tuneLength, 
                              trControl = trControl, 
                              direction = "both", trace = 0, k = 2,scope = Y.f ~ .^2
                              )
    pred <- predict(fit.train, newdata=newX, type = "prob")[, 
                                                            2]
  }
  fit <- list(object = fit.train)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.caret.step.interaction")
  return(out)
}

predict.SL.caret.step.interaction <- function(object, newdata, family, ...) {
  #newdata <- model.matrix(~. -1, newdata)
  newdata <- as.data.frame(newdata)
  
  if (family$family == "gaussian") {
    pred <- predict(object$object, newdata = newdata, type = "raw")
  } else if (family$family == "binomial") {
    pred <- predict(object$object, newdata = newdata, type = "prob")[, 2]
  } else {
    stop("Unsupported family")
  }
  
  return(pred)
}


SL.caret.ksvm <- function (Y, X, newX, family, obsWeights, method = "svmRadial", tuneLength = 10, 
                             trControl = caret::trainControl(method = "adaptive_cv", 
                                                             search = "random",
                                                             summaryFunction = mnLogLoss,
                                                             number = 10, 
                                                             verboseIter = TRUE,
                                                             classProbs = TRUE,
                                                             adaptive = list(min = 2, alpha = 0.05, method = "gls", complete = TRUE)
                             ), 
                             metric = ifelse(family$family == "gaussian", "RMSE", "logLoss"), ...) 
{
  
  X <- model.matrix(~. -1, X)
  newX <- model.matrix(~. -1, newX)
  X <- as.data.frame(X)
  newX <- as.data.frame(newX)
  
  if (family$family == "gaussian") {
    fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                              metric = metric, method = method, tuneLength = tuneLength, 
                              trControl = trControl, 
                              verbose=FALSE, verbosity=0)
    pred <- predict(fit.train, newdata = newX, type = "raw")
  }
  if (family$family == "binomial") {
    Y.f <- as.factor(Y)
    levels(Y.f) <- c("A0", "A1")
    data <- cbind(Y.f,X)
    fit.train <- caret::train(#x = X, y = Y.f, 
                              Y.f ~ ., data=data,
                              weights = obsWeights, 
                              metric=metric, method = method, tuneLength = tuneLength, 
                              trControl = trControl,
                              verbose=FALSE, verbosity=0)
    pred <- predict(fit.train, newdata=newX, type = "prob")[, 
                                                            2]
  }
  fit <- list(object = fit.train)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.caret.ksvm")
  return(out)
}

predict.SL.caret.ksvm<- function(object, newdata, family, ...) {
  newdata <- model.matrix(~. -1, newdata)
  newdata <- as.data.frame(newdata)
  
  if (family$family == "gaussian") {
    pred <- predict(object$object, newdata = newdata, type = "raw")
  } else if (family$family == "binomial") {
    pred <- predict(object$object, newdata = newdata, type = "prob")[, 2]
  } else {
    stop("Unsupported family")
  }
  
  return(pred)
}

SL.caret.glmnet <- function (Y, X, newX, family, obsWeights, method = "glmnet", tuneLength = 10, 
          trControl = caret::trainControl(method = "adaptive_cv", 
                                          search = "random",
                                          summaryFunction = mnLogLoss,
                                          number = 10, 
                                          verboseIter = TRUE,
                                          classProbs = TRUE,
                                          adaptive = list(min = 2, alpha = 0.05, method = "gls", complete = TRUE)
                                          ), 
          metric = ifelse(family$family == "gaussian", "RMSE", "logLoss"), ...) 
{
  
  X <- model.matrix(~. -1, X)
  newX <- model.matrix(~. -1, newX)
  X <- as.data.frame(X)
  newX <- as.data.frame(newX)

    if (family$family == "gaussian") {
    fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                              metric = metric, method = method, tuneLength = tuneLength, 
                              trControl = trControl, 
                              verbose=FALSE, verbosity=0)
    pred <- predict(fit.train, newdata = newX, type = "raw")
  }
  if (family$family == "binomial") {
    Y.f <- as.factor(Y)
    levels(Y.f) <- c("A0", "A1")
    fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                              metric=metric, method = method, tuneLength = tuneLength, 
                              trControl = trControl,
                              verbose=FALSE, verbosity=0)
    pred <- predict(fit.train, newdata=newX, type = "prob")[, 
                                                              2]
  }
  fit <- list(object = fit.train)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.caret.glmnet")
  return(out)
}

predict.SL.caret.glmnet <- function(object, newdata, family, ...) {
  newdata <- model.matrix(~. -1, newdata)
  newdata <- as.data.frame(newdata)
  
  if (family$family == "gaussian") {
    pred <- predict(object$object, newdata = newdata, type = "raw")
  } else if (family$family == "binomial") {
    pred <- predict(object$object, newdata = newdata, type = "prob")[, 2]
  } else {
    stop("Unsupported family")
  }
  
  return(pred)
}

SL.caret.xgboost <- function (Y, X, newX, family, obsWeights, method = "xgbTree", tuneLength = 10, 
                              trControl = caret::trainControl(method = "adaptive_cv", 
                                                              search = "random",
                                                              summaryFunction = mnLogLoss,
                                                              number = 10, 
                                                              verboseIter = TRUE,
                                                              classProbs = TRUE,
                                                              seeds = seeds,
                                                              adaptive = list(min = 2, alpha = 0.05, method = "gls", complete = TRUE)
                              ), 
                              metric = ifelse(family$family == "gaussian", "RMSE", "logLoss"), ...) 
{
  
  X <- model.matrix(~. -1, X)
  newX <- model.matrix(~. -1, newX)
  X <- as.data.frame(X)
  newX <- as.data.frame(newX)
  
  if (family$family == "gaussian") {
    fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                              metric = metric, method = method, tuneLength = tuneLength, 
                              trControl = trControl, 
                              verbose=FALSE, verbosity=0)
    pred <- predict(fit.train, newdata = newX, type = "raw")
  }
  if (family$family == "binomial") {
    Y.f <- as.factor(Y)
    levels(Y.f) <- c("A0", "A1")
    fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                              metric = metric, method = method, tuneLength = tuneLength, 
                              trControl = trControl,
                              verbose=FALSE, verbosity=0)
    pred <- predict(fit.train, newdata = newX, type = "prob")[, 
                                                              2]
  }
  fit <- list(object = fit.train)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.caret.xgboost")
  return(out)
}

predict.SL.caret.xgboost <- function(object, newdata, family, ...) {
  newdata <- model.matrix(~. -1, newdata)
  newdata <- as.data.frame(newdata)
  
  if (family$family == "gaussian") {
    pred <- predict(object$object, newdata = newdata, type = "raw")
  } else if (family$family == "binomial") {
    pred <- predict(object$object, newdata = newdata, type = "prob")[, 2]
  } else {
    stop("Unsupported family")
  }
  
  return(pred)
}




SL.caret.ranger <- function (Y, X, newX, family, obsWeights, method = "ranger", tuneLength = 10, 
                             trControl = caret::trainControl(method = "adaptive_cv", 
                                                             search = "random",
                                                             summaryFunction = mnLogLoss,
                                                             number = 10, 
                                                             verboseIter = TRUE,
                                                             seeds = seeds,
                                                             classProbs = TRUE,
                                                             adaptive = list(min = 2, alpha = 0.05, method = "gls", complete = TRUE)
                             ), 
                             metric = ifelse(family$family == "gaussian", "RMSE", "logLoss"), ...) 
{
  X <- as.data.frame(X)
  newX <- as.data.frame(newX)
  
  
  if (family$family == "gaussian") {
    fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                              metric = metric, method = method, tuneLength = tuneLength,
                              trControl = trControl,
                              importance="permutation"
                              )
    pred <- predict(fit.train, newdata = newX, type = "raw")
  }
  if (family$family == "binomial") {
    Y.f <- as.factor(Y)
    levels(Y.f) <- c("A0", "A1")
    fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                              metric = metric, method = method, tuneLength = tuneLength, 
                              trControl = trControl,
                              importance="permutation"
    )
    pred <- predict(fit.train,newX, type = "prob")[, 2]
  }
  fit <- list(object = fit.train)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.caret.ranger")
  return(out)
}

predict.SL.caret.ranger <- function(object, newdata, family, ...) {

    newdata <- as.data.frame(newdata)
  
  if (family$family == "gaussian") {
    pred <- predict(object$object, newdata = newdata, type = "raw")
  } else if (family$family == "binomial") {
    pred <- predict(object$object, newdata = newdata, type = "prob")[, 2]
  } else {
    stop("Unsupported family")
  }
  
  return(pred)
}


SL.caret.earth <- function (Y, X, newX, family, obsWeights, method = "earth", tuneLength = 10, 
                             trControl = caret::trainControl(method = "adaptive_cv",
                                                             search = "random",
                                                             summaryFunction = mnLogLoss,
                                                             number = 10, 
                                                             verboseIter = TRUE,
                                                             classProbs = TRUE,
                                                             adaptive = list(min = 2, alpha = 0.05, method = "gls", complete = TRUE)
                             ), 
                             metric = ifelse(family$family == "gaussian", "RMSE", "logLoss"), ...) 
{
  X <- as.data.frame(X)
  newX <- as.data.frame(newX)
  
  
  if (family$family == "gaussian") {
    fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                              metric = metric, method = method, tuneLength = tuneLength,
                              trControl = trControl    )
    pred <- predict(fit.train, newdata = newX, type = "raw")
  }
  if (family$family == "binomial") {
    Y.f <- as.factor(Y)
    levels(Y.f) <- c("A0", "A1")
    fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                              metric = metric, method = method, tuneLength = tuneLength, 
                              trControl = trControl    )
    pred <- predict(fit.train,newX, type = "prob")[, 2]
  }
  fit <- list(object = fit.train)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.caret.earth")
  return(out)
}

predict.SL.caret.earth <- function(object, newdata, family, ...) {

  
  if (family$family == "gaussian") {
    pred <- predict(object$object, newdata = newdata, type = "raw")
  } else if (family$family == "binomial") {
    pred <- predict(object$object, newdata = newdata, type = "prob")[, 2]
  } else {
    stop("Unsupported family")
  }
  
  return(pred)
}


SL.caret.rpart <- function (Y, X, newX, family, obsWeights, method = "rpart", tuneLength = 10, 
                            trControl = caret::trainControl(method = "adaptive_cv",
                                                            search = "random",
                                                            summaryFunction = mnLogLoss,
                                                            number = 10, 
                                                            verboseIter = TRUE,
                                                            classProbs = TRUE,
                                                            adaptive = list(min = 2, alpha = 0.05, method = "gls", complete = TRUE)
                            ), 
                            metric = ifelse(family$family == "gaussian", "RMSE", "logLoss"), ...) 
{
  
  X <- model.matrix(~. -1, X)
  newX <- model.matrix(~. -1, newX)
  X <- as.data.frame(X)
  newX <- as.data.frame(newX)
  
  
  if (family$family == "gaussian") {
    fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                              metric = metric, method = method, tuneLength = tuneLength,
                              trControl = trControl    )
    pred <- predict(fit.train, newdata = newX, type = "raw")
  }
  if (family$family == "binomial") {
    Y.f <- as.factor(Y)
    levels(Y.f) <- c("A0", "A1")
    fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                              metric = metric, method = method, tuneLength = tuneLength, 
                              trControl = trControl    )
    pred <- predict(fit.train,newX, type = "prob")[, 2]
  }
  fit <- list(object = fit.train)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.caret.rpart")
  return(out)
}

predict.SL.caret.rpart <- function(object, newdata, family, ...) {
  newdata <- model.matrix(~. -1, newdata)
  newdata <- as.data.frame(newdata)
  
  if (family$family == "gaussian") {
    pred <- predict(object$object, newdata = newdata, type = "raw")
  } else if (family$family == "binomial") {
    pred <- predict(object$object, newdata = newdata, type = "prob")[, 2]
  } else {
    stop("Unsupported family")
  }
  
  return(pred)
}

SL.caret.bartMachine <- function (Y, X, newX, family, obsWeights, method = "bartMachine",
                            trControl = caret::trainControl(method = "repeatedcv",
                                                            summaryFunction = mnLogLoss,
                                                            number = 10, 
                                                            search="grid",
                                                            verboseIter = TRUE,
                                                            classProbs = TRUE,
                                                            ), 
                            metric = ifelse(family$family == "gaussian", "RMSE", "logLoss"), ...) 
{
  
  tuneGrid <- expand.grid(
    num_trees = c(50, 100, 200, 300),
    k = c(2, 3, 4, 5),
    alpha = c(0.95, 0.97, 0.99),
    beta = c(1, 2, 3),
    nu = c(1, 3, 5, 10)
  )
  tuneGrid <- expand.grid(
    num_trees = 50,
    k = 2,
    alpha = 0.95,
    beta = 1,
    nu = 1
  )
  
  X <- as.data.frame(X)
  newX <- as.data.frame(newX)
  
  
  if (family$family == "gaussian") {
    fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                              metric = metric, method = method, tuneGrid = tuneGrid, 
                              trControl = trControl, verbose=FALSE)
    pred <- predict(fit.train, newdata = newX, type = "raw")
  }
  if (family$family == "binomial") {
    Y.f <- as.factor(Y)
    levels(Y.f) <- c("A0", "A1")
    fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                              metric = metric, method = method, tuneGrid = tuneGrid, 
                              trControl = trControl, verbose=FALSE)
    pred <- predict(fit.train,newX, type = "prob")[, 2]
  }
  fit <- list(object = fit.train)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.caret.bartMachine")
  return(out)
}

predict.SL.caret.bartMachine <- function(object, newdata, family, ...) {
  
  
  if (family$family == "gaussian") {
    pred <- predict(object$object, newdata = newdata, type = "raw")
  } else if (family$family == "binomial") {
    pred <- predict(object$object, newdata = newdata, type = "prob")[, 2]
  } else {
    stop("Unsupported family")
  }
  
  return(pred)
}


SL.caret.naive_bayes <- function (Y, X, newX, family, obsWeights, method = "naive_bayes", tuneLength = 10, 
                            trControl = caret::trainControl(method = "adaptive_cv", 
                                                            search = "random",
                                                            summaryFunction = mnLogLoss,
                                                            number = 10, 
                                                            verboseIter = TRUE,
                                                            classProbs = TRUE,
                                                            adaptive = list(min = 2, alpha = 0.05, method = "gls", complete = TRUE)
                            ), 
                            metric = ifelse(family$family == "gaussian", "RMSE", "logLoss"), ...) 
{
  X <- as.data.frame(X)
  newX <- as.data.frame(newX)
  
  
  if (family$family == "gaussian") {
    fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                              metric = metric, method = method, tuneLength = tuneLength,
                              trControl = trControl    )
    pred <- predict(fit.train, newdata = newX, type = "raw")
  }
  if (family$family == "binomial") {
    Y.f <- as.factor(Y)
    levels(Y.f) <- c("A0", "A1")
    fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                              metric = metric, method = method, tuneLength = tuneLength, 
                              trControl = trControl    )
    pred <- predict(fit.train,newX, type = "prob")[, 2]
  }
  fit <- list(object = fit.train)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.caret.naive_bayes")
  return(out)
}

predict.SL.caret.naive_bayes <- function(object, newdata, family, ...) {
  
  
  if (family$family == "gaussian") {
    pred <- predict(object$object, newdata = newdata, type = "raw")
  } else if (family$family == "binomial") {
    pred <- predict(object$object, newdata = newdata, type = "prob")[, 2]
  } else {
    stop("Unsupported family")
  }
  
  return(pred)
}

SL.glmnet_1 <- function (Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10, 
          nlambda = 100, useMin = TRUE, loss = "deviance", ...) 
{
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
    newX <- model.matrix(~-1 + ., newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, 
                             lambda = NULL, type.measure = loss, nfolds = nfolds, 
                             family = family$family, alpha = alpha, nlambda = nlambda, 
                             ...)
  pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, 
                                                                    "lambda.min", "lambda.1se"))
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet_1"
  out <- list(pred = pred, fit = fit)
  return(out)
}

predict.SL.glmnet_1 <- function (object, newdata, remove_extra_cols = T, add_missing_cols = T, 
                                 ...) 
{
  if (!is.matrix(newdata)) {
    newdata <- model.matrix(~-1 + ., newdata)
  }
  original_cols = rownames(object$object$glmnet.fit$beta)
  if (remove_extra_cols) {
    extra_cols = setdiff(colnames(newdata), original_cols)
    if (length(extra_cols) > 0) {
      warning(paste("Removing extra columns in prediction data:", 
                    paste(extra_cols, collapse = ", ")))
      newdata = newdata[, !colnames(newdata) %in% extra_cols, 
                        drop = FALSE]
    }
  }
  if (add_missing_cols) {
    missing_cols = setdiff(original_cols, colnames(newdata))
    if (length(missing_cols) > 0) {
      warning(paste("Adding missing columns in prediction data:", 
                    paste(missing_cols, collapse = ", ")))
      new_cols = matrix(0, nrow = nrow(newdata), ncol = length(missing_cols))
      colnames(new_cols) = missing_cols
      newdata = cbind(newdata, new_cols)
      newdata = newdata[, original_cols, drop = FALSE]
    }
  }
  pred <- predict(object$object, newx = newdata, type = "response", 
                  s = ifelse(object$useMin, "lambda.min", "lambda.1se"))
  return(pred)
}



SL.glmnet_0 <- function (Y, X, newX, family, obsWeights, id, alpha = 0, nfolds = 10, 
                         nlambda = 100, useMin = TRUE, loss = "deviance", ...) 
{
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
    newX <- model.matrix(~-1 + ., newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, 
                             lambda = NULL, type.measure = loss, nfolds = nfolds, 
                             family = family$family, alpha = alpha, nlambda = nlambda, 
                             ...)
  pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, 
                                                                    "lambda.min", "lambda.1se"))
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet_0"
  out <- list(pred = pred, fit = fit)
  return(out)
}

predict.SL.glmnet_0 <- function (object, newdata, remove_extra_cols = T, add_missing_cols = T, 
                                 ...) 
{
  if (!is.matrix(newdata)) {
    newdata <- model.matrix(~-1 + ., newdata)
  }
  original_cols = rownames(object$object$glmnet.fit$beta)
  if (remove_extra_cols) {
    extra_cols = setdiff(colnames(newdata), original_cols)
    if (length(extra_cols) > 0) {
      warning(paste("Removing extra columns in prediction data:", 
                    paste(extra_cols, collapse = ", ")))
      newdata = newdata[, !colnames(newdata) %in% extra_cols, 
                        drop = FALSE]
    }
  }
  if (add_missing_cols) {
    missing_cols = setdiff(original_cols, colnames(newdata))
    if (length(missing_cols) > 0) {
      warning(paste("Adding missing columns in prediction data:", 
                    paste(missing_cols, collapse = ", ")))
      new_cols = matrix(0, nrow = nrow(newdata), ncol = length(missing_cols))
      colnames(new_cols) = missing_cols
      newdata = cbind(newdata, new_cols)
      newdata = newdata[, original_cols, drop = FALSE]
    }
  }
  pred <- predict(object$object, newx = newdata, type = "response", 
                  s = ifelse(object$useMin, "lambda.min", "lambda.1se"))
  return(pred)
}



SL.glmnet_0.5 <- function (Y, X, newX, family, obsWeights, id, alpha = 0.5, nfolds = 10, 
                         nlambda = 100, useMin = TRUE, loss = "deviance", ...) 
{
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
    newX <- model.matrix(~-1 + ., newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, 
                             lambda = NULL, type.measure = loss, nfolds = nfolds, 
                             family = family$family, alpha = alpha, nlambda = nlambda, 
                             ...)
  pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, 
                                                                    "lambda.min", "lambda.1se"))
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet_0.5"
  out <- list(pred = pred, fit = fit)
  return(out)
}

predict.SL.glmnet_0.5 <- function (object, newdata, remove_extra_cols = T, add_missing_cols = T, 
                                   ...) 
{
  if (!is.matrix(newdata)) {
    newdata <- model.matrix(~-1 + ., newdata)
  }
  original_cols = rownames(object$object$glmnet.fit$beta)
  if (remove_extra_cols) {
    extra_cols = setdiff(colnames(newdata), original_cols)
    if (length(extra_cols) > 0) {
      warning(paste("Removing extra columns in prediction data:", 
                    paste(extra_cols, collapse = ", ")))
      newdata = newdata[, !colnames(newdata) %in% extra_cols, 
                        drop = FALSE]
    }
  }
  if (add_missing_cols) {
    missing_cols = setdiff(original_cols, colnames(newdata))
    if (length(missing_cols) > 0) {
      warning(paste("Adding missing columns in prediction data:", 
                    paste(missing_cols, collapse = ", ")))
      new_cols = matrix(0, nrow = nrow(newdata), ncol = length(missing_cols))
      colnames(new_cols) = missing_cols
      newdata = cbind(newdata, new_cols)
      newdata = newdata[, original_cols, drop = FALSE]
    }
  }
  pred <- predict(object$object, newx = newdata, type = "response", 
                  s = ifelse(object$useMin, "lambda.min", "lambda.1se"))
  return(pred)
}

