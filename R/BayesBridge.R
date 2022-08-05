#' Settings for BayesBridge Logistic Regression
#' @export
#' 
setBayesBridge <- function(seed = NULL,
                           n_iter = 250,
                           bridge_exponent = 1,
                           regularizing_slab_size = 1,
                           thin = 1,
                           n_status_update = 10,
                           global_scale = 0.1,
                           coef_sampler_type = "cholesky"){
  if(is.null(seed[1])){
    seed <- as.integer(sample(100000000,1))
  }
  
  param <- list()
  
  attr(param, 'settings') <- list(
    priorfunction = 'Cyclops::createPrior',
    modelType = 'BayesBridge logistic',
    seed = seed[1],
    name = "BayesBridge Logistic Regression",
    n_iter = n_iter,
    bridge_exponent = bridge_exponent,
    regularizing_slab_size = regularizing_slab_size,
    thin = thin,
    n_status_update = n_status_update,
    global_scale = global_scale,
    coef_sampler_type = coef_sampler_type
  )
  
  attr(param, 'modelType') <- 'binary' 
  attr(param, 'saveType') <- 'RtoJson'
  
  result <- list(fitFunction = "fitBayesBridge",
                 param = param)
  class(result) <- "modelSettings"
  
  return(result)
}

#' Predict BayesBridge Logistic Regression
#' @export
#'
# predictBayesBridge <- function(plpModel, data, cohort){
#   start <- Sys.time()
#   #do this for each sample then median
#   prediction <- predictCyclopsType(
#     plpModel$model$postMeans,
#     cohort,
#     data$covariateData,
#     modelType = "logistic"
#   )
#   prediction$value <- 
#   delta <- Sys.time() - start
#   ParallelLogger::logInfo("Prediction took ", signif(delta, 3), " ", attr(delta, "units"))
#   return(prediction)
# }

predictBayesBridge <- function(plpModel, data, cohort, train = FALSE){
  ParallelLogger::logInfo('Starting prediction for test set')
  start <- Sys.time()
  
  if(class(data) == 'plpData'){
    # convert
    matrixObjects <- toSparseM(
      plpData = data, 
      cohort = cohort,
      map = plpModel$covariateImportance %>% 
        dplyr::select(.data$columnId, .data$covariateId)
    )
    newData <- matrixObjects$dataMatrix
    cohort <- matrixObjects$labels
    covariateRef <- matrixObjects$covariateRef
  } else{
    newData <- data
  }
  
  if(class(plpModel) == 'plpModel'){
    model <- plpModel$model
  } else{
    model <- plpModel
  }
  
  out <- c()
  #get posterior medians and predicted probabilities
  for(i in 1:ncol(plpModel$model$coefRaw)){
    coefficients <- plpModel$model$coefRaw[,i]
    names(coefficients) <- c("(Intercept)", covariateRef$covariateId)
    intercept <- coefficients[names(coefficients)%in%'(Intercept)']
    if(length(intercept)==0) intercept <- 0
    coefficients <- coefficients[!names(coefficients)%in%'(Intercept)']
    beta <- as.numeric(coefficients)
    value <- newData %*% beta
    value[is.na(value)] <- 0
    value <- value + intercept
    link <- function(x){
      return(1/(1 + exp(0 - x)))
      }
     value <- link(value)
     out <- cbind(out, value)
  }
  valueMed <- apply(out, 1, median)
  
  #output
  cohort$value <- valueMed
  prediction <- cohort
  
  attr(prediction, "metaData")$modelType <- 'binary'
  prediction$evaluationType <- 'Train'
  
  #system time + wrapup
  delta <- Sys.time() - start
  ParallelLogger::logInfo("Prediction took ", signif(delta, 3), " ", attr(delta, "units"))
  return(prediction)
}


#' Fit BayesBridge Logistic Regression
#' @export
#' 
fitBayesBridge <- function(trainData, param, analysisId, ...){
  settings <- attr(param, 'settings')
  
  start <- Sys.time() 
  
  # check plpData is coo format:
  if (!FeatureExtraction::isCovariateData(trainData$covariateData)){
    stop("Needs correct covariateData")
  }

  
  #remap covariates
  mappedData <- toSparseM(trainData)
  matrixData <- mappedData$dataMatrix
  labels <- mappedData$labels
  covariateRef <- mappedData$covariateRef
  
  ParallelLogger::logInfo('Running BayesBridge')
  #run BayesBridge
  model <- create_model(labels$outcomeCount, matrixData)
  prior <- create_prior(bridge_exponent = settings$bridge_exponent, #param
                        regularizing_slab_size = settings$regularizing_slab_size) #param
  bridge <- instantiate_bayesbridge(model, prior)
  
  n_iter <- settings$n_iter #param
  gibbs_output <- gibbs(bridge, 
                        n_iter = as.integer(settings$n_iter), 
                        init = list(global_scale = settings$global_scale),
                        thin = settings$thin,
                        seed = settings$seed,
                        coef_sampler_type = settings$coef_sampler_type,
                        n_status_update = settings$n_status_update)
  mcmc_samples <- gibbs_output$samples #return all posterior samples
  
  #output modelTrained
  modelTrained <- list()
  modelTrained$coefRaw <- mcmc_samples$coef
  modelTrained$coefficients <- apply(mcmc_samples$coef, 1, median)
  names(modelTrained$coefficients) <- c("(Intercept)", covariateRef$covariateId)
  modelTrained$log_likelihood <- tail(mcmc_samples$logp, 1)
  modelTrained$modelType <- "BayesBridge logistic"
  modelTrained$modelStatus <- "OK"
  attr(modelTrained, "class") <- "plpModel"
  
  #output prediction on train set
  ParallelLogger::logTrace('Getting predictions on train set')
  out <- c()
  #get posterior median predictive values
  for(i in 1:ncol(modelTrained$coefRaw)){
    coefficients <- modelTrained$coefRaw[,i]
    names(coefficients) <- c("(Intercept)", covariateRef$covariateId)
    intercept <- coefficients[names(coefficients)%in%'(Intercept)']
    if(length(intercept)==0) intercept <- 0
    coefficients <- coefficients[!names(coefficients)%in%'(Intercept)']
    beta <- as.numeric(coefficients)
    value <- matrixData %*% beta
    value[is.na(value)] <- 0
    value <- value + intercept
    link <- function(x){
      return(1/(1 + exp(0 - x)))
    }
    value <- link(value)
    out <- cbind(out, value)
  }
  valueMed <- apply(out, 1, median)
  labels$value <- valueMed
  prediction <- labels
  
  attr(prediction, "metaData")$modelType <- 'binary'
  prediction$evaluationType <- 'Train'
  
  #variable importance
  comp <- Sys.time() - start
  
  ParallelLogger::logTrace('Getting variable importance')
  varImp <- data.frame(
    covariateId = as.double(names(modelTrained$coefficients)[names(modelTrained$coefficients)!='(Intercept)']),
    value = modelTrained$coefficients[names(modelTrained$coefficients)!='(Intercept)']
  )
  
  if(sum(abs(varImp$value)>0)==0){
    ParallelLogger::logWarn('No non-zero coefficients')
    varImp <- NULL
  } else {
    ParallelLogger::logInfo('Creating variable importance data frame')
    
    variableImportance <- covariateRef %>% 
      dplyr::collect()  %>%
      dplyr::left_join(varImp, by = 'covariateId') %>%
      dplyr::mutate(covariateValue = ifelse(is.na(.data$value), 0, .data$value)) %>%
      dplyr::select(-.data$value) %>%
      dplyr::arrange(-abs(.data$covariateValue)) %>%
      dplyr::collect()
  } 
  
  variableImportance[is.na(variableImportance)] <- 0
  covariateRef$covariateValue <- variableImportance %>% 
    arrange(covariateId) %>% 
    select(covariateValue)
  covariateRef <- covariateRef %>% arrange(-abs(.data$covariateValue))
  
  result <- list(
    model = modelTrained,
    prediction = prediction,
    settings = list(
      plpDataSettings = attr(trainData, "metaData")$plpDataSettings,
      covariateSettings = attr(trainData, "metaData")$covariateSettings,
      featureEngineering = attr(trainData$covariateData, "metaData")$featureEngineering,
      tidyCovariates = attr(trainData$covariateData, "metaData")$tidyCovariateDataSettings, 
      covariateMap = NULL,
      requireDenseMatrix = F,
      populationSettings = attr(trainData, "metaData")$populationSettings,
      modelSettings = list(
        model = settings$modelType, 
        param = param,
        finalModelParameters = list(
          variance = modelTrained$priorVariance,
          log_likelihood = modelTrained$log_likelihood
        ),
        extraSettings = attr(param, 'settings')
      ),
      splitSettings = attr(trainData, "metaData")$splitSettings,
      sampleSettings = attr(trainData, "metaData")$sampleSettings
    ),
    trainDetails = list(
      analysisId = analysisId,
      cdmDatabaseSchema = attr(trainData, "metaData")$cdmDatabaseSchema,
      outcomeId = attr(trainData, "metaData")$outcomeId,
      cohortId = attr(trainData, "metaData")$cohortId,
      attrition = attr(trainData, "metaData")$attrition, 
      trainingTime = comp,
      trainingDate = Sys.Date()
    ),
    covariateImportance = covariateRef
  )
  
  
  class(result) <- 'plpModel'
  attr(result, 'predictionFunction') <- 'predictBayesBridge'
  attr(result, 'modelType') <- attr(param, 'modelType')
  attr(result, 'saveType') <- attr(param, 'saveType')
  return(result)
}