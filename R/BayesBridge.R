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

predictBayesBridge <- function(plpModel, data, cohort){
  start <- Sys.time()
  
  population <- cohort
  covariateData <- data$covariateData
  
  value <- c()
  #do this for each sample then median
  for(i in 1:ncol(plpModel$model$coefRaw)){
    coefficients <- plpModel$model$coefRaw[,i]
    names(coefficients) <- plpModel$model$coefNames
    intercept <- coefficients[names(coefficients)%in%'(Intercept)']
    if(length(intercept)==0) intercept <- 0
    coefficients <- coefficients[!names(coefficients)%in%'(Intercept)']
    coefficients <- data.frame(beta = as.numeric(coefficients),
                               covariateId = as.numeric(names(coefficients)) #!@ modified 
    )
    coefficients <- coefficients[coefficients$beta != 0, ]
    if(sum(coefficients$beta != 0)>0){
      covariateData$coefficients <- coefficients
      on.exit(covariateData$coefficients <- NULL, add = TRUE)
 
      prediction <- covariateData$covariates %>% 
        dplyr::inner_join(covariateData$coefficients, by= 'covariateId') %>% 
        dplyr::mutate(values = .data$covariateValue*.data$beta) %>%
        dplyr::group_by(.data$rowId) %>%
        dplyr::summarise(value = sum(.data$values, na.rm = TRUE)) %>%
        dplyr::select(.data$rowId, .data$value)

      prediction <- as.data.frame(prediction)
      prediction$value[is.na(prediction$value)] <- 0
      prediction$value <- prediction$value + intercept
        } else{
            warning('Model had no non-zero coefficients so predicted same for all population...')
            prediction <- population
            prediction$value <- rep(0, nrow(population)) + intercept
            }
      link <- function(x) {
        return(1/(1 + exp(0 - x)))
        }
      prediction$value <- link(prediction$value)
      value <- rbind(value, prediction$value)
  }
  
  valueMed <- apply(value, 2, median)
  prediction <- merge(population, prediction, by ="rowId", all.x = TRUE, fill = 0)
  prediction$value <- valueMed
  attr(prediction, "metaData")$modelType <- 'binary'
  
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
  
  ## Format data into sparse matrix and outcome vector
  trainData$covariateData$labels <- trainData$labels %>% 
    dplyr::mutate(
      y = sapply(.data$outcomeCount, function(x) min(1,x)),
      time = .data$survivalTime
    )
  covariates <- filterCovariateIds(param, trainData$covariateData) %>% arrange(rowId) %>% collect() #time?
  
  #remap covariates
  covIds <- unique(covariates$covariateId)
  rowIds <- unique(covariates$rowId)
  map <- tibble(covariateId = covIds,
                mappedCovs = c(1:length(covIds)))
  mapRows <- tibble(rowId = rowIds,
                    mappedRows = c(1:length(rowIds)))
  covariates <- left_join(covariates, map, by = "covariateId", keep = FALSE) %>%
    left_join(mapRows, by = "rowId", keep = FALSE)
  
  #create sparse matrix
  i <- covariates$mappedRows
  j <- covariates$mappedCovs
  maxi <- length(rowIds)
  maxj <- length(covIds)
  data <- sparseMatrix(i, j, x = covariates$covariateValue, dims = c(maxi, maxj)) #what does THIS function do? #time?
  y <- trainData$covariateData$labels %>% select(outcomeCount) %>% collect()
  
 # ParallelLogger::logInfo('Running BayesBridge')
  #run BayesBridge
  model <- create_model(y$outcomeCount, data)
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
  modelTrained$coefNames <- c("(Intercept)", unique(covariates$covariateId))
  modelTrained$coefRaw <- mcmc_samples$coef
  modelTrained$coefficients <- apply(mcmc_samples$coef, 1, median)
  names(modelTrained$coefficients) <- modelTrained$coefNames
  
  #modelTrained$priorVariance <- NULL
  
  modelTrained$log_likelihood <- tail(mcmc_samples$logp, 1)
  modelTrained$modelType <- "BayesBridge logistic"
  modelTrained$modelStatus <- "OK"
  attr(modelTrained, "class") <- "plpModel"
  
  #output prediction on train set
  ParallelLogger::logTrace('Getting predictions on train set')
  tempModel <- list(model = modelTrained)
  attr(tempModel, "modelType") <- attr(param, 'modelType')
  prediction <- predictBayesBridge(
    plpModel = tempModel,
    cohort = trainData$labels, 
    data = trainData
  )
  prediction$evaluationType <- 'Train'
  
  #variable importance, time, cv
  comp <- Sys.time() - start
  
  ParallelLogger::logTrace('Getting variable importance')
  variableImportance <- getVariableImportance(modelTrained, trainData)
  
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
    covariateImportance = variableImportance
  )
  
  
  class(result) <- 'plpModel'
  attr(result, 'predictionFunction') <- 'predictBayesBridge'
  attr(result, 'modelType') <- attr(param, 'modelType')
  attr(result, 'saveType') <- attr(param, 'saveType')
  return(result)
}