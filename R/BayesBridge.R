#' Settings for BayesBridge Logistic Regression
#' @export
#' 
setBayesBridge <- function(seed = NULL){
  if(is.null(seed[1])){
    seed <- as.integer(sample(100000000,1))
  }
  
  param <- list()
  
  attr(param, 'settings') <- list(
    priorfunction = 'Cyclops::createPrior',
    modelType = 'BayesBridge logistic',
    seed = seed[1],
    name = "BayesBridge Logistic Regression"
  )
  
  attr(param, 'modelType') <- 'binary' 
  attr(param, 'saveType') <- 'RtoJson'
  
  result <- list(fitFunction = "fitBayesBridge",
                 param = param)
  class(result) <- "modelSettings"
  
  return(result)
}

#' Fit BayesBridge Logistic Regression
#' @export
#' 
fitBayesBridge <- function(trainData, param, analysisId, ...){
  settings <- attr(param, 'settings')
  
  ## Format data into sparse matrix and outcome vector
  trainData$covariateData$labels <- trainData$labels %>% 
    dplyr::mutate(
      y = sapply(.data$outcomeCount, function(x) min(1,x)),
      time = .data$survivalTime
    )
  covariates <- filterCovariateIds(param, trainData$covariateData) %>% arrange(rowId) %>% collect()
  
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
  data <- sparseMatrix(i, j, x = covariates$covariateValue, dims = c(maxi, maxj))
  y <- trainData$covariateData$labels %>% select(outcomeCount) %>% collect()
  
  #run BayesBridge
  model <- create_model(y$outcomeCount, data)
  prior <- create_prior(bridge_exponent=1, #param
                        regularizing_slab_size = 1.) #param
  bridge <- instantiate_bayesbridge(model, prior)
  
  n_iter <- 250L #param
  gibbs_output <- gibbs(bridge, n_iter, thin = 1,
                        seed = 2021L,
                        coef_sampler_type = "cholesky",
                        n_status_update = 10L)
  mcmc_samples <- gibbs_output$samples
  
  #output modelTrained
  modelTrained <- list()
  coefNames <- c("(Intercept)", unique(covariates$covariateId))
  modelTrained$coefficients <- rowMeans(mcmc_samples$coef)
  names(modelTrained$coefficients) <- coefNames
  modelTrained$priorVariance <- NULL
  modelTrained$log_likelihood <- tail(mcmc_samples$logp, 1)
  modelTrained$modelType <- "BayesBridge logistic"
  modelTrained$modelStatus <- NULL
  attr(modelTrained, "class") <- "plpModel"
  
  #output prediction
  prediction <- NULL
  
  #variable importance, time, cv
  comp <- 0
  cvPerFold <- 0
  variableImportance <- NULL
  
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
      trainingDate = Sys.Date(),
      hyperParamSearch = cvPerFold
    ),
    covariateImportance = variableImportance
  )
  
  
  class(result) <- 'plpModel'
  attr(result, 'predictionFunction') <- 'predictCyclops'
  attr(result, 'modelType') <- attr(param, 'modelType')
  attr(result, 'saveType') <- attr(param, 'saveType')
  return(result)
}
