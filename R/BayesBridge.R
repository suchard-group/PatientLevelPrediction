#' Settings for BayesBridge Logistic Regression
#' @param seed       An option to add a seed when training the model
#' @param n_iter Numeric: total number of MCMC iterations i.e. burn-ins + saved posterior draws
#' @param thin Numeric: Number of iterations per saved samples for “thinning” MCMC to reduce the output size.
#' @param bridge_exponent Numeric: exponent for bridge prior on regression coefficients (<2). The value of 2 corresponds to the Gaussian prior and 1 corresponds to the double-exponential (Bayesian Lasso) prior.
#' @param regularizing_slab_size Numeric: Standard deviation of the Gaussian tail-regularizer on the bridge prior.
#' @param n_status_update Numeric: Number of updates to print during the sampler running
#' @param global_scale Numeric: Reference prior for a scale parameter
#' @param coef_sampler_type An object to specify the sampling method to update regression coefficients:
#' @param params_to_fix A vector to specify parameters to be fixed during MCMC sampling. Currently only supports "global_scale".
#' \itemize{
#'                                         \item{None}{ Chooses a method via a crude heuristic based on model type}
#'                                         \item{cholesky}{ Cholesky decomposition based sampler}
#'                                         \item{cg}{ Conjugate gradient sampler preferred over Cholesky for linear and logistic models with large+sparse design matrix}
#'                                         \item{hmc}{ Hamilton Monte Carlo for other models}
#'                                         } 
#' @export
#' 
setBayesBridge <- function(seed = NULL,
                           n_iter = 250,
                           n_burnin = 100,
                           bridge_exponent = 1,
                           regularizing_slab_size = 1,
                           thin = 1,
                           n_status_update = 10,
                           init = list(global_scale = 0.1),
                           coef_sampler_type = "cholesky",
                           params_to_fix = c(),
                           local_scale_sampler_type = "all",
                           params_to_save = c('coef', 'global_scale', 'logp'),
                           fixed_effects = NULL,
                           mixture = NULL){
  if(is.null(seed[1])){
    seed <- as.integer(sample(100000000,1))
  }
  
  if(is.data.frame(fixed_effects)){
    n_fixed_effects <- nrow(fixed_effects)
  } else{
    n_fixed_effects <- 0
  }
  
  if(is.data.frame(mixture)){
    n_mixture <- nrow(mixture)
    params_to_save <- c(params_to_save, 'gamma')
  } else{
    n_mixture <- 0
  }
  
  param <- list()
  
  attr(param, 'settings') <- list(
    priorfunction = 'Cyclops::createPrior',
    modelType = 'BayesBridge logistic',
    seed = seed[1],
    name = "BayesBridge Logistic Regression",
    n_burnin = n_burnin,
    n_iter = n_iter,
    bridge_exponent = bridge_exponent,
    regularizing_slab_size = regularizing_slab_size,
    thin = thin,
    n_status_update = n_status_update,
    init = init,
    coef_sampler_type = coef_sampler_type,
    params_to_fix = params_to_fix,
    local_scale_sampler_type = local_scale_sampler_type,
    params_to_save = params_to_save,
    fixed_effects = fixed_effects,
    n_fixed_effects = n_fixed_effects,
    mixture = mixture,
    n_mixture = n_mixture
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
    colnames(newData) <- matrixObjects$covariateMap %>% arrange(columnId) %>% pull(covariateId)
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
  
  #Get posterior median
  valueMed <- getLinkPostMedian(x = newData, betas = plpModel$model$samples$coef)
  
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
fitBayesBridge <- function(trainData, modelSettings, analysisId, ...){
  param <- modelSettings$param
  settings <- attr(param, 'settings')
  
  start <- Sys.time() 
  
  # check plpData is coo format:
  if (!FeatureExtraction::isCovariateData(trainData$covariateData)){
    stop("Needs correct covariateData")
  }
  
  #remap covariates
  mappedData <- toSparseM(trainData)
  matrixData <- mappedData$dataMatrix
  colnames(matrixData) <- mappedData$covariateMap %>% arrange(columnId) %>% pull(covariateId)
  labels <- mappedData$labels
  covariateRef <- mappedData$covariateRef
  
  #check for fixed effects
  sd_for_mixture <- Inf
  mean_for_mixture <- Inf
  sd_for_fixed_effect <- Inf
  mean_for_fixed_effect <- Inf
  
  order <- colnames(matrixData)
  
  if(settings$n_mixture > 0){
    if(length(mixture$covariateId) != length(unique(mixture$covariateId))){
      ParallelLogger::logError("Multiple means and standard deviations specified for one covariate.")
    }
    mixture <- settings$mixture %>% filter(covariateId %in% colnames(matrixData)) %>%
      mutate(covariateId = as.character(covariateId))
    matrixDataMixture <- matrixData[,mixture$covariateId]
    matrixData <- matrixData[,!colnames(matrixData) %in% mixture$covariateId]
    sd_for_mixture <- mixture$sd
    mean_for_mixture <- mixture$mean
  } 
  if(settings$n_fixed_effect > 0){
    if(length(fixed_effects$covariateId) != length(unique(fixed_effects$covariateId))){
      ParallelLogger::logError("Multiple means and standard deviations specified for one covariate.")
    }
    fixed_effects <- settings$fixed_effects %>% filter(covariateId %in% colnames(matrixData)) %>%
      mutate(covariateId = as.character(covariateId))
    matrixDataFixed <- matrixData[,fixed_effects$covariateId]
    matrixData <- matrixData[,!colnames(matrixData) %in% fixed_effects$covariateId]
    sd_for_fixed_effect <- fixed_effects$sd
    mean_for_fixed_effect <- fixed_effects$mean
  }
  
  if(exists("matrixDataFixed") && exists("matrixDataMixture")){
    matrixData <- cbind(matrixDataFixed, matrixDataMixture, matrixData)
    n_mixture <- ncol(matrixDataMixture)
  } else if (exists("matrixDataFixed") && !exists("matrixDataMixture")){
    matrixData <- cbind(matrixDataFixed, matrixData)
    n_mixture <- 0
  } else if (!exists("matrixDataFixed") && exists("matrixDataMixture")){
    matrixData <- cbind(matrixDataMixture, matrixData)
    n_mixture <- ncol(matrixDataMixture)
  } else{
    n_mixture <- 0
  }
  
  ParallelLogger::logInfo('Running BayesBridge')
  #run BayesBridge
  start1 <- Sys.time()
  model <- create_model(labels$outcomeCount, matrixData)
  prior <- create_prior(bridge_exponent = settings$bridge_exponent, #param
                        regularizing_slab_size = settings$regularizing_slab_size,
                        n_fixed_effect = as.integer(settings$n_fixed_effect),
                        sd_for_fixed_effect = sd_for_fixed_effect,
                        mean_for_fixed_effect = mean_for_fixed_effect,
                        n_mixture = n_mixture,
                        sd_for_mixture = sd_for_mixture,
                        mean_for_mixture = mean_for_mixture) #param
  bridge <- instantiate_bayesbridge(model, prior)
  
  n_iter <- settings$n_iter #param
  
  # Options for global scale and local scale sampling:
  options <- list()
  if ("global_scale" %in% settings$params_to_fix){
    options <- list("global_scale_update" = NULL)
    } 
  if (settings$local_scale_sampler_type != "all"){
    options <- c(options, list("local_scale_update" = settings$local_scale_sampler_type))
  } else{
    options <- NULL
  }
  
  gibbs_output <- gibbs(bridge, 
                        n_burnin = as.integer(settings$n_burnin),
                        n_iter = as.integer(settings$n_iter), 
                        init = settings$init,
                        thin = settings$thin,
                        seed = settings$seed,
                        coef_sampler_type = settings$coef_sampler_type,
                        n_status_update = settings$n_status_update,
                        params_to_save = settings$params_to_save,
                        options = options)
  
  comp1 <- Sys.time() - start1
  ParallelLogger::logInfo("MCMC took ", signif(comp1, 3), " ", attr(comp1, "units"))
  
  #output modelTrained
  modelTrained <- list()
  modelTrained$samples <- gibbs_output$samples
  rownames(modelTrained$samples$coef) <- c("(Intercept)", colnames(matrixData))
  modelTrained$coefficients <- tibble(betas = apply(modelTrained$samples$coef, 1, median),
                                      covariateIds = rownames(modelTrained$samples$coef))
  modelTrained$modelType <- "BayesBridge logistic"
  modelTrained$modelStatus <- "OK"
  attr(modelTrained, "class") <- "plpModel"
  
  #output prediction on train set
  ParallelLogger::logTrace('Getting predictions on train set')
  start2 <- Sys.time()
  valueMed <- getLinkPostMedian(x = matrixData, betas = modelTrained$samples$coef)
  labels$value <- valueMed
  prediction <- labels
  comp2 <- Sys.time() - start2
  ParallelLogger::logInfo("Prediction for training set took ", signif(comp2, 3), " ", attr(comp2, "units"))
  
  attr(prediction, "metaData")$modelType <- 'binary'
  prediction$evaluationType <- 'Train'
  
  #variable importance
  ParallelLogger::logTrace('Getting variable importance')
  varImp <- data.frame(
    covariateId = as.double(modelTrained$coefficients$covariateIds[modelTrained$coefficients$covariateIds!='(Intercept)']),
    value = modelTrained$coefficients$betas[modelTrained$coefficients$covariateIds!='(Intercept)']
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
    select(covariateValue) %>%
    t() %>%
    as.vector()
  covariateRef <- covariateRef %>% arrange(-abs(.data$covariateValue))
  
  #re-order results for Test set validation
  if("(Intercept)" %in% rownames(modelTrained$samples$coef)){
    order <- c(("Intercept"), order)
  }
  map <- match(rownames(modelTrained$samples$coef), order)
  
  if(!is.null(modelTrained$samples$gamma)){
    gammaMap <- map[-1] - 1
    names(modelTrained$samples$gamma) <- order[-1]
    modelTrained$samples$gamma <- modelTrained$samples$gamma[order(gammaMap),]
  }
  if(!is.null(modelTrained$samples$local_scale)){
    lscaleMap <- map[-1] - 1
    names(modelTrained$samples$local_scale) <- order[-1]
    modelTrained$samples$local_scale <- modelTrained$samples$local_scale[order(lscaleMap),]
  }
  
  
  modelTrained$samples$coef <- modelTrained$samples$coef[order(map),]
  
  #output result
  comp <- Sys.time() - start
  result <- list(
    model = modelTrained,
    prediction = prediction,
    
    preprocessing = list(
      featureEngineering = attr(trainData, "metaData")$featureEngineering,#learned mapping
      tidyCovariates = attr(trainData$covariateData, "metaData")$tidyCovariateDataSettings,  #learned mapping
      requireDenseMatrix = F
    ),
    
    modelDesign = PatientLevelPrediction::createModelDesign(
      targetId = attr(trainData, "metaData")$targetId, 
      outcomeId = attr(trainData, "metaData")$outcomeId, 
      restrictPlpDataSettings = attr(trainData, "metaData")$restrictPlpDataSettings, 
      covariateSettings = attr(trainData, "metaData")$covariateSettings,
      populationSettings = attr(trainData, "metaData")$populationSettings, 
      featureEngineeringSettings = attr(trainData, "metaData")$featureEngineeringSettings,
      preprocessSettings = attr(trainData$covariateData, "metaData")$preprocessSettings,
      modelSettings = modelSettings, 
      splitSettings = attr(trainData, "metaData")$splitSettings,
      sampleSettings = attr(trainData, "metaData")$sampleSettings
    ),
    
    trainDetails = list(
      analysisId = analysisId,
      developmentDatabase = attr(trainData, "metaData")$cdmDatabaseSchema,
      developmentDatabaseId = attr(trainData, "metaData")$cdmDatabaseId,
      attrition = attr(trainData, "metaData")$attrition, 
      trainingTime =  paste(as.character(abs(comp)), attr(comp,'units')),
      trainingDate = Sys.Date(),
      modelName = settings$modelType
    ),
    covariateImportance = covariateRef
  )
  
  
  class(result) <- 'plpModel'
  attr(result, 'predictionFunction') <- 'predictBayesBridge'
  attr(result, 'modelType') <- attr(param, 'modelType')
  attr(result, 'saveType') <- attr(param, 'saveType')
  return(result)
}

getLinkPostMedian <- function(x, betas){
  out <- matrix(NA, nrow = nrow(x), ncol = ncol(betas))
  #get posterior median predictive values
  for(i in 1:ncol(betas)){
    coefficients <- betas[,i]
    intercept <- coefficients[names(coefficients)%in%'(Intercept)']
    if(length(intercept)==0) intercept <- 0
    coefficients <- coefficients[!names(coefficients)%in%'(Intercept)']
    beta <- as.numeric(coefficients)
    value <- x %*% beta
    value[is.na(value)] <- 0
    value <- value + intercept
    link <- function(x){
      return(1/(1 + exp(0 - x)))
    }
    value <- link(value) %>% as.vector()
    out[,i] <- value
  }
  valueMed <- apply(out, 1, median)
  return(valueMed)
}