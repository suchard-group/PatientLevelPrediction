library(Eunomia)
library(bayesbridger)
library(Matrix)
library(FeatureExtraction)
library(PatientLevelPrediction)
library(reticulate)
use_condaenv("bayesbridgemac")
#Learning PLP logistic regression with Eunomia
connectionDetails <- getEunomiaConnectionDetails()
Eunomia::createCohorts(connectionDetails)
cs <- createCovariateSettings(useDemographicsGender = TRUE,
                              useDemographicsAgeGroup = TRUE,
                              useConditionGroupEraLongTerm = TRUE,
                              useConditionGroupEraAnyTimePrior = TRUE,
                              useDrugGroupEraLongTerm = TRUE,
                              useDrugGroupEraAnyTimePrior = TRUE,
                              useVisitConceptCountLongTerm = TRUE,
                              longTermStartDays = -365,
                              endDays = -1)
databaseDetails <- createDatabaseDetails(connectionDetails = connectionDetails,
                                         cdmDatabaseSchema = "main",
                                         cdmDatabaseName = "bayesbridgeTest",
                                         cohortDatabaseSchema = "main",
                                         cohortTable = "cohort",
                                         targetId = 1,
                                         outcomeDatabaseSchema = "main",
                                         outcomeTable = "cohort",
                                         outcomeIds = 3)
restrictPlpDataSettings <- createRestrictPlpDataSettings(sampleSize = 10000)
plpData <- getPlpData(covariateSettings = cs,
                      databaseDetails = databaseDetails,
                      restrictPlpDataSettings = restrictPlpDataSettings)

populationSettings <- createStudyPopulationSettings(binary = T,
                                                    includeAllOutcomes = T,
                                                    firstExposureOnly = FALSE,
                                                    washoutPeriod = 364,
                                                    removeSubjectsWithPriorOutcome = TRUE,
                                                    priorOutcomeLookback = 9999,
                                                    riskWindowStart = 1,
                                                    minTimeAtRisk = 364,
                                                    riskWindowEnd = 9999)
studyPop <- createStudyPopulation(plpData = plpData,
                                  outcomeId = 3,
                                  populationSettings = populationSettings)

splitSettings <- createDefaultSplitSetting(testFraction = 0.25,
                                           splitSeed = 1234,
                                           type = "stratified")
lr <- setLassoLogisticRegression()
lrResults <- runPlp(plpData = plpData,
                    outcomeId = 3,
                    populationSettings = populationSettings,
                    modelSettings = lr,
                    splitSettings = splitSettings)
coeflr <- lrResults$model$model$coefficients

# Artifacts for model fitting
# train <- readRDS("./artifacts/train.rds")
# test <- readRDS("./artifacts/test.rds")
# covariates <- readRDS("./artifacts/covariates.rds")

## configure python for bayesbridge
bayesBridge <- setBayesBridge(n_iter = 1000,
                              n_burnin = 100,
                              bridge_exponent = 1,
                              init = list(global_scale = 0.1),
                              options = list(local_scale_update = "shrunk_only",
                                             global_scale_update = "None",
                                             q_update = "simple"),
                              coef_sampler_type = "cg")
bayesResults <- runPlp(plpData = plpData,
                       outcomeId = 3,
                       populationSettings = populationSettings,
                       modelSettings = bayesBridge,
                       splitSettings = splitSettings)
coef1 <- bayesResults$model$model$coefficients

fixed_effects <- tibble(covariateId = c(8532001, 78272209, 8003),
                        mean = c(10, 20, 30),
                        sd = rep(0.0001, 3))

bayesBridgeFixed <- setBayesBridge(n_iter = 2000,
                              n_burnin = 500,
                              bridge_exponent = 0.5,
                              fixed_effects = fixed_effects,
                              coef_sampler_type = "cg")
bayesResultsFixed <- runPlp(plpData = plpData,
                       outcomeId = 3,
                       populationSettings = populationSettings,
                       modelSettings = bayesBridgeFixed,
                       splitSettings = splitSettings)
coef2 <- bayesResultsFixed$model$model$coefficients

mixture <- tibble(covariateId = c(8532001, 78272209, 8003),
                  mean = c(10, 20, 30),
                  sd = rep(0.0001, 3))
bayesBridgeMixture <- setBayesBridge(n_iter = 5000,
                                   n_burnin = 500,
                                   bridge_exponent = 0.5,
                                   mixture = mixture,
                                   coef_sampler_type = "cg",
                                   local_scale_sampler_type = "shrunk_only",
                                   params_to_save = c('coef', 'global_scale', 'logp', 'local_scale'))
bayesResultsMixture <- runPlp(plpData = plpData,
                              outcomeId = 3,
                              populationSettings = populationSettings,
                              modelSettings = bayesBridgeMixture,
                              splitSettings = splitSettings)
coef3 <- bayesResultsMixture$model$model$coefficients

mixture2 <- tibble(covariateId = c(4027663209, 1728416409, 8003),
                  mean = c(1.67, -0.106, 30),
                  sd = c(0.5, 0.03, 1))
bayesBridgeMixture2 <- setBayesBridge(n_iter = 5000,
                                     n_burnin = 500,
                                     bridge_exponent = 0.5,
                                     mixture = mixture2,
                                     coef_sampler_type = "cg",
                                     local_scale_sampler_type = "shrunk_only",
                                     params_to_save = c('coef', 'global_scale', 'logp', 'local_scale'))
bayesResultsMixture2 <- runPlp(plpData = plpData,
                              outcomeId = 3,
                              populationSettings = populationSettings,
                              modelSettings = bayesBridgeMixture2,
                              splitSettings = splitSettings)
coef4 <- bayesResultsMixture2$model$model$coefficients

plotPlp(bayesResults, saveLocation = "./bayesSmallPlots")

model1 <- createModelDesign(targetId = 1,
                            outcomeId = 3,
                            restrictPlpDataSettings = restrictPlpDataSettings,
                            populationSettings = populationSettings,
                            covariateSettings = cs,
                            splitSettings = splitSettings,
                            modelSettings = bayesBridge)
models <- runMultiplePlp(databaseDetails = databaseDetails,
                         modelDesignList = list(model1),
                         saveDirectory = "./smallMultiple1")

viewMultiplePlp(analysesLocation = "./smallMultiple1")

res <- loadPlpResult("./2022-10-06-/plpResult")
viewPlp(res)

lrRes <- loadPlpResult("./2022-10-17-/plpResult")
viewPlp(lrRes)
