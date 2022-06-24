library(Eunomia)
library(bayesbridger)
library(Matrix)
#Learning PLP logistic regression with Eunomia
connectionDetails <- getEunomiaConnectionDetails()
Eunomia::createCohorts(connectionDetails)
cs <- createCovariateSettings(useDemographicsGender = TRUE,
                              useDemographicsAge = TRUE,
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
                                         cohortId = 1,
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

# Artifacts for model fitting
# train <- readRDS("./artifacts/train.rds")
# test <- readRDS("./artifacts/test.rds")
# covariates <- readRDS("./artifacts/covariates.rds")

## configure python for bayesbridge
configure_python(envname = "bayesbridge")
bayesBridge <- setBayesBridge(n_iter = 500)
bayesResults <- runPlp(plpData = plpData,
                       outcomeId = 3,
                       populationSettings = populationSettings,
                       modelSettings = bayesBridge,
                       splitSettings = splitSettings)
