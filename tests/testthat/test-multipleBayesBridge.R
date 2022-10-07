library(PatientLevelPrediction)
library(bayesbridger)
library(Matrix)
library(reticulate)
library(FeatureExtraction)

#Run 1: sample size 100,000, without Jenna sparse matrix + MKL + renv... results very slow
#Run 2: sample size 100,000, without Jenna sparse matrix
#Run 3: sample size 50,000 (memory issues with larger samples)
#Run 4: new cohorts with above changes

## larger database test
Sys.setenv(DATABASECONNECTOR_JAR_FOLDER = "D:/Drivers")
# DatabaseConnector::downloadJdbcDrivers(dbms = "redshift")
connectionDetailsMDCR <- DatabaseConnector::createConnectionDetails(
  dbms = "redshift",
  server = "ohda-prod-1.cldcoxyrkflo.us-east-1.redshift.amazonaws.com/truven_mdcr",
  port = 5439,
  user = Sys.getenv("username"),
  password = Sys.getenv("password"),
  extraSettings = "ssl=true&sslfactory=com.amazon.redshift.ssl.NonValidatingFactory")
connection <- connect(connectionDetailsMDCR)

# Add the database containing the OMOP CDM data
cdmDatabaseSchema <- 'cdm_truven_mdcr_v1297.dbo'
# Add a shareable name for the database containing the OMOP CDM data
cdmDatabaseName <- 'bayesBridgeTest'
# Add a database with read/write access as this is where the cohorts will be generated
cohortDatabaseSchema <- 'scratch_kli69.dbo'
cohortTable <- "multipleBayesCohorts"

oracleTempSchema <- NULL

databaseDetails <- createDatabaseDetails(
  connectionDetails = connectionDetailsMDCR,
  cdmDatabaseSchema = cdmDatabaseSchema,
  cdmDatabaseName = cdmDatabaseName,
  cohortDatabaseSchema = cohortDatabaseSchema,
  cohortTable = cohortTable,
  outcomeDatabaseSchema = cohortDatabaseSchema,
  outcomeTable = cohortTable,
  cdmVersion = 5)


restrictPlpDataSettings <- createRestrictPlpDataSettings(sampleSize = 50000)
populationSettings <- createStudyPopulationSettings(binary = T,
                                                    includeAllOutcomes = T,
                                                    firstExposureOnly = FALSE,
                                                    washoutPeriod = 364,
                                                    removeSubjectsWithPriorOutcome = TRUE,
                                                    priorOutcomeLookback = 9999,
                                                    riskWindowStart = 1,
                                                    minTimeAtRisk = 364,
                                                    riskWindowEnd = 9999)
cs <- createCovariateSettings(useDemographicsGender = TRUE,
                              useDemographicsAge = TRUE,
                              useConditionGroupEraLongTerm = TRUE,
                              useConditionGroupEraAnyTimePrior = TRUE,
                              useDrugGroupEraLongTerm = TRUE,
                              useDrugGroupEraAnyTimePrior = TRUE,
                              useVisitConceptCountLongTerm = TRUE,
                              longTermStartDays = -365,
                              endDays = -1)

bayesBridge <- setBayesBridge(n_iter = 10000,
                              coef_sampler_type = "cholesky",
                              thin = 10)

configure_python("bayesbridge")

model1 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5412,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model2 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5413,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model3 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5414,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model4 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5415,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model5 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5418,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model6 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5419,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model7 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5420,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model8 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5421,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model9 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5422,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model10 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5423,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model11 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5424,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model12 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5425,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model13 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5426,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model14 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5427,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model15 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5428,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model16 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5429,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model7 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5431,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model8 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5432,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model19 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5433,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model20 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5434,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

model21 <- createModelDesign(
  targetId = 5430,
  outcomeId = 5435,
  restrictPlpDataSettings = restrictPlpDataSettings,
  populationSettings = populationSettings,
  covariateSettings = cs,
  modelSettings = bayesBridge
)

models <- list(model1, model2, model3, model4, model5,
               model6, model7, model8, model9, model10,
               model11, model12, model13, model14, model15,
               model16, model17, model18, model19, model20,
               model21)
results <- runMultiplePlp(
  databaseDetails = databaseDetails,
  modelDesignList = models,
  onlyFetchData = F
)
