library(PatientLevelPrediction)
library(bayesbridger)
library(Matrix)

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
cohortTable <- "bayesBridgeCohorts"

oracleTempSchema <- NULL


#sql
sql <- SqlRender::loadRenderTranslateSql(sqlFilename = "CreateCohortTable.sql",
                                         packageName = "PatientLevelPrediction",
                                         dbms = attr(connection, "dbms"),
                                         oracleTempSchema = oracleTempSchema,
                                         cohort_database_schema = cohortDatabaseSchema,
                                         cohort_table = cohortTable)
DatabaseConnector::executeSql(connection, sql, progressBar = FALSE, reportOverallTime = FALSE)

sql <- SqlRender::loadRenderTranslateSql(sqlFilename = "OHDSI_EU_2019_Angioedema_events.sql",
                                         packageName = "PatientLevelPrediction",
                                         dbms = attr(connection, "dbms"),
                                         oracleTempSchema = oracleTempSchema,
                                         cdm_database_schema = cdmDatabaseSchema,
                                         vocabulary_database_schema = cdmDatabaseSchema,
                                         target_database_schema = cohortDatabaseSchema,
                                         target_cohort_table = cohortTable,
                                         target_cohort_id = 3046)
DatabaseConnector::executeSql(connection, sql)

sql <- SqlRender::loadRenderTranslateSql(sqlFilename = "PvInAction_New_users_of_ACE_inhibitors_as_firstline_monotherapy_for_hypertension.sql",
                                         packageName = "PatientLevelPrediction",
                                         dbms = attr(connection, "dbms"),
                                         oracleTempSchema = oracleTempSchema,
                                         cdm_database_schema = cdmDatabaseSchema,
                                         vocabulary_database_schema = cdmDatabaseSchema,
                                         target_database_schema = cohortDatabaseSchema,
                                         target_cohort_table = cohortTable,
                                         target_cohort_id = 3043)
DatabaseConnector::executeSql(connection, sql)

sql <- SqlRender::loadRenderTranslateSql(sqlFilename = "BookOfOHDSI_Angioedema_events.sql",
                                         packageName = "PatientLevelPrediction",
                                         dbms = attr(connection, "dbms"),
                                         oracleTempSchema = oracleTempSchema,
                                         cdm_database_schema = cdmDatabaseSchema,
                                         vocabulary_database_schema = cdmDatabaseSchema,
                                         target_database_schema = cohortDatabaseSchema,
                                         target_cohort_table = cohortTable,
                                         target_cohort_id = 1770673)
DatabaseConnector::executeSql(connection, sql)

sql <- SqlRender::loadRenderTranslateSql(sqlFilename = "UNLEY_New_ACE_Users.sql",
                                         packageName = "PatientLevelPrediction",
                                         dbms = attr(connection, "dbms"),
                                         oracleTempSchema = oracleTempSchema,
                                         cdm_database_schema = cdmDatabaseSchema,
                                         vocabulary_database_schema = cdmDatabaseSchema,
                                         target_database_schema = cohortDatabaseSchema,
                                         target_cohort_table = cohortTable,
                                         target_cohort_id = 1770675)
DatabaseConnector::executeSql(connection, sql)

sql <- SqlRender::loadRenderTranslateSql(sqlFilename = "PatientLevelPrediction_vignette__O_Ischemic_stroke_events.sql",
                                         packageName = "PatientLevelPrediction",
                                         dbms = attr(connection, "dbms"),
                                         oracleTempSchema = oracleTempSchema,
                                         cdm_database_schema = cdmDatabaseSchema,
                                         vocabulary_database_schema = cdmDatabaseSchema,
                                         target_database_schema = cohortDatabaseSchema,
                                         target_cohort_table = cohortTable,
                                         target_cohort_id = 1769448)
DatabaseConnector::executeSql(connection, sql)

sql <- SqlRender::loadRenderTranslateSql(sqlFilename = "vignette_AF.sql",
                                         packageName = "PatientLevelPrediction",
                                         dbms = attr(connection, "dbms"),
                                         oracleTempSchema = oracleTempSchema,
                                         cdm_database_schema = cdmDatabaseSchema,
                                         vocabulary_database_schema = cdmDatabaseSchema,
                                         target_database_schema = cohortDatabaseSchema,
                                         target_cohort_table = cohortTable,
                                         target_cohort_id = 1769447)
DatabaseConnector::executeSql(connection, sql)


cs <- createCovariateSettings(useDemographicsGender = TRUE,
                              useDemographicsAge = TRUE,
                              useConditionGroupEraLongTerm = TRUE,
                              useConditionGroupEraAnyTimePrior = TRUE,
                              useDrugGroupEraLongTerm = TRUE,
                              useDrugGroupEraAnyTimePrior = TRUE,
                              useVisitConceptCountLongTerm = TRUE,
                              longTermStartDays = -365,
                              endDays = -1)

#Old cohorts
# databaseDetails <- createDatabaseDetails(connectionDetails = connectionDetailsMDCR,
#                                          cdmDatabaseSchema = cdmDatabaseSchema,
#                                          cdmDatabaseName = cdmDatabaseName,
#                                          cohortDatabaseSchema = cohortDatabaseSchema,
#                                          cohortTable = cohortTable,
#                                          cohortId = 3043,
#                                          outcomeDatabaseSchema = cohortDatabaseSchema,
#                                          outcomeTable = cohortTable,
#                                          outcomeIds = 3046)

databaseDetails <- createDatabaseDetails(connectionDetails = connectionDetailsMDCR,
                                         cdmDatabaseSchema = cdmDatabaseSchema,
                                         cdmDatabaseName = cdmDatabaseName,
                                         cohortDatabaseSchema = cohortDatabaseSchema,
                                         cohortTable = cohortTable,
                                         cohortId = 1769447,
                                         outcomeDatabaseSchema = cohortDatabaseSchema,
                                         outcomeTable = cohortTable,
                                         outcomeIds = 1769448)
restrictPlpDataSettings <- createRestrictPlpDataSettings(sampleSize = 50000)
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
                                  outcomeId = 1769448,#3046,
                                  populationSettings = populationSettings)

splitSettings <- createDefaultSplitSetting(testFraction = 0.3,
                                           trainFraction = 0.7,
                                           splitSeed = 1234,
                                           type = "stratified")

lr <- setLassoLogisticRegression()
lrResults <- runPlp(plpData = plpData,
                    outcomeId = 1769448,#3046,
                    populationSettings = populationSettings,
                    modelSettings = lr,
                    splitSettings = splitSettings)

configure_python(envname = "bayesbridge")
bayesBridge <- setBayesBridge(n_iter = 10000,
                              coef_sampler_type = "cg",
                              thin = 10)
bayesResults <- runPlp(plpData = plpData,
                       outcomeId = 1769448,#3046,
                       populationSettings = populationSettings,
                       modelSettings = bayesBridge,
                       splitSettings = splitSettings)

plotPlp(bayesResults, saveLocation = "./largeBayesPlots3")
plotPlp(lrResults, saveLocation = "./largeLrPlots4")

## Plot distribution of regression coefficients

bayesOut <- readRDS("./largeTestBayesFinal/plpResult/runPlp.rds")
lrOut <- readRDS("./largeTestLrFinal/plpResult/runPlp.rds")
coefRaw <- readRds("./9-22 coefs/coef.rds")

coefBayes <- bayesOut$covariateSummary %>% 
  mutate(covariateValue = unname(unlist(covariateValue))) %>%
  na.omit() %>% arrange(-abs(covariateValue))
coefLr <- lrOut$covariateSummary %>%
  na.omit() %>% arrange(-abs(covariateValue)) 

coefBayesPlot <- coefBayes[-1,] 

coefBayes %>% select(covariateId, covariateName, covariateValue) %>% head(10)
coefLr %>% select(covariateId, covariateName, covariateValue) %>% head(10)

ggplot(coefBayesPlot, aes(x = covariateValue)) + 
  geom_histogram(bins=1000) + 
  scale_x_continuous(n.breaks = 20)

ggplot(coefLr, aes(x = covariateValue)) + 
  geom_histogram(bins=500) + 
  scale_x_continuous(n.breaks = 20)
