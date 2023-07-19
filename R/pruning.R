library(Eunomia)
library(bayesbridger)
library(Matrix)
library(FeatureExtraction)
library(PatientLevelPrediction)
library(reticulate)
use_condaenv("bayesbridgemac")
connectionDetails <- getEunomiaConnectionDetails()
Eunomia::createCohorts(connectionDetails)

getSharedCovIds <- function(databaseDetails1, databaseDetails2,
                            covariateSettings1, covariateSettings2,
                            restrictPlpDataSettings){
  databaseDetails <- list(databaseDetails1, databaseDetails2)
  covariateSettings <- list(covariateSettings1, covariateSettings2)
  
  covIds <- list()
  
  for(i in 1:2){
    connection <- DatabaseConnector::connect(databaseDetails[[i]]$connectionDetails)
    dbms <- databaseDetails[[i]]$connectionDetails$dbms
    pathToSql <- system.file(
      paste("sql/", "sql_server", 
            sep = ""),
      'CreateCohorts.sql', 
      package = "PatientLevelPrediction"
    )
    renderedSql <- readChar(pathToSql, file.info(pathToSql)$size) 
    renderedSql <- SqlRender::render(
      sql = renderedSql,
      cdm_database_schema = databaseDetails[[i]]$cdmDatabaseSchema,
      cohort_database_schema = databaseDetails[[i]]$cohortDatabaseSchema,
      cohort_table = databaseDetails[[i]]$cohortTable,
      cdm_version = databaseDetails[[i]]$cdmVersion,
      target_id = databaseDetails[[i]]$targetId,
      study_start_date = restrictPlpDataSettings$studyStartDate,
      study_end_date = restrictPlpDataSettings$studyEndDate,
      first_only = restrictPlpDataSettings$firstExposureOnly,
      washout_period = restrictPlpDataSettings$washoutPeriod,
      use_sample = !is.null(restrictPlpDataSettings$sampleSize),
      sample_number = restrictPlpDataSettings$sampleSize
      )
    renderedSql <- SqlRender::translate(
      sql = renderedSql, 
      targetDialect = dbms, 
      tempEmulationSchema = databaseDetails[[i]]$tempEmulationSchema
      )
    DatabaseConnector::executeSql(connection, renderedSql)
    
    ParallelLogger::logTrace("Fetching cohorts from server")
    start <- Sys.time()
    
    pathToSql <- system.file(
      paste("sql/", "sql_server", 
            sep = ""),
      "GetCohorts.sql", 
      package = "PatientLevelPrediction"
    )
    
    cohortSql <- readChar(pathToSql, file.info(pathToSql)$size) 
    
    cohortSql <- SqlRender::render(
      sql = cohortSql,
      cdm_version = databaseDetails[[i]]$cdmVersion
    )
    
    cohortSql <- SqlRender::translate(
      sql = cohortSql, 
      targetDialect = dbms, 
      tempEmulationSchema = databaseDetails[[i]]$tempEmulationSchema
    )
    cohorts <- DatabaseConnector::querySql(connection, cohortSql)
    colnames(cohorts) <- SqlRender::snakeCaseToCamelCase(colnames(cohorts))
    metaData.cohort <- list(targetId = databaseDetails[[i]]$targetId)
    
    if(nrow(cohorts)==0){
      stop('Target population is empty')
    }
    
    delta <- Sys.time() - start
    ParallelLogger::logTrace(paste("Loading cohorts took", signif(delta, 3), attr(delta, "units")))
    
    covariateData <- FeatureExtraction::getDbCovariateData(
      connection = connection, 
      oracleTempSchema = databaseDetails[[i]]$tempEmulationSchema,
      cdmDatabaseSchema = databaseDetails[[i]]$cdmDatabaseSchema,
      cdmVersion = databaseDetails[[i]]$cdmVersion,
      cohortTable = "#cohort_person",
      cohortTableIsTemp = TRUE,
      rowIdField = "row_id",
      covariateSettings = covariateSettings[[i]]
    )
    
    covIds[[i]] <- covariateData$covariates %>% pull(covariateId) %>% unique()
    }
  out <- intersect(covIds[[1]], covIds[[2]])
  return(out)
}

cs1 <- createCovariateSettings(useDemographicsGender = TRUE,
                               useDemographicsAgeGroup = TRUE,
                               useConditionGroupEraLongTerm = TRUE,
                               useConditionGroupEraAnyTimePrior = TRUE,
                               useDrugGroupEraLongTerm = TRUE,
                               useDrugGroupEraAnyTimePrior = TRUE,
                               useVisitConceptCountLongTerm = TRUE,
                               longTermStartDays = -365,
                               endDays = -1)

covIds <- bayesResults$model$model$coefficients$covariateIds[-1] %>% as.numeric()
cs2 <- createCovariateSettings(
  useDemographicsGender = TRUE,
  useDemographicsAgeGroup = TRUE,
  useConditionGroupEraLongTerm = TRUE,
  useConditionGroupEraAnyTimePrior = TRUE,
  useDrugGroupEraLongTerm = TRUE,
  useDrugGroupEraAnyTimePrior = TRUE,
  useVisitConceptCountLongTerm = TRUE,
  longTermStartDays = -365,
  endDays = -1,
  includedCovariateIds = covIds[1:10])

databaseDetails <- createDatabaseDetails(connectionDetails = connectionDetails,
                                         cdmDatabaseSchema = "main",
                                         cdmDatabaseName = "bayesbridgeTest",
                                         cohortDatabaseSchema = "main",
                                         cohortTable = "cohort",
                                         targetId = 1,
                                         outcomeDatabaseSchema = "main",
                                         outcomeTable = "cohort",
                                         outcomeIds = 3)
test <- getSharedCovIds(databaseDetails, databaseDetails, cs1, cs2,
                        createRestrictPlpDataSettings())


