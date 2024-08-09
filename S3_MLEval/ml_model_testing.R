# Load required libraries
suppressPackageStartupMessages({
  library(MLSeq)
  library(DESeq2)
  library(S4Vectors)
  library(kernlab)
  library(parallel)
  library(doParallel)
  library(foreach)
})

# Print script start message
cat("R script started\n")

# Print R version
cat("R version:", R.version.string, "\n")

# 1 - Function to load data
load_data <- function(filtered_data_path, metadata_path, patients_path) {
  cat("Loading data...\n")
  
  # Load numeric matrix
  PD <- read.table(filtered_data_path, header=TRUE, sep="\t", check.names = FALSE)
  cat("Filtered data dimensions:", nrow(PD), "x", ncol(PD), "\n")
  
  # Set the 'Name' column as row names
  rownames(PD) <- PD$Gene
  # Remove the 'Name' column from the data frame
  PD$Gene <- NULL
  
  # Load metadata
  metadata <- read.csv(metadata_path, header=TRUE, check.names = FALSE)
  cat("Metadata dimensions:", nrow(metadata), "x", ncol(metadata), "\n")
  
  # Load patients
  patients <- read.csv(patients_path, header=TRUE, check.names = FALSE)
  cat("Patients data dimensions:", nrow(patients), "x", ncol(patients), "\n")
  
  return(list(PD = PD, metadata = metadata, patients = patients))
}


# 2 - Filter metadata by patient IDs

filter_and_add_metadata <- function(data_list) {
  cat("Filtering metadata...\n")
    
  # 3.1 - Filter metadata by patient information
  patients_in_matrix_log <- data_list$patients$HudAlphaID %in% colnames(data_list$PD)
  patients_in_matrix <- data_list$patients[patients_in_matrix_log, ]

  cat("The length of patients in matrix:  ", length(patients_in_matrix)) # Here is the problem,
  
  # 3.2 Filter metadata by row information where 'Case Control' is not 'Other', '', or 'POOL'
  metadata_filt <- subset(data_list$metadata, 
                          `Case Control` != "Other" & 
                          `Case Control` != "" & 
                          `Case Control` != "POOL")
    
  cat("Filtering metadata by patient ID...\n") # This is the mistake
    
  metadata_filt <- subset(metadata_filt, HudAlphaSampleName %in% patients_in_matrix)
  metadata_filt <- metadata_filt[, c("HudAlphaSampleName", "Case Control")]
  
  # Create class object
  class <- DataFrame(condition = metadata_filt[, "Case Control"])
  
  # 3.3 - Adding outputs to data_list
  data_list[['metadata_filt']] <- metadata_filt
  data_list[['class']] <- class
  
  return(data_list)
}


# 3 - Function to filter and align data

filter_and_align_data <- function(data_list) {
  cat("Aligning data...\n")
  
  # 4 - Align Numeric Matrix with Metadata
  # 4.1 - Order of HudAlphaSampleName in filtered_metadata
  order_hudalpha <- data_list$metadata_filt$HudAlphaSampleName
  
  # 4.2 Find common column names
  common_columns <- intersect(order_hudalpha, colnames(data_list$PD))
  
  # 4.3 - Check if there are any missing columns
  missing_columns <- setdiff(order_hudalpha, colnames(data_list$PD))
  if (length(missing_columns) > 0) {
    cat("Warning: The following columns are missing in the numeric matrix:\n")
    print(missing_columns)
  }
  
  # Reorder the columns of the numeric matrix to match the order of common HudAlphaSampleName values
  data_list$filtered_PD <- data_list$PD[, common_columns, drop = FALSE]
  
  cat("Filtered PD dimensions:", nrow(data_list$filtered_PD), "x", ncol(data_list$filtered_PD), "\n")
  
  return(data_list)
}

# 4 - Function to prepare data for machine learning

prepare_ml_data <- function(data_list) {
  cat("Preparing data for machine learning...\n")
  set.seed(2128)
  
  # Use data_list$filtered_PD instead of filtered_PD
  PD_rounded <- round(data_list$filtered_PD)
  PD_rounded[] <- lapply(PD_rounded, as.integer)
  
  nTest <- ceiling(ncol(PD_rounded) * 0.3)
  ind <- sample(ncol(PD_rounded), nTest, FALSE)
  
  data.train <- as.matrix(PD_rounded[ ,-ind] + 1)
  data.test <- as.matrix(PD_rounded[ ,ind] + 1)
  
  # Assuming class is also part of data_list
  classtr <- DataFrame(condition = data_list$class[-ind, ])
  classts <- DataFrame(condition = data_list$class[ind, ])
  
  cat(paste("Training data dimensions:", nrow(data.train), "x", ncol(data.train), "\n"))
  cat(paste("Testing data dimensions:", nrow(data.test), "x", ncol(data.test), "\n"))
  
  return(list(data.train = data.train, data.test = data.test, classtr = classtr, classts = classts))
}


# 5 - Function to create DESeq objects

create_deseq_objects <- function(ml_data) {
  cat("Creating DESeq objects...\n")
  
  data.trainS4 <- DESeqDataSetFromMatrix(countData = ml_data$data.train, 
                                         colData = ml_data$classtr,
                                         design = formula(~condition))
  
  data.testS4 <- DESeqDataSetFromMatrix(countData = ml_data$data.test, 
                                        colData = ml_data$classts,
                                        design = formula(~condition))
  
  cat("DESeq objects created successfully\n")
  return(list(data.trainS4 = data.trainS4, data.testS4 = data.testS4))
}


# 6 - Function to run MLSeq models in parallel

run_mlseq_models_parallel <- function(deseq_objects, ref_class = "Case", n_cores = NULL, max_models = NULL) {
  start_time <- Sys.time()
  
  cat("Starting run_mlseq_models_parallel function\n")
  
  # Extract training data from deseq_objects
  train_data <- deseq_objects$data.trainS4
  
  # Handle core specification
  if (is.null(n_cores)) {
    n_cores <- detectCores()
    if (is.na(n_cores) || n_cores < 1) {
      cat("Warning: Could not detect number of cores. Defaulting to 1 core.\n")
      n_cores <- 1
    } else {
      n_cores <- n_cores - 1  # Use one less than detected
    }
  }
  n_cores <- max(1, as.integer(n_cores))  # Ensure at least 1 core and integer value
  cat(paste("Number of cores to be used:", n_cores, "\n"))
  
  cat("Setting up parallel processing...\n")
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  cat("Parallel processing setup complete\n")

  cat("Defining control parameters...\n")
  ctrl <- trainControl(method = "repeatedcv", number = 2, repeats = 2, classProbs = TRUE)
  ctrl_discrete <- discreteControl(method = "repeatedcv", number = 2, repeats = 2)
  ctrl_voom <- voomControl(method = "repeatedcv", number = 2, repeats = 2)
  cat("Control parameters defined\n")

  cat("Defining model configurations...\n")
  all_methods <- availableMethods()
  
  # Separate methods into categories
  continuous_methods <- all_methods[all_methods %in% c("svmLinear", "svmRadial", "svmPoly", "rf", "knn", "nb", "nnet", "rpart", "glmnet", "pam")]
  discrete_methods <- all_methods[all_methods %in% c("PLDA", "PLDA2", "NBLDA")]
  voom_methods <- all_methods[all_methods %in% c("voomNSC", "voomDLDA", "voomDQDA")]
  
  model_configs <- c(
    sapply(continuous_methods, function(method) list(method = method, preProcessing = "deseq-vst", control = ctrl), simplify = FALSE, USE.NAMES = TRUE),
    sapply(discrete_methods, function(method) list(method = method, normalize = "deseq", control = ctrl_discrete), simplify = FALSE, USE.NAMES = TRUE),
    sapply(voom_methods, function(method) list(method = method, normalize = "deseq", control = ctrl_voom), simplify = FALSE, USE.NAMES = TRUE)
  )
  
  # Limit the number of models if specified
  if (!is.null(max_models) && max_models < length(model_configs)) {
    model_configs <- model_configs[1:max_models]
  }
  
  cat(paste("Number of models to be trained:", length(model_configs), "\n"))
  cat("Model configurations defined\n")

  cat("Starting parallel model training. This may take a while...\n")
  models <- foreach(config = model_configs, .packages = c("MLSeq", "caret"), 
                    .errorhandling = "pass") %dopar% {
    model_start_time <- Sys.time()
    set.seed(123)
    tryCatch({
      cat(paste("Training model:", names(config), "\n"))
      if ("preProcessing" %in% names(config)) {
        cat(paste("Using preProcessing:", config$preProcessing, "\n"))
        model <- classify(data = train_data, method = config$method, 
                          preProcessing = config$preProcessing, ref = ref_class, 
                          control = config$control)
      } else {
        cat(paste("Using normalization:", config$normalize, "\n"))
        model <- classify(data = train_data, method = config$method, 
                          normalize = config$normalize, ref = ref_class, 
                          control = config$control)
      }
      model_end_time <- Sys.time()
      cat(paste("Model", names(config), "training complete\n"))
      list(model = model, time = difftime(model_end_time, model_start_time, units = "mins"))
    }, error = function(e) {
      cat(paste("Error in training model", names(config), ":", e$message, "\n"))
      list(error = e$message, time = difftime(Sys.time(), model_start_time, units = "mins"))
    })
  }
  
  cat("Stopping parallel cluster...\n")
  stopCluster(cl)
  cat("Parallel processing complete\n")

  names(models) <- names(model_configs)

  cat("Extracting model information...\n")
  model_info <- lapply(names(models), function(model_name) {
    cat(paste("Processing model:", model_name, "\n"))
    model_result <- models[[model_name]]
    if (!is.null(model_result$error)) {
      cat(paste("Error occurred for model", model_name, "\n"))
      return(data.frame(Model = model_name, Error = model_result$error, Time = model_result$time))
    }
    
    # Extract relevant information from the model
    model <- model_result$model
    info <- tryCatch({
      data.frame(
        Model = model_name,
        Method = ifelse(is.null(model@modelInfo@method), NA, model@modelInfo@method),
        Transformation = ifelse(is.null(model@modelInfo@transformation), NA, model@modelInfo@transformation),
        Normalization = ifelse(is.null(model@modelInfo@normalization), NA, model@modelInfo@normalization),
        PreProcessing = ifelse(is.null(model@modelInfo@preProcessing), NA, model@modelInfo@preProcessing),
        Reference = ifelse(is.null(model@modelInfo@ref), NA, model@modelInfo@ref),
        Accuracy = ifelse(is.null(model@modelInfo@confusionMat$overall["Accuracy"]), NA, model@modelInfo@confusionMat$overall["Accuracy"]),
        Kappa = ifelse(is.null(model@modelInfo@confusionMat$overall["Kappa"]), NA, model@modelInfo@confusionMat$overall["Kappa"]),
        Time = model_result$time
      )
    }, error = function(e) {
      cat(paste("Error extracting information for model", model_name, ":", e$message, "\n"))
      data.frame(
        Model = model_name,
        Error = e$message,
        Time = model_result$time
      )
    })
    
    return(info)
  })
  
  # Combine all model information, filling in missing columns with NA
  all_columns <- unique(unlist(lapply(model_info, names)))
  model_summary <- do.call(rbind, lapply(model_info, function(info) {
    missing_cols <- setdiff(all_columns, names(info))
    info[missing_cols] <- NA
    info[all_columns]
  }))
  
  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "mins")
  
  cat("\nModel summary:\n")
  print(model_summary)
  
  cat(paste("\nTotal execution time:", round(total_time, 2), "minutes\n"))
  cat(paste("Average time per model:", round(mean(model_summary$Time, na.rm = TRUE), 2), "minutes\n"))
  
  return(list(summary = model_summary, models = models, total_time = total_time))
}

# 7 - Main function
main <- function() {
  cat("Entering main function\n")
  
  # Get command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # Print received arguments
  cat("Number of arguments received:", length(args), "\n")
  cat("Arguments received:\n")
  for (i in seq_along(args)) {
    cat("  Argument", i, ":", args[i], "\n")
    if (file.exists(args[i])) {
      cat("    File exists\n")
    } else {
      cat("    File does not exist\n")
    }
  }
  
  # Check if we have the correct number of arguments
  if (length(args) != 3) {
    stop("Error: Three arguments are required: filtered_data_path, metadata_path, patients_path")
  }
  
  # Assign arguments to variables
  filtered_data_path <- args[1]
  metadata_path <- args[2]
  patients_path <- args[3]
  
  cat("\nParsed arguments:\n")
  cat("  filtered_data_path:", filtered_data_path, "\n")
  cat("  metadata_path:", metadata_path, "\n")
  cat("  patients_path:", patients_path, "\n")
  
  # Load data
  data_list <- load_data(filtered_data_path, metadata_path, patients_path)

  # Filter and add metadata
  data_list <- filter_and_add_metadata(data_list)
  
  # Filter and align data
  data_list <- filter_and_align_data(data_list)
  
  # Prepare data for machine learning
  ml_data <- prepare_ml_data(data_list)
  
  # Create DESeq objects
  deseq_objects <- create_deseq_objects(ml_data)
  
  # Run multiple machine learning models
  results <- run_mlseq_models_parallel(deseq_objects, ref_class = "Case", n_cores = NULL, max_models = NULL)
  
  # Save results
  output_dir <- dirname(filtered_data_path)
  output_path <- file.path(output_dir, "ml_model_results.csv")
  write.csv(results$summary, file = output_path, row.names = FALSE)
  
  cat(paste("Results saved to:", output_path, "\n"))
  
  # Return the results (this will be captured in the R script's output)
  return(results$summary)
}

# Run the main function
tryCatch({
  main()
  cat("R script completed successfully\n")
}, error = function(e) {
  cat("Error in R script:", conditionMessage(e), "\n")
  quit(status = 1)
})