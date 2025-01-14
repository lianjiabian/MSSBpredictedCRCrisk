################# 0. Load Required R Packages #################
library(glmnet)
library(pbapply)
library(tidyverse)
library(reshape2)
library(openxlsx)
library(DALEX)
library(readr)
library(gbm)
library(dplyr)
library(caret)
library(ggplot2)
library(pROC)
library(rms)
library(rmda)
library(dcurves)
library(Hmisc)
library(ResourceSelection)
library(DynNom)
library(survey)
library(foreign)
library(plotROC)
library(survival)
library(shapper)
library(iml)
library(e1071)
library(ROCR)
library(corrplot)
library(lattice)
library(Formula)
library(SparseM)
library(riskRegression)
library(pheatmap)
library(fastshap)
library(ingredients)
library(mlr3)
library(table1)
library(tableone)
library(adabag)
library(RColorBrewer)
library(VIM)
library(mice)
library(autoReg)
library(cvms)
library(tibble)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(ROSE)
library(DMwR)
library(scales)
library(kernelshap)
library(shapviz)
library(caret)       
library(rpart)       
library(rpart.plot) 
library(randomForest)  
library(xgboost)         
library(lightgbm)      
library(kknn)        
library(neuralnet)    
library(NeuralNetTools) 
library(gridExtra) 
library(partykit)
library(missForest)
library(MLmetrics)
library(readxl)

################# 1. Preprocessing: Calculate LPKM-SSB #################
# Normalize SSB count data and scale by 1e5
ssb_data <- count_SSB
row_sums <- rowSums(ssb_data)
normalized_ssb <- sweep(ssb_data, 1, row_sums, FUN = "/") * 1e5

# Save the processed data
write.csv(normalized_ssb, file = "W_LPKM_SSB.csv")

################# 2. Calculate F-Score to Differentiate HC, PL, and CRC #################
# Extract training sample information
train_sample_info <- subset(sample_ifo, sample_ID %in% train_W_LPKM_SSB$...1)

# Create a binary result column based on labels
train_sample_info <- train_sample_info %>%
  mutate(result = case_when(
    label == "CRC" ~ 1,
    label == "PL" ~ 1,
    label == "HC" ~ 0,
    TRUE ~ NA_real_  # For unmatched cases
  )) %>%
  select(sample_ID, result)

# Merge with training SSB data
merged_data <- merge(train_sample_info, train_W_LPKM_SSB, by.x = "sample_ID", by.y = "...1")

# Separate label and gene expression data
gene_data <- merged_data %>%
  select(-sample_ID) %>%
  rename(label = result)

# Ensure all columns are numeric
gene_data[] <- lapply(gene_data, function(x) as.numeric(as.character(x)))

# Threshold calculation and binarization
threshold <- mean(as.matrix(gene_data[, -1]), na.rm = TRUE)
binarized_genes <- as.data.frame(lapply(gene_data[, -1], function(x) ifelse(x > threshold, 1, 0)))

# Calculate F1 scores for all genes
f1_scores <- sapply(names(binarized_genes), function(gene) {
  F1_Score(gene_data$label, binarized_genes[[gene]], positive = "1")
})

# Save F1 scores
f1_scores_df <- data.frame(Gene = names(f1_scores), F1_Score = f1_scores)
write.csv(f1_scores_df, file = "train_w_ssb_f_scores.csv")

# Filter genes with F1 score >= 0.85
high_f1_genes <- f1_scores_df %>% filter(F1_Score >= 0.85)
write.csv(high_f1_genes, file = "train_w_ssb_f_scores_085.csv")

################# 3. Feature Selection Using LASSO #################
# Normalize data
raw_data <- w_train %>%
  select(-sample_ID) %>%
  mutate(across(everything(), as.numeric))

# Standardize feature variables
data_scaled <- scale(raw_data[-1])
data_scaled <- data.frame(Result = raw_data$result, data_scaled)

# LASSO regression
lasso_model <- glmnet(as.matrix(data_scaled[-1]), data_scaled$Result, family = "binomial", alpha = 1)
cv_lasso <- cv.glmnet(as.matrix(data_scaled[-1]), data_scaled$Result, family = "binomial", alpha = 1, nfolds = 10)
best_lambda <- cv_lasso$lambda.min

# Extract selected features
selected_features <- coef(cv_lasso, s = best_lambda) %>%
  as.matrix() %>%
  .[. != 0, , drop = FALSE] %>%
  rownames()

# Save selected features
save(selected_features, file = "lasso_feature.Rdata")

################# 4. Binary Classification Modeling #################
# Train and test datasets
load("lasso_feature.Rdata")
load("w_final_data.Rdata")

# Align train and test datasets with selected features
matched_columns <- intersect(selected_features, colnames(w_train))
w_train <- w_train[, matched_columns]
w_test1 <- w_test1[, matched_columns]
w_test2 <- w_test2[, matched_columns]

# Save processed datasets
save(w_train, w_test1, w_test2, file = "w_ssb_lasso_model_data.Rdata")

# Placeholder: Model training and evaluation scripts (logistic regression, KNN, CART, RF, XGBoost, LightGBM, SVM, NN).
# Include detailed sections for each model as needed.

####### 5. Model Calibration and Threshold Performance Analysis#######
# Load necessary libraries
library(runway)
library(readxl)
library(ggplot2)

# Clear the workspace
rm(list = ls())

# Load the dataset
data <- read_excel("test1_pl_crc.xlsx")
colnames(data)

# Define the range of model prediction columns
model_columns <- colnames(data)[c(4:36, 38:40)]  # Adjust this range as needed

# Create a folder for saving the plots
output_dir <- "./Model_test1_pl_crc_Results_Calibration"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)  # Create the folder if it doesn't exist
}

# Iterate over each model column
for (model in model_columns) {
  tryCatch({
    # Create a dataset for the current model
    single_model_dataset <- data.frame(
      outcomes = data$result,      # Actual outcomes
      predictions = data[[model]]  # Predictions for the current model
    )
    
    # Generate the threshold performance plot
    p1 <- threshperf_plot(
      df = single_model_dataset,
      outcome = "outcomes",
      prediction = "predictions",
      positive = "1"  # Ensure the positive class matches the data
    )
    
    # Save the threshold performance plot
    file <- file.path(output_dir, paste0(model, "_Threshold_Performance.pdf"))
    ggsave(filename = file, plot = p1, width = 7, height = 6)
    
    # Generate the calibration plot
    p2 <- cal_plot(
      single_model_dataset,
      outcome = "outcomes",
      prediction = "predictions",
      positive = "1",  # Ensure the positive class matches the data
      n_bins = 6       # Number of bins for calibration
    )
    
    # Save the calibration plot
    file <- file.path(output_dir, paste0(model, "_Calibration_Curve.pdf"))
    ggsave(filename = file, plot = p2, width = 7, height = 6)
    
  }, error = function(e) {
    # Capture errors and provide warnings
    warning(paste("Error processing model:", model, "-", e$message))
  })
}

# Display completion message
message("All plots have been saved to the folder: ", output_dir)

# Display warnings, if any
warnings()





####### 6. Model Performance Comparison - AUC #####

# Clear the workspace
rm(list = ls())

# Load required packages
if (!require(pROC)) {
  install.packages("pROC")
}
library(pROC)
library(caret)
library(ggplot2)
library(dplyr)
library(readxl)

# Load example data
data <- read_excel("train_pl_crc.xlsx")
print(colnames(data))

# Define columns containing model predictions
model_columns <- colnames(data)[c(12:19, 28:36, 39:40)]  # Adjust to match your column range

# Define true labels
true_labels <- factor(data$result)  # Ensure the 'result' column exists and is binary

# Create a folder for saving results
output_dir <- "./train_pl_crc_Results_ROC"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)  # Create directory if it does not exist
}

# Initialize a list to store AUC results and confidence intervals
auc_list <- list()

# Iterate through each model column
for (model in model_columns) {
  tryCatch({
    # Extract predicted probabilities
    predicted_probs <- as.numeric(data[[model]])
    
    # Ensure there are no missing values
    if (any(is.na(predicted_probs))) {
      stop("Missing values detected in predicted probabilities")
    }
    
    # Calculate ROC and AUC
    roc_obj <- roc(response = true_labels, predictor = predicted_probs, levels = rev(levels(true_labels)))
    auc_value <- auc(roc_obj)
    ci_auc <- ci.auc(roc_obj, conf.level = 0.95)
    
    # Create legend text for AUC
    legend_text <- paste0(model, " (AUC = ", round(auc_value, 3), 
                          ", 95% CI: ", round(ci_auc[1], 3), "-", round(ci_auc[3], 3), ")")
    
    # Save AUC results to the list
    auc_list[[model]] <- data.frame(
      Model = model,
      AUC = round(auc_value, 3),
      CI_Lower = round(ci_auc[1], 3),
      CI_Upper = round(ci_auc[3], 3)
    )
    
    # Save ROC plot as a PDF
    roc_file <- file.path(output_dir, paste0(model, "_ROC.pdf"))
    pdf(roc_file, width = 7, height = 6)
    plot(
      roc_obj,
      col = "red",
      lwd = 2,
      main = "ROC Curve",
      xlab = "1 - Specificity",
      ylab = "Sensitivity",
      legacy.axes = TRUE,
      cex.main = 1.6,
      cex.lab = 1.3,
      cex.axis = 1.2
    )
    legend("bottomright", legend = legend_text, col = "red", lty = 1, lwd = 2, cex = 0.6)
    dev.off()
    message("ROC curve saved to: ", roc_file)
    
  }, error = function(e) {
    warning(paste("Error processing model:", model, "-", e$message))
  })
}

# Combine all AUC results
final_auc_df <- do.call(rbind, auc_list)

# View results
print(head(final_auc_df))

# Save AUC results to a CSV file
auc_file <- file.path(output_dir, "AUC_with_CI.csv")
write.csv(final_auc_df, file = auc_file, row.names = FALSE)
message("All AUC results saved to: ", auc_file)

# Visualize AUC and confidence intervals for all models
ggplot(final_auc_df, aes(x = reorder(Model, AUC), y = AUC)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.5, color = "orange") +
  coord_flip() +  # Flip coordinates for better readability
  labs(title = "Model AUC and 95% Confidence Intervals",
       x = "Models",
       y = "AUC") +
  theme_minimal()

# Save AUC comparison plot
auc_plot_file <- file.path(output_dir, "AUC_Comparison_Plot.pdf")
ggsave(filename = auc_plot_file, width = 8, height = 6)
message("AUC comparison plot saved to: ", auc_plot_file)

####### 7. Decision Curve Analysis (DCA) #####

library(dcurves)
rm(list = ls())

# Load example data
data <- read_excel("test1.xlsx")
colnames(data)

# Define model columns for DCA
model_columns <- colnames(data)[c(4:11)]  # Adjust as needed

# Prepare data for DCA
dca_data <- data.frame(y_true = data$result)
for (col in model_columns) {
  dca_data[[col]] <- as.numeric(data[[col]])
}

# Create formula dynamically
dca_formula <- paste("y_true ~", paste(model_columns, collapse = " + "))

# Perform DCA
dca_res <- dca(as.formula(dca_formula), data = dca_data)

# Plot DCA results
colors <- c("#fe9f2a", "#ea4737", "#b14945", "#b67658", "#83bd6c", "#167cb4", "#80bedd", "#683c8c")
plot(
  dca_res,
  main = "Decision Curve Analysis for Multiple Models",
  xlab = "Threshold Probability",
  ylab = "Net Benefit",
  col = colors
)

####### 8.Confusion Matrix and Metrics Export #####

library(caret)
rm(list = ls())

# Load example data
data <- read_excel("test2.xlsx")
colnames(data)

# Define columns for model predictions
model_columns <- colnames(data)[c(4:11, 20:27)]  # Adjust as needed

# Define true labels
true_labels <- factor(data$result)

# Create directory for saving results
output_dir <- "./Confusion_Results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Initialize lists for confusion matrices and metrics
all_confusion_tables <- list()
all_stats <- list()

# Iterate through each model column
for (model in model_columns) {
  tryCatch({
    # Extract predicted labels
    predicted_labels <- factor(data[[model]])
    
    # Compute confusion matrix
    confusion_matrix <- confusionMatrix(predicted_labels, true_labels)
    
    # Save confusion matrix table
    confusion_table <- as.data.frame(as.table(confusion_matrix$table))
    confusion_table$Model <- model
    all_confusion_tables[[model]] <- confusion_table
    
    # Save performance metrics
    stats <- c(confusion_matrix$overall, confusion_matrix$byClass)
    stats <- as.data.frame(t(stats))
    stats$Model <- model
    all_stats[[model]] <- stats
    
  }, error = function(e) {
    warning(paste("Error processing model:", model, "-", e$message))
  })
}

# Combine and save results
final_confusion_table <- do.call(rbind, all_confusion_tables)
final_stats <- do.call(rbind, all_stats)

write.csv(final_confusion_table, file = file.path(output_dir, "Confusion_Tables.csv"), row.names = FALSE)
write.csv(final_stats, file = file.path(output_dir, "Metrics.csv"), row.names = FALSE)

message("Confusion matrices and metrics saved to: ", output_dir)
#######9.Boxplot with Statistical Comparison######
# Load necessary libraries
library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(ggsignif)

# Step 1: Load Excel data
# Adjust file name and path as needed
data <- read_excel("data1.xlsx")

# Step 2: Rename selected columns
# Update column indices to match your dataset
colnames(data)[5:7] <- c("HC", "PL", "CRC")

# Step 3: Transform data to long format
df_long <- data %>%
  select(HC, PL, CRC) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Group",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))  # Remove missing values

# Step 4: Specify group order (HC -> PL -> CRC)
df_long$Group <- factor(df_long$Group, levels = c("HC", "PL", "CRC"))

# Step 5: Define custom colors
my_colors <- c("HC" = "#ef7f51", "PL" = "#78d3ac", "CRC" = "#9355b0")

# Step 6: Create boxplot with overlayed points
p <- ggplot(df_long, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Boxplot without outliers
  scale_fill_manual(values = my_colors) +         # Custom colors
  theme_bw(base_size = 14) +                      # Black-and-white theme
  coord_cartesian(ylim = c(-0.5, 1.5))            # Adjust y-axis limits

# Step 7: Define comparisons and add significance test
comparisons_list <- list(
  c("HC", "PL"),
  c("PL", "CRC"),
  c("HC", "CRC")
)

p <- p + geom_signif(
  comparisons = comparisons_list,
  test = "t.test",           # Change to "wilcox.test" if needed
  map_signif_level = FALSE,  # Show exact p-values
  label = "p.format",        # Display p-value in scientific format
  step_increase = 0.2        # Adjust position of significance annotations
)

# Display the plot
print(p)

# Step 8: Extract and save boxplot data
plot_data <- ggplot_build(p)$data[[1]]

# Handle outliers: Convert list of outliers into comma-separated strings
plot_data$outliers <- sapply(plot_data$outliers, function(x) {
  if (length(x) == 0) {
    return(NA)  # Replace missing outliers with NA
  } else {
    return(paste(x, collapse = ","))
  }
})

# Save plot data to a CSV file
write.csv(plot_data, file = "boxplot_data.csv", row.names = FALSE)
#######10.Binary Classification Using Logistic Regression#######
# Convert the status variable into a binary outcome
pbc$died <- ifelse(pbc$status == 1, 1, 0)  # Status 1 indicates death; convert to binary outcome
pbc$died <- as.factor(pbc$died)  # Convert to a factor type

# Build a logistic regression model
logit_model <- glm(
  died ~ MSSB_KNN + MSSB_NNET + MSSB_LR + Age + Gender,  # Predictor and outcome variables
  data = pbc, 
  family = binomial(link = "logit")  # Specify logistic regression
)

# Generate a nomogram using regplot
library(regplot)  # Install with install.packages("regplot") if not already installed
regplot(
  logit_model,
  observation = pbc[1, ],  # Specify a single observation
  failtime = NULL,         # failtime is not needed for binary classification
  prfail = TRUE,           # Display predicted probabilities
  title = "Nomogram for Binary Classification"  # Add a title
)

# Extract model summary
model_summary <- summary(logit_model)

# Extract the coefficient table
coef_table <- model_summary$coefficients

# Extract p-values for each predictor (including intercept)
p_values <- coef_table[, "Pr(>|z|)"]

# Output p-values
p_values

#########11.Calculate Log-Loss for Models#######
# Load the required library
library(readxl)

# Clear the workspace
rm(list = ls())

# Load the data
data <- read_excel("test2.xlsx")

# Define the Log-Loss function
log_loss <- function(y_true, y_pred) {
  # Ensure predictions are numeric
  y_pred <- as.numeric(y_pred)
  
  # Handle missing or invalid prediction values
  epsilon <- 1e-15
  y_pred <- pmin(pmax(y_pred, epsilon), 1 - epsilon)  # Clip values to avoid 0 or 1
  n <- length(y_true)
  
  # Compute Log-Loss
  loss <- -1/n * sum(y_true * log(y_pred) + (1 - y_true) * log(1 - y_pred), na.rm = TRUE)
  return(loss)
}

# Create an empty data frame to store Log-Loss results
log_loss_results <- data.frame(Model = character(0), Log_Loss = numeric(0))

# Specify columns containing model predictions
model_columns <- colnames(data)[c(4:35)]  # Adjust column indices as necessary

# Compute Log-Loss for each model
for (model in model_columns) {
  model_name <- model  # Use column name as the model name
  
  # Compute Log-Loss using the current model's predictions
  t_log_loss <- log_loss(data$result, data[[model]])
  
  # Append the result to the data frame
  log_loss_results <- rbind(log_loss_results, data.frame(Model = model_name, Log_Loss = t_log_loss))
}

# Save Log-Loss results to CSV files
write.csv(log_loss_results, file = "test2_log_loss.csv", row.names = FALSE)

# Print a message indicating completion
message("Log-Loss results saved to test2_log_loss.csv")

#######12.ROC Curve and DeLong Test for Pairwise AUC Comparison########
# Load required libraries
library(readxl)
library(pROC)

# Clear the workspace
rm(list = ls())

# Load the dataset
data <- read_excel("test2.xlsx")

# Define model columns for ROC calculation
model_columns <- colnames(data)[c(12:19)]  # Adjust indices to match your dataset

# Create a data frame to store AUC and pairwise comparison results
auc_results <- data.frame(
  Model_1 = character(0),
  Model_2 = character(0),
  AUC_1 = numeric(0),
  AUC_2 = numeric(0),
  AUC_diff = numeric(0),
  CI_Lower = numeric(0),
  CI_Upper = numeric(0),
  Z_statistic = numeric(0),
  p_value = numeric(0)
)

# Compute ROC curves and perform pairwise DeLong tests
for (i in 1:(length(model_columns) - 1)) {
  for (j in (i + 1):length(model_columns)) {
    
    # Get model names
    model_1 <- model_columns[i]
    model_2 <- model_columns[j]
    
    # Compute ROC curves
    roc1 <- roc(data$result, data[[model_1]])
    roc2 <- roc(data$result, data[[model_2]])
    
    # Perform DeLong test
    test_result <- roc.test(roc1, roc2, method = "delong")
    
    # Store results in the data frame
    auc_results <- rbind(auc_results, data.frame(
      Model_1 = model_1,
      Model_2 = model_2,
      AUC_1 = auc(roc1),
      AUC_2 = auc(roc2),
      AUC_diff = auc(roc1) - auc(roc2),
      CI_Lower = test_result$conf.int[1],
      CI_Upper = test_result$conf.int[2],
      Z_statistic = test_result$statistic,
      p_value = test_result$p.value
    ))
  }
}

# Print the comparison results
print(auc_results)

# Save the results to a CSV file
write.csv(auc_results, file = "test2_auc_delong.csv", row.names = FALSE)

message("Pairwise AUC comparison results saved to: test2_auc_delong.csv")

########13.KNN Modeling with SHAP Value Analysis#####
# Load necessary libraries
library(caret)
library(kernelshap)
library(shapviz)
library(ggplot2)
library(gridExtra)

# Train a KNN model using the caret package
knn_model <- train(
  Result ~ ., 
  data = traindata_scaled,  
  method = "kknn",  
  tuneGrid = expand.grid(
    kmax = best_k_auc,     # Maximum K value
    distance = 2,          # Distance metric (2 = Euclidean)
    kernel = "triangular"  # Kernel function
  )
)

# Compute SHAP values using kernelshap
explain_kernel_knn <- kernelshap(
  knn_model,
  traindata_scaled[1:n_train, -1],  # Training data without the target column
  bg_X = testdata_scaled[1:n_test, -1],  # Background data from the test set
  pred_fun = function(model, X) {
    probs <- predict(model, X, type = "prob")  # Predict probabilities
    return(as.numeric(probs[, "Yes"]))        # Return probabilities for the "Yes" class
  }
)
# Bar chart for feature importance
sv_importance(explain_kernel_knn, kind = "bar", fill = "#fca50a") +
  theme_bw() +
  ggtitle("Feature Importance - Bar Chart") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Beeswarm plot for feature importance
sv_importance(explain_kernel_knn, kind = "beeswarm") +
  theme_bw() +
  ggtitle("Feature Importance - Beeswarm Plot") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Combined bar and beeswarm plot
sv_importance(explain_kernel_knn, kind = "both") +
  theme_bw() +
  ggtitle("Combined Feature Importance") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# SHAP dependence plot for a specific feature
sv_dependence(
  explain_kernel_knn,
  v = "ANO3",             # Specify the primary feature
  color_var = "CTD-2231H16.1",  # Specify the secondary feature
  size = 4
) +
  theme_bw() +
  ggtitle("Dependence Plot for ANO3") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
#########14.MSSB and CEA Scatter Plot with Septin 9########
# Load necessary packages
library(readxl)   # For reading Excel files
library(ggplot2)  # For plotting
library(dplyr)    # For data manipulation

# Read the Excel file (adjust the file name and path as needed)
a <- read_excel("test1.xlsx")

# Create a data frame for plotting
df <- data.frame(
  mssb = a$MSSB_knn_prob,          # Replace with actual column name
  log5CEA = log(a$CEA_raw) / log(5), # Replace with actual column name
  septin9 = a$septin9_raw          # Replace with actual column name
)

# Add classification labels
df$reault <- ifelse(df$mssb >= 1, "CRC or PL", "HC")
df$mssb_label <- ifelse(df$mssb >= 0.5, "MSSB(+)", "MSSB(-)")
df$cea_label <- ifelse(df$log5CEA >= 1, "CEA(+)", "CEA(-)")
df$septin9_label <- ifelse(df$septin9 == "positive", "Septin9(+)", "Septin9(-)")

# Define thresholds
threshold_cea <- 1
threshold_mssb <- 0.5

# Define colors for Septin 9 categories
my_colors <- c("negative" = "blue", "positive" = "red")

# Create the scatter plot
p <- ggplot(df, aes(x = log5CEA, y = mssb)) +
  geom_point(aes(color = factor(septin9), shape = reault), size = 3, alpha = 0.8) +
  geom_vline(xintercept = threshold_cea, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = threshold_mssb, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = c("CRC or PL" = 17, "HC" = 19)) +
  labs(
    x = "log5CEA",
    y = "MSSB",
    color = "Septin 9",
    shape = "Result",
    title = "Scatter Plot: MSSB vs. log5CEA with Septin 9"
  ) +
  theme_bw(base_size = 14)

# Display the plot
print(p)

# Optional: Restrict axis limits
p + coord_cartesian(xlim = c(-2, 4), ylim = c(-0.5, 1.5))

# Summarize data by CEA and MSSB categories
df_summary <- df %>%
  mutate(
    cea_cat = ifelse(log5CEA >= threshold_cea, "CEA+", "CEA-"),
    mssb_cat = ifelse(mssb >= threshold_mssb, "MSSB+", "MSSB-")
  ) %>%
  group_by(cea_cat, mssb_cat) %>%
  summarise(
    count = n(),
    septin9_pos = sum(septin9 == "positive"),
    septin9_neg = sum(septin9 == "negative"),
    .groups = "drop"
  ) %>%
  mutate(
    total = sum(count),
    pct = paste0(round(100 * count / sum(count), 1), "%")
  )

# View summary table
print(df_summary)

# Annotate plot with summary data
p_annotated <- p +
  annotate("text", x = -1, y = 1.2, label = "CEA- / MSSB+\nn= 63 (%)", size = 5) +
  annotate("text", x = -1, y = -0.2, label = "CEA- / MSSB-\nn= 16 (%)", size = 5) +
  annotate("text", x = 3, y = 1.2, label = "CEA+ / MSSB+\nn= 18.5 (%)", size = 5) +
  annotate("text", x = 3, y = -0.2, label = "CEA+ / MSSB-\nn= 2.5 (%)", size = 5) +
  coord_cartesian(xlim = c(-2, 4), ylim = c(-0.5, 1.5))

# Display annotated plot
print(p_annotated)

#########15.Radar Chart for Multiple Comparisons########

# Load necessary package
library(fmsb)
# Original dataset with comparisons
raw_data <- data.frame(
  Compare = c("HC VS CRC&PL (C1)",
              "HC VS PL (C1)",
              "HC VS CRC (C1)",
              "PL VS CRC (C1)",
              "HC VS PL (C2)"),
  AUC                = c(0.87, 0.87, 0.86, 0.53, 0.92),
  Accuracy           = c(0.87, 0.78, 0.83, 0.62, 0.89),
  Kappa              = c(0.65, 0.57, 0.62, 0.01, 0.73),
  Sensitivity        = c(0.60, 0.60, 0.60, 0.02, 0.71),
  Specificity        = c(0.98, 0.98, 0.99, 0.99, 0.97),
  F1                 = c(0.93, 0.97, 0.97, 0.50, 0.91),
  Balanced_Accuracy  = c(0.79, 0.79, 0.79, 0.50, 0.84)
)

# Set row names to comparison names and remove the first column
df <- raw_data
rownames(df) <- df$Compare
df <- df[ , -1]

# Create max and min rows for radar chart scaling
max_min <- data.frame(
  AUC               = c(1, 0),
  Accuracy          = c(1, 0),
  Kappa             = c(1, 0),
  Sensitivity       = c(1, 0),
  Specificity       = c(1, 0),
  F1                = c(1, 0),
  Balanced_Accuracy = c(1, 0)
)
rownames(max_min) <- c("Max", "Min")

# Combine max/min rows with the data
df_for_radar <- rbind(max_min, df)

# Basic radar chart
radarchart(
  df_for_radar,
  pcol = c("red", "blue", "darkgreen", "brown", "purple"),  # Line colors
  pfcol = scales::alpha(c("red", "blue", "darkgreen", "brown", "purple"), 0.3),  # Fill colors
  plwd = 2,           # Line width
  plty = 1,           # Line type
  cglcol = "grey",    # Gridline color
  cglty = 1,          # Gridline type
  cglwd = 0.8,        # Gridline width
  axislabcol = "grey",# Axis label color
  seg = 5,            # Number of grid segments
  title = "Radar Chart for Multiple Comparisons",
  vlcex = 0.8         # Size of variable labels
)

# Add a legend
legend(
  x = 2, 
  y = 1,
  legend = rownames(df_for_radar)[-c(1, 2)],  # Remove "Max" and "Min" from the legend
  col = c("red", "blue", "darkgreen", "brown", "purple"),
  bty = "n",        # No box around the legend
  pch = 20,         # Points in the legend
  pt.cex = 1.5,     # Point size in the legend
  text.col = "black" # Text color
)

##########16.Performance Indicator Comparison########

# Load necessary libraries
library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)

# Custom colors for the plot
my_colors <- c(
  "#fe9f2a", "#ea4737", "#b14945", "#b67658",
  "#83bd6c", "#167cb4", "#80bedd", "#683c8c"
)

# Clear the workspace
rm(list = ls())

# Read in the data
All_Stats <- read_csv("Confusion_test2_Results/All_Stats.csv")

# Select relevant columns (example: accuracy, F1, precision, recall, etc.)
selected_data <- All_Stats[1:8, c(1:2, 8:14, 18:19)]

# Remove unnecessary columns (e.g., standard deviation)
processed_data <- selected_data[, -11]

# Set model names as rownames and remove the 'model' column
rownames(processed_data) <- processed_data$model
processed_data <- processed_data[, -1]

# Convert rownames into a column for reshaping
df_wide <- processed_data %>% rownames_to_column(var = "Model")

# Reshape data to long format for ggplot
df_long <- df_wide %>%
  pivot_longer(
    cols = !Model,               # All columns except "Model"
    names_to = "Performance_Indicator",
    values_to = "Performance_Value"
  )

# Create the performance indicator comparison plot
p <- ggplot(
  df_long,
  aes(
    x = Performance_Indicator,
    y = Performance_Value,
    group = Model,
    color = Model
  )
) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = my_colors) +  # Custom colors for models
  labs(
    x = "Performance Indicator",
    y = "Performance Value",
    title = "Performance Comparison Across Models"
  ) +
  scale_y_continuous(limits = c(0, 1)) +    # Ensure values stay in [0, 1]
  theme_minimal(base_size = 14) +          # Clean theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Display the plot
print(p)


