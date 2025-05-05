# Load necessary libraries
library(dplyr)
library(data.table)

# Load data using data.table for efficiency
Full.length.assignments <- fread("Full.length.assignment.csv", data.table = FALSE)
ITS.full.assignments <- fread("ITS.full.assignment.csv", data.table = FALSE)
ITS1.assignments <- fread("ITS1.assignment.csv", data.table = FALSE)
ITS2.assignments <- fread("ITS2.assignment.csv", data.table = FALSE)
O18SV9.assignments <- fread("18SV9.assignment.csv", data.table = FALSE)

# Ensure all datasets have the same row order based on 'Full.length.assignments'
datasets <- list(ITSFull = ITS.full.assignments, ITS1 = ITS1.assignments, ITS2 = ITS2.assignments, O18SV9 = O18SV9.assignments)
for (i in seq_along(datasets)) {
  datasets[[i]] <- datasets[[i]][match(row.names(Full.length.assignments), row.names(datasets[[i]])), ]
}

# Function to calculate a score based on non-NA values in a row
calculate_score <- function(row) {
  sum(!is.na(row))
}

# Initialize the final data frame with data from 'Full.length.assignments' and an additional 'source' column
final_data <- Full.length.assignments
final_data$source <- "Full.length"

# Iterate over each row of 'Full.length.assignments'
for (i in 1:nrow(final_data)) {
  base_score <- calculate_score(final_data[i, 1:8])
  for (dataset_name in names(datasets)) {
    dataset <- datasets[[dataset_name]]
    alt_score <- calculate_score(dataset[i, ])
    # Replace row if another dataset provides a row with a higher score
    if (alt_score > base_score) {
      final_data[i, ] <- dataset[i, ]
      final_data$source[i] <- dataset_name
      break  # Stop checking once a better row is found
    }
  }
}

# Calculate scores for the annotations in each row and add as a new column
final_data$score <- apply(final_data[, 2:8], 1, calculate_score)  # excluding 'source' from the score calculation

# Review the final structure of the data
str(final_data)

# Save the final dataset to a CSV file
write.csv(final_data, "final_merged_data.csv", row.names = FALSE)
