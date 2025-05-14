# ------------------------------------------------------------------------------
# Script:       convert_txt_to_csv.R
# Purpose:      Converts data from txt format to csv format
# Author:       Sophia Li
# Affiliation:  Sadreyev lab
# Date:         2025-05-10
#
# Inputs:       
#   - txt_path:           Path to a txt file <input_data>.txt
#
# Outputs:   
#   - <input_data>.csv:   Input data in the form of a csv
#
# Version:      R version 4.4.2 (2024-10-31)
# ------------------------------------------------------------------------------

# ===| A. Convert TXT to CSV |==================================================

convert_txt_to_csv <- function(txt_path) {
  
  # --- Read the delimited text file ---
  data <- read.delim(txt_path, header = TRUE, sep = "\t", row.names = 1,
                     check.names = FALSE)
  
  # --- Write the data as a CSV of the same name ---
  csv_path <- sub("\\.txt$", ".csv", txt_path)
  write.csv(data, file = csv_path, row.names = TRUE)
  message("Wrote CSV to: ", csv_path)
  
  # --- Return invisibly ---
  invisible(NULL)
}