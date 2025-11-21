library(tidyverse)
library(jsonlite)

# Function to fetch disorder data from IEMbase API
fetch_disorder <- function(id) {
  url <- paste0("https://www.iembase.com/api/v2/disorder/", id)

  tryCatch({
    # Add a 1-second sleep to be respectful to the API
    Sys.sleep(1)

    # Fetch and parse JSON
    data <- fromJSON(url)

    # Extract key fields (flatten the nested structure)
    result <- tibble(
      disorder_id = id,
      name = data$name %||% NA,
      name_alt1 = data$name_alt1 %||% NA,
      name_alt2 = data$name_alt2 %||% NA,
      abbr = data$abbr %||% NA,
      gene_sym = data$gene_sym %||% NA,
      hgnc_gene_sym = data$hgnc_gene_sym %||% NA,
      hgnc_gene_name = data$hgnc_gene_name %||% NA,
      ncbi_gene = data$ncbi_gene %||% NA,
      omim_no = data$omim_no %||% NA,
      inheritance = data$inheritance %||% NA,
      treatability = data$treatability %||% NA,
      prevalence = data$prevalence %||% NA,
      error = FALSE
    )

    return(result)

  }, error = function(e) {
    # Return error marker
    return(tibble(
      disorder_id = id,
      name = NA,
      name_alt1 = NA,
      name_alt2 = NA,
      abbr = NA,
      gene_sym = NA,
      hgnc_gene_sym = NA,
      hgnc_gene_name = NA,
      ncbi_gene = NA,
      omim_no = NA,
      inheritance = NA,
      treatability = NA,
      prevalence = NA,
      error = TRUE
    ))
  })
}

# Main download loop
# Check if we have a saved progress file
if (file.exists(here::here("data", "iembase_disorders.RData"))) {
  load(here::here("data", "iembase_disorders.RData"))
  start_id <- max(iembase_disorders$disorder_id) + 1
  cat("Resuming from disorder ID:", start_id, "\n")
} else {
  iembase_disorders <- tibble()
  start_id <- 1
  cat("Starting fresh download from disorder ID 1\n")
}

# Fetch all IDs up to 3000
max_id <- 3000
current_id <- start_id

while (current_id <= max_id) {
  cat("Fetching disorder ID:", current_id, "\n")

  result <- fetch_disorder(current_id)

  # Add to our dataset
  iembase_disorders <- bind_rows(iembase_disorders, result)

  # Track success/error
  if (result$error) {
    cat("  Error encountered (404)\n")
  } else {
    cat("  Success:", result$name, "\n")
  }

  # Save every 10 records
  if (current_id %% 10 == 0) {
    save(iembase_disorders, file = here::here("data", "iembase_disorders.RData"))
    cat("  Saved progress at disorder ID:", current_id, "\n")
  }

  current_id <- current_id + 1
}

# Final save
save(iembase_disorders, file = here::here("data", "iembase_disorders.RData"))
cat("\nDownload complete! Total disorders:", nrow(iembase_disorders), "\n")
cat("Successful entries:", sum(!iembase_disorders$error), "\n")
cat("Final dataset saved to data/iembase_disorders.RData\n")
