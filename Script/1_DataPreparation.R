# 
# Data preprocessing
# 
# Author: Yusuke Toda (yusuke0320.research@gmail.com)
# 

library(openxlsx)

# Create a directory for storing the output files
path_out <- file.path('Result', '1_DataPreparation')
dir.create(path_out, recursive = TRUE, showWarnings = FALSE)



#### 1. Import data ####

path_dat <- file.path('Data', 'CG_database.xlsx')
rawdat <- list()

# Read the column names
rawdat$head <- read.xlsx(path_dat, rows = 1, colNames = FALSE)

# Read the main dataset and assign the column names
rawdat$all <- read.xlsx(path_dat, startRow = 3, colNames = FALSE, na.strings = 'na')
colnames(rawdat$all) <- unlist(rawdat$head)
# View(rawdat$all)



#### 2. Clean data values ####

# Rename columns
dat <- list(df = rawdat$all)
colnames(dat$df)[colnames(dat$df) == 'Pref.']                <- 'Pref'
colnames(dat$df)[colnames(dat$df) == 'Institution/Station']  <- 'InstOrStation'
colnames(dat$df)[colnames(dat$df) == 'Soil Groups ']         <- 'SoilGroups'
colnames(dat$df)[colnames(dat$df) == 'Cultivar name (E)']    <- 'Cultivar'
colnames(dat$df)[colnames(dat$df) == 'DOY']                  <- 'HeadingDOY'
colnames(dat$df)[colnames(dat$df) == 'Head rice']            <- 'HeadRice'
colnames(dat$df)[colnames(dat$df) == 'Milky-white']          <- 'MilkyWhite'
colnames(dat$df)[colnames(dat$df) == 'Basal-white']          <- 'BasalWhite'
colnames(dat$df)[colnames(dat$df) == 'Belly white']          <- 'BellyWhite'
colnames(dat$df)[colnames(dat$df) == 'sieve']                <- 'Sieve'
colnames(dat$df)[colnames(dat$df) == 'Measurement methods '] <- 'MeasurementMethods'
colnames(dat$df)[colnames(dat$df) == 'URL/doi']              <- 'URLorDOI'
colnames(dat$df)[colnames(dat$df) == 'Ref_No.']              <- 'RefNo'

# Extract column metadata
dat$col_info <- data.frame(name = colnames(dat$df))
dat$col_info <- within(dat$col_info, {
  
  # Identify columns containing measured values
  is_val <- seq_along(name) %in% 12:25
  
  # Classify measured values by type
  val_type <- ifelse(
    grepl('^(Tave)|(TaHD..)|(SR)|(RH)$', name), 'Weather', 'Phenotype')
  val_type[!is_val] <- NA
  
})

# Remove trailing spaces from values in the data frame
remove_tail_space <- function(dat) {
  
  for (cn in colnames(dat)) {
    val <- dat[, cn]
    val_with_space <- unique(val[grepl(' $', val)])
    if (length(val_with_space) >= 1) {
      cat(sprintf('Removing unnecessary space in "%s".\n', cn))
      for (txt in val_with_space) {
        cat(sprintf('  "%s"\n', txt))
      }
      dat[, cn] <- gsub(' $', '', val)
    }
  }
  dat
  
}
dat$df <- remove_tail_space(dat$df)

# Detect and remove values outside the expected range
detect_strange_value <- function(df, colinfo) {
  
  # Extract phenotypic variables expressed as percentages
  is_pheno <- rep(FALSE, nrow(colinfo))
  is_pheno[colinfo$val_type == 'Phenotype'] <- TRUE
  is_pheno[colinfo$name == 'HeadingDOY'] <- FALSE
  pheno_df <- df[, is_pheno]
  
  # Identify percentage values outside the range 1–99
  is_strange <- pheno_df >= 100 | pheno_df <= 0
  is_strange[is.na(is_strange)] <- FALSE
  
  # Print the detected values and their locations
  cat('value\n')
  print(pheno_df[rowSums(is_strange) >= 1, ])
  cat('\nis_strange\n')
  print(is_strange[rowSums(is_strange) >= 1, ])
  
  pheno_df[is_strange] <- NA
  df[, is_pheno] <- pheno_df
  df
  
}
dat_clean <- dat
dat_clean$df <- detect_strange_value(dat_clean$df, dat_clean$col_info)

# Correct a cultivar name
dat_clean$df$Cultivar <- gsub('Mutsuhomarre', 'Mutsuhomare', dat_clean$df$Cultivar)



#### 3. Organize data and metadata ####

# Separate measured values from row-level metadata
dat_clean <- within(dat_clean, {
  mat <- as.matrix(df[, col_info$is_val])
  row_info <- df[, !col_info$is_val]
})

# Extract reference and site information
ref_site_info <- dat_clean$row_info[, c('Pref', 'Reference', 'RefNo')]
ref_site_info <- unique(ref_site_info)
ref_site_info <- ref_site_info[order(ref_site_info$RefNo), ]

# Create reference-site identifiers for subsequent analyses
create_ref_site_id <- function(ref_site_info) {
  
  # Create a table of unique references
  ref_info <- unique(ref_site_info[, c('Reference', 'RefNo')])
  
  # Create abbreviated reference names
  ref_info$RefAbb <- ref_info$Reference
  ref_info$RefAbb <- gsub('( et? al\\.?)?, ', '', ref_info$RefAbb)
  ref_info$RefAbb <- gsub('^RCC database- "', 'RCC-', ref_info$RefAbb)
  ref_info$RefAbb <- gsub('"-?$', '', ref_info$RefAbb)
  ref_info$RefAbb <- gsub('The Performance Tests for Recommended Varieties of Rice database', 
                          'RcommendedVarietyTest', ref_info$RefAbb)
  
  # Append lowercase letters to duplicated reference abbreviations
  ref_abb_dup <- unique(ref_info$RefAbb[duplicated(ref_info$RefAbb)])
  for (refabb0 in ref_abb_dup) {
    r <- which(ref_info$RefAbb == refabb0)
    ref_info$RefAbb[r] <- paste0(ref_info$RefAbb[r], letters[seq_along(r)])
  }
  
  # Add the reference abbreviations to the reference-site table
  ref_site_info <- merge(ref_site_info, ref_info, all = TRUE, sort = FALSE)
  
  # Create identifiers by combining reference abbreviations and site information
  ref_site_info$RefSiteID <- with(ref_site_info, paste0(RefAbb, '_', Pref))
  ref_site_info$RefAbb <- NULL
  return(ref_site_info)
  
}
ref_site_info <- create_ref_site_id(ref_site_info)

# Add the generated reference-site identifiers to the dataset
add_column_of_ref_site_id <- function(ref_site_info, info) {
  # Add reference-site identifiers to the data frame `info`.
  # `ref_site_info` is a table of unique reference-site combinations.
  
  # Stop if duplicated identifiers are found in the reference-site table
  if (any(duplicated(ref_site_info$RefSiteID))) {
    stop('Given RefSiteID is duplicated.')
  }
  
  # Separate the identifiers from the remaining reference-site information
  ref_site_id <- ref_site_info$RefSiteID
  ref_site_info <- ref_site_info[, colnames(ref_site_info) != 'RefSiteID']
  
  # Extract the columns required for matching
  info0 <- info[, colnames(ref_site_info)]
  
  # Assign reference-site identifiers
  id <- vector('numeric', nrow(info))
  colsole_width <- getOption('width')
  cat('Searching reference-site ID...\n')
  for (i in seq_along(id)) {
    
    info0_now <- info0[i, ]
    for (j in seq_along(ref_site_id)) {
      
      # Exit the loop when matching reference-site information is found
      if (all(info0_now == ref_site_info[j, ])) break
      
      # Stop if no corresponding entry is found in the reference-site table
      if (j == length(ref_site_id)) {
        cat('No matching information was found for:\n')
        cat(paste(rep('-', colsole_width), collapse = ''))
        cat('\n')
        print(info0_now)
        cat(paste(rep('-', colsole_width), collapse = ''))
        cat('\n')
        stop()
      }
      
    }
    id[i] <- ref_site_id[j]
    
  }
  info$RefSiteID <- id
  cat('Finished.\n')
  return(info)
  
}
dat_clean$row_info <- add_column_of_ref_site_id(ref_site_info, dat_clean$row_info)

# Assign sequential numbers to the reference-site combinations
ref_site_info$RefSiteNo <- seq_len(nrow(ref_site_info))
dat_clean$row_info$RefSiteNo <- with(
  ref_site_info, RefSiteNo[match(dat_clean$row_info$RefSiteID, RefSiteID)])

# Save the processed data and reference-site information
saveRDS(dat_clean, file.path(path_out, 'Data.rds'))
saveRDS(ref_site_info, file.path(path_out, 'Info_Reference-Site.rds'))



#### 4. Inspect the processed values ####

n_refsite <- length(unique(dat_clean$row_info$RefSiteID))

pdf(file.path(path_out, 'Fig_boxplot.pdf'), height = 7, width = 9)
for (i in seq_len(ncol(dat_clean$mat))) {
  
  par(mar = c(8, 4, 2, 2))
  
  plot(NULL, 
       xlim = c(1, n_refsite),
       ylim = range(dat_clean$mat[, i], na.rm = TRUE),
       xaxt = 'n', xlab = '', ylab = colnames(dat_clean$mat)[i])
  
  abline(v = seq_len(n_refsite), col = rgb(0, 0, 0, 0.2))
  
  b <- boxplot(dat_clean$mat[, i] ~ dat_clean$row_info$RefSiteID, 
               xaxt = 'n', add = TRUE)
  
  text(x = seq_along(b$names), 
       y = sum(par('usr')[3:4] * c(1.01, -0.01)),
       labels = b$names, cex = 0.6,
       adj = c(1, 0.5), srt = 90, xpd = TRUE)
  
}
dev.off()
