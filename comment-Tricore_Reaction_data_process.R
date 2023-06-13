# Load required libraries
library('openxlsx')
library(stringr)
library(dplyr)

# Read data from Excel files
df_neg <- read.xlsx('Tricore_Reaction_Reverse_neg_processed_R.xlsx', sheet = 'Know_compound_RT_match')
df_pos <- read.xlsx('Tricore_Reaction_Reverse_pos_processed_R.xlsx', sheet = 'Know_compound_RT_match')

# Read data from CSV file and assign column names
df1 <- read.csv('D:/Hilic_IROA_IS_negative_lib.csv')
colnames(df1) <- c('Metabolite_name', 'Average_Mz', 'Average_Rt_min', 'InChiKey')

# Update column names in df_neg
colnames(df_neg)[colnames(df_neg) %in% c("Average.Rt.min.", "Average.Mz","Metabolite.name")] <- c("Average_Rt_min", "Average_Mz","Metabolite_name")

# Convert column values to appropriate data types
df1['Average_Rt_min'] <- lapply(df1['Average_Rt_min'], as.numeric)
sapply(df1, class)

# Perform data transformations and matching
df_neg$Metabolite_name <- gsub("[.]", "", df_neg$Metabolite_name)
df1$Metabolite_name <- gsub("[.]", "", df1$Metabolite_name)
df_neg$Metabolite_name <- toupper(df_neg$Metabolite_name)
df1$Metabolite_name <- iconv(df1$Metabolite_name, from = 'ISO-8859-1', to = 'utf8')
df1$Metabolite_name <- toupper(df1$Metabolite_name)
df_neg$Metabolite_name <- str_trim(df_neg$Metabolite_name)
df1$Metabolite_name <- str_trim(df1$Metabolite_name)

# Perform matching based on specific conditions
for (i in 1:nrow(df_neg)) {
  for (j in 1:nrow(df1)) {
    if ((df_neg[i,'InChiKey'] == df1[j,'InChiKey']) &&
        (df_neg[i,'Average_Mz'] >= df1[j,'Average_Mz'] - 0.005) &&
        (df_neg[i,'Average_Mz'] <= df1[j,'Average_Mz'] + 0.005) &&
        (df_neg[i,'Average_Rt_min'] >= df1[j,'Average_Rt_min'] - 0.3) &&
        (df_neg[i,'Average_Rt_min'] <= df1[j,'Average_Rt_min'] + 0.3)) {
      
      df_neg[i,'MSI_level'] <- "1"
      df_neg[i,'RT_inlib'] <- df1[j,'Average_Rt_min']
    }
  }
}

# Update column names in df_neg
colnames(df_neg)[colnames(df_neg) %in% c("Average_Rt_min", "Average_Mz","Metabolite_name")] <-c("Average.Rt.min.", "Average.Mz","Metabolite.name")

# Read data from CSV file into df2 and assign column names
df2 <- read.csv('D:/Reverse_IROA_IS_positive_lib.csv')
colnames(df2) <- c('Metabolite_name', 'Average_Mz', 'Average_Rt_min', 'InChiKey')
colnames(df2)

# Set MSI_level and RT_inlib columns in df_pos
df_pos$MSI_level <- ""
df_pos$RT_inlib <-""
df_pos <- df_pos %>% relocate(c('MSI_level','RT_inlib'), .after = Metabolite.name)

# Update column names in df_pos
colnames(df_pos)[colnames(df_pos) %in% c("Average.Rt.min.", "Average_Mz","Metabolite.name")] <- c("Average_Rt_min", "Average_Mz","Metabolite_name")
colnames(df_pos)

# Convert column values in df2 to appropriate data types
df2['Average_Rt_min'] <- lapply(df2['Average_Rt_min'], as.numeric)
sapply(df2, class)

# Perform string transformations on Metabolite_name columns in df_pos and df2
df_pos$Metabolite_name <- gsub("[.]", "", df_pos$Metabolite_name)
df2$Metabolite_name <- gsub("[.]", "", df2$Metabolite_name)
df_pos$Metabolite_name <- toupper(df_pos$Metabolite_name)
df2$Metabolite_name <- iconv(df2$Metabolite_name, from = 'ISO-8859-1', to = 'utf8')
df2$Metabolite_name <- toupper(df2$Metabolite_name)
df_pos$Metabolite_name <- str_trim(df_pos$Metabolite_name)
df2$Metabolite_name <- str_trim(df2$Metabolite_name)

# Perform matching between df_pos and df2 based on specific conditions
for (i in 1:nrow(df_pos)){
  for (j in 1:nrow(df2)){
    if ((df_pos[i,'InChiKey'] == df2[j,'InChiKey']) &&
        (df_pos[i,'Average_Mz'] >= df2[j,'Average_Mz'] - 0.005) &&
        (df_pos[i,'Average_Mz'] <= df2[j,'Average_Mz'] + 0.005) &&
        (df_pos[i,'Average_Rt_min'] >= df2[j,'Average_Rt_min'] - 0.3) &&
        (df_pos[i,'Average_Rt_min'] <= df2[j,'Average_Rt_min'] + 0.3)){
      df_pos[i,'MSI_level'] <- "1"
      df_pos[i,'RT_inlib'] <- df2[j,'Average_Rt_min']
    }
  }   
}

# Update column names in df_pos
colnames(df_pos)[colnames(df_pos) %in% c("Average_Rt_min", "Average_Mz","Metabolite_name")] <-c("Average.Rt.min.", "Average.Mz","Metabolite.name")
colnames(df_pos) <- colnames(df_neg)

# Combine df_neg and df_pos into a single dataframe, df
df <- rbind(df_neg, df_pos)

# Sort df by MSI_level in descending order
df <- df[order(df$MSI_level, decreasing = TRUE),]

# Create MSI_1 dataframe containing rows with MSI_level = '1'
MSI_1 <- subset(df, MSI_level == '1')

# Create df_1 dataframe containing rows with MSI_level not equal to '1'
df_1 <- subset(df, MSI_level != '1')

# Convert Total.score column in df_1 to numeric
df_1['Total.score'] <- lapply(df_1['Total.score'], as.numeric)

# Update MSI_level values in df_1 based on Total.score condition
for (i in 1:nrow(df_1)) {
  if (df_1[i, "Total.score"] >= 95) {
    df_1[i, "MSI_level"] <- "2"
  }
}

# Sort df_1 by MSI_level in descending order
df_1 <- df_1[order(df_1$MSI_level, decreasing = TRUE),]

# Create MSI_2 dataframe containing rows with MSI_level = '2'
MSI_2 <- subset(df_1, MSI_level == '2')

# Create df_2 dataframe containing rows with MSI_level not equal to '2'
df_2 <- subset(df_1, MSI_level != '2')

# Update MSI_level values in df_2 based on Total.score condition
for (i in 1:nrow(df_2)) {
  if (df_2[i, "Total.score"] >= 90) {
    df_2[i, "MSI_level"] <- "3"
  }
}

# Update MSI_level values in df_2 not equal to '3' to '4'
df_2["MSI_level"][df_2["MSI_level"] != "3"] <- '4'

# Combine MSI_1, MSI_2, and df_2 into dt_all dataframe
dt_all <- rbind(MSI_1, MSI_2, df_2)

# Remove trailing spaces in InChiKey column
dt_all$InChiKey <- str_trim(dt_all$InChiKey)

# Separate dt_all into dt_undefine and dt_all based on InChiKey containing "undefined"
dt_undefine <- dt_all[grepl("undefined", dt_all$InChiKey), (invert = TRUE), ]
dt_all <- dt_all[!grepl("undefined", dt_all$InChiKey), (invert = TRUE), ]

# Sort dt_all by InChiKey, MSI_level, and QC_RSD
dt_all <- dt_all[order(dt_all$InChiKey, dt_all$MSI_level, dt_all$QC_RSD),]

# Write dt_all to a CSV file named "remove_rep.csv"
# write.csv(dt_all, "remove_rep.csv", row.names = FALSE)

# Remove duplicate rows based on InChiKey in dt_all and create final dataframe
final <- dt_all[!duplicated(dt_all$InChiKey), ] 

# Perform string transformations and sorting on dt_undefine
dt_undefine$Metabolite.name <- toupper(dt_undefine$Metabolite.name)
dt_undefine$Metabolite.name <- gsub("[.]", "", dt_undefine$Metabolite.name)
dt_undefine$Metabolite.name <- str_trim(dt_undefine$Metabolite.name)
dt_undefine <- dt_undefine[order(dt_undefine$Metabolite.name, dt_undefine$MSI_level, dt_undefine$QC_RSD),]
dt_undefine <- dt_undefine[!duplicated(dt_undefine$Metabolite.name),]

# Define the column range for samples in dt_all
start <- which(colnames(final) == "MS.MS.spectrum") + 1
cols <- grep('Avg', colnames(final))
end <- cols[1] - 1
samples <- names(final)[start:end]
sample_col <- samples[!grepl("QC", samples)]

# Define the columns to select for dt_undefine and data
cols <- c("Metabolite.name", "Average.Rt.min.", "Average.Mz", "Adduct.type", "MSI_level", "InChiKey", "QC_RSD", sample_col)
colnames(dt_undefine)
dt_undefine <- subset(dt_undefine, select = cols)
data <- subset(final, select = cols)

# Read "Chemical Identification.xlsx" file into lb dataframe
lb <- read.xlsx('D:/tricore/Chemical Identification.xlsx', sheet = 'Sheet1')
colnames(lb)

# Select specific columns from lb dataframe
lb_1 <- subset(lb, select = c("InChiKey", "SMILES", "PubChem.CID", "KEGG.ID"))

# Remove duplicate rows based on InChiKey in lb_1
lb_1 <- lb_1[order(lb_1$InChiKey),]
lb_1$InChiKey <- str_trim(lb_1$InChiKey)
data$InChiKey <- str_trim(data$InChiKey)
lb_1 <- lb_1[!duplicated(lb_1$InChiKey), ] 

# Merge data and lb_1 dataframes based on InChiKey, keeping all rows from data
dt_final <- merge(data, lb_1, by = "InChiKey", all.x = TRUE)

# Remove QC_RSD column from dt_final
dt_final$QC_RSD <- NULL

# Select specific columns from lb dataframe
lb_2 <- subset(lb, select = c("Metabolite.name", "SMILES", "PubChem.CID", "KEGG.ID"))

# Perform string transformations and sorting on lb_2 and dt_undefine
lb_2$Metabolite.name <- toupper(lb_2$Metabolite.name)
lb_2$Metabolite.name <- str_trim(lb_2$Metabolite.name)
dt_undefine$Metabolite.name <- str_trim(dt_undefine$Metabolite.name)
lb_2 <- lb_2[!duplicated(lb_2$dt_undefine$Metabolite.name), ] 

# Merge dt_undefine and lb_2 dataframes based on Metabolite.name, keeping all rows from dt_undefine
dt_undefine <- merge(dt_undefine, lb_2, by = "Metabolite.name", all.x = TRUE)

# Remove QC_RSD column from dt_undefine
dt_undefine$QC_RSD <- NULL

# Combine dt_final and dt_undefine into dt_final
dt_final <- rbind(dt_final, dt_undefine)

# Display column names of dt_final
colnames(dt_final)

# Sort dt_final by MSI_level
dt_final <- dt_final[order(dt_final$MSI_level), ]

# Update MSI_level values in dt_final
dt_final['MSI_level'][dt_final["MSI_level"] == "2"] <- '2A'
dt_final['MSI_level'][dt_final["MSI_level"] == "3"] <- '2B'
dt_final['MSI_level'][dt_final["MSI_level"] == "4"] <- '2C'

# Define the column format for dt_final
col_format <- c("Metabolite.name", "Average.Rt.min.", "Average.Mz", "Adduct.type", "MSI_level", "InChiKey", "SMILES", "PubChem.CID", "KEGG.ID", sample_col)
dt_final <- dt_final[, col_format]

# Sort dt_final by MSI_level and InChiKey
dt_final <- dt_final[order(dt_final$MSI_level, dt_final$InChiKey), ]


# Read "Tricore_Reaction_Reverse_neg_processed_R.xlsx" file into unknown1 dataframe
unknown1 <- read.xlsx('Tricore_Reaction_Reverse_neg_processed_R.xlsx', sheet = 'Unknown_compound')

# Read "Tricore_Reaction_Reverse_pos_processed_R.xlsx" file into unknown2 dataframe
unknown2 <- read.xlsx('Tricore_Reaction_Reverse_pos_processed_R.xlsx', sheet = 'Unknown_compound')

# Set column names of unknown2 to match unknown1
colnames(unknown2) <- colnames(unknown1)

# Combine unknown1 and unknown2 into unknow dataframe
unknow <- rbind(unknown1, unknown2)

# Select specific columns from unknow dataframe
unknow <- unknow[, c("Metabolite.name", "Average.Rt.min.", "Average.Mz", "Adduct.type", sample_col)]

# Set the value of "Metabolite.name" column in unknow dataframe to 'Unknown'
unknow$Metabolite.name <- 'Unknown'

# Create a list of datasets with names
dataset_names <- list('Known_compound' = dt_final, 'Unknown_compound' = unknow)

# Write dataset_names to "Tricore_Reaction_combine_data_process.xlsx" file using openxlsx package
openxlsx::write.xlsx(dataset_names, file = 'Tricore_Reaction_combine_data_process.xlsx')
