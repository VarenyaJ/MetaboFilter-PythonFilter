# Load required libraries
library(dplyr)
library(stringr)
library("openxlsx")

# Read the input data from a CSV file
data <- read.csv('Tricore_Reaction_reverse_pos.csv', header = TRUE)
ori_data <- data

# Get column names
colnames(data)

# Find the start index of columns
start <- which(colnames(data) == "MS.MS.spectrum") + 1

# Get the indices of columns containing 'Avg'
cols <- grep('Avg', colnames(data))

# Find the end index
end <- cols[1] - 1

# Replace NA values with 0 for specific columns
for (i in colnames(data)[start:dim(data)[2]]) {
  data[is.na(data[i]), i] <- 0
}

# Filter ISTD and CUDA data
ISTD <- data[grepl("iSTD", data$Metabolite.name), (invert = TRUE), ]
ISTD <- ISTD[!grepl("w/o MS2:", ISTD$Metabolite.name), (invert = TRUE), ]
CUDA <- data[grepl("CUDA", data$Metabolite.name), (invert = TRUE), ]
CUDA <- CUDA[!grepl("w/o MS2", CUDA$Metabolite.name), (invert = TRUE), ]

# Combine ISTD and CUDA data
ISTD <- rbind(CUDA, ISTD)
ls_ISTD <- ISTD$Metabolite.name

# Filter data based on ls_ISTD
data <- subset(data, !(Metabolite.name %in% ls_ISTD))

# Create a function to choose maximum values from a subset
choose_max_of_subset <- function(sub_data) {
  mid = dim(sub_data)[1] %/% 2
  new_row = sub_data[mid,]
  df = sub_data[, start:dim(sub_data)[2]]
  print(dim(sub_data))
  new_row[, start:dim(sub_data)[2]] = apply(df, 2, function(x) max(x))
  return(new_row)
}

# Sort the data
sort_dt <- data[order(data$Average.Mz, data$Average.Rt.min.),]
new_dt <- data.frame()
df_2 <- sort_dt

# Compare and choose maximum values from subsets
i = 1
while (i <= dim(df_2)[1]) {
  a = df_2$Average.Mz[i] + 0.001
  b = df_2$Average.Rt.min.[i] + 0.2
  c = df_2$Average.Rt.min.[i] - 0.2
  d = df_2$Average.Mz[i] - 0.001
  bel = df_2[i:dim(df_2)[1], ]
  sub_set = bel[bel$Average.Mz <= a & bel$Average.Rt.min. >= c & bel$Average.Rt.min. <= b & bel$Average.Mz >= d,]

  if (dim(sub_set)[1] > 1) {
    print(df_2[i, 1:5])
    new_row = choose_max_of_subset(sub_set)
    df_2 <- df_2[setdiff(rownames(df_2), rownames(sub_set)),]
  } else {
    new_row = df_2[i,]
    i = i + 1
  }
  new_dt = rbind(new_dt, new_row)
}

# Combine ISTD and new_dt
data <- rbind(ISTD, new_dt)
#write.csv(data, "data_reaction_pos.csv", row.names = FALSE)
#data<- read.csv('data_neg.csv', header = TRUE )

# Define sample columns
samples <- names(data)[start:end]
samples
sample_col <- samples[!grepl("QC", samples)]
sample_col <- sample_col[!grepl("blank", sample_col)]

sample_dt <- subset(data, select = sample_col)
#cols <- ncol(sample_dt)

QC_cols <- samples[grepl("QC", samples)]
QC_cols <- QC_cols[!grepl("EQC", QC_cols)]
dt_QC <- subset(data, select = QC_cols)

data$RSD <- apply(sample_dt, 1, sd) / rowMeans(sample_dt) * 100
data$fold <- apply(sample_dt, 1, max) / rowMeans(data[, samples[grepl("blank", samples)]])
data$QC_RSD <- apply(dt_QC, 1, sd) / rowMeans(dt_QC) * 100

data <- data %>% relocate(c('fold', 'RSD', 'QC_RSD'), .after = Metabolite.name)

ISTD <- subset(data, (Metabolite.name %in% ls_ISTD))
ori_ISTD <- ISTD
dt_noISTD <- subset(data, !(Metabolite.name %in% ls_ISTD))

sort_data <- dt_noISTD[order(dt_noISTD$fold), ]
sort_data <- subset(sort_data, fold >= 10)
sort_data <- sort_data[order(sort_data$QC_RSD), ]
sort_data <- subset(sort_data, QC_RSD <= 30)

data_unknown <- sort_data[grepl("Unknown|w/o MS2", sort_data$Metabolite.name), (invert = TRUE), ]
data_unknown$Metabolite.name <- 'Unknown'
data_know <- sort_data[!grepl("Unknown|w/o MS2", sort_data$Metabolite.name), (invert = TRUE), ]
data_know["MS.MS.spectrum"][data_know["MS.MS.spectrum"] == ''] <- 'None'
data_unknown <- rbind(data_unknown, data_know[grepl("None", data_know$MS.MS.spectrum), (invert = TRUE), ])
data_know <- data_know[!grepl("None", data_know$MS.MS.spectrum), (invert = TRUE), ]

ISTD <- ISTD[, c("Alignment.ID", "Average.Rt.min.", "Average.Mz", "Metabolite.name", "Adduct.type", samples)]
ISTD$mean <- rowMeans(ISTD[, 6:ncol(ISTD)])
for (i in seq(1, dim(ISTD)[1])) {

  ISTD[i, 6:ncol(ISTD)] = ISTD[i, 6:ncol(ISTD)] / ISTD[i, 'mean']

}
ISTD$mean <- NULL
ISTD_ave <- ISTD[, 6:(ncol(ISTD))]
ave_IS <- colMeans(ISTD_ave)
ISTD <- bind_rows(ISTD, ave_IS)
ISTD[nrow(ISTD), 5] <- "Average"

###
data_know[c('Metabolite_name', 'InChiKey')] <- str_split_fixed(data_know$Metabolite.name, '_', 2)
data_know <- data_know %>% relocate(c('Metabolite_name', 'InChiKey'), .after = Metabolite.name)
data_know$Metabolite.name <- NULL
names(data_know)[names(data_know) == "Metabolite_name"] <- "Metabolite.name"
data_know <- data_know[order(data_know$InChiKey, decreasing = TRUE), ]

data_know["InChiKey"][data_know["InChiKey"] == ''] <- 'None'
data_inchikey <- data_know[!grepl("None", data_know$InChiKey), (invert = TRUE), ]
data_woinchikey <- data_know[grepl("None", data_know$InChiKey), (invert = TRUE), ]

lb <- read.xlsx('D:/tricore/Chemical Identification.xlsx', sheet = 'Sheet1')
names(lb)[names(lb) == "InChIKey"] <- "InChiKey"

lb_1 <- lb[, c('Metabolite.name', 'InChiKey')]
data_woinchikey$InChiKey <- NULL
head(lb_1)

lb_1$Metabolite.name <- gsub("[.]", "", lb_1$Metabolite.name)
data_woinchikey$Metabolite.name <- gsub("[.]", "", data_woinchikey$Metabolite.name)
lb_1$Metabolite.name <- str_trim(lb_1$Metabolite.name)
data_woinchikey$Metabolite.name <- str_trim(data_woinchikey$Metabolite.name)
lb_1$Metabolite.name <- toupper(lb_1$Metabolite.name)
data_woinchikey$Metabolite.name <- toupper(data_woinchikey$Metabolite.name)
lb_1 <- lb_1[!duplicated(lb_1$Metabolite.name), ]

data_woinchikey <- merge(data_woinchikey, lb_1, by = "Metabolite.name", all.x = TRUE)
data_woinchikey <- data_woinchikey %>% relocate("Metabolite.name", .after = "Average.Mz")
data_woinchikey <- data_woinchikey %>% relocate("InChiKey", .after = "Metabolite.name")

#data_woinchikey["InChiKey"][is.na(data_woinchikey["InChiKey"])] <- "None"
#data_getinchikey <- data_woinchikey[!grepl("None", data_woinchikey$InChiKey), (invert = TRUE), ]
#data_noinchikey <- data_woinchikey[grepl("None", data_woinchikey$InChiKey), (invert = TRUE), ]
#write.csv(data_noinchikey, "df_noinchikey_neg.csv", row.names = FALSE)

data_know_f <- rbind(data_inchikey, data_woinchikey)

final_data <- bind_rows(ori_ISTD, data_know_f, data_unknown)

dataset_names <- list('Unprocess_data' = ori_data, 'Processed_data' = final_data, 'Internal_STD' = ISTD, 'Know_compound_RT_match' = data_know_f, 'Unknown_compound' = data_unknown)
openxlsx::write.xlsx(dataset_names, file = 'Tricore_Reaction_Reverse_pos_processed_R.xlsx')
