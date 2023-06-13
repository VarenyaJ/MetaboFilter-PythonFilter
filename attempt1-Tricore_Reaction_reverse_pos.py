'''
Tricore Reaction Reverse steps do?:
    1. Load the input data from a CSV file
    2. Set column names and remove rows with missing values in the 'MS.MS.spectrum' column
    3. Find rows containing 'CUDA' or 'w/o MS2' in the 'Metabolite.name' column
    4. Calculate RSD, fold, and QC_RSD values
    5. Find rows with the maximum fold value
    6. Remove rows with 'None' or 'w/o MS2' in the 'MS.MS.spectrum' column
    7. Add the maximum fold value to the original data using a merge operation?

So basically I need to:
    Reading the CSV file into a pandas DataFrame
    Finding the index of the 'MS.MS.spectrum' column and the end index of the column containing the 'Avg.' value
    Initializing data and removing rows with NaN values in the 'Metabolite.name' column
    Filtering rows where the 'Metabolite.name' is in 'ISTD', 'CUDA', 'iSTD', or 'Unknown'
    Removing rows with 'Unknown' or 'w/o MS2' in the 'Metabolite.name' column
    Initializing data for sorting
    Filtering rows where the 'Metabolite.name' is in 'Unknown' or 'None'
    Calculating RSD, fold, and QC_RSD values
    Replacing missing values in fold and RSD columns with NaN
    Removing rows with NaN values in fold, RSD, and QC_RSD columns
    Merging data_known and data_unknown, keeping rows with unique Metabolite.name values
    Removing rows with 'None' values in the 'MS.MS.spectrum' column

'''

import pandas as pd
from openxlsx import load_workbook

# Load the input data from a CSV file
data = pd.read_csv('Tricore_Reaction_reverse_pos.csv', header=True)

# Set column names and remove rows with missing values in the 'MS.MS.spectrum' column
start = which(data.columns == "MS.MS.spectrum") + 1
start
cols = [col for col in data.columns if col.startswith('Avg')]
end = cols[0] - 1

for i in range(start, data.shape[0]):
    data[cols[i]] = data[cols[i]].fillna(0)

# Remove rows containing 'CUDA' or 'w/o MS2' in the 'Metabolite.name' column
ISTD = data[~data['Metabolite.name'].str.contains('iSTD|CUDA|w/o MS2')]
ISTD = ISTD.set_index('Metabolite.name')

# Remove rows containing 'CUDA' or 'w/o MS2' in the 'Metabolite.name' column, except for the rows with 'iSTD'
CUDA = data[data['Metabolite.name'].str.contains('CUDA') & (~data['Metabolite.name'].str.contains('w/o MS2'))]
CUDA = CUDA.set_index('Metabolite.name')

# Combine ISTD and CUDA DataFrames
ISTD = rbind(ISTD, CUDA)

# Function to choose the maximum value in a subset of rows
def choose_max_of_subset(sub_data):
    mid = len(sub_data) // 2
    new_row = sub_data[mid,]
    df = sub_data.drop(cols, axis=1)
    new_row = apply(df, 1, lambda x: max(x), axis=1)[mid,]
    return new_row

# Find rows with the maximum fold value
data_unknown = data[(data['Average.Mz'] + 0.001) & (data['Average.Rt.min'] - 0.2)]
data_unknown = data_unknown[data_unknown['Average.Mz'] <= (data_unknown['Average.Mz'] - 0.001)]
data_unknown = data_unknown[data_unknown['Average.Rt.min'] <= (data_unknown['Average.Rt.min'] - 0.2)]
data_unknown = data_unknown[data_unknown['Average.Mz'] >= (data_unknown['Average.Mz'] - 0.001)]
data_unknown = data_unknown[data_unknown['Average.Mz'].apply(lambda x: x != 'None')]

# Calculate RSD, fold, and QC_RSD
data['RSD'] = data['Average.Mz'].apply(lambda x: (data[x]['Average.Mz'] - data[x]['Average.Mz'].mean()) / 100 * 100, axis=1)
data['fold'] = data['Average.Mz'].apply(lambda x: data[x].mean() / data[x]['Average.Mz'].mean(), axis=1)
data['QC_RSD'] = data['Average.Rt.min'].apply(lambda x: (data[x]['Average.Rt.min'] - data[x]['Average.Rt.min'].mean()) / 30 * 100, axis=1)

# Find rows with the maximum fold value
data_unknown = data_unknown[(data_unknown['fold'] >= 10) | (data_unknown['QC_RSD'] <= 30)]
data_unknown = data_unknown[(data_unknown['MS.MS.spectrum'] != 'None') | (data_unknown['MS.MS.spectrum'] != 'Unknown')]

# Remove rows with 'None' or 'w/o MS2' in the 'MS.MS.spectrum' column
data_unknown = data_unknown[~data_unknown['MS.MS.spectrum'].str.contains('w/o MS2') | (data_unknown['MS.MS.spectrum'] == 'None')]

# Remove rows with the maximum fold value
data_unknown = data_unknown[data_unknown['fold'] >= 10]

# Add the maximum fold value to the original data
data = data.merge(data_unknown, on='Metabolite.name', how='inner')

################################################################################################################################
#Note-93
################################################################################################################################