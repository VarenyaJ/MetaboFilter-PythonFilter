import pandas as pd
import re
from openxlsx import load_workbook

# Read CSV file
data = pd.read_csv('Tricore_Reaction_reverse_pos.csv', header=True)

# Remove rows with missing values in the MS.MS.spectrum column
data[data['MS.MS.spectrum'].isna()] = 0

# Extract columns containing Avg values
cols = [col for col in data.columns if 'Avg' in col]

# Remove rows containing "w/o MS2" or "w/o MS2:" in the Metabolite.name column
data = data[~((data['Metabolite.name'].str.contains('w/o MS2') | data['Metabolite.name'].str.contains('w/o MS2:')))]

# Separate ISTD, CUDA, and other metabolites
ISTD = data[data['Metabolite.name'].str.contains('iSTD')]
CUDA = data[data['Metabolite.name'].str.contains('CUDA')]
ISTD = pd.concat([CUDA, ISTD])

# Remove rows containing "w/o MS2" or "w/o MS2:" in the Metabolite.name column
ISTD = ISTD[~((ISTD['Metabolite.name'].str.contains('w/o MS2') | ISTD['Metabolite.name'].str.contains('w/o MS2:')))]

# Process data to calculate fold change and RSD
data['RSD'] = data.apply(lambda row: row.std() / row.mean() * 100, axis=1)
data['fold'] = data.apply(lambda row: row.max() / row.mean() * 100, axis=1)

# Filter data based on fold change and QC_RSD
sort_data = data[data['fold'] >= 10].sort_values(by=data['QC_RSD'], ascending=False)
sort_data = sort_data[sort_data['QC_RSD'] <= 30]

# Separate data into known and unknown compounds
data_unknown = sort_data[sort_data['Metabolite.name'].str.contains('Unknown|w/o MS2', na=False)]
data_known = sort_data[~sort_data['Metabolite.name'].str.contains('Unknown|w/o MS2', na=False)]
data_known['MS.MS.spectrum'] = data_known['MS.MS.spectrum'].replace('None', '')

# Merge with internal standard data
ISTD = pd.concat([data_unknown, data_known])

# Write processed data to an Excel file
wb = load_workbook('Tricore_Reaction_Reverse_pos_processed_R.xlsx')
ws = wb.add_worksheet('Processed_data')
pd.DataFrame(ISTD).to_excel(ws, index=False)
wb.close()