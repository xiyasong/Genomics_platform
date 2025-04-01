import pandas as pd
from scipy.stats import spearmanr

# Load the first Excel file
df1 = pd.read_excel('file1.xlsx')

# Load the second Excel file
df2 = pd.read_excel('file2.xlsx')

# Merge the two dataframes based on the common column
merged_df = pd.merge(df1, df2, on='Patient ID')

# Calculate the Spearman's rank correlation coefficient between the two rankings
corr, _ = spearmanr(merged_df['Ranking_x'], merged_df['Ranking_y'])


print('The similarity score (Spearman\'s rank correlation coefficient) is:', corr)
