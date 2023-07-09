import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

pd.set_option('display.max_rows', 10)

# Datensatz einlesen
df = pd.read_csv('/Users/tianxinangelama/Documents/Studium/4. FS/DMS/DMS_data/GFP_AEQVI_Sarkisyan_2016.csv')

# split into multiple columns, 1 column for every mutation
df_ind_col = df

df_ind_col[['m1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11', 'm12', 'm13', 'm14', 'm15']] = \
    df_ind_col['mutant'].str.split(':', 14, expand=True)

# reduce df to fscore, m1 - m15
df_only_fscore_mutations = df_ind_col

df_only_fscore_mutations = df_only_fscore_mutations.drop(['mutant', 'mutated_sequence', 'DMS_score_bin'], axis=1)

df_only_fscore_mutations.head()