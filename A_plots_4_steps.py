import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# dataframe print option: max. 10 rows
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


# 1x mutation
# df mutants with 1x mutation
# copy df
df_1_mutation = df_only_fscore_mutations.copy()

# remain rows with empty values in m2
df_1_mutation = df_1_mutation[df_1_mutation['m2'].isna()]

# list possible mutations in mutants with 1x mutation
# create list with values of m1
possible_mut_1x_mut = df_1_mutation['m1'].tolist()

# delete duplicates
possible_mut_1x_mut = list(set(possible_mut_1x_mut))


# 2x mutations
# df mutants with 2x mutations
# copy df
df_2_mutations = df_only_fscore_mutations.copy()

# delete mutants with 1x mutation
df_2_mutations = df_2_mutations.dropna(subset=['m2'])

# remain rows with empty values in m3
df_2_mutations = df_2_mutations[df_2_mutations['m3'].isna()]

# filter df_2_mutations, only remain mutants with mutations from possible_mut_1x_mut
# copy df
df_2_mutations_filtered = df_2_mutations.copy()

# create a new df for each mut in possible_mut_1x_mut
df_2_mutations_filtered_by = {}
possible_mut_2x_mut_with = {}

for mut in possible_mut_1x_mut:
    # copy df
    df_2_mutations_filtered = df_2_mutations.copy()
    # df remaining rows with mut
    df_2_mutations_filtered = df_2_mutations_filtered[df_2_mutations_filtered.apply(lambda row: mut in row.values, axis=1)]
    # assign to mut
    df_2_mutations_filtered_by[mut] = df_2_mutations_filtered

    # create list with possible mutations of the new df
    possible_mut_2x_mut = list(set(df_2_mutations_filtered.values.flatten()))
    # remove duplicates
    possible_mut_2x_mut = list(set(possible_mut_2x_mut))
    # remove mut if it exists in the list (if, because some mut are not available in df_2_mutations
    if mut in possible_mut_2x_mut:
        possible_mut_2x_mut.remove(mut)
    # remove None
    possible_mut_2x_mut = [value for value in possible_mut_2x_mut if value is not None]
    # remove fscore from list
    possible_mut_2x_mut = [value for value in possible_mut_2x_mut if not isinstance(value, (int, float))]
    # assign to mut
    possible_mut_2x_mut_with[mut] = possible_mut_2x_mut

# merge dfs
# merge all dfs from df_2_mutations_filtered_by
df_2_mutations_filtered_merged = pd.concat(df_2_mutations_filtered_by.values(), ignore_index=False)


# 3x mutations
# df mutants with 3x mutations
# copy df
df_3_mutations = df_only_fscore_mutations.copy()

# delete mutants with < 3 mutations
df_3_mutations = df_3_mutations.dropna(subset=['m3'])

# remain rows with empty values in m4
df_3_mutations = df_3_mutations[df_3_mutations['m4'].isna()]

# filter df_3_mutations with possible_mut_1x_mut
# copy df
df_3_mutations_filtered = df_3_mutations.copy()

# create a new df for each mut in possible_mut_1x_mut
df_3_mutations_filtered_by = {}

for mut in possible_mut_1x_mut:
    # copy df
    df_3_mutations_filtered = df_3_mutations.copy()
    # df remaining rows with mut
    df_3_mutations_filtered = df_3_mutations_filtered[df_3_mutations_filtered.apply(lambda row: mut in row.values, axis=1)]
    # assign to mut
    df_3_mutations_filtered_by[mut] = df_3_mutations_filtered

# filter each df_3_mutations_filtered with respective possible_mut_2x_mut
# create a new df for every mut2 in possible_mut_2x_mut
df_3_mutations_filtered_2_by = {}

# iterate over each df_3_mutations_filtered (each of them was assigned to a mut)
for mut, df_3_mutations_filtered in df_3_mutations_filtered_by.items():
    # iterate over every mut2 in each possible_mut_2x_mut_with[mut]
    for mut2 in possible_mut_2x_mut_with[mut]:
        # filter df_3_mutations_filtered with mut2
        df_3_mutations_filtered_2 = df_3_mutations_filtered_by[mut][df_3_mutations_filtered_by[mut].apply(lambda row: mut2 in row.values, axis=1)]
        # assign every df_3_mutations_filtered_2 to df_3_mutations_filtered_2_by
        key_2 = f"{mut}_{mut2}"
        df_3_mutations_filtered_2_by[key_2] = df_3_mutations_filtered_2

# create a new dictionary df_3_mutations_filtered_2_non_empty_by to store non-empty dfs from df_3_mutations_filtered_2_by
df_3_mutations_filtered_2_non_empty_by = {}

# iterate over mut-mut2 pairs (key_2) of df_3_mutations_filtered_2_by
for key_2, df_3_mutations_filtered_2 in df_3_mutations_filtered_2_by.items():
    if not df_3_mutations_filtered_2.empty:
        df_3_mutations_filtered_2_non_empty_by[key_2] = df_3_mutations_filtered_2

# merge all dfs from df_3_mutations_filtered_2_non_empty_by
df_3_mutations_filtered_2_merged = pd.concat(df_3_mutations_filtered_2_non_empty_by.values(), ignore_index=False)

# create a list for every dataframe of df_3_mutations_filtered_2_non_empty_by, containing the possible mutations, without mut and mut2
# create new empty dictionary
possible_mut_3x_mut_with = {}

# iterate over mut-mut2 pairs (key_2) of df_3_mutations_filtered_2_non_empty_by
for key_2, df_3_mutations_filtered_2 in df_3_mutations_filtered_2_non_empty_by.items():
    # copy df
    df_3_mutations_filtered_2_copy = df_3_mutations_filtered_2.copy()
    # create list with possible mutations of the new df
    possible_mut_3x_mut = list(set(df_3_mutations_filtered_2_copy.values.flatten()))
    # remove duplicates
    possible_mut_3x_mut = list(set(possible_mut_3x_mut))
    # remove mut if it exists in the list (if, because some mut are not available in df_2_mutations)
    if key_2.split('_')[0] in possible_mut_3x_mut:
        possible_mut_3x_mut.remove(key_2.split('_')[0])
    # remove mut2 if it exists in the list (if, because some mut are not available in df_2_mutations)
    if key_2.split('_')[1] in possible_mut_3x_mut:
        possible_mut_3x_mut.remove(key_2.split('_')[1])
    # remove None
    possible_mut_3x_mut = [value for value in possible_mut_3x_mut if value is not None]
    # remove fscore from list
    possible_mut_3x_mut = [value for value in possible_mut_3x_mut if not isinstance(value, (int, float))]
    # assign to key_2
    possible_mut_3x_mut_with[key_2] = possible_mut_3x_mut


# 4x mutations
# df mutants with 4 mutations
# copy df
df_4_mutations = df_only_fscore_mutations.copy()

# delete mutants with < 4 mutations
df_4_mutations = df_4_mutations.dropna(subset=['m4'])

# remain rows with empty values in m5
df_4_mutations = df_4_mutations[df_4_mutations['m5'].isna()]

# filter df_4_mutations with possible_mut_1x_mut
# copy df
df_4_mutations_filtered = df_4_mutations.copy()

# create a new df for each mut in possible_mut_1x_mut
df_4_mutations_filtered_by = {}

for mut in possible_mut_1x_mut:
    # copy df
    df_4_mutations_filtered = df_4_mutations.copy()
    # df remaining rows with mut
    df_4_mutations_filtered = df_4_mutations_filtered[df_4_mutations_filtered.apply(lambda row: mut in row.values, axis=1)]
    # assign to mut
    df_4_mutations_filtered_by[mut] = df_4_mutations_filtered

# filter each df_4_mutations_filtered with respective possible_mut_2x_mut
# create a new df for every mut-mut2 combination, mut2 in possible_mut_2x_mut
df_4_mutations_filtered_2_by = {}

# iterate over each df_4_mutations_filtered (each of them was assigned to a mut from possible_mut_1x_mut)
for mut, df_4_mutations_filtered in df_4_mutations_filtered_by.items():
    # iterate over every mut2 in each possible_mut_2x_mut_with[mut]
    for mut2 in possible_mut_2x_mut_with[mut]:
        # filter df_4_mutations_filtered with mut2
        df_4_mutations_filtered_2 = df_4_mutations_filtered_by[mut][df_4_mutations_filtered_by[mut].apply(lambda row: mut2 in row.values, axis=1)]
        # assign every df_4_mutations_filtered_2 to df_4_mutations_filtered_2_by
        key_2 = f"{mut}_{mut2}"
        df_4_mutations_filtered_2_by[key_2] = df_4_mutations_filtered_2

# create a new dictionary df_4_mutations_filtered_2_non_empty_by to store non-empty dfs from df_4_mutations_filtered_2_by
df_4_mutations_filtered_2_non_empty_by = {}

# iterate over mut-mut2 pairs (key_2) of df_3_mutations_filtered_2_by
for key_2, df_4_mutations_filtered_2 in df_4_mutations_filtered_2_by.items():
    if not df_4_mutations_filtered_2.empty:
        df_4_mutations_filtered_2_non_empty_by[key_2] = df_4_mutations_filtered_2

# filter each df_4_mutations_filtered_2_non_empty_by with respective possible_mut_3x_mut
# create a new df for every mut-mut2-mut3 combination, mut3 in possible_mut_3x_mut
df_4_mutations_filtered_3_by = {}

# iterate over each df_4_mutations_filtered_2_non_empty (each of them was assigned to a mut-mut2 pair)
for key_2, df_4_mutations_filtered_2 in df_4_mutations_filtered_2_non_empty_by.items():
    # check if key_2 is in possible_mut_3x_mut_with (nessesary because some key_2 produced empty df_4_mutations_filtered_2 and therefore key_2 is not in possible_mut_3x_mut_with in these cases)
    if key_2 in possible_mut_3x_mut_with:
        # iterate over every mut3 in each possible_mut_3x_mut_with[key_2]
        for mut3 in possible_mut_3x_mut_with[key_2]:
            # if key_2 is in df_4_mutations_filtered_2_non_empty_by
            if key_2 in df_4_mutations_filtered_2_non_empty_by:
                # filter df_4_mutations_filtered_2 with mut3
                df_4_mutations_filtered_3 = df_4_mutations_filtered_2_non_empty_by[key_2][df_4_mutations_filtered_2_non_empty_by[key_2].apply(lambda row: mut3 in row.values, axis=1)]
                # assign every df_4_mutations_filtered_3 to df_4_mutations_filtered_3_by
                key_3 = f"{key_2}_{mut3}"
                df_4_mutations_filtered_3_by[key_3] = df_4_mutations_filtered_3
            else:
                print(f"Key '{key_2}' does not exist in df_4_mutations_filtered_2_non_empty_by.")
    else:
        print(f"Key '{key_2}' does not exist in possible_mut_3x_mut_with.")

# create a new dictionary df_4_mutations_filtered_3_non_empty_by to store non-empty dfs from df_4_mutations_filtered_3_by
df_4_mutations_filtered_3_non_empty_by = {}

# iterate over mut-mut2-mut3 combinations (key_3) of df_4_mutations_filtered_3_by
for key_3, df_4_mutations_filtered_3 in df_4_mutations_filtered_3_by.items():
    if not df_4_mutations_filtered_3.empty:
        df_4_mutations_filtered_3_non_empty_by[key_3] = df_4_mutations_filtered_3

# merge all df_4_mutations_filtered_3_non_empty_by into one df
df_4_mutations_filtered_3_merged = pd.concat(df_4_mutations_filtered_3_non_empty_by.values())