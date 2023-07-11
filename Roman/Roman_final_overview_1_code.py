# Import the GFP dataset

import numpy as np
import pandas as pd
import csv
with open(r"C:\Users\roman\Desktop\DMS_data\DMS_data\GFP_AEQVI_Sarkisyan_2016.csv") as dms_GFP_datei:
    dms_GFP_datei_object = csv.reader(dms_GFP_datei, delimiter=',')
    print(dms_GFP_datei_object)
    for row in dms_GFP_datei_object:
        print(row)
GFP_dataset = pd.read_csv(r"C:\Users\roman\Desktop\DMS_data\DMS_data\GFP_AEQVI_Sarkisyan_2016.csv")
print(GFP_dataset)

# Clean the Dataset and create a dataframe to work with:

new_AA_of_mutation = []
for index, row in GFP_dataset.iterrows():
    last_character_of_mutation = row["mutant"][-1]
    new_AA_of_mutation.append(last_character_of_mutation)
new_AA_of_mutation_df = pd.DataFrame(new_AA_of_mutation, columns=["New_AS"])
# Last letter in each row

number_mutations_in_mutant = GFP_dataset["mutant"].str.count(":") + 1
number_mutations_singles = number_mutations_in_mutant == 1
# True are all the rows (=mutants) that only have 1 mutation

single_mutants_only_df = new_AA_of_mutation_df[number_mutations_singles]
# This filters all the rows with a "true" value from both dataframes and creates a new dataframe while preserving the filtering.

dms_score_df = []
for index, row in GFP_dataset.iterrows():
    dms_score = row["DMS_score"]
    dms_score_df.append(dms_score)
dms_score_df_alle = pd.DataFrame(dms_score_df, columns=["Fitness_Score"])
# Creates a dataframe with all the dms-scores and the according mutant-number (= experiment number).

dms_score_filtered = dms_score_df_alle[number_mutations_singles]
# Creates a dataframe will the dms-scores of only the single mutants

dms_score_filtered_newAS = single_mutants_only_df.join(dms_score_filtered)
# Combines and creates a new dataframe with the new AA and the corresponding dms-score

import matplotlib.pyplot as plt
dot_size = 3
dms_score_filtered_newAS.plot(x="New_AS", y="Fitness_Score", kind="scatter", s=dot_size)
plt.gcf().set_size_inches(10, 6)
plt.title("Fitness scores of each single mutation based on the new AS")
threshold = 2.5
plt.axhline(threshold, color='red', linestyle='--', label='Threshold')
# Plot shows all the single mutants (new AA) with the corresponding dms-score.
# There is a clear cut off at 2.5 visible. Everything underneath is completely unfunctional.

# Create a new dataframe that is more suitable to work with later down the line.
mutations_pos = []
for index, row in GFP_dataset.iterrows():
    mutations_pos_nummer = row["mutant"][1:-1]
    mutations_pos.append(mutations_pos_nummer)
mutations_pos_df = pd.DataFrame(mutations_pos, columns=["Position"])
# Only removes the first and last character
number_mutations = GFP_dataset["mutant"].str.count(":") + 1
number_mutations_singles = number_mutations == 1
# True are all the rows containing only 1 mutation

single_mutants_only_df_pos = mutations_pos_df[number_mutations_singles]
# Creates a new dataframe that shows the positions within the protein from all the single mutants.

mutations_pos_df_mit_scores = single_mutants_only_df_pos.join(dms_score_filtered)
# Combines and creates a new dataframe showing the position of each single mutation and the respective dms-score (does not show the new amino acid).
new_column = mutations_pos_df_mit_scores["Position"]
Roman_1 = dms_score_filtered_newAS.join(new_column)
Roman_1 = Roman_1[["Position", "New_AS", "Fitness_Score"]]

# Create a heatmap depicting all the available mutations in each position, to gain a better overview.
# However, one heat map would be way too large and incomprehensible, so a subset of 4 heatmaps is created.
# These 4 heatmaps show all the single mutants separated into 4 pieces.

import seaborn as sns
from matplotlib.colors import ListedColormap

mutation_subsets = np.array_split(Roman_1['Position'].unique(), 4)  # Split into 4 subsets

# Creates some subplots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Creates a custom colormap
custom_colors = ["lightgray", "blue"]
cmap = sns.color_palette(custom_colors)

# Creates a shared colorbar with two ticks for the two possibilities yes/no
cbar_ticks = [0.25, 0.75]
cbar_ticklabels = ['Mutation does not exist', 'Mutation exists']


# Defines a custom sorting function for the y-axis values
def position_sort(position):
    if position.isdigit():
        return int(position)
    else:
        return position

# Iterates over subsets and creates 4 heatmaps
for i, subset in enumerate(mutation_subsets):
    row = i // 2
    col = i % 2

    # Filters the data based on each subset
    subset_df = Roman_1[Roman_1['Position'].isin(subset)]

    # Sorts the data in ascending order
    subset_df = subset_df.sort_values('Position')

    # Creates cross-tabulation
    cross_tab = pd.crosstab(subset_df['Position'], subset_df['New_AS'])

    # Sorts the columns alphabetically for better visualization
    cross_tab_sorted = cross_tab.reindex(sorted(cross_tab.columns), axis=1)

    # Sorts the rows based on the custom y-axis sorting function
    cross_tab_sorted = cross_tab_sorted.iloc[sorted(range(len(cross_tab_sorted)), key=lambda x: position_sort(cross_tab_sorted.index[x]))]

    # Creates each heatmap in the corresponding subplot
    ax = sns.heatmap(cross_tab_sorted, cmap=cmap, annot=False, fmt='d', ax=axes[row, col],
                     cbar_kws={'ticks': cbar_ticks, 'drawedges': True}, linewidths=0.5)

    axes[row, col].set_xlabel('New Amino Acid')
    axes[row, col].set_ylabel('Mutation Position')
    axes[row, col].set_title(f'Subset {i+1}')

    # Sets colorbar tick labels
    cbar = ax.collections[0].colorbar
    cbar.set_ticklabels(cbar_ticklabels)

plt.tight_layout()

# A heatmap to see which single mutations are available in the dataset and which are not.
# Only focuses on the new amino acid that occurs, not the old one before the mutation.

# This plot is created to show a global overview of all the positions with their respective dms-scores
# Plot shows all the single mutations that occured at a specific position (x-axis) against the dms-score (y-axis)

mutations_pos_df_mit_scores.plot(x="Position", y="Fitness_Score", kind="scatter")
plt.title("Fitness scores of each single mutation based on the position")

plt.xticks(rotation='vertical')
# Rotates the x-axis labels vertically
plt.gcf().set_size_inches(60, 25)
# Changes the width and height in inches

a_threshold = 2.5
plt.axhline(threshold, color='red', linestyle='--', label='Threshold')
a_threshold_65 = 62
plt.axvline(a_threshold_65, color='red', linestyle='--', label='Threshold_65')
a_threshold_67 = 64
plt.axvline(a_threshold_67, color='red', linestyle='--', label='Threshold_67')
# Inserts some threshold-lines to indicate important positions (chromophor).
# It's interesting to see here, that there are mutations occuring in the chromophore (positions 65-67) that still have a good dms-score.
# Position 66 is the fluorescing amino acid, and position 67 also participates in the cyclization reaction, whereas position 65 has a smaller impact on the fluorescence.
# Therefore, it makes sense that position 66 and 67 show lowered dms-scores, whilst position 65 is more stable.


