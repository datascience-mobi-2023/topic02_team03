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

# Clean the Dataset and create a dataframe to work with:

new_AA_of_mutation = []
for index, row in GFP_dataset.iterrows():
    last_character_of_mutation = row["mutant"][-1]
    new_AA_of_mutation.append(last_character_of_mutation)
new_AA_of_mutation_df = pd.DataFrame(new_AA_of_mutation, columns=["new_AA"])
# Last letter in each row

number_mutations_in_mutant = GFP_dataset["mutant"].str.count(":") + 1
number_mutations_singles = number_mutations_in_mutant == 1
# True are all the rows (=mutants) that only have 1 mutation

single_mutants_only_df = new_AA_of_mutation_df[number_mutations_singles]
# This filters all the rows with a "true" value from both dataframes and creates a new dataframe while preserving the filtering.

dms_score_list = []
for index, row in GFP_dataset.iterrows():
    dms_score = row["DMS_score"]
    dms_score_list.append(dms_score)
dms_score_list_all_mutants = pd.DataFrame(dms_score_list, columns=["DMS-score"])
# Creates a dataframe with all the dms-scores and the according mutant-number (= experiment number).

all_single_dms_scores = dms_score_list_all_mutants[number_mutations_singles]
# Creates a dataframe with all the dms-scores of only the single mutants

all_single_dms_scores_with_new_AA = single_mutants_only_df.join(all_single_dms_scores)
# Combines and creates a new dataframe with the new AA and the corresponding dms-score

import matplotlib.pyplot as plt
plot_dms_scores_singles_new_AA = all_single_dms_scores_with_new_AA.plot(x="new_AA", y="DMS-score", kind="scatter", s=3)
plt.gcf().set_size_inches(10, 6)
plt.title("DMS-scores of each single mutation based on the new AA")
plot_dms_scores_singles_new_AA.set_xlabel("New amino acid")
plot_dms_scores_singles_new_AA.set_ylabel("DMS-score")
plt.axhline(2.5, color='red', linestyle='--', label='Threshold')
plt.show()
# Plot shows all the single mutants (new AA) with the corresponding dms-score.
# There is a clear cut off at 2.5 visible. Everything underneath is completely unfunctional.

# Create a new dataframe that is more suitable to work with later down the line.
mutations_pos_list = []
for index, row in GFP_dataset.iterrows():
    mutations_pos_list_number = row["mutant"][1:-1]
    mutations_pos_list.append(mutations_pos_list_number)
mutations_pos_df = pd.DataFrame(mutations_pos_list, columns=["Position"])
# Only removes the first and last character
number_mutations = GFP_dataset["mutant"].str.count(":") + 1
number_mutations_singles = number_mutations == 1
# True are all the rows containing only 1 mutation

single_mutants_with_pos = mutations_pos_df[number_mutations_singles]
# Creates a new dataframe that shows the positions within the protein from all the single mutants.

mutations_singles_pos_dms_without_new_AA = single_mutants_with_pos.join(all_single_dms_scores)
# Combines and creates a new dataframe showing the position of each single mutation and the respective dms-score (does not show the new amino acid).

new_column_for_pos = mutations_singles_pos_dms_without_new_AA["Position"]
work_with_df_pos_new_AA_dms_score = all_single_dms_scores_with_new_AA.join(new_column_for_pos)
work_with_df_pos_new_AA_dms_score = work_with_df_pos_new_AA_dms_score[["Position", "new_AA", "DMS-score"]]

# Create a heatmap depicting all the available mutations in each position, to gain a better overview.
# However, one heat map would be way too large and incomprehensible, so a subset of 4 heatmaps is created.
# These 4 heatmaps show all the single mutants separated into 4 pieces.

import seaborn as sns

heat_map_subsets_pos = np.array_split(work_with_df_pos_new_AA_dms_score['Position'].unique(), 4)  # Split into 4 subsets

# Creates some subplots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Creates a custom colormap
custom_colors_heatmap = ["lightgray", "blue"]
cmap_heatmap = sns.color_palette(custom_colors_heatmap)

# Creates a shared colorbar with two ticks for the two possibilities yes/no
cbar_ticks_heatmap = [0.25, 0.75]
cbar_ticklabels_heatmap = ['Mutation does not exist', 'Mutation exists']


# Defines a custom sorting function for the y-axis values
def heatmap_sort_by_pos(position):
    if position.isdigit():
        return int(position)
    else:
        return position

# Iterates over subsets and creates 4 heatmaps
for i, subset in enumerate(heat_map_subsets_pos):
    row_heatmap = i // 2
    col_heatmap = i % 2

    # Filters the data based on each subset
    heatmap_subset_pos_new_AA_dms_score = work_with_df_pos_new_AA_dms_score[work_with_df_pos_new_AA_dms_score['Position'].isin(subset)]

    # Sorts the data in ascending order
    heatmap_subset_pos_new_AA_dms_score = heatmap_subset_pos_new_AA_dms_score.sort_values('Position')

    # Creates cross-tabulation
    heatmap_cross_tab = pd.crosstab(heatmap_subset_pos_new_AA_dms_score['Position'], heatmap_subset_pos_new_AA_dms_score['new_AA'])

    # Sorts the columns alphabetically for better visualization
    heatmap_cross_tab_sorted = heatmap_cross_tab.reindex(sorted(heatmap_cross_tab.columns), axis=1)

    # Sorts the rows based on the custom y-axis sorting function
    heatmap_cross_tab_sorted = heatmap_cross_tab_sorted.iloc[sorted(range(len(heatmap_cross_tab_sorted)), key=lambda x: heatmap_sort_by_pos(heatmap_cross_tab_sorted.index[x]))]

    # Creates each heatmap in the corresponding subplot
    ax_heatmap = sns.heatmap(heatmap_cross_tab_sorted, cmap=cmap_heatmap, annot=False, fmt='d', ax=axes[row_heatmap, col_heatmap],
                     cbar_kws={'ticks': cbar_ticks_heatmap, 'drawedges': True}, linewidths=0.5)

    axes[row_heatmap, col_heatmap].set_xlabel('New Amino Acid')
    axes[row_heatmap, col_heatmap].set_ylabel('Mutation Position')
    axes[row_heatmap, col_heatmap].set_title(f'Subset {i+1}')

    # Sets colorbar tick labels
    cbar = ax_heatmap.collections[0].colorbar
    cbar.set_ticklabels(cbar_ticklabels_heatmap)

plt.tight_layout()

# A heatmap to see which single mutations are available in the dataset and which are not.
# Only focuses on the new amino acid that occurs, not the rest one before the mutation.

# This plot is created to show a global overview of all the positions with their respective dms-scores
# Plot shows all the single mutations that occured at a specific position (x-axis) against the dms-score (y-axis)

# Defines the x-axis range for each section in tupels
x_axis_ranges_overview_plot = [(0, 61), (61, 121), (121, 181), (181, 239)]

# Converts the "Position" column to numeric type
mutations_singles_pos_dms_without_new_AA["Position"] = pd.to_numeric(mutations_singles_pos_dms_without_new_AA["Position"])

# Creates subplots with 2 rows and 2 columns
fig, axes = plt.subplots(nrows=2, ncols=2)

# Iterates over the subplots and x-axis ranges
for i, ax in enumerate(axes.flat):
    # Slices the data based on the x-axis range
    x_start, x_end = x_axis_ranges_overview_plot[i]
    overview_plot_pos_dms = mutations_singles_pos_dms_without_new_AA[(mutations_singles_pos_dms_without_new_AA["Position"] >= x_start) & (mutations_singles_pos_dms_without_new_AA["Position"] < x_end)]

    # Plots the scatter for each section
    ax.scatter(overview_plot_pos_dms["Position"], overview_plot_pos_dms["DMS-score"], s=1)
    ax.set_title("DMS-scores of each single mutation based on the position")
    # Sets x-axis ticks to show all data points
    ax.set_xticks(overview_plot_pos_dms["Position"].unique())
    ax.set_xticklabels(overview_plot_pos_dms["Position"].unique(), rotation=90)  # Rotate x-axis labels by 90 degrees
    ax.set_xlabel("Position")  # Set x-axis label
    ax.set_ylabel("DMS-score")  # Set y-axis label
# Customizes individual subplots as needed
# Checks if the specific x-axis labels are present in the current subplot
    if 65 in overview_plot_pos_dms["Position"].unique():
        ax.axvline(x=65, color='red', linestyle='--')
    if 67 in overview_plot_pos_dms["Position"].unique():
        ax.axvline(x=67, color='red', linestyle='--')

    ax.axhline(y=2.5, color='red', linestyle='--')
# Inserts some threshold-lines to indicate important positions (chromophor).

# Show the plot
plt.gcf().set_size_inches(17, 10)  # Adjust the size of the overall figure
plt.tight_layout()  # Ensures proper spacing between subplots
plt.show()

# It's interesting to see here, that there are mutations occuring in the chromophore (positions 65-67) that still have a good dms-score.
# Position 66 is the fluorescing amino acid, and position 67 also participates in the cyclization reaction, whereas position 65 has a smaller impact on the fluorescence.
# Therefore, it makes sense that position 66 and 67 show lowered dms-scores, whilst position 65 is more stable.


