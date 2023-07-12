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


# Creates a new dataframe that is more suitable to work with later down the line.
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

# The columns were grouped by amino acid and by position in new separate dataframes

grouped_by_amino_acid = work_with_df_pos_new_AA_dms_score.groupby('new_AA')
grouped_by_position = work_with_df_pos_new_AA_dms_score.groupby("Position")
mean_dms_by_amino_acid = grouped_by_amino_acid["DMS-score"].mean()
median_dms_by_amino_acid = grouped_by_amino_acid["DMS-score"].median()
mean_dms_by_position = grouped_by_position["DMS-score"].mean()
median_dms_by_position = grouped_by_position["DMS-score"].median()
results_mean_median_AA_df = pd.DataFrame({
    "Mean DMS by Amino Acid": mean_dms_by_amino_acid,
    "Median DMS by Amino Acid": median_dms_by_amino_acid})
results_mean_median_pos_df = pd.DataFrame({
    "Mean DMS by Position": mean_dms_by_position,
    "Median DMS by Position": median_dms_by_position})

# Sorts the DataFrame ascending position wise
results_mean_median_pos_df.sort_values("Position", inplace=True)

# Resets the index
results_mean_median_pos_df.reset_index(drop=True, inplace=True)

print(results_mean_median_pos_df)
print(results_mean_median_AA_df)
# General analysis of the mean and median of the dms-score of both factors: new amino acid and position of the mutation.
# Note: The amount of mutations for each new amino acid and each mutation has to be taken into account.
# If there are few mutations, the values are not very reliable.


# Now we analyse specific amino acid mutations (new amino acids) by searching for the last character of each mutation name.

wanted_specific_new_AA = "R"
# This is the specific amino acid I am looking for
mask_wanted_specific_new_AA = work_with_df_pos_new_AA_dms_score["new_AA"] == wanted_specific_new_AA
# This creates a boolean mask
wanted_specific_new_AA_df = work_with_df_pos_new_AA_dms_score[mask_wanted_specific_new_AA]
# All rows that fulfill the condition are indexed


# Plotting of the results, where everything is dependent on wanted_specific_new_AA
import matplotlib.pyplot as plt
plt.scatter(wanted_specific_new_AA_df["Position"], wanted_specific_new_AA_df["DMS-score"], color='b')
plt.xlabel("Position")
plt.ylabel("DMS-score")
plt.title(f"DMS-score vs. Position for the new amino acid = {wanted_specific_new_AA}")
plt.xticks(rotation='vertical')
# Rotate the x-axis labels vertically
plt.gcf().set_size_inches(15, 6)
plt.show()
# Plots all the found positions of the previous AA-search with their corresponding dms-scores.


import seaborn as sns

# Using the grouped data
grouped_by_amino_acid = work_with_df_pos_new_AA_dms_score.groupby('new_AA')

# Creates an empty list to store the DMS scores for each group
dms_scores_per_group_new_AS = []

# Iterates over each group and extracts the DMS scores
for position, group in grouped_by_amino_acid:
    dms_scores_per_group_new_AS.append(group['DMS-score'])

# Creates a box plot to depict the results
sns.boxplot(x='new_AA', y='DMS-score', data=work_with_df_pos_new_AA_dms_score)
plt.xlabel('new_AA')
plt.ylabel('DMS-score')
plt.title('Distribution of DMS Scores by New_AS')
plt.gcf().set_size_inches(20, 10)
plt.show()

# Box has 50% of all values, from the first quantile to the 3rd (25% to 75%)
# The whiskers reach from 0%-25% and 75%-100%. There are outliers visible.
# Arginine (R) and Proline (P) display big boxes, which means they have a high variance within their group. Therefore, mutations involving a Proline or Arginine have a broad range of dms-scores --> Interesting to analyse.
# It seems, that this effect of dms-score distribution is dependent of the position where the mutation occurs.
# This does not show, how many values each box contains, so if an amino acid only has few values, the boxplot does not convey reliable information.

# Now a violin plot was created to illustrate the distribution of data in an additional dimension (width of the violins)

sns.violinplot(x='new_AA', y='DMS-score', data=work_with_df_pos_new_AA_dms_score)
plt.xlabel('new_AA')
plt.ylabel('DMS-score')
plt.title('Distribution of DMS Scores by New_AA')
plt.gcf().set_size_inches(20, 10)
plt.show()
# A violin plot displays the entire distribution of the data, including information about the kernel density estimation.
# The width of the violin at a particular point represents the density or probability of observing data points at that value. The plot may also include lines or markers to indicate summary statistics such as the median or quartiles.
# Black dot = median
# P and R show a large interquantile range, which means the data is broadly distributed, hence the narrow and elongated shape of the plot. The other amino acids show clear boxplots inside the violin plot and also have more compact distributions.
# This overlaps with the boxplots up above


# Examining the Chromophore

wanted_specific_pos_65 = "65"
# This represents the specific position inside the chromophore that I am looking for
maske_wanted_specific_pos_65 = work_with_df_pos_new_AA_dms_score["Position"] == wanted_specific_pos_65
# Creates a boolean mask
wanted_specific_pos_65_df = work_with_df_pos_new_AA_dms_score[maske_wanted_specific_pos_65]
# All the rows that fulfill this condition will be indexed
print(wanted_specific_pos_65_df)
#--------------------------------------------------------------------------------------------------------------------------
wanted_specific_pos_66 = "66"
maske_wanted_specific_pos_66 = work_with_df_pos_new_AA_dms_score["Position"] == wanted_specific_pos_66
wanted_specific_pos_66_df = work_with_df_pos_new_AA_dms_score[maske_wanted_specific_pos_66]
print(wanted_specific_pos_66_df)
#----------------------------------------------------------------------------------------------------------------------------
wanted_specific_pos_67 = "67"
maske_wanted_specific_pos_67 = work_with_df_pos_new_AA_dms_score["Position"] == wanted_specific_pos_67
wanted_specific_pos_67_df = work_with_df_pos_new_AA_dms_score[maske_wanted_specific_pos_67]
print(wanted_specific_pos_67_df)
# Out of the mutations we have available in our dataset, position 66 and 67 are completely deleterious.
# Position 65 experiences loss of function when mutated to proline, but maintains some activity when mutated to Leucine or Alanine, whilst a mutation to Threonine has the highest remaining activity.
# This is due to the necessary folding of the chromophore to achieve fluorescence.