#import
import pandas as pd
import numpy as np
# read dataset
original_dms_data = pd.read_csv('/Users/liza/Documents/Bioinfo Project/DMS_data/AAAA_GFP_dms_data_original_komplett.csv')
# split first column of df into multiple columns
original_dms_data_col = original_dms_data
only_mutants = original_dms_data["mutant"].to_frame()
original_dms_data_col[['m1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11', 'm12', 'm13', 'm14', 'm15']] = original_dms_data_col['mutant'].str.split(':', 15, expand=True)

# count how many mutations each sequence has
list_mut_count_in_progress = []
for i in range(len(original_dms_data['mutant'])):
    list_mut_count_in_progress.append(original_dms_data['mutant'].iloc[i].count(':'))
list_mut_count_prae = np.array(list_mut_count_in_progress)
list_mut_count = (list_mut_count_prae + 1)
df_mutation_counts = pd.DataFrame(list_mut_count)

#generate more convenient dataframe
working_dataframe_prae = pd.concat([original_dms_data_col, df_mutation_counts], axis="columns")
#drop columns we don´t need at the moment -> working_dataframe
working_dataframe = working_dataframe_prae.drop(['mutant', 'mutated_sequence', 'DMS_score_bin'], axis=1)
working_dataframe.rename(columns={working_dataframe.columns[16]: 'mut_count'}, inplace=True)
#another dataframe for easy access -> nur_fscore_mut_count
nur_fscore_mut_count = working_dataframe.loc[:, ["DMS_score", "mut_count"]]

#goal: one list with all existing mutations in the dataset
#-> all_possible_mutations
working_dataframe_only_ms = working_dataframe.loc[:, ["m1", "m2", "m3", 'm4', 'm5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11', 'm12', 'm13', 'm14', 'm15']]
all_possible_mutations = working_dataframe_only_ms.values.flatten().tolist()
all_possible_mutations = list(set(all_possible_mutations))
#automatically: one "none" value: drop
while None in all_possible_mutations:
    all_possible_mutations.remove(None)
only_mutants_list = only_mutants['mutant']

#checkpoint 1
print('all_possible_mutations finished (1)')

#again: generate more convenient dataframe, add mutation_count to original one
working_dataframe_prae = pd.concat([original_dms_data_col, df_mutation_counts], axis="columns")
#drop other columns
working_dataframe = working_dataframe_prae.drop(['mutant', 'mutated_sequence', 'DMS_score_bin'], axis=1)
working_dataframe.rename(columns={working_dataframe.columns[16]: 'mut_count'}, inplace=True)

#goal: dataframe that contains boolians if the mutations from all_possible_mutations exist in the mutants (Kreuztabelle)
list_of_dfs = []

for i in all_possible_mutations:
    new_column_name = f'{i}'
    new_column_values = [only_mutants_list.str.contains(i, regex= False)]
    new_df = pd.DataFrame({new_column_name: new_column_values})
    new_df_exploded = new_df.explode(new_column_name)
    list_of_dfs.append(new_df_exploded)
result_how_often = pd.concat(list_of_dfs, axis=1)
result_how_often = result_how_often.reset_index(drop=True)

## --> result_how_often (all_possible_mutations (columns), mutants (rows))

# dataframe generated from original df: only fscore (= DMS_score) and mut_count
count_fscore_frame = working_dataframe[['DMS_score', 'mut_count']]

#goal: calculate variance and the values used for variance calculation for each mutation and each mutocunt
#normally: present it in graph (#) -> stored
import matplotlib.pyplot as plt

fig, axes = plt.subplots(nrows=5, ncols=3, figsize=(19, 12))  # Abbildung und Achsenobjekte erstellen
plt.subplots_adjust(wspace=0.4, hspace=0.6)

scatter_plots = []
for j, ax in zip(range(2, 16), axes.flatten()):
    variance_per_mutant_list = []

    for i in all_possible_mutations:
        mut_count_fscore = count_fscore_frame.loc[result_how_often[i] == True]
        fscore_mut = mut_count_fscore['DMS_score'].loc[mut_count_fscore['mut_count'] == j]
        varianz_mut = fscore_mut.var()
        variance_per_mutant_list.append(varianz_mut)

    variance_per_mutant_series = pd.Series(variance_per_mutant_list, index=all_possible_mutations)
    variance_per_mutant_df = variance_per_mutant_series.to_frame()

#-> variance_per_mutant

    how_many_for_variance = []

    for i in all_possible_mutations:
        mut_count_fscore = count_fscore_frame.loc[result_how_often[i] == True]
        fscore_mut = mut_count_fscore['DMS_score'].loc[mut_count_fscore['mut_count'] == j]
        wie_viel_jeweils = len(fscore_mut)
        how_many_for_variance.append(wie_viel_jeweils)

    how_many_for_variance = pd.Series(how_many_for_variance, index=all_possible_mutations)
    how_many_for_variance_df = how_many_for_variance.to_frame()
# -> how_many_for_variance
    how_many_AND_variance_df = pd.concat([how_many_for_variance_df, variance_per_mutant_df], axis = 1)
    how_many_AND_variance_df.columns = ['Anzahl benutzter Werte', 'Varianz']
    how_many_AND_variance_df = how_many_AND_variance_df.dropna()

#scatter plot
    ax.scatter(how_many_AND_variance_df['Anzahl benutzter Werte'],how_many_AND_variance_df['Varianz'], s = j )
    ax.set_xlabel('Anzahl benutzter Werte')
    ax.set_ylabel('Varianz')

    if "V163A" in how_many_AND_variance_df.index:
        ax.scatter(how_many_AND_variance_df['Anzahl benutzter Werte']['V163A'],how_many_AND_variance_df['Varianz']['V163A'], c='red')
    ax.set_title(f'für {j} Mutationen')

    scatter_plots.append(scatter_plot)
saved_plots = []

for scatter_plot in scatter_plots:
    fig.canvas.draw()
    plot_image = np.array(fig.canvas.renderer.buffer_rgba())
    saved_plots.append(plot_image)

plt.close(fig)
#plot saved in variable

#checkpoint 2
print('variances plots finished (2)')

#goal: variances per mutation (mean over all mutcounts)
frame_zum_mitteln_variance = pd.DataFrame(index = all_possible_mutations)
variance_per_mutant_count_list = []

# IMPORTANT: only mutants with a mutcount from 2 to 8 are considered for calculation -> boxplot (see output_variance_analysis.ipynb)
for j, ax in zip(range(2, 8), axes.flatten()):
    variance_per_mutant_list = []

    for i in all_possible_mutations:
        mut_count_fscore = count_fscore_frame.loc[result_how_often[i] == True]
        fscore_mut = mut_count_fscore['DMS_score'].loc[mut_count_fscore['mut_count'] == j]
        varianz_mut = fscore_mut.var()  #variance per mutation per mutcount
        variance_per_mutant_list.append(varianz_mut) #list of variances per mutation per mutcount of all mutations

    variance_per_mutant_df = pd.DataFrame(variance_per_mutant_list, index=all_possible_mutations)
    variance_per_mutant_count_list.append(variance_per_mutant_df)
variance_per_mutant_count_df = pd.concat(variance_per_mutant_count_list, axis=1)
variance_per_mutant_count_df.set_axis(range(2,8), axis=1, inplace=True)
# -> variance_per_mutant_count_df (per mutation per mutant count)
mean_variances_per_mutations = pd.DataFrame(variance_per_mutant_count_df.mean(axis=1, skipna=True), columns=['Mean'])
# mean_variances_per_mutations (all variances per all mutations (rows) per all counts (columns))

#goal: dataframe with the mean fitness scores for each mutcount (-> weights)
mean_fitness_scores = pd.DataFrame(index = range(2,16), columns = ["mean_fitness_score"])
for i in range(2,16):
    fscore_mutcount_mean = count_fscore_frame["DMS_score"].loc[count_fscore_frame["mut_count"] == i].mean()
    mean_fitness_scores.loc[i, "mean_fitness_score"] = fscore_mutcount_mean
#->mean_fitness_scores (weights)

#goal: dataframe that contains the amount of mutants a mutation is part of (how often does a mutation X appear in the dataset = Vorkommen)
how_many_per_mutant_count_list = []

for j, ax in zip(range(2, 8), axes.flatten()):
    how_many_for_variance = []

    for i in all_possible_mutations:
        mut_count_fscore = count_fscore_frame.loc[result_how_often[i] == True]
        fscore_mut = mut_count_fscore['DMS_score'].loc[mut_count_fscore['mut_count'] == j]
        wie_viel_jeweils = len(fscore_mut)
        how_many_for_variance.append(wie_viel_jeweils)

    how_many_per_mutant_df = pd.DataFrame(how_many_for_variance, index=all_possible_mutations)
    how_many_per_mutant_count_list.append(how_many_per_mutant_df)
how_many_per_mutant_count_df = pd.concat(how_many_per_mutant_count_list, axis=1)
how_many_per_mutant_count_df.set_axis(range(2,8), axis=1, inplace=True)

mean_how_many_per_mutations = pd.DataFrame(how_many_per_mutant_count_df.mean(axis=1, skipna=True), columns=['Mean'])

#goal: WEIGHTED DIFFERENCE: weighted_fscore_mean from the mutants WITH the mutation - weighted_fscore_mean from the mutants WITHOUT the mutation

#first: WITH

# Create an empty dataframe to store the results
mean_for_differences_with_neu = pd.DataFrame(index=all_possible_mutations, columns=range(2, 16))

for mutation in all_possible_mutations:
    for mutation_count in range(2, 16):

        index_when_mut_present_weighted = result_how_often.loc[result_how_often[mutation] == True].index
        only_rows_with_mut_weighted = nur_fscore_mut_count[(nur_fscore_mut_count['mut_count'] == mutation_count) & (nur_fscore_mut_count.index.isin(index_when_mut_present_weighted))]
        # Calculate the mean of fitness score for the filtered rows
        mean_fitness_with_mut_score = only_rows_with_mut_weighted['DMS_score'].mean()
        mean_for_differences_with_neu.loc[mutation, mutation_count] = mean_fitness_with_mut_score

weighted_means_WITH = pd.Series(index=mean_for_differences_with_neu.index)

for mutation in mean_for_differences_with_neu.index:
    row_WITH = mean_for_differences_with_neu.loc[mutation]
    non_nan_values_WITH = row_WITH.dropna()

    if len(non_nan_values_WITH) > 0:
        non_nan_weights_WITH = mean_fitness_scores.loc[non_nan_values_WITH.index]['mean_fitness_score']
        weighted_means_WITH[mutation] = np.average(non_nan_values_WITH, weights=non_nan_weights_WITH)

weighted_means_df_WITH = pd.DataFrame({'Weighted Mean WITH': weighted_means_WITH})

#NEXT: sam for WITHOUT


