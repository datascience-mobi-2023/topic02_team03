#import
import pandas as pd
import numpy as np
pd.set_option('display.max_columns', 50)

#read dataset
original_dms_data = pd.read_csv('/Users/liza/Documents/Bioinfo Project/DMS_data/AAAA_GFP_dms_data_original_komplett.csv')
# split first column of df into multiple columns
original_dms_data_col = original_dms_data
only_mutants = original_dms_data["mutant"].to_frame()
original_dms_data_col[['m1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11', 'm12', 'm13', 'm14', 'm15']] = original_dms_data_col['mutant'].str.split(':', 15, expand=True)

#PREPARE DATASET FOR RANKING CALCULATION:
#1. how many mutations does each mutant have? -> mut_count (df_mutation_count)
#2. find out all possible mutations -> all_possible_mutations
#3. df mutations x mutants -> which mutants contain which mutations -> result_how_often
#4. plots + calculation variance x how many values used for variance per mutation count
#-> variance of all fscores of all mutants containing mutation X shows how constant the effect of the mutation is

# count how many mutations each mutant has
list_mut_count_in_progress = []
for i in range(len(original_dms_data['mutant'])):
    list_mut_count_in_progress.append(original_dms_data['mutant'].iloc[i].count(':'))
list_mut_count_prae = np.array(list_mut_count_in_progress)
list_mut_count = (list_mut_count_prae + 1)
df_mutation_counts = pd.DataFrame(list_mut_count)

#concat mutation_count to original df
working_dataframe_prae = pd.concat([original_dms_data_col, df_mutation_counts], axis="columns")

working_dataframe = working_dataframe_prae.drop(['mutant', 'mutated_sequence', 'DMS_score_bin'], axis=1)
working_dataframe.rename(columns={working_dataframe.columns[16]: 'mut_count'}, inplace=True)

#all existing mutations in one list: calculation of all_possible_mutations

working_dataframe_only_ms = working_dataframe.loc[:, ["m1", "m2", "m3", 'm4', 'm5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11', 'm12', 'm13', 'm14', 'm15']]
all_possible_mutations = working_dataframe_only_ms.values.flatten().tolist()
all_possible_mutations = list(set(all_possible_mutations))
#"none" is also a value in the list -> remove
while None in all_possible_mutations:
    all_possible_mutations.remove(None)
only_mutants_list = only_mutants['mutant']

#remodel df so it is more convenient
working_dataframe_prae = pd.concat([original_dms_data_col, df_mutation_counts], axis="columns")
working_dataframe = working_dataframe_prae.drop(['mutant', 'mutated_sequence', 'DMS_score_bin'], axis=1)
working_dataframe.rename(columns={working_dataframe.columns[16]: 'mut_count'}, inplace=True)

#df which mutant from only_mutants_list contains the mutations from all_possible_mutations
list_of_dfs = []
for i in all_possible_mutations:
    new_column_name = f'{i}'
    new_column_values = [only_mutants_list.str.contains(i, regex= False)]
    new_df = pd.DataFrame({new_column_name: new_column_values})
    new_df_exploded = new_df.explode(new_column_name)
    list_of_dfs.append(new_df_exploded)

# concat dfs from the list
result_how_often = pd.concat(list_of_dfs, axis=1)
result_how_often = result_how_often.reset_index(drop=True)

# more convenient df (just mutcount und fscore)
count_fscore_frame = working_dataframe[['DMS_score', 'mut_count']]

#plot: variance of fscores per mutation count and number of values uesed for variance calculation
#-> bottom right corner is the best, because the varaince is the most reliable

import matplotlib.pyplot as plt


fig, axes = plt.subplots(nrows=5, ncols=3, figsize=(19, 12))  # Abbildung und Achsenobjekte erstellen
plt.subplots_adjust(wspace=0.4, hspace=0.6)

for j, ax in zip(range(2, 16), axes.flatten()):
    variance_per_mutant_list = []

    for i in all_possible_mutations:
        mut_count_fscore = count_fscore_frame.loc[result_how_often[i] == True]
        fscore_mut = mut_count_fscore['DMS_score'].loc[mut_count_fscore['mut_count'] == j]
        varianz_mut = fscore_mut.var()
        variance_per_mutant_list.append(varianz_mut)

    variance_per_mutant_series = pd.Series(variance_per_mutant_list, index=all_possible_mutations)
    variance_per_mutant_df = variance_per_mutant_series.to_frame()


#how reliable is the calculated variance? -> how many values are available for the calculation

    how_many_for_variance = []

    for i in all_possible_mutations:
        mut_count_fscore = count_fscore_frame.loc[result_how_often[i] == True]
        fscore_mut = mut_count_fscore['DMS_score'].loc[mut_count_fscore['mut_count'] == j]
        wie_viel_jeweils = len(fscore_mut)
        how_many_for_variance.append(wie_viel_jeweils)

    how_many_for_variance = pd.Series(how_many_for_variance, index=all_possible_mutations)
    how_many_for_variance_df = how_many_for_variance.to_frame()


    how_many_AND_variance_df = pd.concat([how_many_for_variance_df, variance_per_mutant_df], axis = 1)
    how_many_AND_variance_df.columns = ['number of used values', 'Variance']
    how_many_AND_variance_df = how_many_AND_variance_df.dropna()

#scatter plot variance x count
    ax.scatter(how_many_AND_variance_df['number of used values'],how_many_AND_variance_df['Variance'], s = j )
    ax.set_xlabel('number of used values')
    ax.set_ylabel('Variance')

    if "V163A" in how_many_AND_variance_df.index:
        ax.scatter(how_many_AND_variance_df['number of used values']['V163A'],how_many_AND_variance_df['Variance']['V163A'], c='red')
    ax.set_title(f'mutation count = {j}')
                 
#FURTHER PREPARATION FOR RANKINGS:
#1. variance calculation per mutant (not-weighted)
#2. how many values got used for the calculation, same as above but for every mutcount combined -> both in one dataframe
#3. calculate df that contains the difference between the fscore-means of all mutants WITH mutation X and WITHOUT mutation X
# -> how big is the effect of the mutation on existing in general
# -> not weighted


#goal: variance means, calculate variances per mutation per mutation count -> mean of all per mutation
#ONLY MUTCOUNT FROM 2 TO 7 !! -> 5a
#-> not yet weighted
frame_zum_mitteln_variance = pd.DataFrame(index = all_possible_mutations)
variance_per_mutant_count_list = []

# IMPORTANT!! -> only values from mutants with up to 7 mutations (mut_count < 7), because the fscores from mutants with a mut_counts above 7 are generally low (-> boxplot from statistical_test_fscores)
for j, ax in zip(range(2, 8), axes.flatten()):
    variance_per_mutant_list = []

    for i in all_possible_mutations:
        mut_count_fscore = count_fscore_frame.loc[result_how_often[i] == True]
        fscore_mut = mut_count_fscore['DMS_score'].loc[mut_count_fscore['mut_count'] == j]
        varianz_mut = fscore_mut.var()  #die varianz je mutation je anzahl
        variance_per_mutant_list.append(varianz_mut) #liste der Varianzen ALLER Mutationen je anzahl

    variance_per_mutant_df = pd.DataFrame(variance_per_mutant_list, index=all_possible_mutations)
    variance_per_mutant_count_list.append(variance_per_mutant_df)
variance_per_mutant_count_df = pd.concat(variance_per_mutant_count_list, axis=1)
variance_per_mutant_count_df.set_axis(range(2,8), axis=1, inplace=True)

# df variance per mutatuion (rows) per mutation_count (columns)
mean_variances_per_mutations = pd.DataFrame(variance_per_mutant_count_df.mean(axis=1, skipna=True), columns=['Mean'])

#very similar, again
#variance means, calculate variances per mutation per mutation count -> mean of all per mutation
#ALL MUTCOUNTS INVOLVED -> 5c
#-> not yet weighted
frame_zum_mitteln_variance_16 = pd.DataFrame(index=all_possible_mutations)
variance_per_mutant_count_list_16 = []

# da die Fitness-Scores von Mutanten mit einer mut_count über 7 im Allgemeinen niedrig sind.
for j in range(2, 16):
    variance_per_mutant_list_16 = []

    for i in all_possible_mutations:
        mut_count_fscore_16 = count_fscore_frame.loc[result_how_often[i] == True]
        fscore_mut_16 = mut_count_fscore_16['DMS_score'].loc[mut_count_fscore_16['mut_count'] == j]
        varianz_mut_16 = fscore_mut_16.var()  # Varianz für die aktuelle Mutation und Mutation_count
        variance_per_mutant_list_16.append(varianz_mut_16)

    variance_per_mutant_df_16 = pd.DataFrame(variance_per_mutant_list_16, index=all_possible_mutations)
    variance_per_mutant_count_list_16.append(variance_per_mutant_df_16)

variance_per_mutant_count_df_16 = pd.concat(variance_per_mutant_count_list_16, axis=1)
variance_per_mutant_count_df_16.columns = range(2, 16)

mean_variances_per_mutations_16 = pd.DataFrame(variance_per_mutant_count_df_16.mean(axis=1, skipna=True), columns=['Mean'])

#goal: weights for weighted variance in ranking 5d
#variance gets weighted with how many values excist for calculation
how_many_per_mutant_count_list = []

for j, ax in zip(range(2, 16), axes.flatten()):
    how_many_for_variance = []

    for i in all_possible_mutations:
        mut_count_fscore = count_fscore_frame.loc[result_how_often[i] == True]
        fscore_mut = mut_count_fscore['DMS_score'].loc[mut_count_fscore['mut_count'] == j]
        wie_viel_jeweils = len(fscore_mut)
        how_many_for_variance.append(wie_viel_jeweils)

    how_many_per_mutant_df = pd.DataFrame(how_many_for_variance, index=all_possible_mutations)
    how_many_per_mutant_count_list.append(how_many_per_mutant_df)
how_many_per_mutant_count_df = pd.concat(how_many_per_mutant_count_list, axis=1)
how_many_per_mutant_count_df.set_axis(range(2,16), axis=1, inplace=True)

mean_how_many_per_mutations = pd.DataFrame(how_many_per_mutant_count_df.mean(axis=1, skipna=True), columns=['Mean'])

combined_means_variance_how_many = pd.concat([mean_variances_per_mutations, mean_how_many_per_mutations], axis=1)
combined_means_variance_how_many.columns = ['mean_variances_per_mutations', 'mean_how_many_per_mutations']


#goal: differences calculated (effect of mutation X on the fscore)
nur_fscore_mut_count = working_dataframe.loc[:, ["DMS_score", "mut_count"]]
differences_list = []

for i in all_possible_mutations:
    #for all mutants WITH the mutation
    #filter the mutants for the existing of mutation X (-> result_how_often)
    index_when_mut_present = result_how_often.loc[result_how_often[i] == True].index
    #fscores from all mutants that match
    only_rows_with_mut = nur_fscore_mut_count [(nur_fscore_mut_count ['mut_count'] >2) & (nur_fscore_mut_count .index.isin(index_when_mut_present))]

# Calculate the mean of DMS_score for the filtered rows
    mean_dms_score_only_mut = only_rows_with_mut['DMS_score'].mean()
#-------------
    #for all mutants WITHOUT the mutation (same calculation)
    index_when_not_mut_present = result_how_often.loc[result_how_often[i] == False].index

    only_rows_withOUT_mut = working_dataframe[(working_dataframe['mut_count'] >2) & (working_dataframe.index.isin(index_when_not_mut_present))]

# Calculate the mean of DMS_score for the filtered rows
    mean_dms_score_every_but_mut = only_rows_withOUT_mut['DMS_score'].mean()
#----------------
    difference_means = mean_dms_score_only_mut - mean_dms_score_every_but_mut
    differences_list.append(difference_means)

all_differences_means = pd.DataFrame({'Difference': differences_list}, index=all_possible_mutations)

#—————-RANKING 0:
# -> ranked by variance (without taking into account how many values got used for the calculation

sorted_Ranking0 = combined_means_variance_how_many.sort_values(by='mean_variances_per_mutations')

#—————-RANKING 1: (that was the try if the ranking-function works)
# -> ranked by the available values for each mutation (how often does mutation X appear in total)

sorted_Ranking1 = combined_means_variance_how_many.sort_values(by='mean_how_many_per_mutations', ascending=False)

#-------RANKING 1a:
#-> ranked by a rank_score build from the variance and the count of available mutations
#->combination of ranking 0 and 1

combined_means_variance_how_many['Rank'] = combined_means_variance_how_many['mean_variances_per_mutations'].rank(
        ascending=False) - combined_means_variance_how_many['mean_how_many_per_mutations'].rank()

sorted_Ranking1a = combined_means_variance_how_many.sort_values(by='Rank')

# same ranking just with only the stabilizing (Difference > 0)
condition = all_differences_means['Difference'] > 0
sorted_only_stab_Ranking1a = sorted_Ranking1a.drop(all_differences_means.loc[condition].index)

#------RANKING 2:
#-> only ranked by Difference
#-> how big is the effect the mutation has on other excisting mutations (->fscores)

# how often does the mutation X appear in all mutants
list_wie_oft_mut = []
for j in all_possible_mutations:
    matching_indexes = result_how_often.loc[result_how_often[j] == True].index
    wie_oft = len(matching_indexes)
    list_wie_oft_mut.append(wie_oft)
df_wie_oft_muts_insg = pd.DataFrame(list_wie_oft_mut, index=all_possible_mutations)

# difference and count in one dataframe
combined_differenz_wie_oft_mut = pd.concat([all_differences_means, df_wie_oft_muts_insg], axis=1)
combined_differenz_wie_oft_mut.columns = ['Difference', 'wie oft kommt mut insg vor']

ranking2 = combined_differenz_wie_oft_mut.sort_values(by='Difference', ascending=False)

#-----RANKING 3:
#-> ranked by difference and count

list_ranking3 = []
for i in all_possible_mutations:
    score_ranking3 = all_differences_means.loc[i].values[0] * df_wie_oft_muts_insg.loc[i].values[0]
    list_ranking3.append(score_ranking3)
ranking3_unsorted = pd.DataFrame(list_ranking3, index=all_possible_mutations, columns=['ranking3_score'])
ranking3 = ranking3_unsorted.sort_values(by='ranking3_score', ascending=False)

#-----RANKING 4:
#-> calculate own score1
#->score1 = Difference * 1/Variance * how_often
#->variance here: just the mean of all variances per mutcount per mutation (2 to 7)

list_ranking4 = []
for i in all_possible_mutations:
    score_ranking4 = all_differences_means.loc[i].values[0] * df_wie_oft_muts_insg.loc[i].values[0] * (
                    1 / mean_variances_per_mutations.loc[i].values[0])
    list_ranking4.append(score_ranking4)
ranking4_unsorted = pd.DataFrame(list_ranking4, index=all_possible_mutations, columns=['ranking4_score'])
ranking4 = ranking4_unsorted.sort_values(by='ranking4_score', ascending=False)

#------RANKING 5:
#-> calculate own score, just like ranking 4 just different kind of variances

#PREP for 5b: variances weighted with mean_fitness_scores as weights
#dataframe with the mean fscores of each mutation count
mean_fitness_scores = pd.DataFrame(index = range(2,16), columns = ["mean_fitness_score"])
for i in range(2,16):
    fscore_mutcount_mean = count_fscore_frame["DMS_score"].loc[count_fscore_frame["mut_count"] == i].mean()
    mean_fitness_scores.loc[i, "mean_fitness_score"] = fscore_mutcount_mean

weighted_variances_try = pd.Series(index=variance_per_mutant_count_df_16.index)

for mutation in all_possible_mutations:
    row_variance = variance_per_mutant_count_df_16.loc[mutation]
    non_nan_values_variance = row_variance.dropna()

    if len(non_nan_values_variance) > 0:
        non_nan_weights_variance = mean_fitness_scores.loc[non_nan_values_variance.index]['mean_fitness_score']
        weighted_variances_try[mutation] = np.average(non_nan_values_variance, weights=non_nan_weights_variance)

weighted_variances_try_df = pd.DataFrame({'Weighted Variances': weighted_variances_try})

#PREP for 5d: variances weighted with wieviel as weights
#dataframe with the mean fscores of each mutation count

for mutation in all_possible_mutations:
    row_variance = variance_per_mutant_count_df_16.loc[mutation]
    non_nan_values_variance = row_variance.dropna()
    row_weight_how_many = how_many_per_mutant_count_df.loc[mutation].to_frame(name = 'how_many')

    if len(non_nan_values_variance) > 0:
        non_nan_weights_variance = row_weight_how_many.loc[non_nan_values_variance.index]['how_many']
        weighted_variances_try[mutation] = np.average(non_nan_values_variance, weights=non_nan_weights_variance)

weighted_variances_d_df = pd.DataFrame({'Weighted Variances': weighted_variances_try})

#RANKING 5_0: wrong calculation of the variance (only the mutations that have complete data for all mutcounts
#not mean of variance but sum! -> how "complete" is dataset

list_ranking5_0 = []
for i in all_possible_mutations:

    score_ranking5_0 =score_ranking5_0 = (all_differences_means.loc[i].values[0] * df_wie_oft_muts_insg.loc[i].values[0]) / np.nansum(variance_per_mutant_count_df.loc[i].values)

    list_ranking5_0.append(score_ranking5_0)
ranking5_0_unsorted = pd.DataFrame(list_ranking5_0, index=all_possible_mutations, columns=['ranking5_0_score'])
ranking5_0 = ranking5_0_unsorted.sort_values(by='ranking5_0_score', ascending= False)

#RANKING 5a:
# variance calculated with the mean over all mut count varaiances
list_ranking5_a = []
for i in all_possible_mutations:
    score_ranking5_a = (all_differences_means.loc[i].values[0]) * df_wie_oft_muts_insg.loc[i].values[0] / (mean_variances_per_mutations_16.loc[i].values[0])
    list_ranking5_a.append(score_ranking5_a)
ranking5_unsorted_a = pd.DataFrame(list_ranking5_a, index=all_possible_mutations, columns=['ranking5_a_score'])
ranking5_a = ranking5_unsorted_a.sort_values(by='ranking5_a_score', ascending=False)

#RANKING 5b:
#variance calculated with the mean fscores of the corresponding mutcount group
list_ranking5_b = []
for i in all_possible_mutations:

    score_ranking5_b = (all_differences_means.loc[i].values[0]) * df_wie_oft_muts_insg.loc[i].values[0] / (weighted_variances_try_df.loc[i].values[0])
    list_ranking5_b.append(score_ranking5_b)
ranking5_unsorted_b = pd.DataFrame(list_ranking5_b, index=all_possible_mutations, columns=['ranking5_b_score'])
ranking5_b = ranking5_unsorted_b.sort_values(by='ranking5_b_score', ascending= False)

#RANKING 5c:
#variance mean from the variances from the mutcount groups 2 to 7
list_ranking5_c = []
for i in all_possible_mutations:

    score_ranking5_c = (all_differences_means.loc[i].values[0]) * df_wie_oft_muts_insg.loc[i].values[0] / (mean_variances_per_mutations.loc[i].values[0])
    list_ranking5_c.append(score_ranking5_c)
ranking5_unsorted_c = pd.DataFrame(list_ranking5_c, index=all_possible_mutations, columns=['ranking5_c_score'])
ranking5_c = ranking5_unsorted_c.sort_values(by='ranking5_c_score', ascending= False)

#RANKING 5d:
#variance mean weighted with how many values there are for each group
list_ranking5_d = []
for i in all_possible_mutations:

    score_ranking5_d = (all_differences_means.loc[i].values[0]) * df_wie_oft_muts_insg.loc[i].values[0] / (mean_variances_per_mutations.loc[i].values[0])
    list_ranking5_d.append(score_ranking5_d)
ranking5_unsorted_d = pd.DataFrame(list_ranking5_d, index=all_possible_mutations, columns=['ranking5_d_score'])
ranking5_d = ranking5_unsorted_d.sort_values(by='ranking5_d_score', ascending= False)

#------RANKING 6:
#->delta G values ranked on their own
#-> but calculated as differences as fscores

delta_G_data = pd.read_csv('/Users/liza/Downloads/df_ddG.csv')
count_fscore_frame['delta G'] = delta_G_data['Score']

#goal: differences from delta G values
#same calculation as for the fscores above to get the general impact of one mutation on the delta G values
differences_delta_G_list = []

for i in all_possible_mutations:
    index_when_mut_present = result_how_often.loc[result_how_often[i] == True].index

    only_rows_with_mut = count_fscore_frame[(count_fscore_frame['mut_count'] >2) & (count_fscore_frame .index.isin(index_when_mut_present))]

# Calculate the mean of DMS_score for the filtered rows
    mean_delta_G_only_mut = only_rows_with_mut['delta G'].mean()
#-------------
    index_when_not_mut_present = result_how_often.loc[result_how_often[i] == False].index

    only_rows_withOUT_mut = count_fscore_frame[(count_fscore_frame['mut_count'] >2) & (count_fscore_frame.index.isin(index_when_not_mut_present))]

# Calculate the mean of DMS_score for the filtered rows
    mean_delta_G_every_but_mut = only_rows_withOUT_mut['delta G'].mean()
#----------------
    difference_means_delta_G = mean_delta_G_every_but_mut - mean_delta_G_only_mut
    differences_delta_G_list.append(difference_means_delta_G)

all_differences_delta_G_means = pd.DataFrame({'Difference dG': differences_delta_G_list}, index=all_possible_mutations)

# the better the stability of the mutation
#-> difference: WITHOUT - WITH, better mutations have higher difference scores

#ranking delta G
combined_difference_dG_wie_oft_mut = pd.concat([all_differences_delta_G_means, df_wie_oft_muts_insg], axis=1)
combined_difference_dG_wie_oft_mut.columns = ['Difference dG', 'wie oft kommt mut insg vor']

ranking6 = combined_difference_dG_wie_oft_mut.sort_values(by='Difference dG', ascending=False)

ranking6 = ranking6.drop(ranking6[~(ranking6['wie oft kommt mut insg vor'] >= 20)].index)

#-----RANKING 7:
#-> delta G differences combined with ranking 5b (weighted variances with fscore_means)

list_ranking7 = []
for i in all_possible_mutations:

    score_ranking7 = (all_differences_means.loc[i].values[0]) * df_wie_oft_muts_insg.loc[i].values[0] * all_differences_delta_G_means.loc[i].values[0] / (weighted_variances_try_df.loc[i].values[0] )
    list_ranking7.append(score_ranking7)
ranking7_unsorted = pd.DataFrame(list_ranking7, index=all_possible_mutations, columns=['ranking7_score'])
ranking7 = ranking7_unsorted.sort_values(by='ranking7_score', ascending= False)

