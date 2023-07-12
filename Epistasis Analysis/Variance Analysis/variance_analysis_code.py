import numpy as np
import pandas as pd

#%%
original_dms_data = pd.read_csv("/Users/liza/Desktop/Bioinfo Project/DMS_data/GFP_AEQVI_Sarkisyan_2016.csv")
only_scores_column = original_dms_data['DMS_score']
df_only_scores_column = pd.DataFrame(only_scores_column)
only_mutant_column = original_dms_data['mutant']
original_dms_data_col = original_dms_data
only_mutants = original_dms_data["mutant"].to_frame()
original_dms_data_col[['m1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11', 'm12', 'm13', 'm14', 'm15']] = original_dms_data_col['mutant'].str.split(':', 15, expand=True)

#%%
list_mut_count_in_progress = []
for i in range(len(only_mutant_column)):
    list_mut_count_in_progress.append(only_mutant_column.iloc[i].count(':'))
list_mut_count_prae = np.array(list_mut_count_in_progress)
list_mut_count = (list_mut_count_prae + 1)
df_mutation_counts = pd.DataFrame(list_mut_count)
#%%

scores_plus_counts = df_only_scores_column.join(df_mutation_counts)
scores_plus_counts.columns = ['fscores', 'mutation_count']
#%%
#the fscores sorted after mut_count
one_mutation_fscores = scores_plus_counts[(scores_plus_counts['mutation_count'] == 1)]
two_mutation_fscores = scores_plus_counts[(scores_plus_counts['mutation_count'] == 2)]
three_mutation_fscores = scores_plus_counts[(scores_plus_counts['mutation_count'] == 3)]
four_mutation_fscores = scores_plus_counts[(scores_plus_counts['mutation_count'] == 4)]
five_mutation_fscores = scores_plus_counts[(scores_plus_counts['mutation_count'] == 5)]
six_mutation_fscores = scores_plus_counts[(scores_plus_counts['mutation_count'] == 6)]
seven_mutation_fscores = scores_plus_counts[(scores_plus_counts['mutation_count'] == 7)]
eight_mutation_fscores = scores_plus_counts[(scores_plus_counts['mutation_count'] == 8)]
nine_mutation_fscores = scores_plus_counts[(scores_plus_counts['mutation_count'] == 9)]
ten_mutation_fscores = scores_plus_counts[(scores_plus_counts['mutation_count'] == 10)]
eleven_mutation_fscores = scores_plus_counts[(scores_plus_counts['mutation_count'] == 11)]

# nur die functional scores zu den jeweiligen mutations counts
#%% md

#%%
#example variances
#print(np.var(one_mutation_fscores))
#print(np.var(two_mutation_fscores))
#print(np.var(three_mutation_fscores))
#%%
import matplotlib.pyplot as plt

#boxplots to show the distribution of the fscores of each group
#plt.boxplot([one_mutation_fscores['fscores'], two_mutation_fscores['fscores'],three_mutation_fscores['fscores'], four_mutation_fscores['fscores'], five_mutation_fscores['fscores'],six_mutation_fscores['fscores'], seven_mutation_fscores['fscores'], eight_mutation_fscores['fscores'], nine_mutation_fscores['fscores'], ten_mutation_fscores['fscores'], eleven_mutation_fscores['fscores']])

#plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], ['1','2','3','4','5','6','7','8','9','10','11'])
#plt.ylabel('functional scores')

#%%
# mann- whitney-u-test
from scipy.stats import mannwhitneyu
statistic_mwu_78, p_value_mwu_78 = mannwhitneyu(eight_mutation_fscores, seven_mutation_fscores)
#if p_value[0] > 0.01:
    #print('kein signifikanter Unterschied')
#else:
    #print('signifikanter Unterschied')

#%%
#generate a more convenient dataframe
working_dataframe_prae = pd.concat([original_dms_data_col, df_mutation_counts], axis="columns")
working_dataframe = working_dataframe_prae.drop(['mutant', 'mutated_sequence', 'DMS_score_bin'], axis=1)
working_dataframe.rename(columns={working_dataframe.columns[16]: 'mut_count'}, inplace=True)

#%%
#-> all_possible_mutations
from pandas import unique
working_dataframe_only_ms = working_dataframe.loc[:, ["m1", "m2", "m3", 'm4', 'm5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11', 'm12', 'm13', 'm14', 'm15']]
#%%
all_possible_mutations = working_dataframe_only_ms.values.flatten().tolist()
all_possible_mutations = list(set(all_possible_mutations))

while None in all_possible_mutations:
    all_possible_mutations.remove(None)

#%%
#existence of mutations (columns) in mutants (rows), boolians
list_of_dfs = []


for i in all_possible_mutations:
    new_column_name = f'{i}'
    new_column_values = [only_mutants_list.str.contains(i, regex= False)]
    new_df = pd.DataFrame({new_column_name: new_column_values})
    new_df_exploded = new_df.explode(new_column_name)
    list_of_dfs.append(new_df_exploded)

result_how_often = pd.concat(list_of_dfs, axis=1)
result_how_often = result_how_often.reset_index(drop=True)

## result_how_often.to_csv('dataframe_mutanten_Mutationen.csv', index=True)
#%%
#more convenient dataframe with only fscore and mutcount
count_fscore_frame = working_dataframe[['DMS_score', 'mut_count']]

#%%
#for all mutations, for all mutcounts
import matplotlib.pyplot as plt


fig, axes = plt.subplots(nrows=5, ncols=3, figsize=(19, 12))  #plots
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


    how_many_for_variance = []

    for i in all_possible_mutations:
        mut_count_fscore = count_fscore_frame.loc[result_how_often[i] == True]
        fscore_mut = mut_count_fscore['DMS_score'].loc[mut_count_fscore['mut_count'] == j]
        wie_viel_jeweils = len(fscore_mut)
        how_many_for_variance.append(wie_viel_jeweils)

    how_many_for_variance = pd.Series(how_many_for_variance, index=all_possible_mutations)
    how_many_for_variance_df = how_many_for_variance.to_frame()


    how_many_AND_variance_df = pd.concat([how_many_for_variance_df, variance_per_mutant_df], axis = 1)
    how_many_AND_variance_df.columns = ['Anzahl benutzter Werte', 'Varianz']
    how_many_AND_variance_df = how_many_AND_variance_df.dropna()

#scatter plot
    #ax.scatter(how_many_AND_variance_df['Anzahl benutzter Werte'],how_many_AND_variance_df['Varianz'], s = j )
    #ax.set_xlabel('Anzahl benutzter Werte')
    #ax.set_ylabel('Varianz')

    #if "V163A" in how_many_AND_variance_df.index:
    #    ax.scatter(how_many_AND_variance_df['Anzahl benutzter Werte']['V163A'],how_many_AND_variance_df['Varianz']['V163A'], c='red')
    #ax.set_title(f'für {j} Mutationen')



#%% md
#1. Werte pro Mutation mitteln (mittel über Anzahl der Daten und Mittel über Varianz) und die dann ranken
#%%
#build the mean of the varainces over all mutcount per mutations
frame_zum_mitteln_variance = pd.DataFrame(index = all_possible_mutations)
variance_per_mutant_count_list = []

# ATTENTION!! only mutcount from 2 to 8 got used, because everything above has a highly decreasing fscore mean (-> boxplot)
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



# dataframe mit allen varianzen (Zellen) pro alle mutationen (rows) pro alle counts (columns)
#%%
mean_variances_per_mutations = pd.DataFrame(variance_per_mutant_count_df.mean(axis=1, skipna=True), columns=['Mean'])

#%%
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

#%%
mean_how_many_per_mutations = pd.DataFrame(how_many_per_mutant_count_df.mean(axis=1, skipna=True), columns=['Mean'])

#%% md
#TRYING ON ONE MUTATION: calculating the fscore_differences
# (V163A)
#1. fscore mean from all mutants X = mutcount, die V163A beinhalten
#2. fscore mean from all mutants X = mutcount, NICHT V163A beinhalten
#%%
#V163A: mean for dms_score for all mutations containing V163A and a mutcount over 2
index_wann_V163A = result_how_often.loc[result_how_often["V163A"] == True].index

only_V163A_muts = working_dataframe[(working_dataframe['mut_count'] > 2) & (working_dataframe.index.isin(index_wann_V163A))]

# Calculate the mean of DMS_score for the filtered rows
mean_dms_score_only_V163A = only_V163A_muts['DMS_score'].mean()

#%%
#V163A: mean for dms_score for all mutations NOT containing V163A and a mutcount over 2
#-> should be lower than the one containing the mutation -> positive effect on other mutations
index_wann_nicht_V163A = result_how_often.loc[result_how_often["V163A"] == False].index

not_V163A_muts = working_dataframe[(working_dataframe['mut_count'] >2) & (working_dataframe.index.isin(index_wann_nicht_V163A))]

# Calculate the mean of DMS_score for the filtered rows
mean_dms_score_not_V163A = not_V163A_muts['DMS_score'].mean()

#%% md
#--> fscores of mutations WITH V163A are in genereal higher than the ones WITHOUT V163A -> good
#--> comparisment to paper results -> matches

#--> verify the idea with deleterious mutation (Y66C, should be negative according to paper)
#%%
index_wann_Y66C = result_how_often.loc[result_how_often["Y66C"] == True].index

only_Y66C_muts = working_dataframe[(working_dataframe['mut_count'] >2) & (working_dataframe.index.isin(index_wann_Y66C))]

# Calculate the mean of DMS_score for the filtered rows
mean_dms_score_only_Y66C = only_Y66C_muts['DMS_score'].mean()


#%%
index_wann_nicht_Y66C = result_how_often.loc[result_how_often["Y66C"] == False].index

not_Y66C_muts = working_dataframe[(working_dataframe['mut_count'] >2) & (working_dataframe.index.isin(index_wann_nicht_Y66C))]

# Calculate the mean of DMS_score for the filtered rows
mean_dms_score_not_Y66C = not_Y66C_muts['DMS_score'].mean()

#%% md
#--> matches perfectly with data from paper

#-> NOW:
#1. build the difference as a measurement-value for the effect of a mutation


#%%
nur_fscore_mut_count = working_dataframe.loc[:, ["DMS_score", "mut_count"]]

#%%
#takes a very long time to load
#result: difference for each mutation
differences_list = []

for i in all_possible_mutations:
    index_when_mut_present = result_how_often.loc[result_how_often[i] == True].index

    only_rows_with_mut = nur_fscore_mut_count [(nur_fscore_mut_count ['mut_count'] >2) & (nur_fscore_mut_count .index.isin(index_when_mut_present))]

# Calculate the mean of DMS_score for the filtered rows
    mean_dms_score_only_mut = only_rows_with_mut['DMS_score'].mean()
#-------------
    index_when_not_mut_present = result_how_often.loc[result_how_often[i] == False].index

    only_rows_withOUT_mut = working_dataframe[(working_dataframe['mut_count'] >2) & (working_dataframe.index.isin(index_when_not_mut_present))]

# Calculate the mean of DMS_score for the filtered rows
    mean_dms_score_every_but_mut = only_rows_withOUT_mut['DMS_score'].mean()
#----------------
    difference_means = mean_dms_score_only_mut - mean_dms_score_every_but_mut
    differences_list.append(difference_means)

all_differences_means = pd.DataFrame({'Difference': differences_list}, index=all_possible_mutations)

#%%
condition = all_differences_means['Difference'] > 0

# count of rows that have a difference > 0 -> positive effect
num_rows_matching_condition = len(all_differences_means[condition])

#%%
#again variance + count calculation
frame_zum_mitteln_variance = pd.DataFrame(index = all_possible_mutations)
variance_per_mutant_count_list = []

# ACHTUNG: es werden nur counts von 2 bis 7 einbezogen weil die mit mehr sowieso "kaputt" sind!!!
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
mean_variances_per_mutations = pd.DataFrame(variance_per_mutant_count_df.mean(axis=1, skipna=True), columns=['Mean'])

#----------------------

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

#%%
list_wie_oft_mut = []
for j in all_possible_mutations:
    matching_indexes = result_how_often.loc[result_how_often[j] == True].index
    wie_oft = len(matching_indexes)
    list_wie_oft_mut.append(wie_oft)
df_wie_oft_muts_insg = pd.DataFrame(list_wie_oft_mut, index=all_possible_mutations)


#%%
combined_differenz_wie_oft_mut= pd.concat([all_differences_means, df_wie_oft_muts_insg], axis=1)
combined_differenz_wie_oft_mut.columns = ['Difference', 'wie oft kommt mut insg vor']

threshold = 200

condition = combined_differenz_wie_oft_mut['wie oft kommt mut insg vor'] < threshold
sorted_df_with_treshold = combined_differenz_wie_oft_mut.drop(combined_differenz_wie_oft_mut.loc[condition].index)

sorted_combined_differenz_wie_oft_mut = sorted_df_with_treshold.sort_values(by='Difference', ascending= False)

destab_rausnehmen_indeces = sorted_combined_differenz_wie_oft_mut.loc[all_differences_means['Difference'] < 0].index
ranking_without_destab = sorted_combined_differenz_wie_oft_mut.drop(destab_rausnehmen_indeces)

#mutations from paper (TOP 15)
TOP_MUTANTS = ['V163A', 'K166Q', 'I171V', 'K113R', 'K214E', 'K156R']

def highlight_top_mutants(row):
    color = 'red' if row.name in TOP_MUTANTS else 'black'
    return ['color: {}'.format(color)] * len(row)

styled_ranking_without_destab = ranking_without_destab.style.apply(highlight_top_mutants, axis=1)

#with open('formatted_ranking_without_destab.html', 'w') as file:
    #file.write(styled_ranking_without_destab.render())

#%%
plot_differences = plt.plot(combined_differenz_wie_oft_mut['Difference'],combined_differenz_wie_oft_mut['wie oft kommt mut insg vor'], 'o')
plt.xlabel('fscore Differenz')
plt.ylabel('wie viele Werte existieren')
plt.title('fscore_mean gegen Anzahl der Werte')

