import pandas as pd
import numpy as np
#%%
# read dataset
original_dms_data = pd.read_csv('/Users/liza/Documents/Bioinfo Project/DMS_data/AAAA_GFP_dms_data_original_komplett.csv')
# split first column of df into multiple columns
original_dms_data_col = original_dms_data
only_mutants = original_dms_data["mutant"].to_frame()
original_dms_data_col[['m1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11', 'm12', 'm13', 'm14', 'm15']] = original_dms_data_col['mutant'].str.split(':', 15, expand=True)
#%%
# count how many mutations each sequence has
list_mut_count_in_progress = []
for i in range(len(original_dms_data['mutant'])):
    list_mut_count_in_progress.append(original_dms_data['mutant'].iloc[i].count(':'))
list_mut_count_prae = np.array(list_mut_count_in_progress)
list_mut_count = (list_mut_count_prae + 1)
df_mutation_counts = pd.DataFrame(list_mut_count)
#%%
#merge
working_dataframe_prae = pd.concat([original_dms_data_col, df_mutation_counts], axis="columns")
#drop unimportant columns
working_dataframe = working_dataframe_prae.drop(['mutant', 'mutated_sequence', 'DMS_score_bin'], axis=1)
working_dataframe.rename(columns={working_dataframe.columns[16]: 'mut_count'}, inplace=True)

#%%
#all_possible_mutations

from pandas import unique
working_dataframe_only_ms = working_dataframe.loc[:, ["m1", "m2", "m3", 'm4', 'm5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11', 'm12', 'm13', 'm14', 'm15']]
all_possible_mutations = working_dataframe_only_ms.values.flatten().tolist()
all_possible_mutations = list(set(all_possible_mutations))

while None in all_possible_mutations:
    all_possible_mutations.remove(None)
only_mutants_list = only_mutants['mutant']

#%%
#merge
working_dataframe_prae = pd.concat([original_dms_data_col, df_mutation_counts], axis="columns")

working_dataframe = working_dataframe_prae.drop(['mutant', 'mutated_sequence', 'DMS_score_bin'], axis=1)
working_dataframe.rename(columns={working_dataframe.columns[16]: 'mut_count'}, inplace=True)

#%%
##if alla_possible_mutations exist in mutants -> boolians-table

list_of_dfs = []

for i in all_possible_mutations:
    new_column_name = f'{i}'
    new_column_values = [only_mutants_list.str.contains(i, regex= False)]
    new_df = pd.DataFrame({new_column_name: new_column_values})
    new_df_exploded = new_df.explode(new_column_name)
    list_of_dfs.append(new_df_exploded)


result_how_often = pd.concat(list_of_dfs, axis=1)
result_how_often = result_how_often.reset_index(drop=True)

print("result_how_often finsished (1)")
## result_how_often.to_csv('dataframe_mutanten_Mutationen.csv', index=True)
#%%
# dataframe aus original machen der nur mutcount und fscore hat
count_fscore_frame = working_dataframe[['DMS_score', 'mut_count']]
#%%
#plots für jeden mutcount varianz gegen how_many

import matplotlib.pyplot as plt


#fig, axes = plt.subplots(nrows=5, ncols=3, figsize=(19, 12))  # Abbildung und Achsenobjekte erstellen
#plt.subplots_adjust(wspace=0.4, hspace=0.6)

for j, ax in zip(range(2, 16), axes.flatten()):
    variance_per_mutant_list = []

    for i in all_possible_mutations:
        mut_count_fscore = count_fscore_frame.loc[result_how_often[i] == True]
        fscore_mut = mut_count_fscore['DMS_score'].loc[mut_count_fscore['mut_count'] == j]
        varianz_mut = fscore_mut.var()
        variance_per_mutant_list.append(varianz_mut)

    variance_per_mutant_series = pd.Series(variance_per_mutant_list, index=all_possible_mutations)
    variance_per_mutant_df = variance_per_mutant_series.to_frame()


#rausfinden wie viele Daten wir jeweils haben zum berechnen
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

#scatter plot erstellen, mit benennungen
    #ax.scatter(how_many_AND_variance_df['Anzahl benutzter Werte'],how_many_AND_variance_df['Varianz'], s = j )
    #ax.set_xlabel('Anzahl benutzter Werte')
    #ax.set_ylabel('Varianz')

    #if "V163A" in how_many_AND_variance_df.index:
        #ax.scatter(how_many_AND_variance_df['Anzahl benutzter Werte']['V163A'],how_many_AND_variance_df['Varianz']['V163A'], c='red')
    #ax.set_title(f'für {j} Mutationen')

print("variances and how_many finished (2)")
#%% md
#VORBEREITUNG RANKINGS: Definitionen etc

#%%
#Varianzen mitteln für jede Mutation, über die Mutationscounts hinweg
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

# Varianz je Mutation je Mutationanzahl

# dataframe mit allen varianzen (Zellen) pro alle mutationen (rows) pro alle counts (columns)
mean_variances_per_mutations = pd.DataFrame(variance_per_mutant_count_df.mean(axis=1, skipna=True), columns=['Mean'])

#%%
print("variances_per_mutant_count_df finished (3)")
#%%

#%%
##mean of the variance weighted -> sum of all mean_variance_permutcount * mean_fitness_score_per_mutcount divided by the sum of all used mean_fitness_score_per_mutcount, because apparently that´s necessary?

weighted_mean_per_mutant_df = pd.DataFrame(index=all_possible_mutations)

for mutation in all_possible_mutations:
    weighted_mean = 0
    total_weight = 0

    for mutation_count in range(2, 15):
        variance = variance_per_mutant_count_df.loc[mutation, mutation_count]

        if not np.isnan(variance):
            mean_fitness_score = mean_fitness_scores.at[mutation_count, 'mean_fitness_score']
            weighted_mean += variance * mean_fitness_score
            total_weight += mean_fitness_score

    if total_weight != 0:
        weighted_mean /= total_weight

    weighted_mean_per_mutant_df.at[mutation, 'weighted_mean'] = weighted_mean
print("weighted variance finished (4)")
#%% md
--> 384 von 1810 sind gleich 0 weil die Mutation nur einmal pro mutcount vorkommt, müssen bei den rankings nicht beachtet werden, zu wenig Werte

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

mean_how_many_per_mutations = pd.DataFrame(how_many_per_mutant_count_df.mean(axis=1, skipna=True), columns=['Mean'])

#%%
combined_means_variance_how_many = pd.concat([mean_variances_per_mutations, mean_how_many_per_mutations], axis=1)
combined_means_variance_how_many.columns = ['mean_variances_per_mutations', 'mean_how_many_per_mutations']

print("ranking prep finished (5)")
#%% md
#----------------RANKING 0 weighted: nur nach Varianz der fscores
#%%
sorted_Ranking0_weighted = weighted_mean_per_mutant_df.sort_values(by='weighted_mean')

#%%
TOP_MUTANTS = ['V163A', 'K166Q', 'V68M', 'E172A', 'A206V', 'T43N', 'H25Q', 'S205T', 'E6K', 'T62S', 'I171V', 'T203I', 'Y39N', 'E111V', 'E32A']

# Funktion zum Formatieren der Zeilen und Hervorheben der Werte in TOP_MUTANTS
#def highlight_top_mutants(row):
 #   color = 'red' if row.name in TOP_MUTANTS else 'black'
    #return ['color: {}'.format(color)] * len(row)

# Anwendung der Formatierungsfunktion auf das gesamte DataFrame
#styled_ranking0_weighted= sorted_Ranking0_weighted.style.apply(highlight_top_mutants, axis=1)

# Den formatierten DataFrame als HTML-Datei speichern
#with open('formatted_ranking0_weighted.html', 'w') as file:
 #   file.write(styled_ranking0_weighted.render())
#%% md




----------------RANKING 2 weighted: nur nach fscore_mean Differenz
#%%
# Create an empty dataframe to store the results
mean_for_differences_with = pd.DataFrame(index=all_possible_mutations, columns=range(2, 16))

for mutation in all_possible_mutations:
    for mutation_count in range(2, 16):

        index_when_mut_present_weighted = result_how_often.loc[result_how_often[mutation] == True].index

        only_rows_with_mut_weighted = nur_fscore_mut_count_weighted[(nur_fscore_mut_count_weighted['mut_count'] == mutation_count) & (nur_fscore_mut_count_weighted.index.isin(index_when_mut_present_weighted))]

        # Calculate the mean of fitness score for the filtered rows
        mean_fitness_withOUT_mut_score = only_rows_with_mut_weighted['DMS_score'].mean()

        # Store the mean fitness score in the result dataframe
        mean_for_differences_with.loc[mutation, mutation_count] = mean_fitness_withOUT_mut_score

#%%
# Assuming your dataframe is named 'data' and 'mean_fitness_scores' is another dataframe with weights

weights = mean_fitness_scores['mean_fitness_score'].values

# Calculate weighted means per row, ignoring NaN values
weighted_means = np.zeros(len(mean_for_differences_with))

for i, row in enumerate(mean_for_differences_with.values):
    non_nan_values = row[~pd.isna(row)]
    non_nan_weights = weights[~pd.isna(row)]

    if len(non_nan_values) > 0:
        weighted_means[i] = np.average(non_nan_values, weights=non_nan_weights)

# Create a new dataframe to store the weighted means
weighted_means_with_df = pd.DataFrame({'Weighted Mean': weighted_means}, index=mean_for_differences_with.index)

print("weighted means WITH finished (6)")
#%%
# Create an empty dataframe to store the results
mean_for_differences_withOUT = pd.DataFrame(index=all_possible_mutations, columns=range(2, 16))

for mutation in all_possible_mutations:
    for mutation_count in range(2, 16):

        index_when_mut_NOT_present_weighted = result_how_often.loc[result_how_often[mutation] == False].index

        only_rows_withOUT_mut_weighted = nur_fscore_mut_count_weighted[(nur_fscore_mut_count_weighted['mut_count'] == mutation_count) & (nur_fscore_mut_count_weighted.index.isin(index_when_mut_NOT_present_weighted))]

        # Calculate the mean of fitness score for the filtered rows
        mean_fitness_withOUT_mut_score = only_rows_withOUT_mut_weighted['DMS_score'].mean()

        # Store the mean fitness score in the result dataframe
        mean_for_differences_withOUT.loc[mutation, mutation_count] = mean_fitness_withOUT_mut_score

#%%
weights = mean_fitness_scores.values.ravel()

weighted_means = np.average(mean_for_differences_withOUT, axis=1, weights= weights)

# Create a new dataframe to store the weighted means
weighted_means_withOUT_df = pd.DataFrame({'Weighted Mean': weighted_means}, index=mean_for_differences_withOUT.index)

print("weighted means withOUT finished (7)")
#%%
all_differences_means_weighted = weighted_means_with_df -weighted_means_withOUT_df

#%% md

#%%
#VORBEREITUNG copy&paste
list_wie_oft_mut = []
for j in all_possible_mutations:
    matching_indexes = result_how_often.loc[result_how_often[j] == True].index
    wie_oft = len(matching_indexes)
    list_wie_oft_mut.append(wie_oft)
df_wie_oft_muts_insg = pd.DataFrame(list_wie_oft_mut, index=all_possible_mutations)

#%%
#VORBEREITUNG copy&paste
#code für ranking aus anderem dokument aber mit den sachen von oben berücksichtigt, alle destab raus
combined_differenz_wie_oft_mut_weighted= pd.concat([all_differences_means_weighted, df_wie_oft_muts_insg], axis=1)
combined_differenz_wie_oft_mut_weighted.columns = ['Difference_weighted', 'wie oft kommt mut insg vor']

#%%
ranking2_weighted = combined_differenz_wie_oft_mut_weighted.sort_values(by='Difference_weighted', ascending= False)

#%%
#TOP_MUTANTS = ['V163A', 'K166Q', 'V68M', 'E172A', 'A206V', 'T43N', 'H25Q', 'S205T', 'E6K', 'T62S', 'I171V', 'T203I', 'Y39N', 'E111V', 'E32A']

# Funktion zum Formatieren der Zeilen und Hervorheben der Werte in TOP_MUTANTS
#def highlight_top_mutants(row):
    #color = 'red' if row.name in TOP_MUTANTS else 'black'
    #return ['color: {}'.format(color)] * len(row)

# Anwendung der Formatierungsfunktion auf das gesamte DataFrame
#styled_ranking2_weighted= ranking2_weighted.style.apply(highlight_top_mutants, axis=1)

# Den formatierten DataFrame als HTML-Datei speichern
#with open('formatted_ranking2_weighted.html', 'w') as file:
    #file.write(styled_ranking2_weighted.render())
print("ranking 2 finished (8)")
#%% md




#----------------RANKING 3 weighted: nach Differenz gewichtet nach Anzahl
#%%
list_ranking3_weighted = []
for i in all_possible_mutations:
    score_ranking3_weighted = all_differences_means_weighted.loc[i].values[0] * df_wie_oft_muts_insg.loc[i].values[0]
    list_ranking3_weighted.append(score_ranking3_weighted)

ranking3_unsorted_weighted = pd.DataFrame(list_ranking3_weighted, index=all_possible_mutations, columns=['ranking3_score_weighted'])
ranking3_weighted = ranking3_unsorted_weighted.sort_values(by='ranking3_score_weighted', ascending= False)

#%%
#TOP_MUTANTS = ['V163A', 'K166Q', 'V68M', 'E172A', 'A206V', 'T43N', 'H25Q', 'S205T', 'E6K', 'T62S', 'I171V', 'T203I', 'Y39N', 'E111V', 'E32A']

# Funktion zum Formatieren der Zeilen und Hervorheben der Werte in TOP_MUTANTS
#def highlight_top_mutants(row):
 #   color = 'red' if row.name in TOP_MUTANTS else 'black'
  #  return ['color: {}'.format(color)] * len(row)

# Anwendung der Formatierungsfunktion auf das gesamte DataFrame
#styled_ranking3_weighted= ranking3_weighted.style.apply(highlight_top_mutants, axis=1)

# Den formatierten DataFrame als HTML-Datei speichern
#with open('formatted_ranking3_weighted.html', 'w') as file:
    #file.write(styled_ranking3_weighted.render())
print("ranking3 weighted finished (9)")
#%% md

#--> bei den Werten wo die Varianz wegen zu wenigen Werten gleich 0 ist, entstehen durch den Bruch INF werte, müssen ignoriert werden
#---> theoretisch später raussortierbar
#---> kein if clause, weil sonst unterschiedliche Berechnung ("unfair")


#----------------RANKING 4 weighted: nach eigenem score1:
#-> score1 = Differenz * 1/Varianz * Anzahl muts
#%%
list_ranking4_weighted  = []
for i in all_possible_mutations:

    score_ranking4_weighted = all_differences_means_weighted .loc[i].values[0] * df_wie_oft_muts_insg.loc[i].values[0] * (1/mean_variances_per_mutations.loc[i].values[0])
    list_ranking4_weighted .append(score_ranking4_weighted )
ranking4_unsorted_weighted  = pd.DataFrame(list_ranking4_weighted , index=all_possible_mutations, columns=['ranking4_score_weighted '])
ranking4_weighted  = ranking4_unsorted_weighted .sort_values(by='ranking4_score_weighted ', ascending= False)

#%%
#TOP_MUTANTS = ['V163A', 'K166Q', 'V68M', 'E172A', 'A206V', 'T43N', 'H25Q', 'S205T', 'E6K', 'T62S', 'I171V', 'T203I', 'Y39N', 'E111V', 'E32A']

# Funktion zum Formatieren der Zeilen und Hervorheben der Werte in TOP_MUTANTS
#def highlight_top_mutants(row):
   # color = 'red' if row.name in TOP_MUTANTS else 'black'
    #return ['color: {}'.format(color)] * len(row)

# Anwendung der Formatierungsfunktion auf das gesamte DataFrame
#styled_ranking4_weighted = ranking4_weighted .style.apply(highlight_top_mutants, axis=1)

# Den formatierten DataFrame als HTML-Datei speichern
#with open('formatted_ranking4_weighted .html', 'w') as file:
    #file.write(styled_ranking4_weighted .render())

#%% md




#----------------RANKING 5 weighted: nach eigenem score2 :
#-> score2 = Differenz * aggregierte Varianz
#-> Varianz: gewichten mit mutcount (aggregierte Varianz)
#(Summe aller Varianzen*1/Anzahl muts)/(Gesamtzahl Mutationen)
###macht keinen sinn, aggregierte Varianz nochmal ? , ist auch nicht gut

#aggregierte Varianz: bezieht streuung über ganzes Datenset mi ein
#gewichtete Varianz: Verhalten innerhalb der Gruppe

#--> NUR DIFFERENCE ist gewichtet miteinbezogen!! Varianz nicht, nicht geeigenet
#%%
list_ranking5_weighted1 = []
for i in all_possible_mutations:

    score_ranking5_weighted1 = (all_differences_means_weighted.loc[i].values[0] * df_wie_oft_muts_insg.loc[i].values[0]* 51714) / np.sum(variance_per_mutant_count_df.loc[i].values)
    list_ranking5_weighted1.append(score_ranking5_weighted1)
ranking5_unsorted_weighted1 = pd.DataFrame(list_ranking5_weighted1, index=all_possible_mutations, columns=['ranking5_score_weighted1'])
ranking5_weighted1 = ranking5_unsorted_weighted1.sort_values(by='ranking5_score_weighted1', ascending= False)
print(ranking5_weighted1)
#%%
#TOP_MUTANTS = ['V163A', 'K166Q', 'V68M', 'E172A', 'A206V', 'T43N', 'H25Q', 'S205T', 'E6K', 'T62S', 'I171V', 'T203I', 'Y39N', 'E111V', 'E32A']

# Funktion zum Formatieren der Zeilen und Hervorheben der Werte in TOP_MUTANTS
#def highlight_top_mutants(row):
 #   color = 'red' if row.name in TOP_MUTANTS else 'black'
   # return ['color: {}'.format(color)] * len(row)

# Anwendung der Formatierungsfunktion auf das gesamte DataFrame
#styled_ranking5_weighted1= ranking5_weighted1.style.apply(highlight_top_mutants, axis=1)

# Den formatierten DataFrame als HTML-Datei speichern
#with open('formatted_ranking5_weighted.html', 'w') as file:
   # file.write(styled_ranking5_weighted1.render())
print("ranking 5 finished (10)")
#%% md





#RANKING 6 weighted:
#-> delta G values ranked on their own

#%%
delta_G_data = pd.read_csv('/Users/liza/Downloads/df_ddG.csv')
#%%
count_fscore_frame['delta G'] = delta_G_data['Score']

#%%

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
print("delta G differences (11)")
#je niedriger dG desto stabiler ist Protein
#Difference "ohne - mit": kleiner besser, wenn mit besser ist als ohne ist difference positiv
#%% md


#RANKING 7 weighted:
#-> delta G Differenzwerte mit Ranking 5 verrechnet
#%%
list_ranking7_weighted = []
for i in all_possible_mutations:

    score_ranking7_weighted = (all_differences_means_weighted.loc[i].values[0] * df_wie_oft_muts_insg.loc[i].values[0]* 51714) / (np.sum(variance_per_mutant_count_df.loc[i].values) * all_differences_delta_G_means.loc[i].values[0])
    list_ranking7_weighted.append(score_ranking7_weighted)
ranking7_unsorted_weighted = pd.DataFrame(list_ranking7_weighted, index=all_possible_mutations, columns=['ranking7_weighted_score'])
ranking7_weighted = ranking7_unsorted_weighted.sort_values(by='ranking7_weighted_score', ascending= False)

#%%
#TOP_MUTANTS = ['V163A', 'K166Q', 'V68M', 'E172A', 'A206V', 'T43N', 'H25Q', 'S205T', 'E6K', 'T62S', 'I171V', 'T203I', 'Y39N', 'E111V', 'E32A']

# Funktion zum Formatieren der Zeilen und Hervorheben der Werte in TOP_MUTANTS
#def highlight_top_mutants(row):
    #color = 'red' if row.name in TOP_MUTANTS else 'black'
    #return ['color: {}'.format(color)] * len(row)

# Anwendung der Formatierungsfunktion auf das gesamte DataFrame
#styled_ranking7_weighted= ranking7_weighted.style.apply(highlight_top_mutants, axis=1)

# Den formatierten DataFrame als HTML-Datei speichern
#with open('formatted_ranking7_weighted.html', 'w') as file:
    #file.write(styled_ranking7_weighted.render())
#%%
