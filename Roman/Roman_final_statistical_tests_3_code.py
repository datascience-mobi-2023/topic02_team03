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

# Making boxplots to depict the distribution of DMS-scores according to the position of the mutation.

import matplotlib.pyplot as plt
import seaborn as sns

grouped_by_position = work_with_df_pos_new_AA_dms_score.groupby('Position')
dms_scores_per_group_pos_list = []

# Iterates over each group and extracts the DMS scores
for position, group in grouped_by_position:
    dms_scores_per_group_pos_list.append(group['DMS-score'])

sns.boxplot(x='Position', y='Fitness_Score', data=work_with_df_pos_new_AA_dms_score)
plt.xlabel('Position')
plt.ylabel('Fitness Score')
plt.title('Distribution of DMS Scores by Mutation Position')
plt.xticks(range(0, len(grouped_by_position), 10))
plt.gcf().set_size_inches(25, 12)
plt.show()

# This plot is very incomprehensible and confusing, but is meant to give a small overview before starting with statistical tests.

# Now the goal is to understand which parameter (position or new amino acid) contributes the most to the effect on the dms-score.
# Firstly, an analysis of variance (ANOVA) was performed.

from scipy.stats import f_oneway

# Performs the ANOVA
statistic, p_value = f_oneway(*dms_scores_per_group_pos_list)
# dms_scores_per_group was defined up above as the data grouped by position and extracting the corresponding dms-scores

print("ANOVA")
print(f"Statistic: {statistic}")
print(f"P-value: {p_value}")

# In the case of a one-way ANOVA, the test statistic follows an F-distribution.
# The F-statistic represents the ratio of the between-group variability to the within-group variability.
# A larger F-statistic indicates a greater difference between the groups' means relative to the variability within each group.
# Extension of t-test for independent samples and checks if there are statistically significant differences between >2 groups
# ANOVA without repeated measurements because the groups are independent
# Independent variable = position
# Dependent variable = dms-score
# H0 hypothesis = No difference between the means of the individual groups
# H1 hypothesis = There is a difference between at least 2 groups

# The same thing was done for the new amino acid

# Grouping the data by new amino acid
grouped_by_amino_acid = work_with_df_pos_new_AA_dms_score.groupby('new_AA')

# Creates an empty list to store the DMS scores for each group
dms_scores_per_group_new_AA_list = []

# Iterates over each group and extracts the DMS scores
for position, group in grouped_by_amino_acid:
    dms_scores_per_group_new_AA_list.append(group['DMS-score'])

# Performs the ANOVA
statistic, p_value = f_oneway(*dms_scores_per_group_new_AA_list)

# Print the test results
print("ANOVA")
print(f"Statistic: {statistic}")
print(f"P-value: {p_value}")

# The following code shows more information about the results of the ANOVA
# This incorporates all 3 factors (new amino acid, position and dms-score)
import statsmodels.api as sm
from statsmodels.formula.api import ols


# Performs the ANOVA
model_ANOVA_3_components = ols('DMS-score ~ Position + new_AA', data=work_with_df_pos_new_AA_dms_score).fit()
anova_table_3_components = sm.stats.anova_lm(model_ANOVA_3_components)

# Print the ANOVA table
print('ANOVA:')
print(anova_table_3_components)

# df = degrees of freedom
# sum_sq = The sum of squares quantifies the variation explained by each factor. Higher values indicate more significant effects.
# mean_sq = It is obtained by dividing the sum of squares by the respective degrees of freedom.
# The mean square represents the average variation explained by each factor.
# F = It is the ratio of the mean square of each factor to the mean square of the residual.
# The F-value indicates the significance of the factor in explaining the variation in the DMS score.
# Higher F-values suggest a more significant effect.
# PR(>F) = A smaller p-value suggests stronger evidence against the null hypothesis and indicates a significant effect.

# The factor "position" seems to have significant influence. Therefor, there seem to be significant differences in the dms-score between different mutated positions.
# The new amino acid also seems to have a significant influence.
# The residual value represents the unexplained variance of the data set, that cannot be described by position or new amino acid.
# There seems to be a correlation between the dms-score and the position or new amino acid.
# However, it is very complicated to analyse the correlation, as the three components cannot directly be compared, because the new amino acid is non-numerical.
# This difficulty is further increased by not having all mutational combinations in our data.


# Next, an eta-squared test was performed to analyze the effect size of each parameter:

# Results from the previous ANOVA
sum_sq_position = anova_table_3_components.loc['Position', 'sum_sq']
SST = anova_table_3_components['sum_sq'].sum()

sum_sq_new_as = anova_table_3_components.loc['new_AA', 'sum_sq']
SST = anova_table_3_components['sum_sq'].sum()

eta_squared_position = sum_sq_position / SST
print(eta_squared_position)

eta_squared_new_AA = sum_sq_new_AA / SST
print(eta_squared_new_AA)

# Effect Size Calculation: Common effect size measures for ANOVA include eta-squared (η²) or partial eta-squared (η²p).
# These measures indicate the proportion of variance in the DMS score that can be explained by the mutated position or the new amino acid.
# The higher the eta-squared value, the more influence the mutated position has on the DMS score.

# 58% of the variance can be explained by the position of the mutation.
# There seems to be a strong relation between mutated position and dms-score.
# Only 7% of the variance can be explained by the new amino acid, so it seems to not have that big of an impact.
# The rest of the variance is in the residual.


from statsmodels.stats.multitest import multipletests

# Extracts the p-values from the ANOVA results
p_values_from_ANOVA = anova_table_3_components['PR(>F)']

# Applies the Bonferroni correction
corrected_p_values_from_ANOVA = multipletests(p_values_from_ANOVA, method='bonferroni')[1]

# Adds the corrected p-values to the ANOVA table
anova_table_3_components['Corrected P-value'] = corrected_p_values_from_ANOVA

print(anova_table_3_components)

# The correction did not have a lot of influence, as the results are still in the same order of e-82 (position) and e-23 (new_AA).
# The Bonferroni correction is a conservative correction method that controls the family-wise error rate (FWER) by dividing the desired significance level (e.g., 0.05) by the number of comparisons.

# It is important to mention here, that a requirement for ANOVA and thus the reliability of the results regarding the eta-squared test is a normally distributed dataset.
# Therefore, now it will be tested whether our data follows a normal distribution (hint, it does not).


from scipy.stats import shapiro

# Shapiro-Wilk-Test
stat, p = shapiro(work_with_df_pos_new_AA_dms_score["DMS-score"])
alpha = 0.05
# The alpha value was set to 5% as a threshold value

if p > alpha:
    print("Data is normally distributed (Do not discard H0")
    # H0 = the sample came from a normally distributed population
else:
    print("Data is not normally distributed (Discard H0)")
print()
print(p)

# Furthermore, a graphical test in the form of a Q-Q-plot was performed to check for normal distribution.
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import probplot

# Extracts the column with the dms-scores
q_q_dms_scores = work_with_df_pos_new_AA_dms_score['DMS-score'].values

# Creates the Q-Q-plot
probplot(q_q_dms_scores, dist="norm", plot=plt)


plt.title("Q-Q-plot of the DMS-scores")
plt.xlabel("Theoretical quantile")
plt.ylabel("Data quantile")

# Plot anzeigen
plt.show()

#-------------------------------------------------------------------------------------------------
# Next, a Mann-Whitney-U-test was performed
# This test is a nonparametric test.
# The null hypothesis is that, for randomly selected values X and Y from two populations, the probability of X being greater than Y is equal to the probability of Y being greater than X.
# Therefore it checks, whether the probabilities of X or Y to occur, are the same or close to each other.
# It applies rank sums (aka Wilcoxon rank-sum test)

import scipy.stats as stats

grouped_position_MWU = work_with_df_pos_new_AA_dms_score.groupby('Position')
grouped_AA_MWU = work_with_df_pos_new_AA_dms_score.groupby("new_AA")


for position, group_position_loop_MWU in grouped_position_MWU:
    for AA, group_AA_loop_MWU in grouped_AA_MWU:
        dms_scores_position = group_position_loop_MWU['DMS-score'].values
        dms_scores_new_AA = group_AA_loop_MWU['DMS-score'].values

        statistic, p_value = stats.mannwhitneyu(dms_scores_position, dms_scores_new_AA, alternative='two-sided')

        alpha = 0.05

        print(f"Position: {position}, amino acid: {AA}")
        print(f"Mann-Whitney-U-Test")
        print(f"Test statistic: {statistic}")
        print(f"P-value: {p_value}")

        if p_value < alpha:
            print("There is a significant difference.")
        else:
            print("There is no significant difference.")

        print()

        # This code compares all the data from the parameters individually and one by one for now.

# MUSS NOCH KORRIGIERT WERDEN BZW DIE VARIABLEN ANDERS UND AUCH ALLES MUSS NOCH DURCHGELAUFEN WERDEN
import pandas as pd
import scipy.stats as stats


#Ich will meine Werte sortieren
def position_sort(position):
    if position.isdigit():
        return int(position)
    else:
        return position


#Leere Liste
results = []
#Gruppierung basierend auf Position und Aminosäure
grouped_position = Roman_1.groupby('Position')
grouped_AS = Roman_1.groupby("New_AS")

# Konvertiert die Positionsspalte in numerische Werte
Roman_1['Position'] = pd.to_numeric(Roman_1['Position'])

#Whitney-Test für alle Kombinationen von Position und AS durchgeführt --> Signifikanter Unterschied in den Fitness-Scores zwischen verschiedenen Positionen und AS?
for position, group_position in grouped_position:
    for AS, group_AS in grouped_AS:
        fitness_scores_position = group_position['Fitness_Score'].values
        fitness_scores_AS = group_AS['Fitness_Score'].values

        statistic, p_value = stats.mannwhitneyu(fitness_scores_position, fitness_scores_AS, alternative='two-sided')
        #Alpha-Wert selbst gesetzt
        alpha = 0.05

        result = {'Position': position, 'Aminosäure': AS, 'Teststatistik': statistic, 'P-Wert': p_value,
                  'Signifikanz': "Ja" if p_value < alpha else "Nein"}

        results.append(result)
#Ergebnisse ausspucken lassen
summary_df = pd.DataFrame(results)
#Sortieren nach aufsteigender Positionsnummer
summary_df = summary_df.sort_values('Position', ascending=True)
print(summary_df)

#Wenn p-Wert <0,05 (alpha), dann signifikant --> Fitness-Score zwischen den Gruppen nicht gleich

##Wenn der p-Wert klein ist (typischerweise kleiner als 0,05), deutet dies darauf hin, dass es einen statistisch signifikanten Unterschied in den Fitness-Scores zwischen den Positionen und Aminosäuren gibt.

##Wenn der p-Wert größer als 0,05 ist, wird die Nullhypothese beibehalten, was bedeutet, dass kein statistisch signifikanter Unterschied in den Fitness-Scores zwischen den Positionen und Aminosäuren vorliegt.

#Verwendet für große Stichproben und nicht normalverteilt

#Der Mann-Whitney-U-Test wird verwendet, um zu prüfen, ob es statistisch signifikante Unterschiede in den Verteilungen der Fitness-Scores zwischen den Positionen und Aminosäuren gibt.

#Nullhypothese: Es gibt keinen signifikanten Unterschied in den Fitness-Scores zwischen den verschiedenen Positionen und Aminosäuren.
