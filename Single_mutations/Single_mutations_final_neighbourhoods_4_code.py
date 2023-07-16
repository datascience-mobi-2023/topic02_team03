# Import GFP dataset
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

last_character_list = []
for index, row in GFP_dataset.iterrows():
    last_character = row["mutant"][-1]
    last_character_list.append(last_character)
last_character_df = pd.DataFrame(last_character_list, columns=["new_AA"])
# Last character of each row = new amino acid

number_mutations = GFP_dataset["mutant"].str.count(":") + 1
number_mutations_Single = number_mutations == 1
# True are all the rows that only have 1 mutation

single_mutants_df = last_character_df[number_mutations_Single]
# Filters all true rows from both dataframes and creates a new one. The filtering of true values is automatically being transferred.

dms_score_df = []
for index, row in GFP_dataset.iterrows():
    dms_score = row["DMS_score"]
    dms_score_df.append(dms_score)
dms_score_df_all = pd.DataFrame(dms_score_df, columns=["DMS-score"])
# Creates a dataframe with all DMS-scores and the respective experiment number

dms_score_filtered = dms_score_df_all[number_mutations_Single]
# Dataframe with all DMS-scores of the single mutants

dms_score_filtered_new_AA = single_mutants_df.join(dms_score_filtered)
# Combines and creates a new dataframe with the new amino acid and the dms-score

mutations_pos = []
for index, row in GFP_dataset.iterrows():
    mutations_pos_number = row["mutant"][1:-1]
    mutations_pos.append(mutations_pos_number)
mutations_pos_df = pd.DataFrame(mutations_pos, columns=["Position"])
# Removes the first and last character = position

number_mutations = GFP_dataset["mutant"].str.count(":") + 1
number_mutations_Single = number_mutations == 1
# True are all the rows that only contain 1 mutation.

single_mutants_df_pos = mutations_pos_df[number_mutations_Single]
# Creates a new dataframe that has the positions of all single mutants.

mutations_pos_df_with_scores = single_mutants_df_pos.join(dms_score_filtered)
# Combines and creates a new dataframe with the position and the DMS-score (without the new amino acid)

new_column = mutations_pos_df_with_scores["Position"]
Neighbourhoods_1 = dms_score_filtered_new_AA.join(new_column)
Neighbourhoods_1 = Neighbourhoods_1[["Position", "new_AA", "DMS-score"]]
# Creates a new dataframe with position, new amino acid and DMS-score

GFP_dataset["MutationCount"] = GFP_dataset["mutant"].str.count(":") + 1
# Counts the number of mutations in each row

single_mutations_df = GFP_dataset[GFP_dataset["MutationCount"] == 1].copy()
# Filters the DataFrame to select rows with a single mutation

single_mutations_df.drop("MutationCount", axis=1, inplace=True)
# Drops the "MutationCount" column if not needed anymore

sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

AA_properties = pd.read_csv(r"C:\Users\roman\Desktop\AS Eigenschaften\aminoacids.csv")
# Import the amino acid property dataset

columns_to_delete = ['Name', 'Abbr', 'Molecular Formula', 'Residue Formula']
AA_properties.drop(columns=columns_to_delete, inplace=True)
delete_columns = ["carbon", "hydrogen", "nitrogen", "oxygen", "sulfur"]
AA_properties.drop(columns=delete_columns, inplace=True)
# Cleans the dataframe from unnecessary properties

row_index_to_delete = 12
AA_properties_to_work_with = AA_properties.drop(row_index_to_delete)
# Deletes a special amino acid that does not occur in our GFP dataset


# Generate mutated sequences and calculate properties for each mutant
import pandas as pd
import re

property_maps = {}
for property_name in AA_properties_to_work_with.columns[1:]:
    property_maps[property_name] = dict(zip(AA_properties_to_work_with["Letter"].dropna(), AA_properties_to_work_with[property_name].dropna()))
# Creates a dictionary of property maps for each property

results = []
# Creates a list to store the results

for i in range(len(sequence) - 6):
    neighbourhood = sequence[i:i + 7]  # Gets the 7-amino-acid neighbourhood sequence
    result = {"Neighbourhood": neighbourhood}
    for property_name, property_map in property_maps.items():
        score = sum(property_map.get(AA, 0) for AA in
                    neighbourhood)  # Calculates the sum of property values in the neighbourhood
        result[property_name] = score
    results.append(result)
# Calculates the neighbourhood score for every 7-amino-acid window for each property

unmutated_all_neighbourhoods = pd.DataFrame(results)
# This acts as the unmutated reference!

unmutated_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
single_mutations_df["position"] = single_mutations_df["mutant"].apply(lambda x: int(re.search(r'\d+', x).group()))
results = []

for _, row in single_mutations_df.iterrows():
    mutation = row["mutant"]
    position = row["position"]
    mutated_sequence = unmutated_sequence[:position - 1] + mutation[-1] + unmutated_sequence[position:]

    result = {"Mutation": mutation}


    for i in range(len(mutated_sequence) - 6):
        neighbourhood = mutated_sequence[i:i + 7]
        neighbourhood_properties = {"Neighbourhood": neighbourhood}
        for property_name, property_map in property_maps.items():
            neighbourhood_properties[property_name] = sum(property_map.get(AA, 0) for AA in neighbourhood)
        results.append(neighbourhood_properties)

# Cleaning the dataframe
Neighbourhoods_2 = pd.DataFrame(results)

Neighbourhoods_2 = Neighbourhoods_2.drop_duplicates()

Neighbourhoods_2_filtered = Neighbourhoods_2.drop_duplicates()

Neighbourhoods_2_filtered = pd.concat([unmutated_all_neighbourhoods, Neighbourhoods_2_filtered], axis=0, ignore_index=True)

Neighbourhoods_3 = Neighbourhoods_2_filtered.drop_duplicates()

Neighbourhoods_4 = Neighbourhoods_3.drop(["pKx3"], axis=1)

Neighbourhoods_5 = Neighbourhoods_4.drop_duplicates()
# Final cleaned dataframe with all the single mutants and with all non-duplicated frames.
# Rows 0-231 are the unmutated reference.

#Name: name of the amino acid.
#Abbr: abbreviation of the amino acid.
#Letter: letter of the amino acid.
# Molecular Weight: molecular weight.
# Molecular Formula: molecular formula.
# Residue Formula: residue formula.
# Residue Weight: residue weight (-H20)
# pKa1: the negative of the logarithm of the dissociation constant for the -COOH group.
# pKb2: the negative of the logarithm of the dissociation constant for the -NH3 group.
# pKx3: the negative of the logarithm of the dissociation constant for any other group in the molecule.
# pl4: the pH at the isoelectric point.
# H: hydrophobicity.
# VSC: volumes of side chains amino acids.
# P1: polarity.
# P2: polarizability.
# SASA: solvent accesible surface area.
# NCISC: net charge of side chains.

import pandas as pd

Neighbourhoods_6 = pd.DataFrame(columns=["Neighbourhood"])
columns_to_copy = ["Molecular Weight", "Residue Weight", "pKa1", "pKb2", "pl4", "H", "VSC", "P1", "P2", "SASA", "NCISC"]

unmutated_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

for index, row in Neighbourhoods_5.iterrows():
    neighbourhood = row["Neighbourhood"]  # Gets the neighbourhood value of each row
    neighbourhood_length = len(neighbourhood)

    for i in range(
            len(unmutated_sequence) - 6):
        mismatch_found = False  # To track if a mismatch occurs

        for j in range(neighbourhood_length):
            if j == 3:  # Finds the 4th position of each neighbourhood
                if neighbourhood[j] == unmutated_sequence[
                    i + j]:  # Checks, if the 4th position matches the original sequence
                    mismatch_found = True  #If there is a match, then true
                    break  # Breaks the loop if there is a match (I dont want a match!)
            else:
                if neighbourhood[j] != unmutated_sequence[
                    i + j]:  # Checks, if the not-4th positions match the original sequence
                    mismatch_found = True  # If an amino acid matches outside of the 4th position, then true
                    break  # Breaks the loop if any other position but the 4th matches

        if not mismatch_found:  # If there was no mismatch, except for the 4th position
            Neighbourhoods_6 = pd.concat([Neighbourhoods_6, pd.DataFrame({"Neighbourhood": [neighbourhood]})],
                                ignore_index=True)  # If no mismatch, then append the neighbourhood to the new dataframe
            break  #Breaks the code because there can only be 1 match maximum

Neighbourhoods_6 = Neighbourhoods_6.join(Neighbourhoods_5[columns_to_copy])  # Adds the old columns I want to keep

import pandas as pd

Neighbourhoods_7 = pd.DataFrame(columns=["Neighbourhood"] + columns_to_copy)

unmutated_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

rows = []
for _, row in Neighbourhoods_6.iterrows():
    neighbourhood = row["Neighbourhood"]
    unmutated_row = neighbourhood[:3] + unmutated_sequence[3] + neighbourhood[4:]
    properties = {col: row[col] for col in columns_to_copy}
    mutated_df = pd.DataFrame(data=[[neighbourhood] + list(row[columns_to_copy])],
                              columns=["Neighbourhood"] + columns_to_copy)
    unmutated_df = pd.DataFrame(data=[[unmutated_row] + list(properties.values())],
                                columns=["Neighbourhood"] + columns_to_copy)
    rows.append(mutated_df)
    rows.append(unmutated_df)

Neighbourhoods_7 = pd.concat(rows, ignore_index=True)
# This contains an error, where he takes the 3rd position of the original sequence (absolute) and not relative to the neighbourhood how I need it

property_names = ["Molecular Weight", "Residue Weight", "pKa1", "pKb2", "pl4", "H"]  # List of property names to calculate

Neighbourhoods_8 = Neighbourhoods_7.copy()
for property_name in property_names:
    Neighbourhoods_8[property_name] = Neighbourhoods_8["Neighbourhood"].apply(lambda x: sum(AA_properties_to_work_with[property_name].loc[AA_properties_to_work_with["Letter"] == aa].values[0] for aa in x))

# Lambda quickly defines a function without using def
# This function calculates the sum of specific property values for each letter in the "Neighbourhood" value.
# Loc finds the correct character/letter and takes the values from the columns
# Values 0 takes the first value

sequence = "MSKGEEL"
molecular_weight = sum(AA_properties_to_work_with["Molecular Weight"].loc[AA_properties_to_work_with["Letter"] == aa].values[0] for aa in sequence)
print("Molecular Weight:", molecular_weight)
# This is a test calculation to see if the properties of the neighbourhoods are calculated correctly.
# Comparison to manually calculating this yields the same result!

Neighbourhoods_9 = Neighbourhoods_8.copy()

filtered_mutations = [mutant for mutant in single_mutations_df['mutant'] if not (len(mutant) == 3 and mutant[1] == '3')]
# Filters out mutations where 'n' in 'XnY' mutation name is 3

mutants_cloned = [mutant for mutant in filtered_mutations for _ in range(2)]
# Clones each mutant and write twice consecutively in a new list

mutants_cloned = mutants_cloned[:len(Neighbourhoods_9)]
# Truncates the list of mutants to match the length of the DataFrame

Neighbourhoods_9.insert(Neighbourhoods_9.columns.get_loc('Neighbourhood'), 'Mutation', mutants_cloned)
# Inserts the 'Mutation' column left of the 'Neighbourhood' column

Neighbourhoods_10 = Neighbourhoods_9.copy()

modified_mutations_name = []
modified_neighbourhood = []
# Creates a list to store the modified mutation names

for index, row in Neighbourhoods_9.iterrows():
    mutation = row['Mutation']
    neighbourhood = row['Neighbourhood']

    X = mutation[0]
    # Extracts the X character from the mutation name

    if index % 2 != 0:
        modified_neighbourhood.append(neighbourhood[:3] + X + neighbourhood[4:])
    else:
        modified_neighbourhood.append(neighbourhood)
   # Modifies the 4th position of the neighbourhood sequence based on X for odd rows

    if index % 2 != 0:
        modified_mutations_name.append(mutation + '-unmut')
    else:
        modified_mutations_name.append(mutation)
  # Appends '-unmut' to every other mutation name

Neighbourhoods_10['Mutation'] = modified_mutations_name
Neighbourhoods_10['Neighbourhood'] = modified_neighbourhood
# Updates the 'Mutation' and 'Neighbourhood' columns with the modified values

property_names = ["Molecular Weight", "Residue Weight", "pKa1", "pKb2", "pl4", "H", "VSC", "P1", "P2", "SASA",
                  "NCISC"]  # List of property names to calculate

Neighbourhoods_11 = Neighbourhoods_10.copy()
for property_name in property_names:
    Neighbourhoods_11[property_name] = Neighbourhoods_11["Neighbourhood"].apply(
        lambda x: sum(AA_properties_to_work_with[property_name].loc[AA_properties_to_work_with["Letter"] == aa].values[0] for aa in x))

# This calculates all values for the newly inserted unmutated sequences and fills them in (to update the values)
# The values seem to be correct according to manual checking!
# All the mutations bevor position 4 and after 235 cannot be seen in this dataframe, because they can never be in the middle of a neighbourhood of length 7

# Enter a specific mutation in the little window that pops up!

mutation_name = input("Enter the mutation name: ")
Dataframe_for_specific_mutation = Neighbourhoods_11[Neighbourhoods_11['Mutation'].str.contains(mutation_name)]
print(Dataframe_for_specific_mutation)
# works!


