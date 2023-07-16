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

new_AA_of_mutation = []
for index, row in GFP_dataset.iterrows():
    last_character_of_mutation = row["mutant"][-1]
    new_AA_of_mutation.append(last_character_of_mutation)
new_AA_of_mutation_df = pd.DataFrame(new_AA_of_mutation, columns=["new_AA"])


number_mutations = GFP_dataset["mutant"].str.count(":") + 1
number_mutations_singles = number_mutations == 1

single_mutants_only_df = new_AA_of_mutation_df[number_mutations_singles]

dms_score_list = []
for index, row in GFP_dataset.iterrows():
    dms_score = row["DMS_score"]
    dms_score_list.append(dms_score)
dms_score_list_all_mutants = pd.DataFrame(dms_score_list, columns=["DMS-score"])

all_single_dms_scores = dms_score_list_all_mutants[number_mutations_singles]

all_single_dms_scores_with_new_AA = single_mutants_only_df.join(all_single_dms_scores)

mutations_pos_list = []
for index, row in GFP_dataset.iterrows():
    mutations_pos_list_number = row["mutant"][1:-1]
    mutations_pos_list.append(mutations_pos_list_number)
mutations_pos_df = pd.DataFrame(mutations_pos_list, columns=["Position"])

number_mutations = GFP_dataset["mutant"].str.count(":") + 1
number_mutations_singles = number_mutations == 1

single_mutants_with_pos = mutations_pos_df[number_mutations_singles]
mutations_singles_pos_dms_without_new_AA = single_mutants_with_pos.join(all_single_dms_scores)

new_column_for_pos = mutations_singles_pos_dms_without_new_AA["Position"]
work_with_df_pos_new_AA_dms_score = all_single_dms_scores_with_new_AA.join(new_column_for_pos)
work_with_df_pos_new_AA_dms_score = work_with_df_pos_new_AA_dms_score[["Position", "new_AA", "DMS-score"]]

# Creates a separate dataframe to work with later and compare
GFP_dataset["MutationCount"] = GFP_dataset["mutant"].str.count(":") + 1

work_with_df_pos_bin_mutant_dms_score_mut_sequence = GFP_dataset[GFP_dataset["MutationCount"] == 1].copy()

# Drops the "MutationCount" column because its not needed anymore
work_with_df_pos_bin_mutant_dms_score_mut_sequence.drop("MutationCount", axis=1, inplace=True)

# Now the unmutated sequence of our GFP is determined by manually reversing a single mutation:

unmutated_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

# Import a dataframe with the amino acid porperties:
AA_properties_dataset = pd.read_csv(r"C:\Users\roman\Desktop\AS Eigenschaften\aminoacids.csv")

columns_to_delete_AA_properties = ['Name', 'Abbr', 'Molecular Formula', 'Residue Formula', "carbon", "hydrogen", "nitrogen", "oxygen", "sulfur"]
AA_properties_dataset.drop(columns=columns_to_delete_AA_properties, inplace=True)

# Delete the 12th row because it contains a special amino acid that is not present in GFP
row_index_to_delete_AA_properties = 12
AA_properties_work_with = AA_properties_dataset.drop(row_index_to_delete_AA_properties)

# These are the property abbreviations and explanations:

# Name: name of the amino acid.
# Abbr: abbreviation of the amino acid.
# Letter: letter of the amino acid.
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


# Generate mutated sequences and calculate properties for each mutant
import pandas as pd
import re

# Creates a dictionary of property maps for each property
AA_property_maps = {}
for property_name in AA_properties_work_with.columns[1:]:
    AA_property_maps[property_name] = dict(zip(AA_properties_work_with["Letter"].dropna(), AA_properties_work_with[property_name].dropna()))

# Creates a list to store the Neighbourhood_results_list
Neighbourhood_results_list = []

# This code determines the neighbourhood score for every 7-amino-acid window for each property
for i in range(len(unmutated_sequence) - 6):
    neighbourhood = unmutated_sequence[i:i+7]
    # Picks the 7-amino-acid neighbourhood sequence
    Neighbourhood_result = {"Neighbourhood": neighbourhood}
    for property_name, property_map in AA_property_maps.items():
        score = sum(property_map.get(AA, 0) for AA in neighbourhood)
        # Calculates the sum of property values in the neighbourhood
        Neighbourhood_result[property_name] = score
    Neighbourhood_results_list.append(Neighbourhood_result)

# Create a dataframe from the Neighbourhood_results_list
Neighbourhood_results_df = pd.DataFrame(Neighbourhood_results_list)

# This code above creates the unmutated reference for later analysis

work_with_df_pos_bin_mutant_dms_score_mut_sequence["position"] = work_with_df_pos_bin_mutant_dms_score_mut_sequence["mutant"].apply(lambda x: int(re.search(r'\d+', x).group()))

Neighbourhood_extraction_results_first_list = []
for _, row in work_with_df_pos_bin_mutant_dms_score_mut_sequence.iterrows():
    extract_mutation = row["mutant"]
    extract_position = row["position"]
    determine_mutated_sequence = unmutated_sequence[:extract_position - 1] + extract_mutation[-1] + unmutated_sequence[extract_position:]
    # Generates mutated sequence using string slicing so that the mutation is determined and the reference sequence is changed to the mutated sequence.
    # M[1] checks the second character in the mutant name
    # intM[1] turn the second character into an integer --> Important for extraction of the position so he understands all the numbers as one coherent number
    # #int(mutation[1])-1: Subtracting 1 from the position adjusts it to a zero-based index
    # unmutated_sequence[:int(mutation[1])-1]: This retrieves the subsequence of unmutated_sequence from the beginning up to the position of the mutation (exclusive).
    # It captures the amino acids before the mutation position.
    # M[-1] checks the last character of the mutant name
    # unmutated_sequence[int(mutation[1]):]: This retrieves the subsequence of unmutated_sequence from the position of the mutation to the end.
    # It captures the amino acids after the mutation position.
    # +: This concatenates the three parts together, resulting in the complete mutated sequence.
    Neighbourhood_extraction_result = {"Mutation": extract_mutation, "Neighbourhoods": []}

    # Calculate properties for each neighbourhood sequence
    for i in range(len(determine_mutated_sequence) - 6):
        # Cuts off the end of the general sequence, to always have neighbourhoods of 7 amino acids.
        # Iteration over N has to occur with N-1 iteration length
        neighbourhood = determine_mutated_sequence[i:i+7]
        # Retrieve the 7-amino-acid neighbourhood sequence via slicing
        neighbourhood_properties = {"Neighbourhood": neighbourhood}
        # Creates a dictionary to save properties of each neighbourhood
        for property_name, property_map in AA_property_maps.items():
            neighbourhood_properties[property_name] = sum(property_map.get(AA, 0) for AA in neighbourhood)
            # Calculates the property value for the neighbourhood
        Neighbourhood_extraction_result["Neighbourhoods"].append(neighbourhood_properties)
    Neighbourhood_extraction_results_first_list.append(Neighbourhood_extraction_result)

# Convert the "Neighbourhoods" column to strings
Neighbourhood_extraction_results_first_list_modified = [
    {**entry, "Neighbourhoods": str(entry["Neighbourhoods"])} for entry in Neighbourhood_extraction_results_first_list
]

Neighbourhood_extraction_results_df = pd.DataFrame(Neighbourhood_extraction_results_first_list_modified)
Neighbourhood_extraction_results_df_without_duplicates = Neighbourhood_extraction_results_df.drop_duplicates()
Neighbourhood_extraction_results_df_without_duplicates_concatenated = pd.concat([Neighbourhood_results_df, Neighbourhood_extraction_results_df_without_duplicates], ignore_index=True).drop_duplicates()
Neighbourhoods_without_pKx3 = Neighbourhood_extraction_results_df_without_duplicates_concatenated.drop(["pKx3"], axis=1)

# First preliminary version to work with.
# It contains all single mutants with non-duplicated frame-shifts along the respective mutated sequence.
# Rows 0-231 are the unmutated reference sequence as neighbourhoods.

# Now it was quickly tested, whether the results make sense based on some easy calculations, and they do.
# It was manually tested, whether the results make sense, in particular the amount of rows in the dataframe.
# By calculating what the code does, the same amount of rows was calculated as are present.

# Now the dataframe will undergo some severe changes until its finalized version:


columns_for_neighbourhood_to_add = pd.DataFrame(columns=["Neighbourhood"])
columns_to_copy = ["Molecular Weight", "Residue Weight", "pKa1", "pKb2", "pl4", "H", "VSC", "P1", "P2", "SASA", "NCISC"]

unmutated_sequence = str(unmutated_sequence)
for index, row in Neighbourhoods_without_pKx3.iterrows():
    neighbourhood_rows_values = row["Neighbourhood"]
    neighbourhood_rows_values = str(neighbourhood_rows_values)
    neighbourhood_length = len(neighbourhood_rows_values)

    for i in range(len(unmutated_sequence) - 6):
        mismatch_found = False
        # To track if there is a mismatch

        for j in range(neighbourhood_length):
            if j == 3:
                # 4th position
                if neighbourhood_rows_values[j] == unmutated_sequence[i+j]:
                    # Checks if the 4th amino acid matches the unmutated sequence
                    mismatch_found = True
                    break
                    # Breaks the loop if there is a match
            else:
                if neighbourhood_rows_values[j] != unmutated_sequence[i+j]:
                   mismatch_found = True
                   # If a mismatch occurs outside of the 4th position it breaks the loop
                   break

        if not mismatch_found:
            columns_for_neighbourhood_to_add = pd.concat([columns_for_neighbourhood_to_add, pd.DataFrame({"Neighbourhood": [neighbourhood_rows_values]})], ignore_index=True)
            # If there is not a mismatch, it adds the neighbourhood to the dataframe
            break
            # Code has to break because there can only be one match
columns_for_neighbourhood_to_add = columns_for_neighbourhood_to_add.join(Neighbourhoods_without_pKx3[columns_to_copy])


Neighbourhoods_ready_to_improve = pd.DataFrame(columns=["Neighbourhood"] + columns_to_copy)

Neighbourhoods_ready_to_improve_rows = []  # Initialize an empty list to store the rows

# Iterate over your data and add rows to the list
for i in range(len(determine_mutated_sequence) - 6):
    neighbourhood = determine_mutated_sequence[i:i + 7]
    neighbourhood_properties = {"Neighbourhood": neighbourhood}

    for property_name, property_map in AA_property_maps.items():
        neighbourhood_properties[property_name] = sum(property_map.get(AA, 0) for AA in neighbourhood)

    Neighbourhoods_ready_to_improve_rows.append(pd.Series(neighbourhood_properties))

# Check if the rows list is not empty before concatenating
if Neighbourhoods_ready_to_improve_rows:
    Neighbourhoods_ready_to_improve = pd.concat(Neighbourhoods_ready_to_improve_rows, axis=1).T.reset_index(drop=True)
else:
    Neighbourhoods_ready_to_improve = pd.DataFrame()

property_names_to_calculate = ["Molecular Weight", "Residue Weight", "pKa1", "pKb2", "pl4", "H"]

# Calculates the property values for each row in Neighbourhoods_ready_to_improve and update the corresponding columns
Neighbourhoods_improving_step_by_step = Neighbourhoods_ready_to_improve.copy()
for property_name in property_names_to_calculate:
    Neighbourhoods_improving_step_by_step[property_name] = Neighbourhoods_improving_step_by_step["Neighbourhood"].apply(lambda x: sum(AA_properties_work_with[property_name].loc[AA_properties_work_with["Letter"] == aa].values[0] for aa in x))
    # Lambda quickly defines a function. Lambda argument: expression
    # This function calculates the sum of specific property values for each letter in the "Neighbourhood" value.
    # loc searches for the correct letter (amino acid) and retrieves the values from the columns

# This calculates all the values for the unmutated sequences and adds them to the dataframe.
# All mutations before position 4 and after position 235 cannot be incorporated, because they can never be in the middle of the neighbourhood.

Neighbourhood_test_sequence = "MSKGEEL"
molecular_weight = sum(AA_properties_work_with["Molecular Weight"].loc[AA_properties_work_with["Letter"] == aa].values[0] for aa in Neighbourhood_test_sequence)
print("Molecular Weight:", molecular_weight)
# To confirm the validity. Calculating this by hand yields the same result.


Neighbourhoods_improving_step_by_step_without_edge_mut = Neighbourhoods_improving_step_by_step.copy()
# Filters out mutations where 'n' in 'XnY' mutation name is 3
Neighbourhood_filtered_mutations_without_edge = [mutant for mutant in work_with_df_pos_bin_mutant_dms_score_mut_sequence['mutant'] if not (len(mutant) == 3 and mutant[1] == '3')]

# Clones each mutant and write twice consecutively in a new list
Neighbourhood_mutants_cloned = [mutant for mutant in Neighbourhood_filtered_mutations_without_edge for _ in range(2)]

# Truncates the list of mutants to match the length of the DataFrame
Neighbourhood_mutants_cloned = Neighbourhood_mutants_cloned[:len(Neighbourhoods_improving_step_by_step_without_edge_mut)]

# Inserts the 'Mutation' column left of the 'Neighbourhood' column
Neighbourhoods_improving_step_by_step_without_edge_mut.insert(Neighbourhoods_improving_step_by_step_without_edge_mut.columns.get_loc('Neighbourhood'), 'Mutation', Neighbourhood_mutants_cloned)

# Creates a copy of the original DataFrame
Neighbourhoods_improving_step_by_step_penultimate = Neighbourhoods_improving_step_by_step_without_edge_mut.copy()

# Creates a list to store the modified mutation names
Neighbourhoods_modified_mutations_name_list = []
modified_neighbourhood_list = []

for index, row in Neighbourhoods_improving_step_by_step_without_edge_mut.iterrows():
    mutation_for_neighbourhoods = row['Mutation']
    neighbourhood_rows_values = row['Neighbourhood']

    # Extracts the X character from the mutation name
    X_mutation_name = mutation_for_neighbourhoods[0]

    # Modify the 4th position of the neighbourhood sequence based on X in XnY for odd rows
    if index % 2 != 0:
        modified_neighbourhood_list.append(neighbourhood_rows_values[:3] + X_mutation_name + neighbourhood_rows_values[4:])
    else:
        modified_neighbourhood_list.append(neighbourhood_rows_values)

    # Append '-unmut' to every other mutation name
    if index % 2 != 0:
        Neighbourhoods_modified_mutations_name_list.append(mutation_for_neighbourhoods + '-unmut')
    else:
        Neighbourhoods_modified_mutations_name_list.append(mutation_for_neighbourhoods)

# Update the 'Mutation' and 'Neighbourhood' columns with the modified values
Neighbourhoods_improving_step_by_step_penultimate['Mutation'] = Neighbourhoods_modified_mutations_name_list
Neighbourhoods_improving_step_by_step_penultimate['Neighbourhood'] = modified_neighbourhood_list



property_names_for_advanced_neighbourhoods = ["Molecular Weight", "Residue Weight", "pKa1", "pKb2", "pl4", "H", "VSC", "P1", "P2", "SASA", "NCISC"]

# Calculates the property values for each row in Neighbourhoods_ready_to_improve and updates the corresponding columns
Final_neighbourhoods_df_to_work_with = Neighbourhoods_improving_step_by_step_penultimate.copy()
for property_name in property_names_for_advanced_neighbourhoods:
    Final_neighbourhoods_df_to_work_with[property_name] = Final_neighbourhoods_df_to_work_with["Neighbourhood"].apply(lambda x: sum(AA_properties_work_with[property_names_for_advanced_neighbourhoods].loc[AA_properties_work_with["Letter"] == aa].values[0] for aa in x))
%store Final_neighbourhoods_df_to_work_with
# Final_neighbourhoods_df_to_work_with is the final and corrected version of the dataframe.

# Now I want to build a small interface where I can search for a specific mutation.
# The output should be the mutated and the unmutated rows of the dataframe with all the property values calculated.

mutation_name_in_neighbourhoods_specific = input("Enter the mutation name: ")
# This opens a small window to enter the mutation name.
# Convert 'Mutation' column to string

mutation_name_in_neighbourhoods_specific_result = Final_neighbourhoods_df_to_work_with[Final_neighbourhoods_df_to_work_with['Mutation'].str.contains(mutation_name_in_neighbourhoods_specific)]
print(mutation_name_in_neighbourhoods_specific_result)
# Now the specific effects of every amino acid can be analysed and perhaps the cause of a worse dms-score lies in a change of properties.
# This is a bit speculative however. The next step would have to be to perform 3D structure modelling of GFP...

