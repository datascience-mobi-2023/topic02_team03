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


#---------------------------------------------------------------------------------------------------------------------------------------
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
        property_score_result[property_name] = score
    Neighbourhood_results_list.append(Neighbourhood_result)

# Create a dataframe from the Neighbourhood_results_list
Neighbourhood_results_df = pd.DataFrame(Neighbourhood_results_list)

# This code above creates the unmutated reference for later analysis

work_with_df_pos_bin_mutant_dms_score_mut_sequence["position"] = work_with_df_pos_bin_mutant_dms_score_mut_sequence["mutant"].apply(lambda x: int(re.search(r'\d+', x).group()))
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
    Neighbourhood_extraction_result = {"Mutation": extract_mutation}

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
        Neighbourhood_extraction_result.append(neighbourhood_properties)

Roman_2 = pd.DataFrame(Neighbourhood_extraction_result)
#Roman_2
# Drop duplicate rows with identical values
Roman_2_filtered = Roman_2.drop_duplicates()
# Print the filtered DataFrame
#(Roman_2_filtered)
Roman_2_filtered = pd.concat([Neighbourhood_results_df, Roman_2_filtered], axis=0, ignore_index=True)
#print(Roman_2_filtered)
Roman_3 = Roman_2_filtered.drop_duplicates()
#Roman_3
Roman_4 = Roman_3.drop(["pKx3"], axis=1)
#Roman_4
Roman_5 = Roman_4.drop_duplicates()
print(Roman_5)
#Endgültiges gereinigtes dataframe mit allen single mutanten mit allen nicht-duplikativen frame-Wanderungen. Zeile 0 bis 231 ist di nicht mutierte Referenz