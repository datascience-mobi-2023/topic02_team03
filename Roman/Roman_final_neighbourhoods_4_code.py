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
Roman_2_filtered = Roman_2.drop_duplicates()
Roman_3 = pd.concat([Neighbourhood_results_df, Roman_2_filtered], ignore_index=True).drop_duplicates()
Roman_5 = Roman_3.drop(["pKx3"], axis=1)

# First preliminary version to work with.
# It contains all single mutants with non-duplicated frame-shifts along the respective mutated sequence.
# Rows 0-231 are the unmutated reference sequence as neighbourhoods.

# Now it was quickly tested, whether the results make sense based on some easy calculations, and they do.
# It was manually tested, whether the results make sense, in particular the amount of rows in the dataframe.
# By calculating what the code does, the same amount of rows was calculated as are present.

# Now the dataframe will undergo some severe changes until its finalized version:


Roman_6 = pd.DataFrame(columns=["Neighbourhood"])
columns_to_copy = ["Molecular Weight", "Residue Weight", "pKa1", "pKb2", "pl4", "H", "VSC", "P1", "P2", "SASA", "NCISC"]

for index, row in Roman_5.iterrows():
    neighbourhood = row["Neighbourhood"] #Neighbourhood-Wert jeder Zeile nehmen
    neighbourhood_length = len(neighbourhood)

    for i in range(len(unmutated_sequence) - 6): #entfernt wieder 6 Positionen damit es aufkommt und wandert über alle möglichen Positionen
        mismatch_found = False #um zu tracken, ob es einen mismatch gibt

        for j in range(neighbourhood_length): #läuft durch alle Positionen eines neighbourhoods
            if j == 3: #aktuelle Position soll die vierte sein
                if neighbourhood[j] == unmutated_sequence[i+j]: #Prüft ob die vierte AS die Position der original-Sequenz matched
                    mismatch_found = True #wenn match, dann True
                    break # Bricht den loop, wenn die vierte Position matched
            else:
                if neighbourhood[j] != unmutated_sequence[i+j]: #prüft ob die die NICHT vierten AS die Position der original-Sequenz NICHT matched
                   mismatch_found = True #Wenn eine AS außerhalb der vierten Position nicht matched, dann True
                   break # Bricht den loop, wenn eine andere Position außer der vierten nicht matched

        if not mismatch_found: #wenn kein mismatch außer in der vierten Position
            Roman_6 = pd.concat([Roman_6, pd.DataFrame({"Neighbourhood": [neighbourhood]})], ignore_index=True) #wenn kein mismatch, dann Neighbourhood dem neuen Dataframe hinzufügen
            break #Weil ein match gefunden wurde, abbrechen, weil es keinen zweiten Match geben kann

Roman_6 = Roman_6.join(Roman_5[columns_to_copy]) #gibt mir meine alten Spalten dazu


Roman_7 = pd.DataFrame(columns=["Neighbourhood"] + columns_to_copy)

rows = []
for _, row in Roman_6.iterrows():
    neighbourhood = row["Neighbourhood"]
    unmutated_row = neighbourhood[:3] + unmutated_sequence[3] + neighbourhood[4:]
    properties = {col: row[col] for col in columns_to_copy}
    mutated_df = pd.DataFrame(data=[[neighbourhood] + list(row[columns_to_copy])], columns=["Neighbourhood"] + columns_to_copy)
    unmutated_df = pd.DataFrame(data=[[unmutated_row] + list(properties.values())], columns=["Neighbourhood"] + columns_to_copy)
    rows.append(mutated_df)
    rows.append(unmutated_df)

Roman_7 = pd.concat(rows, ignore_index=True)

#klappt mit fehler aber concat
#Fehler, der nimmt immer die 3. Position der unmutierten Sequenz absolut, aber nicht relativ in dem Bereich, wo die neighbourhood ist

property_names = ["Molecular Weight", "Residue Weight", "pKa1", "pKb2", "pl4", "H"]  # List of property names to calculate

# Calculate the property values for each row in Roman_7 and update the corresponding columns
Roman_8 = Roman_7.copy()  # Create a copy of Roman_7
for property_name in property_names:
    Roman_8[property_name] = Roman_8["Neighbourhood"].apply(lambda x: sum(ASE_clear[property_name].loc[ASE_clear["Letter"] == aa].values[0] for aa in x))
    #lambda definiert schnell eine Funktion ohne es vorher mit "def" machen zu müssen. lambda argument: expression
       #This function calculates the sum of specific property values for each letter in the "Neighbourhood" value.
    #loc sucht den richtigen Buchstaben und nimmt davon die Werte der Spalten
    #values 0 nimmt den ersten Wert

#Berechnet mir alle Werte für die neu eingetragenen unmutierten Sequenzen und trägt sie ein.
#Roman_8 ist die Version, woran man alles sehen kann.
#Alle Mutationen vor AS 4 und nach AS 235 können nicht angezeigt werden, weil sie nie in der Mitte sein können.

sequence = "MSKGEEL"
molecular_weight = sum(ASE_clear["Molecular Weight"].loc[ASE_clear["Letter"] == aa].values[0] for aa in sequence)
print("Molecular Weight:", molecular_weight)
#Zum Überprüfen der Berechnung

# Create a new DataFrame Roman_9 as a copy of Roman_8
Roman_9 = Roman_8.copy()

# Filter out mutations where 'n' in 'XnY' mutation name is 3
filtered_mutations = [mutant for mutant in single_mutations_df['mutant'] if not (len(mutant) == 3 and mutant[1] == '3')]

# Clone each mutant and write twice consecutively in a new list
mutants_cloned = [mutant for mutant in filtered_mutations for _ in range(2)]

# Truncate the list of mutants to match the length of the DataFrame
mutants_cloned = mutants_cloned[:len(Roman_9)]

# Insert the 'Mutation' column left of the 'Neighbourhood' column
Roman_9.insert(Roman_9.columns.get_loc('Neighbourhood'), 'Mutation', mutants_cloned)

# Create a copy of the original DataFrame
Roman_10 = Roman_9.copy()

# Create a list to store the modified mutation names
modified_mutations_name = []
modified_neighbourhood = []

for index, row in Roman_9.iterrows():
    mutation = row['Mutation']
    neighbourhood = row['Neighbourhood']

    # Extract the X character from the mutation name
    X = mutation[0]

    # Modify the 4th position of the neighbourhood sequence based on X for odd rows
    if index % 2 != 0:
        modified_neighbourhood.append(neighbourhood[:3] + X + neighbourhood[4:])
    else:
        modified_neighbourhood.append(neighbourhood)

    # Append '-unmut' to every other mutation name
    if index % 2 != 0:
        modified_mutations_name.append(mutation + '-unmut')
    else:
        modified_mutations_name.append(mutation)

# Update the 'Mutation' and 'Neighbourhood' columns with the modified values
Roman_10['Mutation'] = modified_mutations_name
Roman_10['Neighbourhood'] = modified_neighbourhood



property_names = ["Molecular Weight", "Residue Weight", "pKa1", "pKb2", "pl4", "H", "VSC", "P1", "P2", "SASA", "NCISC"]  # List of property names to calculate

# Calculate the property values for each row in Roman_7 and update the corresponding columns
Roman_11 = Roman_10.copy()  # Create a copy of Roman_7
for property_name in property_names:
    Roman_11[property_name] = Roman_11["Neighbourhood"].apply(lambda x: sum(ASE_clear[property_name].loc[ASE_clear["Letter"] == aa].values[0] for aa in x))
    #lambda definiert schnell eine Funktion ohne es vorher mit "def" machen zu müssen. lambda argument: expression
       #This function calculates the sum of specific property values for each letter in the "Neighbourhood" value.
    #loc sucht den richtigen Buchstaben und nimmt davon die Werte der Spalten
    #values 0 nimmt den ersten Wert

Roman_11
#Berechnet mir alle Werte für die neu eingetragenen unmutierten Sequenzen und trägt sie ein. Macht es nochmal, um sicherzugehen
#RECHENWERTE STIMMEN!!!!!
#Roman_11 ist die Version, woran man alles sehen kann.
#Alle Mutationen vor AS 4 und nach AS 235 können nicht angezeigt werden, weil sie nie in der Mitte sein können.

#Roman_11 nach gewünschter Mutation filtern


mutation_name = input("Enter the mutation name: ")
#öffnet ein kleines Fenster zum Eintippen

what_I_see_df = Roman_11[Roman_11['Mutation'].str.contains(mutation_name)]
print(what_I_see_df)
#Funktioniert, muss nur den Code vorher fixen um richtige Werte zu erzeugen

