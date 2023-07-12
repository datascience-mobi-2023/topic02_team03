import pandas as pd
import seaborn as sns
import numpy as np

original_dms_data = pd.read_csv("/Users/liza/Desktop/Bioinfo Project/DMS_data/GFP_AEQVI_Sarkisyan_2016.csv")
# alle Daten einlesen, komplette Tabelle
only_smal_dms_dat = original_dms_data['mutant']
# nur die Spalte mit den Definitionen der Mutanten
zweite_Zeile = only_smal_dms_dat.iloc[3]
#auf eine Zeile zugreifen in der datei mit nur einer Spalte


list_mut_count_in_progress = []
for i in range(len(only_smal_dms_dat)):
    list_mut_count_in_progress.append(only_smal_dms_dat.iloc[i].count(':'))
list_mut_count = np.array(list_mut_count_in_progress)
print(list_mut_count +1)


sns.histplot(data=list_mut_count + 1, discrete = True)
