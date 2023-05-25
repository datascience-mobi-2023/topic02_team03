import numpy as np
import pandas as pd

original_dms_data = pd.read_csv("/Users/liza/Desktop/Bioinfo Project/DMS_data/GFP_AEQVI_Sarkisyan_2016.csv")
print(original_dms_data)

original_dms_data_mask=original_dms_data['DMS_score']>=2      ###als mask definieren was über dem threshhold 2 liegt
filtered_dms_scores = original_dms_data[original_dms_data_mask]   ###nur die Zeilen die in mask definiert wurden
print(filtered_dms_scores)

plt.savefig  # -> in ython anzeigen
plt.show()

##nur noch 30000 Werte anstatt 51000 -> Werte raus die einen DMS_score unter 2 haben, also wahrscheinlich ohne FUnktion sind