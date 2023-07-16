# team03 of topic02 
DMS analysis of green fluorescent protein (GFP)

This project investigates the impact of singular and multiple mutations on the dms-score of green fluorescent protein (GFP). The analysis begins by examining the effects of mutation position and amino acid substitutions on the dms-scores of single mutants. Additionally, the study explores the physiological properties of neighboring amino acids, to identify causality for a lowered dms-score. Moreover, this project delves into the phenomenon of epistasis in GFP, investigating the interplay of multiple mutations and identifying buried amino acid residues that exhibit enhanced epistatic effects. Furthermore, the influence of the number of mutations on the dms-score is analyzed, leading to the development of weighted rankings to identify mutations with positive epistatic effects. Finally, pyrosetta is utilized to predict the difference in free energy (ΔG) between the unfolded and folded protein. The findings provide valuable insights into the factors influencing the dms-score in GFP mutants and have implications for protein engineering and design.

Single mutation:
 In this part of the project, the analysis focused primarily on single mutations. The main objective was to determine and quantify the relative impact of the mutation position  versus the newly created amino acid on the resulting dms-score. Furthermore, the physiological impact of the change in amino acid to the fluorescence activity was determined   by analyzing amino acid-specific parameters such as the volume of side chains (VSC), net charge index of side chains (NCISC), polarity (P1) and solvent-accessible surface      area (SASA).
 
 How to use the repository:
 1) Go to the directory "single mutations"
 2) 

Epistasis
Additionally, analyses were performed to obtain information about the ability of mutations to stabilize the fitness scores of proteins that already contain other mutations and get quantified data to rank them accordingly. To make the calculations as specific as possible, weighted factors were computed, and the effect of mutation count on protein fitness got included.

Protein stability prediction with PyRosetta 
The stability of proteins is important as it can affect proper folding and functional performance. Protein stability can be characterized by thermodynamic stability, which refers to the resistance against denaturation, as well as kinetic stability, which measures the resistance of a protein against irreversible inactivation (Liu, Xun et al. 2019). Since no experimentally determined parameters capturing stability were available in the dataset used (Sarkisyan, Bolotin et al. 2016), we utilized the software platform PyRosetta to predict a stability parameter, enabling us to investigate the relationship between function and stability in mutated GFP. 
