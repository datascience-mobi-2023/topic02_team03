# team03 of topic02 
##### Roman Kurley, Tianxin Angela Ma, Lisa Duttenhöfer, Rebecca Ress
DMS analysis of green fluorescent protein (GFP)

The final report (GFP_DMS_report_duttenhoefer_kurley_ma_ress) is to be found in the repository, containing all important information and analysis results for the project.

This project investigates the impact of singular and multiple mutations on the dms-score of green fluorescent protein (GFP). The analysis begins by examining the effects of mutation position and amino acid substitutions on the dms-scores of single mutants. Additionally, the study explores the physiological properties of neighboring amino acids, to identify causality for a lowered dms-score. Moreover, this project delves into the phenomenon of epistasis in GFP, investigating the interplay of multiple mutations and identifying buried amino acid residues that exhibit enhanced epistatic effects. Furthermore, the influence of the number of mutations on the dms-score is analyzed, leading to the development of weighted rankings to identify mutations with positive epistatic effects. Finally, pyrosetta is utilized to predict the difference in free energy (ΔG) between the unfolded and folded protein. The findings provide valuable insights into the factors influencing the dms-score in GFP mutants and have implications for protein engineering and design.

# basic analysis
To start of we did general basic plots visualizing the dataset. The code in "basic_plots_proposal" results in the two plots used in the project proposal. 
See the plots in the two png files. Codes for general sorting of the data set are in the directory "general_sorting".

# Single mutation 
##### (Roman Kurley)
In this part of the project, the analysis focused primarily on single mutations. The main objective was to determine and quantify the relative impact of the mutation position versus the newly created amino acid on the resulting dms-score. Furthermore, the physiological impact of the change in amino acid to the fluorescence activity was determined by analyzing amino acid-specific parameters such as the volume of side chains (VSC), net charge index of side chains (NCISC), polarity (P1) and solvent accessible surface area (SASA).
 
## How to use the repository for single mutations:

#### Directory "single mutations" contains everything about single mutations. This is the important directory!
   
-Single_mutants_final_overview_1_code.py: code behind the analysis of the first overview 

-Single_mutants_final_overview_1_plots.ipnyb: Plots and outputs of the code above (DMS scores of each single mutation based on the new AA, heatmaps showing
available mutations, DMS scores of each single mutation based on the position)
     
-Single_mutants_final_basic_analysis_2_code.py: code behind the basic analysis.

-Single_mutants_final_basic_analysis_2_plots.ipnyb: Plots and outputs of the code above (Mean and median by position, mean and median by new amino acid, DMS 
score vs. position for the new amino acid = R scatter, Distribution of DMS scores by new amino acid boxplot, Distribution of DMS scores by new amino acid
violin plot, analysis of chromophore single mutations)
     
-Single_mutants_final_statistical_tests_3_code.py: code behind the statistical tests.

-Single_mutants_final_statistical_tests_3_plots_true.ipnyb: Plots and outputs of the code above (Distribution of DMS scores by mutation position boxplots,
ANOVA results, eta-squared results, Bonferroni correction, Shapiro-wilks results, q-q-plot of the DMS scores, Mann-Whitney-U results individually for each
mutation and summarized, Mann-Whitney-U results binary scatter, Kruskal-Wallis results, Friedman results)
      
-Single_mutants_final_neighbourhoods_4_code.py: code behind the neighbourhood analysis.

-Single_mutants_final_neighbourhoods_4_plots.ipnyb: Plots and outputs of the code above (confirmation that calculations in the dataframe are correct, enter a specific single mutation to see the corresponding mutated and unmutated neighbourhoods + properties)
      
-Single_mutants_final_GFP_variants_5.ipnyb: implemented blast results, blast results from NCBI website, CLustal Omega results, GFP variants + properties, FP colors + properties)

#### Directory "plots_and_outputs_for_presentation_and_report" contains all plots and tables/outputs used in the final presentation and in the report.

#### Directory "uncleaned_work_single_mutations" contains older work-in-progress versions of code.
 

# Epistasis Analysis 
Additionally, analyses were performed to obtain information about the ability of mutations to stabilize the fitness scores of proteins that already contain other mutations and get quantified data to rank them accordingly. To make the calculations as specific as possible, weighted factors were computed, and the effect of mutation count on protein fitness got included.

## Sequential Mutants
##### (Tianxin Angela Ma)
The data set was filtered to find sequential mutants to determine whether epistasis can be observed. All codes for this part of the project can be found in the directory “Sequential_Mutants”. The filtering was done with “sequential_mutants_filtering.ipynb”, and plots were generated with “sequential_mutants_plots.ipynb”. Plots used for the report can also be found as png-files in the subdirectory "plot_paths". The directory “further codes” contains previous code versions and trials.

## Structure Analysis 
##### (Tianxin Angela Ma)
To analyze the impact of amino acid residue orientation on epistasis, double mutation mutants were grouped according to the orientation of both amino acid residues. All files for this part of the project can be found in the directory “Structure_Analysis”. The needed structure information is stored in “AS_residues.xlsx”, where “1” stands for buried amino acid residues and “0” for surface-exposed ones. “AS_residues_in_out” contains the entire code for processing and plot generation, and plots used for the report are saved as png-files.

## Variance Analysis 
##### (Lisa Duttenhöfer)
To analyze the influence the mutation count has on the fitness of the protein, two plots were generated, showing the influence (see "final plots (output) for the png files). 
The directory "further code (trials etc)" contains codes developed in the process of the final plots, showing the development process of the final code.
The code in "variance_analysis_code.py"is the combined code from the notebooks in the directory "jupyter notebooks".
For results see the png files. 

## Rankings 
##### (Lisa Duttenhöfer)
To quantify the effect mutations have on the mutant's fitness due to positive epistatic effects, ranking and different variance calculations got developed and applied in different rankings.  
Again the directory "further code (trials etc)" contains the original code produced while the development of the final idea ("different_ranking_ideas") and a trial code performing cross-validation on ranking 5 and five folds. Due to the missing of a valid testing dataset, it could only be performed on parts of the training set and was therefore not included in the final report.
The directory "python files" contains the python files with the shortened code for both the weighted and unweighted rankings. 
The directory "ranking notebooks" contains the notebooks with the complete code for all rankings, commented for better understanding ("unweighted rankings.ipynb" and "weighted rankings.ipynb").
For results (output) there are html files for all existing rankings ordered in three folders.
  
- "complete unweighted rankings (output)" : the whole tables with all the important data of the mutations for all unweighted rankings (FSD is not calculated weighted)
   
- "complete weighted rankings (output)" : the whole tables with all the important data of the mutations for all weighted rankings (FSD is calculated weighted)
  
- "ranked mutation names (output) : all rankings without the data making up the scores, just the ranked mutation names

To see which ranking is based on which calculation see the juypyter notebook. 
The rankings included in the final presentation were:
1a (a), 3 (b), 3_w (c), 5_a (d), 5_a_w (e), 5_b (f), 5_b_w (g), 5_c (h), 5_c_w (i), 5_d (j), 7_w (k)

# Protein stability prediction with PyRosetta 
##### (Rebecca Ress)
The stability of proteins is important as it can affect proper folding and functional performance. Protein stability can be characterized by thermodynamic stability, which refers to the resistance against denaturation, as well as kinetic stability, which measures the resistance of a protein against irreversible inactivation (Liu, Xun et al. 2019). Since no experimentally determined parameters capturing stability were available in the dataset used (Sarkisyan, Bolotin et al. 2016), we utilized the software platform PyRosetta to predict a stability parameter, enabling us to investigate the relationship between function and stability in mutated GFP. 
## How to use the repository for predicting ddG with PyRosetta
The "PyRosetta" folder contains the "Final Version" and the "Other Analysis" subfolders.
- The "Final Version" subfolder contains only imagine file ("correlation_plot_DMS_score_ddG") and a Jupyter notebook "Final_Version_Prediction_with_Pyrosetta_and_Analysis" which contains everything for the report and presentation i.a. the prediction of ddG with PyRosetta and the analysis
- The "Other Analysis" subfolder contains the "Additional Analysis" folder and the "Tm Values" subfolder.
    - The "Other Analysis" subfolder contains a Jupyter notebook "REB_ddG_additional_analysis_final" which stores all additional analysis and fails, and the Jupyter notebooks "REB_Prediction_ddG_version_1" and "REB_Prediction_ddG_version_2" which contain older versions of ddG's prediction
    - The "Tm Values" ​​subfolder contains a Jupyter notebook with additional predictions for the Tm parameter

## used packages
numpy
pandas
matplotlib
seaborn
scipy.stats
statsmodels.api
ols
re
statsmodels.stats.multitest
Bio.Blast
Bio
Bio.Align
Bio.Seq
Bio.SeqRecord
Bio.Align.Applications
