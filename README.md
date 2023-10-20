# Voronoi
Voronoi entropy scripts

We investigate the correlation between the Voronoi Entropy (VE) of ligand molecules and their affinity to receptors to test the hypothesis that less ordered ligands have higher mobility of molecular groups and therefore higher probability to attach to receptors.
VE of 1144 ligands is calculated using SMILES-based 2D graphs representing molecular structure. The affinity of the ligands with the SARS-CoV-2 main protease is obtained from the BindingDB Database as half-maximal inhibitory concentration (IC50) data. 
The VE distribution is close to the Gaussian, 0.4≤Sv≤1.66 and a strong correlation with IC50 is found, IC50=-275*Sv+613 nM, indicating to the correlation between ligand complexity and affinity. 
On the contrary, the Shannon Entropy (SE) descriptor failed to provide enough evidence to reject the null hypothesis (p-value > 0.05), indicating that the spatial arrangement of atoms is crucial for molecular mobility and binding. 
