import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import csv
import math
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi

# Create lists to store IC50 values, entropy values, SMILES strings, and number of atoms
ic50_values = []
entropy_values = []
smiles_list = []
num_atoms_list = []

# Open the "smiles_ic50.dat" file and read the SMILES strings, IC50 values
with open("smiles_ic50.dat", "r") as file:
    for line in file:
        line = line.strip()  # Remove newline character
        smiles, ic50 = line.split(",")  # Split the line into SMILES and IC50
        ic50 = float(ic50)  # Convert IC50 to float

        # Generate a molecule from the SMILES
        mol = Chem.MolFromSmiles(smiles)

        # Check if the molecule is valid
        if mol is not None:
            # Get the number of atoms
            num_atoms = mol.GetNumAtoms()

            # Generate a molecular graph
            g = Chem.RWMol(mol)
            
            # Generate a 2D conformer
            AllChem.Compute2DCoords(mol)

            # Get 2D atomic coordinates
            coords_2d = np.array(mol.GetConformer().GetPositions()[:, :2])
            
            # Compute Voronoi diagram
            vor = Voronoi(coords_2d)

            # Count occurrences of each polygon class
            polygon_class_counts = {}

            for region_index in vor.point_region:
                region_vertices = vor.regions[region_index]

                # Exclude unbounded regions
                if -1 not in region_vertices:
                    polygon_class = len(region_vertices)

                    if polygon_class in polygon_class_counts:
                        polygon_class_counts[polygon_class] += 1
                    else:
                        polygon_class_counts[polygon_class] = 1

            # Calculate total number of bounded regions and calculate proportions
            total_bounded_regions = sum(polygon_class_counts.values())
            proportions = {polygon_class: count / total_bounded_regions for polygon_class, count in polygon_class_counts.items()}

            # Calculate Voronoi entropy
            voronoi_entropy = -sum(p * math.log(p) for p in proportions.values() if p > 0)

            # Add IC50 value, entropy value, SMILES, and number of atoms to the lists
            ic50_values.append(ic50)
            entropy_values.append(voronoi_entropy)
            smiles_list.append(smiles)
            num_atoms_list.append(num_atoms)

# Save the data to a CSV file
with open("entropy_results.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["SMILES", "IC50", "Number_of_Atoms", "Voronoi_Entropy"])  # Write the header row
    for smiles, ic50, num_atoms, voronoi_entropy in zip(smiles_list, ic50_values, num_atoms_list, entropy_values):
        writer.writerow([smiles, ic50, num_atoms, voronoi_entropy])

print("Calculation and data saving complete.")

