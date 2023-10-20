import csv
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem

# Define the CSV file path
csv_file = "Filtered_IC50_Data.csv"  # Replace with the path to your CSV file

# Create an empty dictionary to store PubChem CID and IC50 values
cid_ic50_data = {}

# Create an empty set to store unique CIDs
unique_cids = set()

# Read the CSV file and extract the CIDs
with open(csv_file, "r") as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header row
    
    for row in reader:
        cid = int(row[0])  # Assuming the CID is in the first column
        ic50 = float(row[1])  # Assuming the IC50 value is in the second column
        
        # Add the CID to the set if it doesn't exist
        unique_cids.add(cid)
        
        # Retrieve compound information using PubChem CID
        compound = pcp.Compound.from_cid(cid)
        
        # Calculate the number of atoms using RDKit
        molecule = Chem.MolFromSmiles(compound.isomeric_smiles)
        num_atoms = molecule.GetNumAtoms() if molecule is not None else 0
        num_torsions = AllChem.CalcNumRotatableBonds(molecule)
        
        # Add PubChem CID, IC50 value, and number of atoms to the dictionary
        cid_ic50_data[cid] = {
            "CID": cid,
            "IC50": ic50,
#            "CompoundName": compound.iupac_name,
            "SMILES": compound.canonical_smiles,
            "MolecularWeight": compound.molecular_weight,
            "NumAtoms": num_atoms,
            "NumTorsions": num_torsions,
            # Add more desired compound properties
        }

# Export the CID and IC50 data to a new CSV file
output_csv_file = "CID_IC50_Complete_Data_DS.csv"  # Replace with the desired output file name

with open(output_csv_file, "w", newline="") as file:
    writer = csv.DictWriter(file, fieldnames=["CID", "IC50", "SMILES", "MolecularWeight", "NumAtoms", "NumTorsions"])
    writer.writeheader()
    writer.writerows(cid_ic50_data.values())

# Convert the set of unique CIDs back to a list
#cid_list = list(unique_cids)
