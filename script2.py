import csv

# Define the CSV file path
csv_file = "CID_IC50_Complete_Data_DS.csv"  # Replace with the path to your CSV file

# Create an empty list to store SMILES strings and IC50 values
smiles_ic50_list = []

# Read the CSV file and extract the SMILES and IC50 values
with open(csv_file, "r") as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header row
    
    for row in reader:
        cid = int(row[0])  # Assuming the CID is in the first column
        smiles = row[2]    # Assuming the SMILES is in the third column
        ic50 = float(row[1])  # Assuming the IC50 value is in the second column
        
        smiles_ic50_list.append((smiles, ic50))

# Save the SMILES and IC50 values to a "smiles_ic50.dat" file
with open("smiles_ic50.dat", "w") as file:
    for smiles, ic50 in smiles_ic50_list:
        file.write(f"{smiles},{ic50}\n")
