import csv
import math

def shannon_entropy(data):
    total_count = len(data)
    frequency = {}

    for item in data:
        if item in frequency:
            frequency[item] += 1
        else:
            frequency[item] = 1

    entropy = -sum((count / total_count) * math.log2(count / total_count) for count in frequency.values())
    return entropy

input_csv = r'C:\Users\aglik-pc\Desktop\entropy_results.csv'
output_csv = r'C:\Users\aglik-pc\Desktop\entropy_results_with_shennon.csv'

with open(input_csv, 'r') as file:
    reader = csv.reader(file)
    # Skip the first line (header)
    next(reader)
    rows = [row for row in reader]

for row in rows:
    entropy = shannon_entropy(row[0])
    row.append(entropy)

with open(output_csv, 'w', newline='') as file:
    writer = csv.writer(file)
    # Write the header row
    writer.writerow(['SMILES', 'IC50', 'Number_of_Atoms','Voronoi_Entropy','Shennon_Enthropy'])
    writer.writerows(rows)