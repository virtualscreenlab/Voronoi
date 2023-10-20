import csv

# Define the CSV file paths
input_csv_file = "Database_corr.csv"
output_csv_file = "Filtered_IC50_Data.csv"

# Create a set to store IC50 values
ic50_values = set()

# Read the input CSV file, filter rows, and write to the output CSV file
filtered_rows = []

try:
    with open(input_csv_file, "r") as input_file:
        reader = csv.reader(input_file)
        header = next(reader)
        filtered_rows.append(header)  # Include the header row

        for row in reader:
            try:
                ic50 = float(row[1])  # Assuming the IC50 value is in the second column

                # Add the IC50 value to the set if it's less than or equal to 1000
                if ic50 <= 1000:
                    if ic50 not in ic50_values:
                        ic50_values.add(ic50)
                    filtered_rows.append(row)
            except (ValueError, IndexError):
                # Handle potential conversion errors or missing columns
                pass

    with open(output_csv_file, "w", newline="") as output_file:
        writer = csv.writer(output_file)
        writer.writerows(filtered_rows)

except FileNotFoundError:
    print("Input CSV file not found.")
except Exception as e:
    print("An error occurred:", str(e))

