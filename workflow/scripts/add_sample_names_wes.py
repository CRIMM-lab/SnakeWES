import csv
import sys

# Check for the correct number of command-line arguments
if len(sys.argv) != 4:
	print("Usage: python script.py input_file.tsv output_file.tsv samples_list.txt")
	sys.exit(1)

# Input and output file paths
input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

# List of samples
samples_file_path = sys.argv[3] 

try:
	with open(samples_file_path, 'r') as file:
		
		# Read the file contents and split them into a list (assuming each line is an item)
		samples_list = file.read().splitlines()
		
		print("File contents as a list:")
		
		print(samples_list)

except FileNotFoundError:
	print(f"File '{samples_file_path}' not found.")

except Exception as e:
	print(f"An error occurred: {str(e)}")

# Function to process a row and return the modified row
def process_row(row):
	# Extract the last 10 columns from the row
	last_10_columns = row[-len(samples_list):]

	# Initialize variables to store samples meeting the condition and their count
	meeting_samples = []
	count = 0

	# Iterate through the samples and corresponding columns
	for sample, column in zip(samples_list, last_10_columns):
		# Check if the column is not equal to "./."
		if column != "./.":
			count += 1
			meeting_samples.append(sample)

	# Add the meeting_samples as a comma-separated string and count after the 4th column
	row.insert(7, ",".join(meeting_samples))
	row.insert(8, count)

	return row


# Read the input file, process each row, and write to the output file

try: 
	
	with open(input_file_path, 'r', newline='') as tsvfile:
		reader = csv.reader(tsvfile, delimiter='\t')
		rows = list(reader)
		
		# Process the data rows
		processed_rows = [process_row(row) for row in rows]

except FileNotFoundError:
	print(f"Error: The file '{input_file_path}' was not found.")

except PermissionError:
	print(f"Error: Permission denied to access '{input_file_path}'.")

except csv.Error as e:
	print(f"CSV error occurred: {e}")

except Exception as e:
	print(f"An unexpected error occurred: {e}")


# Write the modified data to the output file

try:    
	
	with open(output_file_path, 'w', newline='') as outfile:
		writer = csv.writer(outfile, delimiter='\t')
		writer.writerows(processed_rows)

except PermissionError:
	print(f"Error: Permission denied to write to '{output_file_path}'.")

except csv.Error as e:
	print(f"CSV error occurred while writing: {e}")

except Exception as e:
	print(f"An unexpected error occurred while writing: {e}")


print(f"Data processed and written to {output_file_path}")

