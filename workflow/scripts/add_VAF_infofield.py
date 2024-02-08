import vcf

# Input VCF file and output VCF file
input_vcf_file = "input.vcf"
output_vcf_file = "output.vcf"

# Define a function to calculate the new AF value
def calculate_AF(record):
    AO = record.INFO.get('AO', 0)
    RO = record.INFO.get('RO', 0)
    if AO + RO == 0:
        return 0.0
    return AO / (AO + RO)

# Open the input VCF file
vcf_reader = vcf.Reader(filename=input_vcf_file)

# Open the output VCF file for writing
vcf_writer = vcf.Writer(open(output_vcf_file, 'w'), vcf_reader)

# Process each variant in the input VCF file
for record in vcf_reader:
    # Calculate the new AF value and update the INFO field
    record.INFO['AF'] = calculate_AF(record)

    # Write the updated record to the output VCF file
    vcf_writer.write_record(record)

print(f"AF values updated and written to {output_vcf_file}")

#-----------------------------
## withput PyVCF

import vcf

# Input VCF file and output VCF file
input_vcf_file = "input.vcf"
output_vcf_file = "output.vcf"

# Define a function to calculate the new AF value
def calculate_AF(record):
    # Split the AO field by commas and convert to integers
    AO_values = [int(value) for value in record.INFO.get('AO', [0])]

    # If there are two values in AO, use them, otherwise use 0
    if len(AO_values) == 2:
        AO, RO = AO_values
    else:
        AO, RO = 0, record.INFO.get('RO', 0)

    if AO + RO == 0:
        return 0.0
    return AO / (AO + RO)

# Open the input VCF file
vcf_reader = vcf.Reader(filename=input_vcf_file)

# Open the output VCF file for writing
vcf_writer = vcf.Writer(open(output_vcf_file, 'w'), vcf_reader)

# Process each variant in the input VCF file
for record in vcf_reader:
    # Calculate the new AF value and update the INFO field
    record.INFO['AF'] = calculate_AF(record)

    # Write the updated record to the output VCF file
    vcf_writer.write_record(record)

print(f"AF values updated and written to {output_vcf_file}")
