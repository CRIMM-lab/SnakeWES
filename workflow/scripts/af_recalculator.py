import vcf
import argparse

def update_af(info):
    ao = info.get('AO', 0)
    ro = info.get('RO', 0)

    if isinstance(ao, list):
        ao = ao[0] if ao else 0

    if isinstance(ro, list):
        ro = ro[0] if ro else 0

    if ao and ro:
        new_af = ro / (ao + ro)
        info['AF'] = new_af

    return info

def process_vcf(input_vcf, output_vcf):
    vcf_reader = vcf.Reader(filename=input_vcf)
    vcf_writer = vcf.Writer(open(output_vcf, 'w'), vcf_reader)

    for record in vcf_reader:
        # Update the AF value in the INFO field
        record.INFO = update_af(record.INFO)

        # Write the updated record to the output VCF file
        vcf_writer.write_record(record)

    print(f"VCF file '{input_vcf}' processed. Updated file saved to '{output_vcf}'.")

def main():
    parser = argparse.ArgumentParser(description="Update AF values in a VCF file.")
    parser.add_argument("input_vcf", help="Input VCF file")
    parser.add_argument("output_vcf", help="Output VCF file")
    
    args = parser.parse_args()
    
    process_vcf(args.input_vcf, args.output_vcf)

if __name__ == "__main__":
    main()

