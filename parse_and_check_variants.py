import sys
import os

clin_sig_mapping = {
    '0': 'Uncertain significance',
    '1': 'Not provided',
    '2': 'Benign',
    '3': 'Likely benign',
    '4': 'Likely pathogenic',
    '5': 'Pathogenic',
    '6': 'Drug response',
    '7': 'Histocompatibility',
    '255': 'Other'
}

def parse_info(info):
    fields = {}
    for pair in info.split(';'):
        if '=' in pair:
            key, value = pair.split('=', 1)
            fields[key] = value
    return fields

def parse_variant_line(line):
    fields = line.strip().split('\t')
    chrom, pos, var_id, ref, alt, qual, filter_, info, format_, *data = fields

    info_dict = parse_info(info)

    gene_name = "Unknown"
    gene_id = "Unknown"
    impact = "Unknown"
    if "ANN" in info_dict:
        annotations = info_dict["ANN"].split(',')
        if annotations:
            first_annotation = annotations[0].split('|')
            if len(first_annotation) > 3:
                gene_name = first_annotation[3]  
            if len(first_annotation) > 4:
                gene_id = first_annotation[4]  
            if len(first_annotation) > 2:
                impact = first_annotation[2]  

    condition = info_dict.get("CLNDBN", "Unknown").split('|')[0]  
    significance_raw = info_dict.get("CLNSIG", "Unknown")
    significance = "Unknown"
    if significance_raw != "Unknown":
        significance_code = significance_raw.split('|')[0]
        significance = clin_sig_mapping.get(significance_code, "Unknown")

    return chrom, var_id, gene_name, gene_id, condition, impact, significance

def parse_vcf(input_vcf, output_txt):
    if not os.path.exists(input_vcf):
        print(f"Error: Input file '{input_vcf}' not found.")
        sys.exit(2)

    variant_count = 0
    with open(input_vcf, 'r') as infile, open(output_txt, 'w') as outfile:
        outfile.write("Chromosome\tVariant_ID\tGene_Name\tGene_ID\tCondition\tImpact\tClinical_significance\n")

        for line in infile:
            if line.startswith("#"): 
                continue

            parsed_variant = parse_variant_line(line)
            outfile.write('\t'.join(parsed_variant) + '\n')
            variant_count += 1

    if variant_count > 0:
        print(f"Parsing complete. {variant_count} variants written to {output_txt}.")
        sys.exit(0)  
    else:
        print("No variants found.")
        sys.exit(1) 

if __name__ == "__main__":
    if len(sys.argv) != 3:
        #print("Usage: python3 parse_and_check_variants.py <input_vcf> <output_txt>")
        sys.exit(1)

    input_vcf = sys.argv[1] 
    output_txt = sys.argv[2] 
    parse_vcf(input_vcf, output_txt)