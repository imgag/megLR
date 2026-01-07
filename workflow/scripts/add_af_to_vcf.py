#!/usr/bin/env python
"""
Add AF (Allele Frequency) field to VCF FORMAT column.
Calculates AF from AD (Allele Depth) field if present.
"""

import sys
import argparse


def parse_ad_field(ad_str):
    """
    Parse AD field and return list of depths.
    AD format: REF_depth,ALT1_depth,ALT2_depth,...
    """
    # Handle missing or empty AD values
    if not ad_str or ad_str == '.':
        return None
    try:
        depths = [int(x) for x in ad_str.split(',')]
        return depths
    except (ValueError, AttributeError):
        return None


def calculate_af(ad_depths):
    """
    Calculate AF from AD depths.
    AF = ALT_depth / (REF_depth + ALT_depth)
    For multiallelic sites: AF = ALT_i_depth / total_depth
    """
    if not ad_depths or len(ad_depths) < 2:
        return None
    
    ref_depth = ad_depths[0]
    alt_depths = ad_depths[1:]
    total_depth = sum(ad_depths)
    
    if total_depth == 0:
        return None
    
    # Calculate AF for each ALT allele
    afs = [f"{alt_depth / total_depth:.6f}" for alt_depth in alt_depths]
    return ','.join(afs)


def process_vcf(input_file, output_file):
    """
    Process VCF file to add AF field to FORMAT column.
    """
    af_header_added = False
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Handle meta-information header lines
            if line.startswith('##'):
                outfile.write(line)
                continue
            
            # Add AF header before the #CHROM line
            if line.startswith('#CHROM') and not af_header_added:
                outfile.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
                af_header_added = True
                outfile.write(line)
                continue
            
            # Process variant lines
            if not line.startswith('#'):
                fields = line.rstrip('\n').split('\t')
                
                # VCF has: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample1 sample2 ...
                if len(fields) < 10:
                    # Not enough fields, write as-is
                    outfile.write(line)
                    continue
                
                format_field = fields[8]
                format_tags = format_field.split(':')
                
                # Check if AF already exists
                if 'AF' in format_tags:
                    outfile.write(line)
                    continue
                
                # Find AD field index
                ad_index = -1
                if 'AD' in format_tags:
                    ad_index = format_tags.index('AD')
                
                # Add AF to FORMAT field
                format_tags.append('AF')
                fields[8] = ':'.join(format_tags)
                
                # Process each sample column
                for i in range(9, len(fields)):
                    sample_data = fields[i].split(':')
                    
                    # Calculate AF from AD if available
                    af_value = '.'
                    if ad_index >= 0 and ad_index < len(sample_data):
                        ad_str = sample_data[ad_index]
                        ad_depths = parse_ad_field(ad_str)
                        if ad_depths:
                            calculated_af = calculate_af(ad_depths)
                            if calculated_af:
                                af_value = calculated_af
                    
                    # Add AF value to sample data
                    sample_data.append(af_value)
                    fields[i] = ':'.join(sample_data)
                
                outfile.write('\t'.join(fields) + '\n')


def main():
    # Check if running under Snakemake
    if 'snakemake' in globals():
        input_file = snakemake.input.vcf
        output_file = snakemake.output.vcf
    else:
        parser = argparse.ArgumentParser(
            description='Add AF field to VCF FORMAT column'
        )
        parser.add_argument(
            '-i', '--input',
            required=True,
            help='Input VCF file'
        )
        parser.add_argument(
            '-o', '--output',
            required=True,
            help='Output VCF file'
        )
        
        args = parser.parse_args()
        input_file = args.input
        output_file = args.output
    
    try:
        process_vcf(input_file, output_file)
    except Exception as e:
        print(f"Error processing VCF: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
