#!/usr/bin/env python
# coding: utf-8 

'''
Author: Chen Zhaoyang
'''

import os
import pandas as pd
from collections import defaultdict
import sys
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import re
from modules.module1 import LTRModule1
from modules.module2 import Module2
from modules.module3 import Module3
from modules.module4 import DomainAnnotator, LTRClassifier
from modules.module5 import LTRFilter
from modules.module6 import LTRCleaner

def extract_sequences_by_ids(file_content, id_list):
    """
    Extract sequences and their associated information for specified LTR IDs
    
    Parameters:
        file_content: Content containing LTR sequences
        id_list: List of LTR IDs to extract
        
    Returns:
        dict: A dictionary with complete LTR headers as keys and sequences as values
    """
    result_dict = {}
    current_header = None 
    sequence = []
    
    for line in file_content.split('\n'):
        if line.startswith('>'):
            # Save previous sequence if its ID is in the target list
            if current_header and current_header.split()[0][1:] in id_list:
                result_dict[current_header] = ''.join(sequence)
            
            # Get new header
            current_header = line
            sequence = []
            
        elif line.strip() and current_header:
            sequence.append(line.strip())
    
    # Process the last sequence
    if current_header and current_header.split()[0][1:] in id_list:
        result_dict[current_header] = ''.join(sequence)
        
    return result_dict

def write_sequences_to_file(sequences_dict, output_filename):
    """
    Write sequences dictionary to a FASTA format file
    
    Parameters:
        sequences_dict: Dictionary containing headers and sequences
        output_filename: Name of the output file
    """
    with open(output_filename, 'w') as outfile:
        for header, sequence in sequences_dict.items():
            # Write header line
            outfile.write(f"{header}\n")
            # Write sequence with line breaks every 60 characters
            for i in range(0, len(sequence), 60):
                outfile.write(f"{sequence[i:i+60]}\n")

def extract_ltr_sequences(input_fasta, output_fasta):
    """
    Extract 5' and 3' LTR sequences from complete LTR-RT sequences
    """
    # Regular expression for parsing coordinates
    pattern = r'coords=(\d+):(\d+)\s+5ltr=(\d+):(\d+)\s+3ltr=(\d+):(\d+)'
    ltr_sequences = []
    
    for record in SeqIO.parse(input_fasta, "fasta"):
        match = re.search(pattern, record.description)
        if not match:
            continue

        # Get coordinates
        coords_start = int(match.group(1))
        ltr5_start, ltr5_end = int(match.group(3)), int(match.group(4))
        ltr3_start, ltr3_end = int(match.group(5)), int(match.group(6))

        # Calculate relative positions
        full_seq = str(record.seq)
        ltr5_seq = full_seq[ltr5_start-coords_start:ltr5_end-coords_start+1]
        ltr3_seq = full_seq[ltr3_start-coords_start:ltr3_end-coords_start+1]

        # Create new sequence records
        # 5' LTR
        ltr_sequences.append(SeqRecord(
            Seq(ltr5_seq),
            id=f"{record.id}_5ltr_{ltr5_start}_{ltr5_end}",
            description=""
        ))
        
        # 3' LTR
        ltr_sequences.append(SeqRecord(
            Seq(ltr3_seq),
            id=f"{record.id}_3ltr_{ltr3_start}_{ltr3_end}",
            description=""
        ))

    # Write to output file
    SeqIO.write(ltr_sequences, output_fasta, "fasta")
    print(f"Extracted {len(ltr_sequences)} LTR-RT sequences")

    return len(ltr_sequences)
def Args():
    parser = argparse.ArgumentParser(
        prog='fillter',
        description='Filtering of candidate LTR-RTs .',
    )
    parser.add_argument("--sequence", 
                       type=str,
                       help="input LTR-RTs sequences in fasta format [required]")
    parser.add_argument("--genome", 
                       type=str,
                       help="input genome in fasta format [required]")
    parser.add_argument("--output", 
                       type=str,
                       default="output",  # Default output directory
                       help="output directory name [default: output]")
    
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)
          
    return args

if __name__ == '__main__':
    args = Args()
    # argements
    dir = args.output
    input_fasta = args.sequence
    genome_file = args.genome

    ## Module1
    filter_module = LTRModule1(
        max_gap_size=10,              # Maximum allowed gap size
        min_internal_len=100,         # Minimum internal region length
        max_internal_len=15000,       # Maximum internal region length
        min_ltr_internal_ratio=0.05,  # Minimum LTR:internal ratio
        max_ltr_internal_ratio=50.0,  # Maximum LTR:internal ratio
        trf_output_dir="./trf_results"  # TRF results output directory
    )
    dir_1 = dir + "/module1"
    output_1 = dir_1 + "/module1.fasta"
    if not os.path.exists(dir_1):
        os.makedirs(dir_1)
    
    # try:
    #     filter_module.filter_candidates(input_fasta, output_1)
    #     temp_dir = os.getcwd() + "/trf_results"
    #     os.rmdir(temp_dir)  # Clean up TRF results directory
    # except Exception as e:
    #     logging.error(f"Filtering process error: {e}")
    #     raise

    ## Module2
    dir_2 = dir + "/module2"
    if not os.path.exists(dir_2):
        os.makedirs(dir_2)
    input_2_fa = output_1

    processor = Module2(
        genome_file, 
        blast_path="blastn",  # If in PATH, can use program name directly
        makeblastdb_path="makeblastdb"
    )
    processor.process_ltrs(input_2_fa, dir_2)

    ## Module3
    dir_3 = dir + "/module3"
    if not os.path.exists(dir_3):
        os.makedirs(dir_3)

    input_3_fa = dir_2 + "/boundary_modified.fasta"
    
    ltr = Module3(genome_file, input_3_fa, dir_3, r=1.3e-8)
    ltr.process()

    ## Module4
    dir_4 = dir + "/module4"
    if not os.path.exists(dir_4):
        os.makedirs(dir_4)
    input_4_fa = dir_3 + "/passed.fa"
    input_4_tsv = dir_4 + "/ltr_rt_annotation"
    output_4_tsv = dir_4 + "/ltr_rt_annotation_domains.tsv"
    output_4_res = dir_4 + "/domain_analysis_results.tsv"
    if not os.path.exists(input_fasta):
        raise FileNotFoundError(f"Input file not found: {input_fasta}")
        
    print(f"Input file: {input_fasta}")
    print(f"File size: {os.path.getsize(input_fasta)} bytes")

    # Create annotator instance and run annotation
    annotator = DomainAnnotator()
    annotator.annotate_ltr_rt(input_4_fa, input_4_tsv)

    classifier = LTRClassifier()

    try:
        # Read input file
        df = pd.read_csv(output_4_tsv, sep='\t')
        print(f"Successfully read file: {output_4_tsv}")
        print(f"Number of data rows: {len(df)}")

        # Process data
        results = classifier.process_ltr_data(df)

        # Output results
        print("\nLTR-RT classification results:")
        print("-" * 80)
        
        if results:
            # Create output DataFrame (results are already sorted)
            output_df = pd.DataFrame(results)[['Sequence_ID', 'Classification', 'Domain_count', 'Domain_positions']]

            # Print to console
            for result in results:
                print(f"Sequence ID: {result['Sequence_ID']}")
                print(f"Classification: {result['Classification']}")
                print(f"Domain count: {result['Domain_count']}")
                print(f"Domain positions: {result['Domain_positions']}")
                print("-" * 80)

            # Save to file
            if output_4_res:
                output_df.to_csv(output_4_res, sep='\t', index=False)
                print(f"\nResults saved to: {output_4_res}")
            else:
                print("\nOutput table:")
                print(output_df.to_string(index=False))

            # Statistics
            print(f"\nStatistics:")
            print(f"Total sequences: {len(results)}")
            print(f"Copia types: {sum(1 for r in results if r['Classification'] == 'Copia')}")
            print(f"Gypsy types: {sum(1 for r in results if r['Classification'] == 'Gypsy')}")
            print(f"Unknown types: {sum(1 for r in results if r['Classification'] == 'unknown')}")
            print(f"\nNote: Results are sorted in ascending order by LTR number (e.g., LTR_1, LTR_2, LTR_10...)")

        else:
            print("No sequences met the criteria (at least 3 domains required)")
            if output_4_res:
                # Create empty file
                with open(output_4_res, 'w') as f:
                    f.write("Sequence_ID\tClassification\tDomain_count\tDomain_positions\n")
                    
    except FileNotFoundError:
        print(f"Error: File not found {output_4_tsv}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: An error occurred while processing the file - {str(e)}")
        sys.exit(1)


    ## Module5 
    dir_5 = dir + "/module5"
    if not os.path.exists(dir_5):
        os.makedirs(dir_5)
    output_4_fa = dir_4 + '/output_sequences.fasta'
    input_5_fa = output_4_fa
    with open(output_4_res, 'r') as f:
        ltr_ids = f.readlines()[1:]  # begin from second line
        ltr_ids = [line.split('\t')[0] for line in ltr_ids]
    with open(input_fasta, 'r') as f:
        content = f.read()

    sequences = extract_sequences_by_ids(content, ltr_ids)
    write_sequences_to_file(sequences, output_4_fa)
    try:
        filter = LTRFilter(input_5_fa, dir_5)
        filter.run()
    except Exception as e:
        print(f"Error: An error occurred while running the filter - {str(e)}")

    ## Module6
    dir_6 = dir + "/module6"
    if not os.path.exists(dir_6):
        os.makedirs(dir_6)
    input_6_fa = dir_5 + "/filtered.fasta"

    ltr_ids = []
    with open(dir_2 + '/boundary_results.tsv','r') as file:
        # Skip first line
        next(file)

        # Read remaining lines
        for line in file:
            # Split each line into columns
            columns = line.strip().split()  # Default split by whitespace, modify if different delimiter

            # Ensure the line has at least 7 columns
            if len(columns) >= 7:
                # Check if the 7th column is 'pass'
                if columns[6] == 'pass' or columns[6]== "truncated":
                    # Add the 1st column to the list
                    ltr_ids.append(columns[0])
    sequences = extract_sequences_by_ids(content, ltr_ids)
    library_LTR_RT_file = dir_5 + "/LTR_library.fa"
    library_file = dir_5 + "/LTRs_library.fa"
    output_6_pa = dir_6 + '/final_output.fasta'
    write_sequences_to_file(sequences, library_LTR_RT_file)
    extract_ltr_sequences(library_LTR_RT_file, library_file)

    cleaner = LTRCleaner(genome_file=genome_file, 
                         candidates_file=input_6_fa, 
                         library_file=library_file, 
                         output_dir=dir_6)
    cleaner.clean()

    ltr_ids = []
    with open(dir_6 + '/LTR_processing_log.tsv','r') as file:
        # Skip first line
        next(file)

        # Read remaining lines
        for line in file:
            # Split each line into columns
            columns = line.strip().split()  # Default split by whitespace, modify if different delimiter

            # Ensure the line has at least 7 columns
            if len(columns) >= 7:
                # Check if the 7th column is 'pass'
                if columns[0] == 'PASS':
                    # Add the 1st column to the list
                    ltr_ids.append(columns[0])
    sequences = extract_sequences_by_ids(content, ltr_ids)
    write_sequences_to_file(sequences, output_6_pa)