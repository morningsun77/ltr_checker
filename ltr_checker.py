#!/usr/bin/env python
# coding: utf-8 

'''
Author: Zhao-Yang Chen
'''

import torch
from random import seed,randint
import subprocess
from time import sleep
import argparse
from numpy import array
import os
import multiprocessing
from Bio import SeqIO
import time
from tqdm import tqdm
import logging
from Harvest_format import LTRHarvestProcessor, LTRElement
from filter import *
from model import fasta_to_onehot_to_predict, predict, DETECT, fasta2one_hot
from sequence_utils import genome_dict, split_num_l, split_ltr_range, ltr_range
from ltr_annotators import ltr, ltr_detector, run_ltr_annotation, run_ltr_detector
from output_parser import *



parser = argparse.ArgumentParser()

parser.add_argument('--genome',type=str,required=True,help="genome file")
parser.add_argument('--threads',type=int,required=False,default=10,help='number of threads, default=10')
parser.add_argument('--output',type=str,required=True,help="output path")
parser.add_argument('--stride',type=int,required=False,default=10000,help="stride size for sliding windows in LTR-RT identification, in bp, default=10000")
parser.add_argument('--max',type=int,required=False,default=15000,help="max separation distance of two LTR,default=15000")
parser.add_argument('--min',type=int,required=False,default=1000,help="min separation distance of two LTR, default=1000")
parser.add_argument('--max_ltr',type=int,required=False,default=7000,help="max length of LTR, default=7000")
parser.add_argument('--min_ltr',type=int,required=False,default=100,help="min length of LTR, default=100")
parser.add_argument('--tgca',type=str,required=False,default="no",help="whether require the presence of TGCA, default=no")
parser.add_argument('--tsd',type=str,required=False,default="no",help="whether require the presence of TSD, default=no")
parser.add_argument('--model',type=str,required=False,help="path of model")
parser.add_argument('--split',type=int,required=False,default=2,help="how many chromosme segments you want to split into, default=2")
parser.add_argument('--device',type=str,required=False,default="cpu",help="cpu or cuda, default=cpu")
parser.add_argument('--method',type=str,required=False,default="ltr_finder",help="method to use: ltr_finder, ltr_harvest, ltrdetector, or all, default=ltr_finder")
parser.add_argument('--identity', type=int, required=False, default=85, help="minimum identity between 5' and 3' LTRs for LtrDetector, default=85")
parser.add_argument('--filter',type=str,required=False,default="yes",help="whether to run filtering modules (Module 1-6), options: yes/no, default=yes")
parser.add_argument('--output_format',type=str,required=False,default="ltr_harvest",help="unified output format when using --method all: ltr_harvest, ltr_finder, or ltrdetector, default=ltr_harvest")
args = parser.parse_args()

if __name__ == '__main__':
    print("#"*41)
    print("#" + " "*39 + "#")
    print("#" + " "*39 + "#")
    print("#" + " "*14 + "ltr_checker" + " "*14 + "#")
    print("#" + " "*39 + "#")
    print("#" + " "*39 + "#")
    print("#"*41)
    start_time = time.time()
    ## argparse
    file = args.genome
    file = os.path.abspath(file)
    threads = args.threads
    stride = args.stride
    max_len_threshold = args.max
    min_len_threshold = args.min
    max_ltr = args.max_ltr
    min_ltr = args.min_ltr
    tg_ca = args.tgca
    split = args.split
    TSD = args.tsd
    device = args.device
    method = args.method.lower()
    total_win_len = 50000
    outputDir = args.output
    outputDir = os.path.abspath(outputDir)
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    if args.model:
        model_path = os.path.abspath(args.model)
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        model_path = os.path.join(script_dir, 'model', 'model_cnn.pth')
    identity = args.identity
    output_format = args.output_format.lower()
    run_filter = args.filter.lower()
    
    # Validate method selection
    if method not in ["ltr_finder", "ltr_harvest", "ltrdetector", "all"]:
        print(f"Error: Unsupported method '{method}'. Please use 'ltr_finder', 'ltr_harvest', 'ltrdetector', or 'all'.")
        exit(1)
    
    # Validate output format selection
    if output_format not in ["ltr_harvest", "ltr_finder", "ltrdetector"]:
        print(f"Error: Unsupported output format '{output_format}'. Please use 'ltr_harvest', 'ltr_finder', or 'ltrdetector'.")
        exit(1)

    # load model
    model = DETECT()
    model_state_dict = torch.load(model_path, map_location=torch.device(device))
    model.load_state_dict(model_state_dict)
    model.eval()
    genome = genome_dict(file)

    # Detect LTR-RTs by model
    print("Begin to detect LTR-RTs by model.")
    model_res = []
    for k,v in genome.items():
        chr = fasta_to_onehot_to_predict(v,k,split,stride,model)
        model_res.append(chr)
    model = None

    # Determine LTR-RT regions
    region = []
    for Dict in model_res:
        for k,v in Dict.items():
            cc = ltr_range(k,v,len(genome[k]),10000)
            region.append(cc)
    Region = sum(region,[])

    # Define environment variables
    for i in range(1,threads + 1,1):
        locals()["temp" + str(i)] = []

    # Distribute LTR-RT regions to environment variables
    if len(Region) < threads:
        threads = len(Region)
        start = 0
        end = len(Region) // threads
        num = (len(Region) // threads)
        for i in range(1,threads + 1,1):
            if i < threads:
                for index in range(start,end,1):
                    locals()["temp" + str(i)].append(Region[index])
                start += num
                end += num
            elif i == threads:
                for index in range(start,len(Region),1):
                    locals()["temp" + str(i)].append(Region[index])
    else:
        start = 0
        end = len(Region) // threads
        for i in range(1,threads + 1,1):
            num = (len(Region) // threads)
            if i < threads:
                for index in range(start,end,1):
                    locals()["temp" + str(i)].append(Region[index])
                start += num
                end += num
            elif i == threads:
                for index in range(start,len(Region)-1,1):
                    locals()["temp" + str(i)].append(Region[index])

    seq_id_file = []
    for i in range(1,threads + 1,1):
        with open(outputDir + '/' +"temp" +str(i) +'.txt','w') as f:
            seq_id_file.append("temp" +str(i))
            for seq in locals()["temp" + str(i)]:
                seq = seq.split('\t')
                Seq = seq[0] + ":" + str(seq[1]) + "-" + str(seq[2])
                Seq = Seq + '\n'
                f.writelines(Seq)

    seq_id_file = []
    BED = []
    for i in range(1,threads + 1,1):
        seq = "temp" + str(i)
        seq_id_file.append(seq)
        f = open(outputDir + '/' +seq + '.txt','r')
        for txt in f:
            txt = txt.split('\n')[0]
            txt = txt.split(':')
            bed = txt[0] + '\t' + txt[1]
            bed = bed.split('-')
            final_txt = bed[0] + '\t' + bed[1]
            BED.append(final_txt)
    
    ltr_seq = dict()
    for i in BED:
        i = i.split("\t")
        key = i[0] + ":" + str(i[1]) + "-" + str(i[2])
        seq = genome[i[0]][int(i[1]):int(i[2])]
        if len(seq) > 0:  # Only add non-empty sequences
            ltr_seq[key] = seq
    
    # Annotate LTR-RTs with the selected method
    if method == "all":
        # Run all three method in sequence: ltr_finder -> ltr_harvest -> ltrdetector
        print("\n" + "="*50)
        print("Running all annotation methods in sequence...")
        print("="*50)
        
        # 1. Run LTR_FINDER
        print("\n[1/3] Running LTR_FINDER...")
        run_ltr_annotation(threads, max_len_threshold, min_len_threshold, max_ltr, min_ltr, 
                        outputDir, tg_ca, TSD, ltr_seq, seq_id_file, "ltr_finder", output_suffix="all")
        
        # 2. Run LTR_HARVEST
        print("\n[2/3] Running LTR_HARVEST...")
        run_ltr_annotation(threads, max_len_threshold, min_len_threshold, max_ltr, min_ltr, 
                        outputDir, tg_ca, TSD, ltr_seq, seq_id_file, "ltr_harvest", output_suffix="all")
        
        # 3. Run LtrDetector
        print("\n[3/3] Running LtrDetector...")
        run_ltr_detector(outputDir, min_len_threshold, max_len_threshold, min_ltr, max_ltr, 
                        identity, threads, ltr_seq, seq_id_file, output_suffix="all")
        
        print("\n" + "="*50)
        print("All methods completed successfully!")
        print("Output files:")
        print("  - ltr_finder.out")
        print("  - ltr_harvest.out")
        print("  - ltr_detector.out")
        print("="*50)
        
        # Unify outputs to the specified format
        unify_outputs(outputDir, output_format)
        
        # Remove temporary files
        for i in range(1, threads + 1, 1):
            if os.path.exists(outputDir + '/temp' + str(i) + '.txt'):
                os.remove(outputDir + '/temp' + str(i) + '.txt')
            if os.path.exists(outputDir + '/temp' + str(i)):
                import shutil
                shutil.rmtree(outputDir + '/temp' + str(i))
    
    elif method == "ltrdetector":
        # LtrDetector uses its own multithreading and does not require our multiprocessing method
        run_ltr_detector(outputDir, min_len_threshold, max_len_threshold, min_ltr, max_ltr, identity, threads, ltr_seq, seq_id_file)
    else:
        # Use the original ltr_finder or ltr_harvest
        run_ltr_annotation(threads, max_len_threshold, min_len_threshold, max_ltr, min_ltr, 
                        outputDir, tg_ca, TSD, ltr_seq, seq_id_file, method)
    
    # Generate LTR-RT FASTA file regardless of method or output format chosen
    print("\n" + "="*50)
    print("Generating LTR-RT sequences in FASTA format...")
    print("="*50)
    try:
        generate_ltr_fasta(outputDir, file)
        print("LTR-RT FASTA file successfully generated!")
    except Exception as e:
        print(f"Warning: Failed to generate FASTA file: {e}")
        print("You may need to manually process the output files.")
    
    ltr_candidates_fasta = outputDir + '/ltr_candidates.fasta'
    
    

    # argements
    dir = outputDir
    input_fasta = ltr_candidates_fasta
    genome_file = args.genome

    # Check if filtering modules should be run
    if run_filter == "no":
        end_time = time.time()
    else:
        ## Module1
        print("\n" + "="*50)
        print("Module 1: Initial Quality Filtering")
        print("="*50)
        print("Starting initial quality filtering of LTR-RT candidates...")
        print(f"  - Filtering candidates with excessive gaps (N's)")
        print(f"  - Removing tandem repeats")
        print(f"  - Validating internal region lengths")
        print(f"  - Checking LTR:internal region ratios")
        filter_module = LTRModule1(
            max_gap_size=10,              # Max consecutive N length allowed
            min_internal_len=100,         # Min internal region length
            max_internal_len=15000,       # Max internal region length
            min_ltr_internal_ratio=0.05,  # Min LTR:internal ratio
            max_ltr_internal_ratio=50.0,  # Max LTR:internal ratio
            trf_output_dir="./trf_results"  # TRF results output directory
        )
        dir_1 = dir + "/module1"
        output_1 = dir_1 + "/module1.fasta"
        if not os.path.exists(dir_1):
            os.makedirs(dir_1)
        
        try:
            filter_module.filter_candidates(input_fasta, output_1)
            temp_dir = os.getcwd() + "/trf_results"
            os.rmdir(temp_dir)  # Remove TRF results directory
        except Exception as e:
            logging.error(f"Filtering process error: {e}")
            raise

        ## Module2
        print("\n" + "="*50)
        print("Module 2: LTR Boundary Refinement")
        print("="*50)
        print("Refining LTR boundaries using BLAST alignment...")
        print(f"  - Aligning left and right LTRs")
        print(f"  - Adjusting boundaries based on alignment results")
        print(f"  - Validating refined boundaries")
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
        print("\n" + "="*50)
        print("Module 3: Nested Insertion Detection and Filtering")
        print("="*50)
        print("Detecting and filtering nested LTR insertions...")
        print(f"  - Checking for nested insertions")
        print(f"  - Removing redundant elements")
        print(f"  - Validating element integrity")
        dir_3 = dir + "/module3"
        if not os.path.exists(dir_3):
            os.makedirs(dir_3)

        input_3_fa = dir_2 + "/boundary_modified.fasta"
        
        ltr = Module3(genome_file, input_3_fa, dir_3, r=1.3e-8)
        ltr.process()

        ## Module4
        print("\n" + "="*50)
        print("Module 4: Protein Domain Annotation and Classification")
        print("="*50)
        print("Annotating protein domains and classifying LTR-RTs...")
        print(f"  - Searching for RT, INT, GAG, and other domains")
        print(f"  - Classifying elements as Copia, Gypsy, or Unknown")
        print(f"  - Generating domain position maps")
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
            print("\nLTR retrotransposon classification results:")
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
                print(f"Copia: {sum(1 for r in results if r['Classification'] == 'Copia')}")
                print(f"Gypsy: {sum(1 for r in results if r['Classification'] == 'Gypsy')}")
                print(f"Unknown: {sum(1 for r in results if r['Classification'] == 'unknown')}")
                print(f"\nNote: Results are sorted in ascending order by LTR ID (e.g., LTR_1, LTR_2, LTR_10...)")

            else:
                print("No sequences met the criteria (at least 3 domains required)")
                if output_4_res:
                    # Create empty file
                    with open(output_4_res, 'w') as f:
                        f.write("Sequence_ID\tClassification\tDomain_count\tDomain_positions\n")
                        
        except FileNotFoundError:
            print(f"错误：找不到文件 {output_4_tsv}")
            sys.exit(1)
        except Exception as e:
            print(f"错误：处理文件时出错 - {str(e)}")
            sys.exit(1)


        ## Module5 
        print("\n" + "="*50)
        print("Module 5: Sequence Similarity Filtering")
        print("="*50)
        print("Filtering highly similar sequences...")
        print(f"  - Clustering sequences by similarity")
        print(f"  - Removing redundant sequences")
        print(f"  - Keeping representative sequences from each cluster")
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
            print(f"Program execution failed: {str(e)}")

        ## Module6
        print("\n" + "="*50)
        print("Module 6: Final Cleanup and Validation")
        print("="*50)
        print("Performing final cleanup and validation...")
        print(f"  - Removing false positives")
        print(f"  - Final quality check")
        print(f"  - Generating final LTR-RT library")
        dir_6 = dir + "/module6"
        if not os.path.exists(dir_6):
            os.makedirs(dir_6)
        input_6_fa = dir_5 + "/filtered.fasta"

        ltr_ids = []
        with open(dir_2 + '/boundary_results.tsv','r') as file:
            # Skip the first line
            next(file)

            # Read the remaining lines
            for line in file:
                # Split each line into columns
                columns = line.strip().split()  # Split by whitespace by default; modify if using a different delimiter

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
            # Skip the first line
            next(file)

            # Read the remaining lines
            for line in file:
                # Split each line into columns
                columns = line.strip().split()  # Split by whitespace by default; modify if using a different delimiter

                # Ensure the line has at least 7 columns
                if len(columns) >= 7:
                    # Check if the 7th column is 'pass'
                    if columns[0] == 'PASS':
                        # Add the 1st column to the list
                        ltr_ids.append(columns[0])
        sequences = extract_sequences_by_ids(content, ltr_ids)
        write_sequences_to_file(sequences, output_6_pa)
        write_sequences_to_file(sequences, outputDir + '/final_ltr_rt.fasta')
        fasta_to_gff(output_6_pa, dir_6 + '/final_ltr_rt.gff3')
        end_time = time.time()
        print("ltr_checker cost {} s.".format(end_time - start_time))