import os
import subprocess
from time import sleep
import multiprocessing
def ltr(max_len_threshold, min_len_threshold, max_ltr, min_ltr, outputDir, tg_ca, TSD, dict, seq_id_file, software="ltr_finder"):
    """
    Annotate LTR-RTs using LTR_FINDER or LTR_HARVEST
    
    args:
        max_len_threshold: Maximum distance between two LTRs
        min_len_threshold: Minimum distance between two LTRs
        max_ltr: Maximum length of LTR
        min_ltr: Minimum length of LTR
        outputDir: Output directory
        tg_ca: Whether to include TGCA
        TSD: Whether to include TSD
        dict: Sequence dictionary
        seq_id_file: Sequence ID file
        software: Software to use, "ltr_finder" or "ltr_harvest"

    returns:
        Annotation results list
    """
    path = outputDir + '/' + seq_id_file
    f = open(outputDir + '/' + seq_id_file + '.txt','r')
    fasta = open(outputDir + '/' + seq_id_file + '.fasta','w')
    for seq_id in f:
        seq_id = seq_id.split('\n')[0]
        if seq_id in dict and len(dict[seq_id]) > 0:  # Ensure the sequence exists and is not empty
            fasta.write('>' + seq_id+'\n')
            fasta.write(dict[seq_id]+'\n')
    fasta.close()
    sleep(0.01)
    if os.path.exists(path) != True:
        os.mkdir(path)
    sleep(0.01)
    print(seq_id_file + " is start")
    print(path)
    num = int(seq_id_file.split('temp')[1])
    sleep(num)
    
    if software == "ltr_finder":
        if tg_ca:
            finder_filter = '1111'
        else:
            finder_filter = '0000'
        if TSD:
            finder_filter += '1'
        else:
            finder_filter += '0'
        finder_filter += '000000'
        
        command = 'ltr_finder -F ' + finder_filter + ' -D ' + str(max_len_threshold) + ' -d ' + str(min_len_threshold) + ' -w2 -C -p 20 -M 0.85 -L ' + str(max_ltr) + ' -l ' + str(min_ltr) + ' ' + outputDir + '/' + seq_id_file + '.fasta'
        output = subprocess.Popen(command, stdout=subprocess.PIPE, cwd=path, shell=True)
        sleep(0.01)
        res = output.stdout
        res = res.read().decode().split('\n')
        print(seq_id_file + " is finish!")
        
        Res = []
        for i in res:
            if '[' in i:
                Res.append(i)
    
    elif software == "ltr_harvest":
        # Create index
        command1 = 'gt suffixerator -db ' + outputDir + '/' + seq_id_file + '.fasta -indexname ' + path + '/genome.fa -tis -suf -lcp -des -ssp -sds -dna'
        subprocess.run(command1, shell=True, check=True)

        # Run LTR_HARVEST
        motif_option = "-motif TGCA -motifmis 1 " if tg_ca == "yes" else ""
        tsd_option = "-mintsd 4 -maxtsd 6 " if TSD == "yes" else ""
        
        command2 = 'gt ltrharvest -index ' + path + '/genome.fa -minlenltr ' + str(min_ltr) + ' -maxlenltr ' + str(max_ltr) + ' ' + tsd_option + motif_option + '-similar 85 -vic 10 -seed 20 -seqids yes > ' + path + '/genome.fa.harvest.scn'
        subprocess.run(command2, shell=True, check=True)
        
        # Read result file
        with open(path + '/genome.fa.harvest.scn', 'r') as harvest_file:
            lines = harvest_file.readlines()

        # Process results
        Res = []
        for line in lines:
            if line.startswith('#') or not line.strip():
                continue

            # Parse LTR_HARVEST result line
            fields = line.strip().split()
            if len(fields) >= 11:  # Ensure there are enough fields
                s_ret = fields[0]
                e_ret = fields[1]
                l_ret = fields[2]
                s_lLTR = fields[3]
                e_lLTR = fields[4]
                l_lLTR = fields[5]
                s_rLTR = fields[6]
                e_rLTR = fields[7]
                l_rLTR = fields[8]
                sim = fields[9]
                seq_nr = fields[10]

                # Read original sequence ID
                with open(outputDir + '/' + seq_id_file + '.txt', 'r') as f:
                    seq_ids = [line.strip() for line in f if line.strip()]
                    if int(seq_nr) < len(seq_ids):
                        seq_id = seq_ids[int(seq_nr)]
                    else:
                        seq_id = f"unknown_{seq_nr}"
                
                # Format result line
                formatted_line = f"{s_ret} {e_ret} {l_ret} {s_lLTR} {e_lLTR} {l_lLTR} {s_rLTR} {e_rLTR} {l_rLTR} {sim} {seq_id}"
                Res.append(formatted_line)
        
    if os.path.exists(outputDir + '/' + seq_id_file + '.fasta'):
        os.remove(outputDir + '/' + seq_id_file + '.fasta')
    
    return Res


def run_ltr_annotation(threads, max_len_threshold, min_len_threshold, max_ltr, min_ltr, 
                      outputDir, tg_ca, TSD, ltr_seq, seq_id_file, software, output_suffix=""):
    """
    Use the selected software to annotate LTR-RTs and process results
    args:
        threads: Number of threads
        max_len_threshold: Maximum distance between two LTRs
        min_len_threshold: Minimum distance between two LTRs
        max_ltr: Maximum length of LTR
        min_ltr: Minimum length of LTR
        outputDir: Output directory
        tg_ca: Whether to include TGCA
        TSD: Whether to include TSD
        ltr_seq: LTR sequence dictionary
        seq_id_file: List of sequence ID files
        software: Software to use, "ltr_finder" or "ltr_harvest"
        output_suffix: Output file suffix (used to distinguish between different software outputs)
    """
    print(" "*10)
    print(f"Begin to annotate LTR-TRs by {software.upper()}.")
    pool = multiprocessing.Pool(processes=threads)

    # Validate sequence ID files to ensure they contain at least one non-empty sequence
    valid_seq_id_file = []
    for seq_id in seq_id_file:
        has_valid_sequence = False
        with open(outputDir + '/' + seq_id + '.txt', 'r') as f:
            for line in f:
                seq_id_full = line.strip()
                if seq_id_full in ltr_seq and len(ltr_seq[seq_id_full]) > 0:
                    has_valid_sequence = True
                    break
        if has_valid_sequence:
            valid_seq_id_file.append(seq_id)

    # Use valid sequence ID files
    localresults = [pool.apply_async(ltr, args=[max_len_threshold, min_len_threshold, 
                                            max_ltr, min_ltr, outputDir, tg_ca, 
                                            TSD, ltr_seq, seq_id, software]) 
                for seq_id in valid_seq_id_file]


    
    pool.close()
    pool.join()
    localTables = [p.get() for p in localresults]

    total_txt = sum(localTables, [])
    
    if software == "ltr_finder":
        for i in total_txt:
            line = i.split('\t')
            line = line[1:]
            element_1 = line[0].replace(':', '-').split('-')
            element_2 = line[1].split('-')
            element_2 = str(int(element_2[0]) + int(element_1[1]) - 1) + '-' + str(int(element_2[1]) + int(element_1[1]))
            element_1 = element_1[0]
            line[0] = element_1
            line[1] = element_2
            line = '\t'.join(str(i) for i in line)
            Index = total_txt.index(i)
            total_txt[Index] = line

        title = 'SeqID\tLocation\tLTR len\tInserted element len\tTSR\tPBS\tPPT\tRT\tIN (core)\tIN (c-term)\tRH\tStrand\tScore\tSharpness\tSimilarity'
        output_file = 'ltr_finder.out' if output_suffix else 'ltr_check.out'
        
        total_txt.insert(0, title)

        with open(outputDir + '/' + output_file, 'w') as f:
            for i in total_txt:
                i = i + '\n'
                f.writelines(i)
    
    elif software == "ltr_harvest":
        # LTR_HARVEST results are already formatted in the ltr function, use directly here
        output_file = 'ltr_harvest.out'

        # Merge all temporary files' harvest results
        harvest_results = []

        # Check if there are available temporary directories
        temp_dirs_exist = False
        for i in range(1, threads + 1):
            temp_dir = outputDir + '/temp' + str(i)
            if os.path.exists(temp_dir):
                temp_dirs_exist = True
                break
                
        if temp_dirs_exist:
            # Write the final harvest results file
            with open(outputDir + '/' + output_file, 'w') as f:
                f.write("# args=-minlenltr " + str(min_ltr) + " -maxlenltr " + str(max_ltr) + "\n")
                f.write("# predictions are reported in the following way\n")
                f.write("# s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-id\n")
                f.write("# where:\n")
                f.write("# s = starting position\n")
                f.write("# e = ending position\n")
                f.write("# l = length\n")
                f.write("# ret = LTR-retrotransposon\n")
                f.write("# lLTR = left LTR\n")
                f.write("# rLTR = right LTR\n")
                f.write("# sim = similarity\n")
                f.write("# seq-id = sequence identifier\n\n")

                # Collect all harvest results from temporary directories
                for i in range(1, threads + 1):
                    temp_dir = outputDir + '/temp' + str(i)
                    harvest_file_path = temp_dir + '/genome.fa.harvest.scn'
                    if os.path.exists(harvest_file_path):
                        with open(harvest_file_path, 'r') as harvest_file:
                            for line in harvest_file:
                                if not line.startswith('#') and line.strip():
                                    parts = line.strip().split()
                                    if len(parts) >= 11:
                                        # Get sequence ID
                                        seq_nr = parts[10]
                                        seq_id_file_path = outputDir + '/temp' + str(i) + '.txt'
                                        if os.path.exists(seq_id_file_path):
                                            with open(seq_id_file_path, 'r') as seq_file:
                                                seq_ids = [line.strip() for line in seq_file if line.strip()]
                                                if int(seq_nr) < len(seq_ids):
                                                    full_seq_id = seq_ids[int(seq_nr)]
                                                    
                                                    # Extract chromosome name and offset from the full ID
                                                    if ":" in full_seq_id:
                                                        chrom = full_seq_id.split(":")[0]
                                                        offset = int(full_seq_id.split(":")[1].split("-")[0])

                                                        # Adjust position information by adding the offset
                                                        s_ret = str(int(parts[0]) + offset)
                                                        e_ret = str(int(parts[1]) + offset)
                                                        s_lLTR = str(int(parts[3]) + offset)
                                                        e_lLTR = str(int(parts[4]) + offset)
                                                        s_rLTR = str(int(parts[6]) + offset)
                                                        e_rLTR = str(int(parts[7]) + offset)

                                                        # Update result line
                                                        parts[0] = s_ret
                                                        parts[1] = e_ret
                                                        parts[3] = s_lLTR
                                                        parts[4] = e_lLTR
                                                        parts[6] = s_rLTR
                                                        parts[7] = e_rLTR
                                                        parts[10] = chrom  # Retain only the chromosome name
                                                    else:
                                                        parts[10] = full_seq_id
                                                else:
                                                    parts[10] = f"unknown_{seq_nr}"
                                            
                                            f.write(" ".join(parts) + "\n")
        else:
            # If no temporary directory exists, use total_txt directly
            with open(outputDir + '/' + output_file, 'w') as f:
                f.write("# args=-minlenltr " + str(min_ltr) + " -maxlenltr " + str(max_ltr) + "\n")
                f.write("# predictions are reported in the following way\n")
                f.write("# s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-id\n")
                f.write("# where:\n")
                f.write("# s = starting position\n")
                f.write("# e = ending position\n")
                f.write("# l = length\n")
                f.write("# ret = LTR-retrotransposon\n")
                f.write("# lLTR = left LTR\n")
                f.write("# rLTR = right LTR\n")
                f.write("# sim = similarity\n")
                f.write("# seq-id = sequence identifier\n\n")
                
                for result in total_txt:
                    # Parse result line
                    if "\t" in result:
                        parts = result.split("\t")
                    else:
                        parts = result.split()
                    
                    if len(parts) >= 11:
                        # Obtain and handle sequence ID
                        full_seq_id = parts[10]
                        
                        # Extract chromosome name and offset from the full ID
                        if ":" in full_seq_id:
                            chrom = full_seq_id.split(":")[0]
                            offset = int(full_seq_id.split(":")[1].split("-")[0])
                            
                            # Adjust position information by adding the offset
                            s_ret = str(int(parts[0]) + offset)
                            e_ret = str(int(parts[1]) + offset)
                            s_lLTR = str(int(parts[3]) + offset)
                            e_lLTR = str(int(parts[4]) + offset)
                            s_rLTR = str(int(parts[6]) + offset)
                            e_rLTR = str(int(parts[7]) + offset)

                            # Update result line
                            parts[0] = s_ret
                            parts[1] = e_ret
                            parts[3] = s_lLTR
                            parts[4] = e_lLTR
                            parts[6] = s_rLTR
                            parts[7] = e_rLTR
                            parts[10] = chrom  # Retain only the chromosome name
                        
                        f.write(" ".join(parts) + "\n")
                    else:
                        f.write(result + "\n")
        
    # Remove temporary files (cleanup only when output_suffix is not used)
    if not output_suffix:
        for i in range(1, threads + 1, 1):
            if os.path.exists(outputDir + '/temp' + str(i) + '.txt'):
                os.remove(outputDir + '/temp' + str(i) + '.txt')
            if os.path.exists(outputDir + '/temp' + str(i)):
                import shutil
                shutil.rmtree(outputDir + '/temp' + str(i))

def ltr_detector(max_len_threshold, min_len_threshold, max_ltr, min_ltr, outputDir, tg_ca, TSD, dict, seq_id_file, identity=85):
    """
    Annotate LTR-RTs using LtrDetector
    Parameters:
        outputDir: Output directory
        min_len: Minimum length of LTR-RT
        max_len: Maximum length of LTR-RT
        min_ltr: Minimum length of LTR
        max_ltr: Maximum length of LTR
        identity: Minimum identity
        threads: Number of threads
        dict: Sequence dictionary
        seq_id_file: Sequence ID file

    Returns:
        List of annotation results
    """
    # Create temporary directory for fasta files
    fasta_dir = os.path.join(outputDir, 'ltr_detector_input')
    if not os.path.exists(fasta_dir):
        os.makedirs(fasta_dir)

    # Read sequence IDs and create fasta files
    f = open(outputDir + '/' + seq_id_file + '.txt', 'r')
    for seq_id in f:
        seq_id = seq_id.split('\n')[0]
        # Create fasta file, ensure the filename ends with .fa
        fasta_file = os.path.join(fasta_dir, seq_id.replace(':', '_').replace('-', '_') + '.fa')
        with open(fasta_file, 'w') as fasta:
            fasta.write('>' + seq_id + '\n')
            fasta.write(dict[seq_id] + '\n')
    f.close()

    # Run LtrDetector
    command = f'LtrDetector -fasta {fasta_dir} -destDir {outputDir} -minLen {min_len_threshold} -maxLen {max_len_threshold} -minLenLTR {min_ltr} -maxLenLTR {max_ltr} -id {identity} -nThreads {threads}'
    subprocess.run(command, shell=True, check=True)

    # Collect results
    results = []
    for filename in os.listdir(outputDir):
        if filename.endswith('Detector.bed'):
            with open(os.path.join(outputDir, filename), 'r') as bed_file:
                # Read all lines
                lines = bed_file.readlines()
                # Only keep actual data lines
                for line in lines:
                    if not line.strip() or "SeqID" in line or "Start" in line:
                        continue
                    results.append(line.strip())
    
    # Remove temporary directory
    import shutil
    shutil.rmtree(fasta_dir)
    
    return results

def run_ltr_detector(outputDir, min_len_threshold, max_len_threshold, min_ltr, max_ltr, identity, threads, ltr_seq, seq_id_file, output_suffix=""):
    """
    Annotate LTR-RTs using LtrDetector and process results

    Parameters:
        outputDir: Output directory
        min_len_threshold: Minimum length of LTR-RT
        max_len_threshold: Maximum length of LTR-RT
        min_ltr: Minimum length of LTR
        max_ltr: Maximum length of LTR
        identity: Minimum identity
        threads: Number of threads
        ltr_seq: LTR sequence dictionary
        seq_id_file: Sequence ID file list
        output_suffix: Output file suffix (used to distinguish between different software outputs)
    """
    print(" "*10)
    print("Begin to annotate LTR-TRs by LTRDETECTOR.")

    # Create a single fasta file containing all sequences
    fasta_dir = os.path.join(outputDir, 'ltr_detector_input')
    if not os.path.exists(fasta_dir):
        os.makedirs(fasta_dir)

    # Create a merged fasta file for all sequences
    merged_fasta = os.path.join(fasta_dir, 'all_sequences.fa')
    with open(merged_fasta, 'w') as f_out:
        for seq_id in seq_id_file:
            with open(outputDir + '/' + seq_id + '.txt', 'r') as f:
                for line in f:
                    seq_id_full = line.strip()
                    if seq_id_full in ltr_seq:
                        f_out.write('>' + seq_id_full + '\n')
                        f_out.write(ltr_seq[seq_id_full] + '\n')

    # Run LtrDetector
    command = f'LtrDetector -fasta {fasta_dir} -destDir {outputDir} -minLen {min_len_threshold} -maxLen {max_len_threshold} -minLenLTR {min_ltr} -maxLenLTR {max_ltr} -id {identity} -nThreads {threads}'
    subprocess.run(command, shell=True, check=True)

    # Collect and process results
    detector_results = []

    # Add standard LtrDetector header lines (two lines)
    header_line1 = "SeqID\tRetrotransposon\tLeft_LTR\tRight_LTR\t\t\tLeft_TSD\tRight_TSD\tPolypurine Tract\t\t\tTG\tCA"
    header_line2 = "Start\tEnd\tStart\tEnd\tStart\tEnd\tID\tStart\tEnd\tStart\tEnd\tStart\tEnd\tStrand\tPurine%\tStart\tEnd"
    detector_results.append(header_line1)
    detector_results.append(header_line2)

    # Process all output files
    detector_files_to_clean = []  # Keep track of files to clean up
    for filename in os.listdir(outputDir):
        if filename.endswith('Detector.bed'):
            detector_files_to_clean.append(os.path.join(outputDir, filename))
            with open(os.path.join(outputDir, filename), 'r') as bed_file:
                lines = bed_file.readlines()

                # Skip header lines
                data_start = False
                for i, line in enumerate(lines):
                    # Skip header lines
                    if "SeqID" in line or "Start" in line:
                        continue

                    # Process data lines
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 8:
                            seq_id = parts[0]
                            # Check if the format is "chr:start-end"
                            if ":" in seq_id and "-" in seq_id:
                                chrom = seq_id.split(":")[0]
                                offset = int(seq_id.split(":")[1].split("-")[0])

                                # Update position information by adding the offset
                                parts[1] = str(int(parts[1]) + offset)  # Retrotransposon Start
                                parts[2] = str(int(parts[2]) + offset)  # Retrotransposon End
                                parts[3] = str(int(parts[3]) + offset)  # Left_LTR Start
                                parts[4] = str(int(parts[4]) + offset)  # Left_LTR End
                                parts[5] = str(int(parts[5]) + offset)  # Right_LTR Start
                                parts[6] = str(int(parts[6]) + offset)  # Right_LTR End

                                # Update TSD positions (if they exist)
                                if parts[8] != "---":
                                    parts[8] = str(int(parts[8]) + offset)  # Left_TSD Start
                                if parts[9] != "---":
                                    parts[9] = str(int(parts[9]) + offset)  # Left_TSD End
                                if parts[10] != "---":
                                    parts[10] = str(int(parts[10]) + offset)  # Right_TSD Start
                                if parts[11] != "---":
                                    parts[11] = str(int(parts[11]) + offset)  # Right_TSD End

                                # Update PPT positions (if they exist)
                                if parts[12] != "---":
                                    parts[12] = str(int(parts[12]) + offset)  # PPT Start
                                if parts[13] != "---":
                                    parts[13] = str(int(parts[13]) + offset)  # PPT End

                                # Update TG/CA positions (if they exist)
                                if len(parts) > 16 and parts[16] != "---":
                                    parts[16] = str(int(parts[16]) + offset)  # TG Start
                                if len(parts) > 17 and parts[17] != "---":
                                    parts[17] = str(int(parts[17]) + offset)  # TG End

                                # Replace sequence ID with chromosome ID
                                parts[0] = chrom
                            
                            detector_results.append('\t'.join(parts))

    # Write final results to file
    output_file = 'ltr_detector.out'
    with open(outputDir + '/' + output_file, 'w') as f:
        for line in detector_results:
            f.write(line + '\n')

    # Clean up all Detector.bed files
    for bed_file in detector_files_to_clean:
        if os.path.exists(bed_file):
            os.remove(bed_file)
            print(f"Cleaned up: {os.path.basename(bed_file)}")

    # Clean up temporary files
    import shutil
    if os.path.exists(fasta_dir):
        shutil.rmtree(fasta_dir)

    # Delete all_sequencesDetector.bed file
    detector_bed_file = os.path.join(outputDir, 'all_sequencesDetector.bed')
    if os.path.exists(detector_bed_file):
        os.remove(detector_bed_file)
        print(f"Cleaned up: all_sequencesDetector.bed")

    # Clean up temporary files (only when not using output_suffix)
    if not output_suffix:
        for i in range(1, threads + 1, 1):
            if os.path.exists(outputDir + '/temp' + str(i) + '.txt'):
                os.remove(outputDir + '/temp' + str(i) + '.txt')
            if os.path.exists(outputDir + '/temp' + str(i)):
                import shutil
                shutil.rmtree(outputDir + '/temp' + str(i))