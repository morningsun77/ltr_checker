import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import re
from pathlib import Path
import shutil
import csv
from datetime import datetime

class LTRCleaner:
    def __init__(self, genome_file, candidates_file, library_file, output_dir):
        self.genome_file = genome_file
        self.candidates_file = candidates_file
        self.library_file = library_file
        self.output_dir = Path(output_dir)
        self.mask_threshold = 0.8
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Copy library file to output directory (to avoid path issues)
        self.local_library = self.output_dir / "LTR_library.fa"
        shutil.copy2(self.library_file, self.local_library)

        # Initialize processing log
        self.processing_log = []

    def add_log_entry(self, seq_id, status, reason="", mask_ratio=None, internal_length=None):
        """Add processing log entry"""
        entry = {
            'sequence_id': seq_id,
            'status': status,
            'reason': reason,
            'internal_mask_ratio': f"{mask_ratio:.4f}" if mask_ratio is not None else "N/A",
            'internal_length': internal_length if internal_length is not None else "N/A",
            'mask_threshold': self.mask_threshold
        }
        self.processing_log.append(entry)

    def write_processing_log(self):
        """Write processing log to file"""
        log_file = self.output_dir / "LTR_processing_log.tsv"
        
        with open(log_file, 'w', newline='') as f:
            fieldnames = ['sequence_id', 'status', 'reason', 'internal_mask_ratio', 
                         'internal_length', 'mask_threshold']
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')

            # Write header
            writer.writeheader()

            # Write records
            for entry in self.processing_log:
                writer.writerow(entry)

        # Generate summary
        summary_file = self.output_dir / "LTR_processing_summary.txt"
        total = len(self.processing_log)
        passed = sum(1 for e in self.processing_log if e['status'] == 'PASS')
        failed = total - passed

        # Count failure reasons
        fail_reasons = {}
        for entry in self.processing_log:
            if entry['status'] == 'FAIL' and entry['reason']:
                fail_reasons[entry['reason']] = fail_reasons.get(entry['reason'], 0) + 1
        
        with open(summary_file, 'w') as f:
            f.write(f"LTR Processing Summary\n")
            f.write(f"=====================\n")
            f.write(f"Processing Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write(f"Total sequences processed: {total}\n")
            f.write(f"Sequences passed: {passed} ({passed/total*100:.1f}%)\n")
            f.write(f"Sequences failed: {failed} ({failed/total*100:.1f}%)\n")
            f.write(f"Mask threshold used: {self.mask_threshold}\n\n")
            
            if fail_reasons:
                f.write("Failure reasons breakdown:\n")
                for reason, count in sorted(fail_reasons.items(), key=lambda x: x[1], reverse=True):
                    f.write(f"  - {reason}: {count}\n")

        print(f"\nProcessing log has been saved to: {log_file}")
        print(f"Processing summary has been saved to: {summary_file}")

    def create_blast_db(self):
        """Create BLAST database for RepeatMasker"""
        print("Creating BLAST database...")
        
        # Use makeblastdb to create database
        cmd = [
            "makeblastdb",
            "-in", str(self.local_library),
            "-dbtype", "nucl",
            "-out", str(self.local_library)
        ]
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print("BLAST database created successfully")
            return True
        except subprocess.CalledProcessError as e:
            print(f"Failed to create BLAST database: {e}")
            if e.stderr:
                print("Error:", e.stderr)
            return False
        except FileNotFoundError:
            print("Warning: makeblastdb not found, attempting to run RepeatMasker directly")
            return True

    def run_repeatmasker(self):
        """Run RepeatMasker for sequence masking"""
        # Create BLAST database
        self.create_blast_db()

        # Copy candidate sequences to output directory
        local_candidates = self.output_dir / "candidates.fa"
        shutil.copy2(self.candidates_file, local_candidates)
        
        cmd = [
            "RepeatMasker",
            "-engine", "ncbi",  # 使用ncbi而不是rmblast
            "-lib", str(self.local_library),
            "-q",
            "-no_is",
            "-norna",
            "-nolow",
            "-div", "40",
            "-cutoff", "225",
            str(local_candidates)
        ]
        
        try:
            print("Executing RepeatMasker command:", " ".join(cmd))

            # Run in output directory
            result = subprocess.run(
                cmd,
                cwd=str(self.output_dir),
                check=True,
                capture_output=True,
                text=True
            )
            
            print("RepeatMasker ran successfully")
            if result.stdout:
                print("Output:", result.stdout)
            return True
            
        except subprocess.CalledProcessError as e:
            print(f"RepeatMasker failed to run: {e}")
            if e.stdout:
                print("Output:", e.stdout)
            if e.stderr:
                print("Error:", e.stderr)

            # Trying different parameters
            print("\nTrying simplified parameters...")
            return self.run_repeatmasker_simple()
            
        except Exception as e:
            print(f"Error occurred while running RepeatMasker: {e}")
            return False

    def run_repeatmasker_simple(self):
        """Run RepeatMasker with simplified parameters"""
        local_candidates = self.output_dir / "candidates.fa"
        
        cmd = [
            "RepeatMasker",
            "-lib", str(self.local_library),
            "-q",
            "-nolow",
            str(local_candidates)
        ]
        
        try:
            print("Executing simplified RepeatMasker command:", " ".join(cmd))
            result = subprocess.run(
                cmd,
                cwd=str(self.output_dir),
                check=True,
                capture_output=True,
                text=True
            )
            print("RepeatMasker ran successfully (simplified mode)")
            return True
        except subprocess.CalledProcessError as e:
            print(f"RepeatMasker failed to run (simplified mode): {e}")
            return False

    def parse_ltr_positions(self, header):
        """Parse LTR position information from FASTA header"""
        pattern = r'chr=(.*?)\s+coords=(\d+):(\d+)\s+5ltr=(\d+):(\d+)\s+3ltr=(\d+):(\d+)'
        match = re.search(pattern, header)
        if match:
            return {
                'chr': match.group(1),
                'start': int(match.group(2)),
                'end': int(match.group(3)),
                '5ltr_start': int(match.group(4)),
                '5ltr_end': int(match.group(5)), 
                '3ltr_start': int(match.group(6)),
                '3ltr_end': int(match.group(7))
            }
        return None

    def remove_masked_portions(self, sequence):
        """Remove masked portions (lowercase letters and X)"""
        cleaned = ''.join([base for base in sequence if base.isupper() and base not in 'NX'])
        return cleaned

    def process_masked_sequences(self):
        """Process masked sequences, focusing on internal regions"""
        # Initially record the original candidate sequence IDs
        original_ids = set()
        for record in SeqIO.parse(self.candidates_file, "fasta"):
            original_ids.add(record.id)

        # Find masked file
        masked_file = None
        possible_names = [
            self.output_dir / "candidates.fa.masked",
            self.output_dir / "candidates.masked",
            self.output_dir / f"{Path(self.candidates_file).name}.masked"
        ]
        
        for possible_file in possible_names:
            if possible_file.exists():
                masked_file = possible_file
                break
        
        if not masked_file:
            print(f"Error: Masked file not found, tried the following files:")
            for f in possible_names:
                print(f"  - {f}")
            # Record all sequences as failed
            for seq_id in original_ids:
                self.add_log_entry(seq_id, "FAIL", "Masked file not found")
            return False
            
        output_file = self.output_dir / "cleaned_LTRs.fa"
        
        cleaned_records = []
        processed_ids = set()
        
        for record in SeqIO.parse(masked_file, "fasta"):
            processed_ids.add(record.id)
            pos = self.parse_ltr_positions(record.description)
            
            if not pos:
                print(f"Warning: Failed to parse LTR positions for {record.id}")
                self.add_log_entry(record.id, "FAIL", "Failed to parse LTR positions")
                continue

            # Get the full masked sequence
            full_seq = str(record.seq)
            seq_len = len(full_seq)

            # Calculate positions relative to the start of the sequence
            coords_start = pos['start']
            ltr5_start_rel = pos['5ltr_start'] - coords_start
            ltr5_end_rel = pos['5ltr_end'] - coords_start + 1
            ltr3_start_rel = pos['3ltr_start'] - coords_start
            ltr3_end_rel = pos['3ltr_end'] - coords_start

            # Validate positions
            if ltr5_end_rel > seq_len or ltr3_end_rel > seq_len:
                print(f"Warning: {record.id} LTR positions exceed sequence length, skipping")
                self.add_log_entry(record.id, "FAIL", "LTR positions exceed sequence length")
                continue

            # Split into three parts
            ltr_5prime = full_seq[ltr5_start_rel:ltr5_end_rel]
            internal_region = full_seq[ltr5_end_rel:ltr3_start_rel]
            ltr_3prime = full_seq[ltr3_start_rel:ltr3_end_rel]

            # Check internal region
            if len(internal_region) == 0:
                print(f"Warning: {record.id} No internal region, skipping")
                self.add_log_entry(record.id, "FAIL", "No internal region", 
                                 internal_length=0)
                continue

            # Calculate the number of masked bases in the internal region
            internal_masked_count = sum(1 for base in internal_region 
                                      if base.islower() or base in 'NXnx')
            internal_mask_ratio = internal_masked_count / len(internal_region)

            print(f"{record.id} - Internal region mask ratio: {internal_mask_ratio:.2%}")
            
            if internal_mask_ratio < self.mask_threshold:
                # Remove masked portions from internal region
                internal_cleaned = self.remove_masked_portions(internal_region)
                
                # Reconstruct sequence: keep LTRs unchanged, clean internal region only
                cleaned_seq = ltr_5prime.upper() + internal_cleaned + ltr_3prime.upper()

                # Create new sequence record
                new_record = SeqRecord(
                    Seq(cleaned_seq),
                    id=record.id,
                    description=record.description
                )
                cleaned_records.append(new_record)
                self.add_log_entry(record.id, "PASS", "", 
                                 mask_ratio=internal_mask_ratio,
                                 internal_length=len(internal_region))
            else:
                print(f"  -> Excluded: Internal region {internal_mask_ratio:.1%} masked")
                self.add_log_entry(record.id, "FAIL", 
                                 f"Internal region mask ratio ({internal_mask_ratio:.3f}) exceeds threshold ({self.mask_threshold})",
                                 mask_ratio=internal_mask_ratio,
                                 internal_length=len(internal_region))

        # Check for unprocessed sequences (may be lost during RepeatMasker)
        missing_ids = original_ids - processed_ids
        for seq_id in missing_ids:
            self.add_log_entry(seq_id, "FAIL", "Sequence not found in masked file")

        # Write cleaned sequences
        if cleaned_records:
            SeqIO.write(cleaned_records, output_file, "fasta")
            print(f"\nCleaned sequences written to: {output_file}")

        # Write processing log
        self.write_processing_log()
        
        return len(cleaned_records) > 0

    def clean(self):
        """Execute the full cleaning process"""
        print("Starting nested insertion removal process...")
        print(f"Using LTR sequence library: {self.library_file}")
        print(f"Processing candidate sequences file: {self.candidates_file}")
        print(f"Masking threshold: {self.mask_threshold}\n")

        if self.run_repeatmasker():
            if self.process_masked_sequences():
                print("\nNested insertion removal completed successfully!")
                return True
        else:
            # If RepeatMasker fails, log all sequences
            print("RepeatMasker failed to run, logging all sequences...")
            for record in SeqIO.parse(self.candidates_file, "fasta"):
                self.add_log_entry(record.id, "FAIL", "RepeatMasker failed to run")
            self.write_processing_log()

        print("\nNested insertion removal failed!")
        return False

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


if __name__ == "__main__":
    genome_file = "/proj/nobackup/hpc2nstor2024-028/zhychen/bin/software/ltr_checker/test/test.fasta"
    candidates_file = "/proj/nobackup/hpc2nstor2024-028/zhychen/bin/software/ltr_checker/result/module5/filtered.fasta"
    library_file = "/proj/nobackup/hpc2nstor2024-028/zhychen/bin/software/ltr_checker/result/module6/LTR_library.fa"
    output_dir = "/proj/nobackup/hpc2nstor2024-028/zhychen/bin/software/ltr_checker/result/module6"

    dir = "/proj/nobackup/hpc2nstor2024-028/zhychen/bin/software/ltr_checker/result"
    dir_2 = dir + "/module2"
    dir_6 = dir + "/module6"
    dir_5 = dir + "/module5"

    with open("/proj/nobackup/hpc2nstor2024-028/zhychen/bin/software/ltr_checker/result/ltr_candidates.fasta", 'r') as f:
        content = f.read()
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
    
    
    
    
    cleaner = LTRCleaner(
        genome_file=genome_file,
        candidates_file=candidates_file,
        library_file=library_file,
        output_dir=output_dir
    )
    
    cleaner.clean()


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