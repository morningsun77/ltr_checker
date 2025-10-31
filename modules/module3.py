import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
import math

class Module3:
    """Processing LTR Insertion Time and Similarity Calculation"""
    def __init__(self, genome_file, input_file, output_dir, r=1.3e-8):
        self.genome_file = genome_file
        self.input_file = input_file 
        self.output_dir = output_dir
        self.results = []
        self.rate = r  # mutation rate
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
    def parse_ltr_info(self, header):
        """Parse LTR header information"""
        info = {}
        parts = header.strip().split()
        info['id'] = parts[0][1:] # Remove '>'
        
        for part in parts[1:]:
            key, value = part.split('=')
            info[key] = value
            
        return info
        
    def extract_sequences(self):
        """Extract 5'LTR and 3'LTR sequences from the genome"""
        genome = SeqIO.to_dict(SeqIO.parse(self.genome_file, "fasta"))
        
        ltr_pairs = []
        with open(self.input_file) as f:
            for line in f:
                if line.startswith('>'):
                    info = self.parse_ltr_info(line)

                    # Extract 5'LTR position
                    five_start, five_end = map(int, info['5ltr'].split(':'))
                    # Extract 3'LTR position
                    three_start, three_end = map(int, info['3ltr'].split(':'))

                    # Extract sequences from genome
                    chr_id = info['chr']
                    five_seq = str(genome[chr_id].seq[five_start-1:five_end])
                    three_seq = str(genome[chr_id].seq[three_start-1:three_end])
                    
                    ltr_pairs.append({
                        'header': line.strip(),
                        'info': info,
                        '5ltr': five_seq,
                        '3ltr': three_seq
                    })
                    
        return ltr_pairs

    def align_and_calculate(self, five_seq, three_seq):
        """Use MUSCLE for sequence alignment and calculate divergence rate"""
        # Get MUSCLE path
        try:
            MUSCLE_PATH = subprocess.check_output(['which', 'muscle'], 
                                            text=True).strip()
            print(f"Found MUSCLE at: {MUSCLE_PATH}")  # Debug info
        except subprocess.CalledProcessError:
            print("Could not find MUSCLE in PATH")
            return 0, 0
        
        if not os.path.exists(MUSCLE_PATH):
            print(f"MUSCLE executable not found at {MUSCLE_PATH}")
            return 0, 0

        # Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)

        # Write temporary file
        temp_fa = os.path.join(self.output_dir, "temp.fa")
        with open(temp_fa, 'w') as f:
            f.write(f">5LTR\n{five_seq}\n>3LTR\n{three_seq}\n")

        print(f"Input sequences written to: {temp_fa}")  # Debug info

        # Run MUSCLE
        temp_aln = os.path.join(self.output_dir, "temp.aln")
        cmd = [MUSCLE_PATH, '-align', temp_fa, '-output', temp_aln]  # MUSCLE 5.x syntax

        try:
            print(f"Executing command: {' '.join(cmd)}")  # Debug info
            result = subprocess.run(cmd,
                                capture_output=True,
                                text=True,
                                check=True)
            print("MUSCLE alignment completed")  # Debug info
            
        except subprocess.CalledProcessError as e:
            print(f"MUSCLE alignment failed with error:")
            print(f"stdout: {e.stdout}")
            print(f"stderr: {e.stderr}")

            # Try old MUSCLE syntax
            cmd_old = [MUSCLE_PATH, '-in', temp_fa, '-out', temp_aln]
            print(f"Trying old MUSCLE syntax: {' '.join(cmd_old)}")
            try:
                result = subprocess.run(cmd_old,
                                    capture_output=True,
                                    text=True,
                                    check=True)
                print("MUSCLE alignment completed with old syntax")
            except subprocess.CalledProcessError as e:
                print(f"Old syntax also failed:")
                print(f"stdout: {e.stdout}")
                print(f"stderr: {e.stderr}")
                if os.path.exists(temp_fa):
                    os.remove(temp_fa)
                return 0, 0

        # Check if output file exists
        if not os.path.exists(temp_aln):
            print(f"MUSCLE failed to generate alignment file: {temp_aln}")
            if os.path.exists(temp_fa):
                os.remove(temp_fa)
            return 0, 0
            
        try:
            # Read alignment results
            print(f"Reading alignment from: {temp_aln}")  # Debug info
            alignment = AlignIO.read(temp_aln, "fasta")
            seq1, seq2 = str(alignment[0].seq), str(alignment[1].seq)

            # Calculate valid sites and differences
            valid_sites = 0
            diff_sites = 0
            for i in range(len(seq1)):
                if seq1[i].upper() in 'ATGC' and seq2[i].upper() in 'ATGC':
                    valid_sites += 1
                    if seq1[i].upper() != seq2[i].upper():
                        diff_sites += 1

            print(f"Valid sites: {valid_sites}, Different sites: {diff_sites}")  # Debug info

        finally:
            # Clean up temporary files
            if os.path.exists(temp_fa):
                os.remove(temp_fa)
            if os.path.exists(temp_aln):
                os.remove(temp_aln)
            
        if valid_sites == 0:
            return 0, 0
                
        diff_rate = diff_sites / valid_sites
        similarity = 1 - diff_rate

        print(f"Similarity: {similarity:.4f}, Difference rate: {diff_rate:.4f}")  # Debug info

        return similarity, diff_rate
    
    def calculate_time(self, diff_rate):
        """Calculate insertion time"""
        try:
            # JC correction
            K = -3/4 * math.log(1 - 4/3 * diff_rate)
        except ValueError:
            return float('inf')

        # Calculate time, mutation rate is 1.3e-8
        r = self.rate  # Mutation rate
        time_my = K / (2 * r) / 1000000
        
        return time_my

    def process(self):
        """Process all sequences"""
        # Output files
        gff_file = os.path.join(self.output_dir, "results.gff")
        fasta_file = os.path.join(self.output_dir, "passed.fa")
        
        ltr_pairs = self.extract_sequences()
        
        with open(gff_file, 'w') as gff, open(fasta_file, 'w') as fa:
            for ltr in ltr_pairs:
                # Calculate similarity and difference rate
                similarity, diff_rate = self.align_and_calculate(
                    ltr['5ltr'], ltr['3ltr'])

                # Calculate insertion time
                insert_time = self.calculate_time(diff_rate)

                # Check filtering criteria
                status = "pass"
                filter_reason = []
                
                if similarity < 0.85 or similarity > 1.0:
                    status = "no_pass"
                    filter_reason.append("similarity_out_of_range")
                    
                if insert_time > 5:
                    status = "no_pass" 
                    filter_reason.append("too_old")

                # Write to GFF file
                gff.write(f"{ltr['info']['chr']}\tLTR_retriever\tLTR_retrotransposon\t"
                        f"{ltr['info']['coords']}\t.\t.\t.\t.\t"
                        f"ID={ltr['info']['id']};status={status};"
                        f"similarity={similarity:.3f};time={insert_time:.3f}MY;"
                        f"filter_reason={','.join(filter_reason)}\n")

                # If passed filtering, write complete FASTA record
                if status == "pass":
                    # Read full sequence from input file
                    with open(self.input_file) as input_fa:
                        write_seq = False
                        for line in input_fa:
                            if line.startswith('>'):
                                if line.strip().split()[0][1:] == ltr['info']['id']:
                                    fa.write(f"{line.strip()} sim={similarity:.3f}\n")
                                    write_seq = True
                                else:
                                    write_seq = False
                            elif write_seq:
                                fa.write(line)


    

if __name__ == "__main__":

    genome_file = "/proj/nobackup/hpc2nstor2024-028/zhychen/bin/software/ltr_checker/test/test.fasta"
    input_file = "/proj/nobackup/hpc2nstor2024-028/zhychen/bin/software/ltr_checker/result/module2/boundary_modified.fasta"
    output_dir = "/proj/nobackup/hpc2nstor2024-028/zhychen/bin/software/ltr_checker/result/module3/"

    ltr = Module3(genome_file, input_file, output_dir, r=1.3e-8)
    ltr.process()