#!/usr/bin/env python3
import os
import subprocess
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import re
import pandas as pd
from collections import defaultdict
import sys

class DomainAnnotator:
    def __init__(self, min_len=50):
        self.min_len = min_len
        # Get the directory where the current script is located
        current_dir = os.path.dirname(os.path.abspath(__file__))
        # Construct HMM file path (located in the dataset directory of the parent directory)
        self.hmm_db = os.path.abspath(os.path.join(
            current_dir,
            '..',
            'dataset',
            'REXdb_protein_database_viridiplantae_v3.0.hmm'
        ))
        # Store original sequence information
        self.original_seqs = {}
        # Store protein fragment position information
        self.fragment_positions = {}
        # Store LTR-TR position information in the genome
        self.ltr_tr_positions = {}
        
    def check_hmm_database(self):
        """Check if HMM database has been compressed and indexed"""
        required_extensions = ['.h3f', '.h3i', '.h3m', '.h3p']
        missing_files = []
        
        for ext in required_extensions:
            if not os.path.exists(self.hmm_db + ext):
                missing_files.append(self.hmm_db + ext)
                
        if missing_files:
            print("HMM database needs to be compressed and indexed.")
            print("Running hmmpress...")
            try:
                subprocess.run(['hmmpress', self.hmm_db], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error running hmmpress: {e}")
                raise
                
    def parse_ltr_tr_positions(self, record):
        """
        Parse LTR-TR position information from sequence ID or description
        Assume format is similar to: seq_id:start-end or seq_id_start_end
        """
        # Try different formats for parsing
        patterns = [
            r'(\w+):(\d+)-(\d+)',  # format: seq_id:start-end
            r'(\w+)_(\d+)_(\d+)',   # format: seq_id_start_end
            r'(\w+)\|(\d+)\|(\d+)', # format: seq_id|start|end
        ]
        
        for pattern in patterns:
            match = re.search(pattern, record.id)
            if not match:
                match = re.search(pattern, record.description)
            
            if match:
                seq_id = match.group(1)
                start = int(match.group(2))
                end = int(match.group(3))
                return seq_id, start, end

        # If unable to parse position information, return entire sequence as LTR-TR
        return record.id, 1, len(record.seq)
    
    def extract_ltr_tr_sequences(self, input_fasta):
        """
        Extract LTR-TR internal sequences from input sequences
        """
        ltr_tr_seqs = []
        
        for record in SeqIO.parse(input_fasta, "fasta"):
            # Parse LTR-TR positions
            seq_id, ltr_start, ltr_end = self.parse_ltr_tr_positions(record)

            # Store original sequence information
            self.original_seqs[record.id] = {
                'seq': str(record.seq),
                'length': len(record.seq),
                'description': record.description
            }

            # Store LTR-TR position information
            self.ltr_tr_positions[record.id] = {
                'genomic_id': seq_id,
                'ltr_start': ltr_start,
                'ltr_end': ltr_end,
                'ltr_length': ltr_end - ltr_start + 1
            }

            # Extract LTR-TR sequence (if position information is valid)
            if ltr_start > 0 and ltr_end <= len(record.seq):
                ltr_tr_seq = str(record.seq)[ltr_start-1:ltr_end]
            else:
                ltr_tr_seq = str(record.seq)
            
            ltr_tr_seqs.append((record.id, ltr_tr_seq))
            
        return ltr_tr_seqs
    
    def write_ltr_tr_sequences(self, ltr_tr_seqs, output_file):
        """Write LTR-TR sequences to file"""
        with open(output_file, 'w') as f:
            for seq_id, seq in ltr_tr_seqs:
                f.write(f">{seq_id}\n{seq}\n")
                
    def run_transeq(self, input_fasta, output_file):
        """Run EMBOSS transeq for 6-frame translation"""
        cmd = [
            'transeq',
            '-sequence', input_fasta,
            '-outseq', output_file,
            '-frame', '6',
            '-table', '1'
        ]
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            if result.stderr:
                print(f"Transeq stderr: {result.stderr}")
        except subprocess.CalledProcessError as e:
            print(f"Error running transeq: {e}")
            raise
            
    def run_hmmscan(self, protein_file, output_file):
        """Run hmmscan for domain annotation"""
        # First check HMM database
        self.check_hmm_database()

        # Check input file
        if not os.path.exists(protein_file):
            raise FileNotFoundError(f"Protein file not found: {protein_file}")

        # Check if input file is empty
        if os.path.getsize(protein_file) == 0:
            raise ValueError("Input protein file is empty")
            
        cmd = [
            'hmmscan',
            '--domtblout', output_file,
            '--cpu', '4',
            '-E', '1e-3',
            '--domE', '1e-3',
            '--noali',
            self.hmm_db,
            protein_file
        ]
        
        print(f"\nRunning command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )
            if result.stderr:
                print(f"Hmmscan stderr: {result.stderr}")
        except subprocess.CalledProcessError as e:
            print(f"Error running hmmscan: {e}")
            print(f"STDOUT: {e.stdout}")
            print(f"STDERR: {e.stderr}")
            raise
            
    def calculate_nucleotide_position(self, protein_pos, frame, ltr_length):
        """
        Calculate nucleotide position based on protein position and reading frame
        
        Args:
            protein_pos: Position in protein sequence (0-based)
            frame: Reading frame (1-6)
            ltr_length: Length of LTR-TR sequence
        
        Returns:
            (start, end, strand): Start and end positions in nucleotide sequence, and strand direction
        """
        # Frame 1-3: Positive strand; Frame 4-6: Negative strand
        if frame <= 3:
            # Positive strand
            strand = '+'
            frame_offset = frame - 1
            nucl_start = frame_offset + (protein_pos * 3)
            nucl_end = nucl_start + 3
        else:
            # Negative strand (needs reverse calculation)
            strand = '-'
            frame_offset = frame - 4
            # For negative strand, position needs to be calculated from the end of the sequence
            nucl_end = ltr_length - frame_offset - (protein_pos * 3)
            nucl_start = nucl_end - 3
            
        return nucl_start, nucl_end, strand
            
    def process_protein_sequences(self, trans_file):
        """Process translated protein sequences, ignoring stop codons"""
        processed_seqs = []
        seq_count = 0
        
        for record in SeqIO.parse(trans_file, "fasta"):
            seq_count += 1

            # Parse sequence ID and frame information
            try:
                original_id, frame = record.id.rsplit('_', 1)
                frame = int(frame)
            except ValueError:
                original_id = record.id
                frame = 1

            # Get LTR-TR length
            if original_id in self.ltr_tr_positions:
                ltr_length = self.ltr_tr_positions[original_id]['ltr_length']
            else:
                print(f"Warning: LTR-TR position info not found for {original_id}")
                continue
            
            # Keep the complete protein sequence intact without splitting (including stop codons)
            seq_str = str(record.seq)

            # If sequence length meets minimum length requirement
            if len(seq_str) >= self.min_len:
                new_id = f"{original_id}_frame{frame}"

                # Calculate nucleotide position in the genomic sequence
                nucl_start, _, strand = self.calculate_nucleotide_position(
                    0, frame, ltr_length
                )
                nucl_end_pos = len(seq_str) - 1
                _, nucl_end, _ = self.calculate_nucleotide_position(
                    nucl_end_pos, frame, ltr_length
                )

                # Ensure positions are correct (for negative strand, start should be less than end)
                if strand == '-':
                    nucl_start, nucl_end = min(nucl_start, nucl_end), max(nucl_start, nucl_end)

                # Store position information
                self.fragment_positions[new_id] = {
                    'original_id': original_id,
                    'frame': frame,
                    'strand': strand,
                    'protein_start': 0,
                    'protein_end': len(seq_str),
                    'nucl_start': nucl_start,
                    'nucl_end': nucl_end,
                    'fragment_length': len(seq_str)
                }
                
                processed_seqs.append((new_id, seq_str))
                        
        print(f"Processed {seq_count} sequences, generated {len(processed_seqs)} sequences")
        return processed_seqs
        
    def write_processed_seqs(self, processed_seqs, output_file):
        """Write processed sequences to file"""
        with open(output_file, 'w') as f:
            for seq_id, seq in processed_seqs:
                f.write(f">{seq_id}\n{seq}\n")

    def parse_hmmscan_output(self, hmmscan_output):
        """Parse hmmscan output file"""
        domain_hits = defaultdict(list)
        
        try:
            with open(hmmscan_output) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    fields = line.strip().split()
                    if len(fields) < 23:
                        continue
                        
                    hit = {
                        'query': fields[3],
                        'domain': fields[0],
                        'e_value': float(fields[11]),
                        'score': float(fields[13]),
                        'domain_start': int(fields[15]),
                        'domain_end': int(fields[16]),
                        'ali_start': int(fields[17]),
                        'ali_end': int(fields[18])
                    }
                    
                    domain_hits[fields[3]].append(hit)
        except Exception as e:
            print(f"Error parsing hmmscan output: {e}")
            raise
            
        return domain_hits

    def calculate_domain_positions(self, fragment_id, ali_start, ali_end):
        """
        Calculate domain positions in the genomic and LTR-TR sequences
        
        Args:
            fragment_id: Protein fragment ID
            ali_start: Domain start position in the protein sequence (1-based)
            ali_end: Domain end position in the protein sequence (1-based)
        
        Returns:
            dict: A dictionary containing genomic and LTR-TR position information
        """
        if fragment_id not in self.fragment_positions:
            return None
            
        frag_info = self.fragment_positions[fragment_id]
        original_id = frag_info['original_id']
        
        if original_id not in self.ltr_tr_positions:
            return None
            
        ltr_info = self.ltr_tr_positions[original_id]

        # Convert to 0-based
        ali_start_0 = ali_start - 1
        ali_end_0 = ali_end - 1

        # Calculate domain position in LTR-TR sequence
        if frag_info['strand'] == '+':
            ltr_domain_start = frag_info['nucl_start'] + (ali_start_0 * 3)
            ltr_domain_end = frag_info['nucl_start'] + ((ali_end_0 + 1) * 3) - 1
        else:
            # Negative strand requires reverse calculation
            ltr_domain_end = frag_info['nucl_end'] - (ali_start_0 * 3)
            ltr_domain_start = frag_info['nucl_end'] - ((ali_end_0 + 1) * 3) + 1
            ltr_domain_start, ltr_domain_end = min(ltr_domain_start, ltr_domain_end), max(ltr_domain_start, ltr_domain_end)

        # Calculate domain position in the genomic sequence
        genomic_domain_start = ltr_info['ltr_start'] + ltr_domain_start
        genomic_domain_end = ltr_info['ltr_start'] + ltr_domain_end
        
        return {
            'genomic_start': genomic_domain_start,  # 1-based
            'genomic_end': genomic_domain_end,      # 1-based
            'ltr_start': ltr_domain_start + 1,      # 1-based
            'ltr_end': ltr_domain_end + 1,          # 1-based
            'strand': frag_info['strand'],
            'frame': frag_info['frame']
        }

    def annotate_ltr_rt(self, input_fasta, output_prefix):
        """Main annotation process"""
        # Check input file
        if not os.path.exists(input_fasta):
            raise FileNotFoundError(f"Input file not found: {input_fasta}")
            
        print(f"Starting annotation process...")
        print(f"Input file size: {os.path.getsize(input_fasta)} bytes")

        # Create temporary files
        with tempfile.TemporaryDirectory() as tmpdir:
            ltr_tr_file = os.path.join(tmpdir, "ltr_tr_sequences.fasta")
            trans_file = os.path.join(tmpdir, "translated.fasta")
            processed_file = os.path.join(tmpdir, "processed.fasta")
            hmmscan_out = os.path.join(tmpdir, "hmmscan.out")

            # 1. Extract LTR-TR sequences
            print("\nExtracting LTR-TR sequences...")
            ltr_tr_seqs = self.extract_ltr_tr_sequences(input_fasta)
            self.write_ltr_tr_sequences(ltr_tr_seqs, ltr_tr_file)
            print(f"Extracted {len(ltr_tr_seqs)} LTR-TR sequences")

            # 2. Run 6-frame translation on LTR-TR sequences
            print("\nRunning 6-frame translation on LTR-TR sequences...")
            self.run_transeq(ltr_tr_file, trans_file)

            # 3. Process protein sequences (without stop codons)
            print("\nProcessing protein sequences...")
            processed_seqs = self.process_protein_sequences(trans_file)
            
            if not processed_seqs:
                print("Warning: No valid sequences found after processing")
                return
                
            print(f"\nWriting {len(processed_seqs)} processed sequences...")
            self.write_processed_seqs(processed_seqs, processed_file)

            # 4. Run hmmscan
            print("\nRunning hmmscan...")
            self.run_hmmscan(processed_file, hmmscan_out)

            # 5. Parse results
            print("\nParsing results...")
            domain_hits = self.parse_hmmscan_output(hmmscan_out)

            # 6. Write results
            output_file = f"{output_prefix}_domains.tsv"
            self.write_results(domain_hits, output_file)
            print(f"\nResults written to: {output_file}")

    def write_results(self, domain_hits, output_file):
        """Write annotation results to a file, including genomic and LTR-TR position information"""
        with open(output_file, 'w') as f:
            # Write detailed header
            headers = [
                "Sequence_ID", "Domain", "E-value", "Score",
                "Domain_start", "Domain_end", "Ali_start", "Ali_end",
                "Strand", "Frame", 
                "Genomic_start", "Genomic_end",
                "LTR_TR_start", "LTR_TR_end",
                "Fragment_ID"
            ]
            f.write("\t".join(headers) + "\n")

            # Write each hit's information
            hit_count = 0
            for query in domain_hits:
                for hit in domain_hits[query]:
                    # Get original sequence ID
                    fragment_id = hit['query']
                    if fragment_id in self.fragment_positions:
                        original_id = self.fragment_positions[fragment_id]['original_id']
                    else:
                        original_id = fragment_id.split('_frame')[0]

                    # Calculate position information
                    pos_info = self.calculate_domain_positions(
                        fragment_id, hit['ali_start'], hit['ali_end']
                    )
                    
                    if pos_info:
                        f.write(f"{original_id}\t{hit['domain']}\t{hit['e_value']}\t"
                               f"{hit['score']}\t{hit['domain_start']}\t{hit['domain_end']}\t"
                               f"{hit['ali_start']}\t{hit['ali_end']}\t"
                               f"{pos_info['strand']}\t{pos_info['frame']}\t"
                               f"{pos_info['genomic_start']}\t{pos_info['genomic_end']}\t"
                               f"{pos_info['ltr_start']}\t{pos_info['ltr_end']}\t"
                               f"{fragment_id}\n")
                    else:
                        # If position cannot be calculated, still output basic information
                        f.write(f"{original_id}\t{hit['domain']}\t{hit['e_value']}\t"
                               f"{hit['score']}\t{hit['domain_start']}\t{hit['domain_end']}\t"
                               f"{hit['ali_start']}\t{hit['ali_end']}\t"
                               f"NA\tNA\tNA\tNA\tNA\tNA\t{fragment_id}\n")
                    hit_count += 1
                    
            print(f"Wrote {hit_count} domain hits to {output_file}")

class LTRClassifier:
    """LTR retrotransposon classifier for processing and classifying LTR retrotransposon data"""
    
    def __init__(self):
        """Initialize LTR classifier"""
        pass
    
    def extract_domain_type(self, domain_name):
        """Extract domain type from domain name (GAG, PROT, INT, RT, RH)"""
        domain_name_upper = domain_name.upper()
        
        # Classify the type according to domain name characteristics
        if 'GAG' in domain_name_upper:
            return 'GAG'
        elif 'PROT' in domain_name_upper:
            return 'PROT'
        elif 'INT' in domain_name_upper:
            return 'INT'
        elif 'RH' in domain_name_upper:
            return 'RH'
        elif 'RT' in domain_name_upper:
            return 'RT'
        else:
            return None
    
    def check_overlap(self, start1, end1, start2, end2):
        """Check if two intervals overlap"""
        return not (end1 < start2 or end2 < start1)
    
    def classify_ltr(self, domains):
        """Classify LTR based on the order of domains"""
        # Get positions of INT, RT, RH
        int_pos = []
        rt_pos = []
        rh_pos = []
        
        for d in domains:
            if d['type'] == 'INT':
                int_pos.append(d['genomic_start'])
            elif d['type'] == 'RT':
                rt_pos.append(d['genomic_start'])
            elif d['type'] == 'RH':
                rh_pos.append(d['genomic_start'])

        # If none of these key domains are present, return unknown
        if not (int_pos and rt_pos):
            return 'unknown'

        # Get representative positions (take minimum)
        int_min = min(int_pos) if int_pos else float('inf')
        rt_min = min(rt_pos) if rt_pos else float('inf')
        rh_min = min(rh_pos) if rh_pos else float('inf')

        # Determine order
        # Copia: INT before RT
        # Gypsy: INT after RT
        if rt_min < int_min:
            return 'Gypsy'
        elif int_min < rt_min:
            return 'Copia'
        else:
            return 'unknown'
    
    def extract_ltr_number(self, seq_id):
        """Extract LTR number from Sequence_ID"""
        match = re.search(r'LTR_(\d+)', seq_id)
        if match:
            return int(match.group(1))
        return float('inf')  # If no match is found, place it at the end

    def process_ltr_data(self, data):
        """Process LTR data and classify"""
        results = []

        # Group by Sequence_ID
        grouped = data.groupby('Sequence_ID')
        
        for seq_id, group in grouped:
            # Extract all domain information
            domains = []
            domain_dict = defaultdict(list)
            
            for _, row in group.iterrows():
                domain_type = self.extract_domain_type(row['Domain'])
                if domain_type:
                    domain_info = {
                        'type': domain_type,
                        'e_value': float(row['E-value']),
                        'genomic_start': int(row['Genomic_start']),
                        'genomic_end': int(row['Genomic_end']),
                        'ltr_start': int(row['LTR_TR_start']),
                        'ltr_end': int(row['LTR_TR_end'])
                    }
                    domain_dict[domain_type].append(domain_info)

            # For each domain type, keep the one with the smallest e-value
            selected_domains = []
            for domain_type, domain_list in domain_dict.items():
                # Sort by e-value and select the smallest
                domain_list.sort(key=lambda x: x['e_value'])
                best_domain = domain_list[0]

                # Check for overlap with already selected domains
                overlap_with_rt = False
                overlap_with_others = []
                
                for selected in selected_domains:
                    if self.check_overlap(best_domain['genomic_start'], best_domain['genomic_end'],
                                   selected['genomic_start'], selected['genomic_end']):
                        if selected['type'] == 'RT':
                            overlap_with_rt = True
                        else:
                            overlap_with_others.append(selected)

                # If overlap with RT, prioritize keeping RT
                if overlap_with_rt and domain_type != 'RT':
                    continue

                # If current is RT and overlaps with others, remove others
                if domain_type == 'RT' and overlap_with_others:
                    for other in overlap_with_others:
                        selected_domains.remove(other)
                
                selected_domains.append(best_domain)

            # Check if the minimum requirement of 3 domains is met
            if len(selected_domains) < 3:
                continue

            # Sort by genomic position
            selected_domains.sort(key=lambda x: x['genomic_start'])

            # Classify
            classification = self.classify_ltr(selected_domains)

            # Build result
            result = {
                'Sequence_ID': seq_id,
                'Classification': classification,
                'Domain_count': len(selected_domains),
                'Domains': []
            }

            # Record the position information for each domain
            for domain in selected_domains:
                domain_info = f"{domain['type']}:{domain['ltr_start']}-{domain['ltr_end']}"
                result['Domains'].append(domain_info)
            
            result['Domain_positions'] = '; '.join(result['Domains'])
            results.append(result)

        # Sort by LTR number
        results.sort(key=lambda x: self.extract_ltr_number(x['Sequence_ID']))
        
        return results
    
    def process_file(self, input_file, output_file=None):
        """Process input file and generate results"""
        try:
            # Read input file
            df = pd.read_csv(input_file, sep='\t')
            print(f"Successfully read file: {input_file}")
            print(f"Number of rows: {len(df)}")

            # Process data (results are already sorted by LTR number)
            results = self.process_ltr_data(df)

            # Create output DataFrame
            if results:
                output_df = pd.DataFrame(results)[['Sequence_ID', 'Classification', 'Domain_count', 'Domain_positions']]

                # Save to file
                if output_file:
                    output_df.to_csv(output_file, sep='\t', index=False)
                    print(f"Results saved to: {output_file}")

                return output_df
            else:
                print("No sequences met the criteria (at least 3 domains required)")
                return None
                
        except Exception as e:
            print(f"Error: {str(e)}")
            return None



def main():
    ## Module4-1
    # Set input and output files
    input_fasta = "/proj/nobackup/hpc2nstor2024-028/zhychen/bin/software/ltr_checker/result/module3/passed.fa"
    output_prefix = "/proj/nobackup/hpc2nstor2024-028/zhychen/bin/software/ltr_checker/result/module4/ltr_rt_annotation"

    # Check input file
    if not os.path.exists(input_fasta):
        raise FileNotFoundError(f"Input file not found: {input_fasta}")
        
    print(f"Input file: {input_fasta}")
    print(f"File size: {os.path.getsize(input_fasta)} bytes")

    # Create annotator instance and run annotation
    annotator = DomainAnnotator()
    annotator.annotate_ltr_rt(input_fasta, output_prefix)

    
    ## Module4-2
    # Set input and output files
    # Create classifier instance
    classifier = LTRClassifier()
    
    input_file = '/proj/nobackup/hpc2nstor2024-028/zhychen/bin/software/ltr_checker/result/module4/ltr_rt_annotation_domains.tsv'
    output_file = '/proj/nobackup/hpc2nstor2024-028/zhychen/bin/software/ltr_checker/result/module4/domain_analysis_results.tsv'
    
    try:
        # Read input file
        df = pd.read_csv(input_file, sep='\t')
        print(f"Successfully read file: {input_file}")
        print(f"Number of rows: {len(df)}")

        # Process data (results are already sorted by LTR number)
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
            if output_file:
                output_df.to_csv(output_file, sep='\t', index=False)
                print(f"\nResults saved to: {output_file}")
            else:
                print("\nOutput table:")
                print(output_df.to_string(index=False))

            # Statistics
            print(f"\nStatistics:")
            print(f"Total sequences: {len(results)}")
            print(f"Copia: {sum(1 for r in results if r['Classification'] == 'Copia')}")
            print(f"Gypsy: {sum(1 for r in results if r['Classification'] == 'Gypsy')}")
            print(f"Unknown: {sum(1 for r in results if r['Classification'] == 'unknown')}")
            print(f"\nNote: Results are sorted in ascending order by LTR number (e.g., LTR_1, LTR_2, LTR_10...)")
            
        else:
            print("No sequences met the criteria (at least 3 domains required)")
            if output_file:
                # Create empty file
                with open(output_file, 'w') as f:
                    f.write("Sequence_ID\tClassification\tDomain_count\tDomain_positions\n")
                    
    except FileNotFoundError:
        print(f"Error: File not found {input_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: An error occurred while processing the file - {str(e)}")
        sys.exit(1)
    
if __name__ == "__main__":
    main()