#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
LTR Filtering Tool
Used for filtering LTR sequences, based on BLAST alignment results for analysis
"""

import os
import subprocess
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline

class LTRFilter:
    def __init__(self, input_file, output_dir):
        """
        Initialize LTR Filter

        Args:
            input_file (str): Path to the input FASTA file
            output_dir (str): Path to the output directory
        """
        self.input_file = input_file
        self.output_dir = output_dir

        # Get the current script directory
        current_dir = os.path.dirname(os.path.abspath(__file__))
        # Database file directory
        dataset_dir = os.path.join(current_dir, '..', 'dataset')

        # Set database file paths
        self.db_paths = {
            'DNA': os.path.join(dataset_dir, 'Tpases020812DNA.txt'),
            'LINE': os.path.join(dataset_dir, 'Tpases020812LINE.txt'),
            'MAKER': os.path.join(dataset_dir, 'alluniRefprexp082813.txt')
        }
        
        self.results = {}
        self.filtered_sequences = []
        self.sequence_lengths = self._get_sequence_lengths()
        
        # ===== 新增：预加载LTR坐标 =====
        self.ltr_coordinates = self._extract_all_ltr_coordinates()
        print(f"Parsed LTR coordinates for {len(self.ltr_coordinates)} sequences")

    def _extract_all_ltr_coordinates(self):
        """
        Preload LTR coordinate information for all sequences

        Returns:
            dict: {seq_id: {'ltr5': (start, end), 'ltr3': (start, end)}}
        """
        all_coords = {}
        
        try:
            import re
            for record in SeqIO.parse(self.input_file, "fasta"):
                desc = record.description

                # Parse location: 5ltr=106522:109606 3ltr=114996:118080
                ltr5_match = re.search(r'5ltr=(\d+):(\d+)', desc)
                ltr3_match = re.search(r'3ltr=(\d+):(\d+)', desc)
                
                if ltr5_match and ltr3_match:
                    all_coords[record.id] = {
                        'ltr5': (int(ltr5_match.group(1)) - 1, int(ltr5_match.group(2))),
                        'ltr3': (int(ltr3_match.group(1)) - 1, int(ltr3_match.group(2)))
                    }
        except Exception as e:
            print(f"Error occurred while extracting coordinates: {str(e)}")

        return all_coords

    def get_ltr_coordinates(self, seq_id):
        """
        Get LTR Coordinates for Specified Sequence
        """
        return self.ltr_coordinates.get(seq_id, None)

    def _get_sequence_lengths(self):
        """
        Get Length Information for Input Sequences

        Returns:
            dict: Mapping of sequence ID to sequence length
        """
        lengths = {}
        try:
            for record in SeqIO.parse(self.input_file, "fasta"):
                lengths[record.id] = len(record.seq)
        except Exception as e:
            print(f"Error occurred while reading sequence lengths: {str(e)}")
        return lengths
    
    def make_blast_db(self):
        """
        Create BLAST databases for each database
        """
        for db_type, db_path in self.db_paths.items():
            if not os.path.exists(db_path):
                raise FileNotFoundError(f"Database file not found: {db_path}")

            db_name = os.path.join(self.output_dir, f"{db_type}_db")
            cmd = f"makeblastdb -in {db_path} -dbtype prot -out {db_name}"
            print(f"Running command: {cmd}")

            try:
                result = subprocess.run(cmd, shell=True, check=True,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     text=True)
                print(f"Database {db_type} created successfully")
            except subprocess.CalledProcessError as e:
                print(f"Failed to create database: {e.stderr}")
                raise
    
    def run_blast(self, query, db_type, output):
        """
        Run BLAST search (updated version, including location information)
        
        Args:
            query (str): Query Sequence File Path
            db_type (str): Database Type
            output (str): Output File Path
        """
        db_path = os.path.join(self.output_dir, f"{db_type}_db")
        
        if not os.path.exists(f"{db_path}.phr"):
            raise FileNotFoundError(f"BLAST database file not found: {db_path}.phr")

        print(f"Running BLAST search: {db_type}")

    def process_blast_results(self, blast_output):
        """
        Process BLAST output file, including location information

        Args:
            blast_output (str): BLAST output file path

        Returns:
            dict: Processed BLAST results, including location information
        """
        hits = {}
        
        # 首先更新BLAST命令的输出格式
        # outfmt="6 qseqid sseqid pident length qstart qend evalue"
        
        with open(blast_output) as f:
            for line in f:
                if not line.strip():
                    continue
                    
                parts = line.strip().split('\t')
                
                # 根据您提供的样例数据，格式应该是：
                # qseqid sseqid pident length qstart qend evalue
                if len(parts) < 7:
                    print(f"警告：行格式不完整: {line}")
                    continue
                    
                qseqid = parts[0]
                sseqid = parts[1]
                pident = float(parts[2])
                length = int(parts[3])
                qstart = int(parts[4])
                qend = int(parts[5])
                evalue = float(parts[6])
                
                # Filter Unqualified Alignments
                if pident < 30 or length < 30:  # 30 amino acids = 90bp
                    continue
                    
                if evalue > 0.001:
                    continue
                
                if qseqid not in hits:
                    hits[qseqid] = []
                    
                # Handle Strand Orientation
                # Ensure start < end, record strand direction
                if qstart <= qend:
                    actual_start = qstart
                    actual_end = qend
                    strand = '+'
                else:
                    actual_start = qend
                    actual_end = qstart
                    strand = '-'
                
                hits[qseqid].append({
                    'subject': sseqid,
                    'length': length,
                    'pident': pident,
                    'evalue': evalue,
                    'q_start': actual_start,
                    'q_end': actual_end,
                    'q_start_raw': qstart,  # Keep Original Value
                    'q_end_raw': qend,      # Keep Original Value
                    'strand': strand,
                    'align_length_nt': length * 3  # Nucleotide Length
                })
                
        return hits

    def calculate_coverage(self, hits):
        """
        Calculate the non-overlapping coverage of the sequence

        Args:
            hits (list): BLAST alignment results list

        Returns:
            tuple: (Total Coverage Length, Coverage Regions List)
        """
        if not hits:
            return 0, []
        
        # Extract all alignment regions
        regions = []
        for hit in hits:
            regions.append({
                'start': hit['q_start'],
                'end': hit['q_end'],
                'subject': hit.get('subject', 'unknown'),
                'identity': hit['pident']
            })
        
        # Merge Overlapping Regions
        merged_regions = self.merge_overlapping_regions(regions)
        
        # Calculate Total Coverage Length
        total_coverage = sum(region['end'] - region['start'] for region in merged_regions)
        
        return total_coverage, merged_regions
    
    def merge_overlapping_regions(self, regions):
        """
        Merge Overlapping Regions, ensuring each position is only counted once

        Args:
            regions (list): List of regions, each containing start, end, etc.

        Returns:
            list: Merged regions list
        """
        if not regions:
            return []
        
        # Sort by start position
        sorted_regions = sorted(regions, key=lambda x: x['start'])
        
        merged = []
        current = sorted_regions[0].copy()
        
        for region in sorted_regions[1:]:
            if region['start'] <= current['end']:
                # Regions overlap, merge them
                current['end'] = max(current['end'], region['end'])
                # Record all involved subjects
                if 'merged_from' not in current:
                    current['merged_from'] = [current['subject']]
                current['merged_from'].append(region['subject'])
            else:
                # No overlap, save current region and start a new one
                merged.append(current)
                current = region.copy()

        # Add the last region
        merged.append(current)
        
        return merged
    
    def check_region_coverage(self, seq_id, seq_length, all_hits, ltr_coords):
        """
        Check protein-coding coverage in different regions of LTR-RT

        Args:
            seq_id: Sequence ID
            seq_length: Total length of the sequence
            all_hits: All BLAST hit results
            ltr_coords: LTR coordinates

        Returns:
            (bool, str): (Whether it passed the check, Failure reason)
        """
        if not ltr_coords:
            return True, "No LTR coordinates available"
        
        # Define three regions
        regions = {
            '5_LTR': (ltr_coords['ltr5'][0] - ltr_coords['ltr5'][0] + 1, ltr_coords['ltr5'][1] - ltr_coords['ltr5'][0]),
            'internal': (ltr_coords['ltr5'][1] - ltr_coords['ltr5'][0] + 1, ltr_coords['ltr3'][0] - ltr_coords['ltr5'][0]),
            '3_LTR': (ltr_coords['ltr3'][0] - ltr_coords['ltr5'][0] + 1, ltr_coords['ltr3'][1] - ltr_coords['ltr5'][0])
        }
        
        # Check each region
        for region_name, (region_start, region_end) in regions.items():
            region_length = region_end - region_start
            if region_length <= 0:
                continue
            # Collect hits that fall within this region
            region_hits = []
            for hit in all_hits:
                hit_start = hit['q_start']
                hit_end = hit['q_end']

                # Calculate overlap
                overlap_start = max(hit_start, region_start)
                overlap_end = min(hit_end, region_end)
                
                if overlap_start < overlap_end:
                    region_hits.append({
                        'start': overlap_start,
                        'end': overlap_end,
                        'subject': hit.get('subject', 'unknown'),
                        'identity': hit.get('pident', 0)
                    })

            # Calculate coverage for this region
            if region_hits:
                merged = self.merge_overlapping_regions(region_hits)
                region_coverage = sum(r['end'] - r['start'] for r in merged)
                coverage_ratio = region_coverage / region_length
                
                print(f"    {region_name}: {region_coverage}bp / {region_length}bp = {coverage_ratio*100:.1f}%")

                # Apply 70% rule
                if coverage_ratio > 0.7:
                    return False, f"{region_name} is {coverage_ratio*100:.0f}% protein-coding"
        
        return True, "Pass region check"
    
    def write_annotation(self, output_file):
        """
        Write annotation file 

        Args:
            output_file (str): Output file path
        """
        with open(output_file, 'w') as f:
            f.write("#Sequence_ID\tStatus\tReason\tLength\tCoverage\tHas_LTR_coords\n")
            for seq_id in self.sequence_lengths.keys():
                status = "pass" if seq_id in self.filtered_sequences else "no_pass"
                reason = self.results.get(seq_id, {}).get('reason', 'NA')
                length = self.sequence_lengths[seq_id]
                coverage = self.results.get(seq_id, {}).get('coverage', 0)
                has_coords = "Yes" if seq_id in self.ltr_coordinates else "No"
                f.write(f"{seq_id}\t{status}\t{reason}\t{length}\t{coverage}\t{has_coords}\n")

    def write_fasta(self, output_file):
        """
        Write filtered FASTA file

        Args:
            output_file (str): Output file path
        """
        with open(output_file, 'w') as out_f:
            for record in SeqIO.parse(self.input_file, "fasta"):
                if record.id in self.filtered_sequences:
                    SeqIO.write(record, out_f, "fasta")

    def run(self):
        """
        Run main process (enhanced version, including region check)
        """
        print("Starting LTR filtering process...")
        
        # 1. Create BLAST database
        print("\nCreating BLAST database...")
        self.make_blast_db()

        # Store all hits from all databases for comprehensive analysis
        all_seq_hits = {}

        # 2. Run BLAST for each database
        for db_type in self.db_paths.keys():
            blast_output = os.path.join(self.output_dir, f"blast_{db_type}.txt")
            print(f"\nProcessing database: {db_type}")

            try:
                self.run_blast(self.input_file, db_type, blast_output)
                
                if os.path.exists(blast_output):
                    hits = self.process_blast_results(blast_output)

                    # Collect all hits for subsequent region analysis
                    for seq_id, seq_hits in hits.items():
                        if seq_id not in all_seq_hits:
                            all_seq_hits[seq_id] = []
                        all_seq_hits[seq_id].extend(seq_hits)

                    # Show statistics for individual database
                    for seq_id, seq_hits in hits.items():
                        seq_length = self.sequence_lengths.get(seq_id, 0)
                        coverage, merged_regions = self.calculate_coverage(seq_hits)

                        print(f"\n  Sequence {seq_id} in {db_type} database:")
                        print(f"    Number of matching regions: {len(seq_hits)}")
                        print(f"    Coverage length: {coverage} bp ({coverage/seq_length*100:.1f}%)")

            except Exception as e:
                print(f"Error processing database {db_type}: {str(e)}")
                import traceback
                traceback.print_exc()
                continue

        # 3. Comprehensive analysis for each sequence (including region check)
        print("\n\n================== Comprehensive Analysis ==================")
        for seq_id in self.sequence_lengths.keys():
            seq_length = self.sequence_lengths[seq_id]
            print(f"\nAnalyzing sequence {seq_id} (length: {seq_length} bp):")

            if seq_id in all_seq_hits:
                all_hits = all_seq_hits[seq_id]
                total_coverage, merged_regions = self.calculate_coverage(all_hits)

                print(f"  Total coverage: {total_coverage} bp ({total_coverage/seq_length*100:.1f}%)")

                # Check 1: Total coverage length > 1000bp
                if total_coverage > 1000:
                    self.results[seq_id] = {
                        'pass': False,
                        'reason': f">{1000}bp ({total_coverage}bp) matched to non-LTR databases",
                        'coverage': total_coverage
                    }
                    print(f"  ✗ Failed: Exceeded 1000bp ({total_coverage}bp) matched to non-LTR databases")
                    continue

                # Check 2: Total coverage ratio > 70%
                if seq_length > 0 and total_coverage > seq_length * 0.7:
                    self.results[seq_id] = {
                        'pass': False,
                        'reason': f">70% ({total_coverage/seq_length*100:.1f}%) matched to non-LTR databases",
                        'coverage': total_coverage
                    }
                    print(f"  ✗ Failed: Exceeded 70% ({total_coverage/seq_length*100:.1f}%) matched to non-LTR databases")
                    continue

                # Check 3: Region-specific check (new)
                ltr_coords = self.get_ltr_coordinates(seq_id)
                if ltr_coords:
                    print(f"  Performing region check...")
                    passed, reason = self.check_region_coverage(seq_id, seq_length, all_hits, ltr_coords)
                    if not passed:
                        self.results[seq_id] = {
                            'pass': False,
                            'reason': reason,
                            'coverage': total_coverage
                        }
                        print(f"  ✗ Failed: {reason}")
                        continue
                    else:
                        print(f"  ✓ Passed region check")
                else:
                    print(f"  ! Warning: Unable to retrieve LTR coordinates, skipping region check")

                # Passed all checks
                print(f"  ✓ Passed all checks")
            else:
                print(f"  ✓ No matches (passed)")

            # If passed all checks, add to filtered sequences list
            if seq_id not in self.results:
                self.filtered_sequences.append(seq_id)
        
        # 4. Output Results
        print("\n\n================== Output Results ==================")
        annotation_file = os.path.join(self.output_dir, "annotation.txt")
        fasta_file = os.path.join(self.output_dir, "filtered.fasta")
        
        self.write_annotation(annotation_file)
        self.write_fasta(fasta_file)

        # Output statistics
        print(f"\nProcessing complete!")
        print(f"Annotation file: {annotation_file}")
        print(f"FASTA file: {fasta_file}")
        print(f"\nStatistics:")
        print(f"  Input sequences: {len(self.sequence_lengths)}")
        print(f"  Passed sequences: {len(self.filtered_sequences)}")
        print(f"  Failed sequences: {len(self.results)}")

        # Output failure reason statistics
        if self.results:
            print("\nFailure reason classification:")
            reason_counts = {}
            for seq_id, info in self.results.items():
                reason = info['reason']
                # Simplified reason classification
                if '>1000bp' in reason:
                    key = '>1000bp matched'
                elif '>70%' in reason:
                    key = '>70% matched (whole sequence)'
                elif 'LTR' in reason and 'protein-coding' in reason:
                    key = 'LTR region is protein-coding'
                elif 'internal' in reason:
                    key = 'Internal region is protein-coding'
                else:
                    key = 'Other'
                reason_counts[key] = reason_counts.get(key, 0) + 1
            
            for reason, count in sorted(reason_counts.items()):
                print(f"  - {reason}: {count}")

def main():
    """
    Main function
    """
    input_file = "/proj/nobackup/hpc2nstor2024-028/zhychen/bin/software/ltr_checker/result/module4/output_sequences.fasta"
    output_dir = "/proj/nobackup/hpc2nstor2024-028/zhychen/bin/software/ltr_checker/result/module5"
    
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        filter = LTRFilter(input_file, output_dir)
        filter.run()
    except Exception as e:
        print(f"Error occurred while running the program: {str(e)}")

if __name__ == "__main__":
    main()