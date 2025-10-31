#!/usr/bin/env python3

import os
import tempfile
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

class Module2:
    def __init__(self, genome_file, search_range=15, blast_path="blastn", makeblastdb_path="makeblastdb"):
        """Initialization"""
        self.genome_dict = self._load_genome(genome_file)
        self.search_range = search_range
        self.blast_path = blast_path
        self.makeblastdb_path = makeblastdb_path
        
        # Define motif priority - TGCA has highest priority
        self.motif_priority = ["TGCA", "TGCT", "TACA", "TACT", "TGGA", "TATA", "TGTA", "TCCA"]

        # Check BLAST tool availability
        self._check_blast_availability()

    def _check_blast_availability(self):
        """Check BLAST tool availability"""
        try:
            subprocess.run([self.blast_path, "-version"], 
                         capture_output=True, check=True)
            subprocess.run([self.makeblastdb_path, "-version"], 
                         capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            raise RuntimeError("BLAST+ tools not found. Please install BLAST+ or provide correct paths.")

    def _load_genome(self, genome_file):
        """Load genome sequences"""
        genome_dict = {}
        for record in SeqIO.parse(genome_file, "fasta"):
            genome_dict[record.id] = str(record.seq)
        return genome_dict

    def _parse_fasta_header(self, header):
        """Parse FASTA header information"""
        info = {}
        parts = header.split()
        info['id'] = parts[0].replace('>', '')
        
        for part in parts:
            if 'chr=' in part:
                info['chr'] = part.split('=')[1]
            elif 'coords=' in part:
                coords = part.split('=')[1].split(':')
                info['start'] = int(coords[0])
                info['end'] = int(coords[1])
            elif '5ltr=' in part:
                coords = part.split('=')[1].split(':') 
                info['5ltr_start'] = int(coords[0])
                info['5ltr_end'] = int(coords[1])
            elif '3ltr=' in part:
                coords = part.split('=')[1].split(':')
                info['3ltr_start'] = int(coords[0])
                info['3ltr_end'] = int(coords[1])
        return info
    
    ##  ========== Module 2-1 ==========

    def _enhanced_blast_ltr_alignment(self, ltr5_seq, ltr3_seq):
        """Use real BLAST to align 5' and 3' LTR sequences - Core functionality of Module 2"""
        if len(ltr5_seq) == 0 or len(ltr3_seq) == 0:
            return {
                'has_alignment': False,
                'identity': 0.0,
                'alignment_pairs': 0,
                'boundary_shift': 0,
                'alignments': []
            }
        
        # Check if sequence lengths are sufficient for alignment
        if len(ltr5_seq) < 50 or len(ltr3_seq) < 50:
            return {
                'has_alignment': False,
                'identity': 0.0,
                'alignment_pairs': 1,
                'boundary_shift': 0,
                'alignments': []
            }
        
        # Create temporary files
        query_file = None
        subject_file = None
        
        try:
            query_file = self._create_temp_fasta([ltr5_seq], "ltr5_query")
            subject_file = self._create_temp_fasta([ltr3_seq], "ltr3_subject")

            # Run BLAST alignment
            blast_output = self._run_blastn(
                query_file, subject_file, 
                evalue=1e-5, word_size=7, max_target_seqs=50
            )

            # Parse BLAST results
            alignments = self._parse_blast_output(blast_output)
            
            if not alignments:
                return {
                    'has_alignment': False,
                    'identity': 0.0,
                    'alignment_pairs': 0,
                    'boundary_shift': 0,
                    'alignments': []
                }
            
            # Analyze alignment results
            significant_alignments = 0
            best_length = 0 
            
            for alignment in alignments:  
                significant_alignments += 1
                
                if alignment['length'] > best_length:
                    best_identity = alignment['identity']
                    best_length = alignment['length']
                    # Calculate boundary offset
                    best_shift = max(abs(alignment['q_start'] - 1), abs(alignment['s_start'] - 1))
                    
            return {
                'has_alignment': significant_alignments > 0,
                'identity': best_identity / 100.0,  # Convert to 0-1 range
                'alignment_pairs': significant_alignments,
                'boundary_shift': best_shift,
                'alignments': alignments
            }
            
        finally:
            # Clean up temporary files
            for temp_file in [query_file, subject_file]:
                if temp_file and os.path.exists(temp_file):
                    os.unlink(temp_file)
    
    def _create_temp_fasta(self, sequences, prefix="temp_seq"):
        """Create temporary FASTA file"""
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', 
                                              prefix=prefix, delete=False)
        for i, seq in enumerate(sequences):
            if seq and len(seq) > 0:
                temp_file.write(f">seq_{i}\n{seq}\n")
        temp_file.close()
        return temp_file.name
    
    def _run_blastn(self, query_file, subject_file, evalue=10, word_size=11, 
                   max_target_seqs=10, outfmt=6):
        """Run blastn alignment"""
        try:
            # Create database for subject sequence
            db_file = subject_file + ".db"
            subprocess.run([
                self.makeblastdb_path, "-in", subject_file, 
                "-dbtype", "nucl", "-out", db_file
            ], capture_output=True, check=True)

            # Run blastn
            result = subprocess.run([
                self.blast_path, "-query", query_file, "-db", db_file,
                "-evalue", str(evalue), "-word_size", str(word_size),
                "-max_target_seqs", str(max_target_seqs),
                "-outfmt", str(outfmt), "-dust", "no"
            ], capture_output=True, text=True, check=True)

            # Clean up database files
            for ext in [".nhr", ".nin", ".nsq"]:
                db_ext_file = db_file + ext
                if os.path.exists(db_ext_file):
                    os.unlink(db_ext_file)
            
            return result.stdout
            
        except subprocess.CalledProcessError as e:
            print(f"BLAST error: {e}")
            return ""
    
    def _parse_blast_output(self, blast_output):
        """Parse BLAST output"""
        alignments = []
        if not blast_output.strip():
            return alignments
            
        for line in blast_output.strip().split('\n'):
            if line.strip():
                fields = line.split('\t')
                if len(fields) >= 12:
                    alignment = {
                        'query_id': fields[0],
                        'subject_id': fields[1],
                        'identity': float(fields[2]),
                        'length': int(fields[3]),
                        'mismatches': int(fields[4]),
                        'gaps': int(fields[5]),
                        'q_start': int(fields[6]),
                        'q_end': int(fields[7]),
                        's_start': int(fields[8]),
                        's_end': int(fields[9]),
                        'evalue': float(fields[10]),
                        'bitscore': float(fields[11])
                    }
                    alignments.append(alignment)
        return alignments
    
    def _find_2bp_palindromic_motif(self, seq, ltr5_start, ltr5_end, ltr3_start, ltr3_end, search_range=100):
        """Search for 2bp palindromic motifs near internal boundaries

        For 5' LTR: Search for ...CA at its 3' end (right boundary)
        For 3' LTR: Search for TG... at its 5' end (left boundary)
        """
        # Search 3' end of 5' LTR (internal boundary)
        end_search_begin = max(0, ltr5_end - search_range)
        end_search_end = min(len(seq), ltr5_end + search_range)
        end_region = seq[end_search_begin:end_search_end]
        
        ca_positions = []
        for i in range(len(end_region) - 1):
            if end_region[i:i+2].upper() == "CA":
                # CA's end position as new end position for 5' LTR
                ca_positions.append(end_search_begin + i + 2)

        # Search 5' end of 3' LTR (internal boundary)
        start_search_begin = max(0, ltr3_start - search_range)
        start_search_end = min(len(seq), ltr3_start + search_range)
        start_region = seq[start_search_begin:start_search_end]
        
        tg_positions = []
        for i in range(len(start_region) - 1):
            if start_region[i:i+2].upper() == "TG":
                # TG's start position as new start position for 3' LTR
                tg_positions.append(start_search_begin + i)

        # Return positions closest to original predictions
        new_5ltr_end = None
        new_3ltr_start = None
        
        if ca_positions:
            new_5ltr_end = min(ca_positions, key=lambda x: abs(x - ltr5_end))
            
        if tg_positions:
            new_3ltr_start = min(tg_positions, key=lambda x: abs(x - ltr3_start))
        
        return new_5ltr_end, new_3ltr_start
    
    def _check_tsd(self, chr_seq, genome_start5, genome_end5, genome_start3, genome_end3):
        """Check TSD - Outside of LTR-RT elements in genomic sequence

        Args:
            chr_seq: Chromosome/genomic sequence
            genome_start5: 5' LTR start position in genome (absolute coordinates)
            genome_end5: 5' LTR end position in genome (absolute coordinates)
            genome_start3: 3' LTR start position in genome (absolute coordinates)
            genome_end3: 3' LTR end position in genome (absolute coordinates)
        """
        # TSD is located outside of the entire LTR-RT
        # Upstream TSD: 8bp before 5' LTR start position
        upstream_start = max(0, genome_start5 - 8)
        upstream = chr_seq[upstream_start:genome_start5]
        
        # Downstream TSD: 8bp after 3' LTR end position
        downstream_end = min(len(chr_seq), genome_end3 + 8)
        downstream = chr_seq[genome_end3:downstream_end]
        
        size_priority = [5, 4, 6]
        
        best_match = None
        best_mismatch = float('inf')
        
        for size in size_priority:
            if size > len(upstream) or size > len(downstream):
                continue

            # Search upstream region from right to left (close to LTR)
            for i in range(len(upstream) - size, -1, -1):
                upstream_tsd = upstream[i:i+size]

                # Search downstream region from left to right (close to LTR)
                for j in range(len(downstream) - size + 1):
                    downstream_tsd = downstream[j:j+size]
                    mismatch_count = self._calculate_mismatch(upstream_tsd, downstream_tsd)
                    
                    if mismatch_count <= 1:
                        if (mismatch_count < best_mismatch or 
                            (mismatch_count == best_mismatch and best_match is None)):
                            best_match = {
                                'upstream_tsd': upstream_tsd,
                                'downstream_tsd': downstream_tsd,
                                'size': size,
                                'mismatch': mismatch_count,
                                'upstream_pos': upstream_start + i,  # Absolute position in genome
                                'downstream_pos': genome_end3 + j    # Absolute position in genome
                            }
                            best_mismatch = mismatch_count
                        
                        if size == 5 and mismatch_count == 0:
                            return True, best_match
            
            if best_match and best_match['size'] == size:
                break
        
        if best_match:
            return True, best_match
        else:
            return False, None

    def find_boundary(self, seq, ltr5_start, ltr5_end, ltr3_start, ltr3_end, chr_seq, chr_start):
        """Enhanced boundary detection using real BLAST alignment - Module 2"""
        
        # MODULE 2 Step 1: Perform inter-LTR alignment using real BLAST
        ltr5_seq = seq[ltr5_start:ltr5_end] if ltr5_start < ltr5_end else ""
        ltr3_seq = seq[ltr3_start:ltr3_end] if ltr3_start < ltr3_end else ""
        
        alignment_result = self._enhanced_blast_ltr_alignment(ltr5_seq, ltr3_seq)

        # Apply Module 2 exclusion criteria
        if not alignment_result['has_alignment']:
            return (ltr5_start, ltr5_end, ltr3_start, ltr3_end), None, "no_self_alignment"
        
        if alignment_result['boundary_shift'] > 100:
            return (ltr5_start, ltr5_end, ltr3_start, ltr3_end), None, "boundary_shift_excessive"
        
        if alignment_result['alignment_pairs'] >= 8:
            return (ltr5_start, ltr5_end, ltr3_start, ltr3_end), None, "heavily_nested_insertions"

        # MODULE 2 Step 2: Search for 2bp palindromic motif at internal boundaries
        new_5ltr_end, new_3ltr_start = self._find_2bp_palindromic_motif(
            seq, ltr5_start, ltr5_end, ltr3_start, ltr3_end, 100)

        # If new internal boundaries are found, use them
        if new_5ltr_end is not None:
            ltr5_end = new_5ltr_end
        if new_3ltr_start is not None:
            ltr3_start = new_3ltr_start
        
        # Calculate genomic absolute coordinates
        genome_start5 = chr_start + ltr5_start
        genome_end5 = chr_start + ltr5_end
        genome_start3 = chr_start + ltr3_start
        genome_end3 = chr_start + ltr3_end

        # Check TSD using genomic sequence
        has_tsd, tsd_info = self._check_tsd(chr_seq, genome_start5, genome_end5, 
                                           genome_start3, genome_end3)
        
        if has_tsd:
            return (ltr5_start, ltr5_end, ltr3_start, ltr3_end), tsd_info, "pass_module2"
        else:
            return (ltr5_start, ltr5_end, ltr3_start, ltr3_end), None, "no_tsd_validation"
    
    ##  ========== Module 2-2 ==========

    def _enhanced_check_flanking_alignment(self, chr_seq, genome_start5, genome_end5, 
                                          genome_start3, genome_end3):
        """使用真实BLAST检测侧翼序列比对 - Module 3核心功能
        使用基因组序列和绝对坐标
        """
        flank_size = 50
        repeat_size = 10

        # Extract 60bp region (50bp flanking + 10bp repeat region)
        upstream_5_start = max(0, genome_start5 - flank_size)
        upstream_5_end = min(len(chr_seq), genome_start5 + repeat_size)
        upstream_5 = chr_seq[upstream_5_start:upstream_5_end]
        
        downstream_5_start = max(0, genome_end5 - repeat_size)
        downstream_5_end = min(len(chr_seq), genome_end5 + flank_size)
        downstream_5 = chr_seq[downstream_5_start:downstream_5_end]
        
        upstream_3_start = max(0, genome_start3 - flank_size)
        upstream_3_end = min(len(chr_seq), genome_start3 + repeat_size)
        upstream_3 = chr_seq[upstream_3_start:upstream_3_end]
        
        downstream_3_start = max(0, genome_end3 - repeat_size)
        downstream_3_end = min(len(chr_seq), genome_end3 + flank_size)
        downstream_3 = chr_seq[downstream_3_start:downstream_3_end]
        
        # Verify flanking sequence alignment with BLAST
        # Check upstream region ("a" and "c" regions)
        if self._blast_sequence_alignment(upstream_5, upstream_3, min_identity=60.0, min_length=36):
            return False  # Found flanking alignment, marked as false positive

        # Check downstream region ("b" and "d" regions)
        if self._blast_sequence_alignment(downstream_5, downstream_3, min_identity=60.0, min_length=36):
            return False  # Found flanking alignment, marked as false positive

        return True  # Flanking sequences are clean, passed check

    def _blast_sequence_alignment(self, seq1, seq2, min_identity=60.0, min_length=36):
        """Use BLAST to verify whether two sequences have significant alignment"""
        if len(seq1) < min_length or len(seq2) < min_length:
            return False
        
        query_file = None
        subject_file = None
        
        try:
            query_file = self._create_temp_fasta([seq1], "flanking_query")
            subject_file = self._create_temp_fasta([seq2], "flanking_subject")
            
            # Run BLAST alignment with relaxed parameters to detect extended alignment
            blast_output = self._run_blastn(
                query_file, subject_file, 
                evalue=10, word_size=7, max_target_seqs=10
            )
            
            alignments = self._parse_blast_output(blast_output)
            
            for alignment in alignments:
                # Check if significant alignment conditions are met
                if (alignment['identity'] >= min_identity and 
                    alignment['length'] >= min_length and
                    alignment['evalue'] <= 1.0):

                    # Calculate coverage - at least 36bp with 60% coverage
                    coverage = alignment['length'] / min(len(seq1), len(seq2))
                    if coverage >= 0.6:  # 60% coverage
                        return True
            
            return False
            
        except Exception as e:
            print(f"BLAST alignment error: {e}")
            return False
            
        finally:
            # Clean up temporary files
            for temp_file in [query_file, subject_file]:
                if temp_file and os.path.exists(temp_file):
                    os.unlink(temp_file)

    def _exhaustive_tsd_motif_search(self, chr_seq, genome_start5, genome_end5, 
                                    genome_start3, genome_end3):
        """TSD-motif joint search - using genomic coordinates"""
        # Extract 11bp region: 8bp flanking + 3bp element terminus
        left_start = max(0, genome_start5 - 8)
        left_end = min(len(chr_seq), genome_start5 + 3)
        left_11bp = chr_seq[left_start:left_end]
        
        right_start = max(0, genome_end3 - 3)
        right_end = min(len(chr_seq), genome_end3 + 8)
        right_11bp = chr_seq[right_start:right_end]
        
        # Search for motifs in priority order: TGCA first
        for motif in self.motif_priority:
            if motif == "TGCA":
                result = self._search_canonical_structure(left_11bp, right_11bp)
            else:
                start_motif = motif[:2]
                end_motif = motif[2:]
                result = self._search_motif_structure(left_11bp, right_11bp, start_motif, end_motif)
            
            if result['found']:
                result['motif'] = motif  # Record the type of motif found
                return result

        # If no known motifs are found, look for the longest matching k-mer as TSD candidates
        kmer_result = self._find_longest_kmer_tsd(left_11bp, right_11bp)
        return kmer_result
    
    def _search_canonical_structure(self, left_11bp, right_11bp):
        """Search for canonical TGCA motif with 5bp TSD"""
        return self._search_motif_structure(left_11bp, right_11bp, "TG", "CA")
    
    def _search_motif_structure(self, left_11bp, right_11bp, start_motif, end_motif):
        """Search for specific motif structure with TSD"""
        # Search for motif in left region
        start_pos = left_11bp.upper().rfind(start_motif)
        if start_pos == -1:
            return {'found': False}
        
        # Search for motif in right region
        end_pos = right_11bp.upper().find(end_motif)
        if end_pos == -1:
            return {'found': False}
        
        # Extract potential TSD (5bp preferred)
        for tsd_size in [5, 4, 6]:
            if start_pos >= tsd_size and len(right_11bp) - end_pos >= tsd_size + 2:
                left_tsd = left_11bp[start_pos - tsd_size:start_pos]
                right_tsd = right_11bp[end_pos + 2:end_pos + 2 + tsd_size]
                
                if len(left_tsd) == tsd_size and len(right_tsd) == tsd_size:
                    mismatch = self._calculate_mismatch(left_tsd, right_tsd)
                    if mismatch <= 1:  # 1bp mismatch allowed
                        return {
                            'found': True,
                            'motif': start_motif + end_motif,
                            'upstream_tsd': left_tsd,
                            'downstream_tsd': right_tsd,
                            'size': tsd_size,
                            'mismatch': mismatch
                        }
        
        return {'found': False}
    
    def _find_longest_kmer_tsd(self, left_11bp, right_11bp):
        """Search for longest matching k-mer as TSD candidates"""
        max_kmer_size = min(8, len(left_11bp), len(right_11bp))

        for kmer_size in range(max_kmer_size, 2, -1):  # Start from longest
            for i in range(len(left_11bp) - kmer_size + 1):
                left_kmer = left_11bp[i:i + kmer_size]
                for j in range(len(right_11bp) - kmer_size + 1):
                    right_kmer = right_11bp[j:j + kmer_size]
                    
                    mismatch = self._calculate_mismatch(left_kmer, right_kmer)
                    if mismatch <= 1:
                        return {
                            'found': True,
                            'motif': 'unknown',
                            'upstream_tsd': left_kmer,
                            'downstream_tsd': right_kmer,
                            'size': kmer_size,
                            'mismatch': mismatch
                        }
        
        return {'found': False}
    
    def _calculate_mismatch(self, seq1, seq2):
        """Calculate the number of mismatches between two sequences of equal length"""
        if len(seq1) != len(seq2):
            return float('inf')
        return sum(1 for a, b in zip(seq1.upper(), seq2.upper()) if a != b)
    
    def _additional_validation_checks(self, seq, ltr5_start, ltr5_end, ltr3_start, ltr3_end):
        """Additional validation checks for LTR structure"""
        reasons = []
        
        ltr5_length = ltr5_end - ltr5_start
        ltr3_length = ltr3_end - ltr3_start
        
        if ltr5_length < 100 or ltr5_length > 7000:
            reasons.append(f"5'LTR length abnormal({ltr5_length}bp)")
        
        if ltr3_length < 100 or ltr3_length > 7000:
            reasons.append(f"3'LTR length abnormal({ltr3_length}bp)")
        
        length_diff = abs(ltr5_length - ltr3_length)
        if length_diff > 0.1 * max(ltr5_length, ltr3_length):
            reasons.append(f"LTR length difference too large({length_diff}bp)")

        if reasons:
            return "truncated", "; ".join(reasons)
        else:
            return "pass", ""
    
    def check_ltr_status(self, seq, ltr5_start, ltr5_end, ltr3_start, ltr3_end, 
                        chr_seq, chr_start):
        """Enhanced LTR status check using real BLAST alignment - Module 2-2"""
        
        # Calculate genomic absolute coordinates
        genome_start5 = chr_start + ltr5_start
        genome_end5 = chr_start + ltr5_end
        genome_start3 = chr_start + ltr3_start
        genome_end3 = chr_start + ltr3_end
        
        # MODULE 3 Step 1: Critical flanking sequence alignment check
        flanking_clean = self._enhanced_check_flanking_alignment(
            chr_seq, genome_start5, genome_end5, genome_start3, genome_end3)
        
        if not flanking_clean:
            return "false", "flanking_alignment_detected", None
        
        # MODULE 3 Step 2: Exhaustive TSD-motif joint search
        tsd_motif_result = self._exhaustive_tsd_motif_search(
            chr_seq, genome_start5, genome_end5, genome_start3, genome_end3)
        
        if tsd_motif_result['found']:
            # Check if one of the 8 known motifs was found
            motif = tsd_motif_result.get('motif', 'unknown')
            if motif == 'unknown' or motif not in self.motif_priority:
                # If no typical 8 motifs were found, mark as truncated
                return "truncated", "Found TSD but missing typical motif", {
                    'upstream_tsd': tsd_motif_result['upstream_tsd'],
                    'downstream_tsd': tsd_motif_result['downstream_tsd'],
                    'size': tsd_motif_result['size'],
                    'mismatch': tsd_motif_result['mismatch'],
                    'motif': motif
                }
            
            # Found complete TSD and motif structure
            enhanced_tsd_info = {
                'upstream_tsd': tsd_motif_result['upstream_tsd'],
                'downstream_tsd': tsd_motif_result['downstream_tsd'],
                'size': tsd_motif_result['size'],
                'mismatch': tsd_motif_result['mismatch'],
                'motif': motif
            }

            # Ensure TSD size is within 4-6bp range
            if tsd_motif_result['size'] < 4 or tsd_motif_result['size'] > 6:
                return "truncated", f"TSD size abnormal({tsd_motif_result['size']}bp)", enhanced_tsd_info

            # Additional validation checks
            status, reasons = self._additional_validation_checks(
                seq, ltr5_start, ltr5_end, ltr3_start, ltr3_end)
            
            if status == "pass":
                return "pass", f"intact LTR-RT motif:{enhanced_tsd_info['motif']}", enhanced_tsd_info
            else:
                return "truncated", reasons, enhanced_tsd_info
        
        # MODULE 3 Step 3: If complete structure is not found, perform basic check
        status = "truncated"
        reasons = []

        # Check LTR length
        ltr5_length = ltr5_end - ltr5_start
        ltr3_length = ltr3_end - ltr3_start
        
        if ltr5_length < 100 or ltr5_length > 7000:
            reasons.append(f"5'LTR length abnormal({ltr5_length}bp)")
        
        if ltr3_length < 100 or ltr3_length > 7000:
            reasons.append(f"3'LTR length abnormal({ltr3_length}bp)")

        length_diff = abs(ltr5_length - ltr3_length)
        if length_diff > 0.1 * max(ltr5_length, ltr3_length):
            reasons.append(f"LTR length difference too large({length_diff}bp)")

        # Check TSD again
        has_tsd, tsd_info = self._check_tsd(chr_seq, genome_start5, genome_end5, 
                                           genome_start3, genome_end3)
        if not has_tsd:
            reasons.append("TSD not detected")
        else:
            # Even if TSD is present but complete motif is absent, still mark as truncated
            reasons.append("Missing complete TSD-motif structure")

        return status, "; ".join(reasons) if reasons else "Structure incomplete", tsd_info
    

    def process_ltrs(self, ltr_file, output_dir):
        """Process LTR sequences and output results (maintain original output format)"""
        os.makedirs(output_dir, exist_ok=True)
        
        fasta_out_path = os.path.join(output_dir, "boundary_modified.fasta")
        results_out_path = os.path.join(output_dir, "boundary_results.tsv")
        
        with open(fasta_out_path, 'w') as fasta_out, open(results_out_path, 'w') as results_out:
            results_out.write("ID\tChromosome\tOriginal_5LTR\tModified_5LTR\tOriginal_3LTR\tModified_3LTR\t"
                            "Status\tReason\tTSD_Size\tTSD_Upstream\tTSD_Downstream\tTSD_Mismatch\n")
            
            for record in SeqIO.parse(ltr_file, "fasta"):
                info = self._parse_fasta_header(record.description)
                
                if info['chr'] not in self.genome_dict:
                    print(f"Warning: Chromosome {info['chr']} not found in genome file, skipping {info['id']}")
                    continue
                
                chr_seq = self.genome_dict[info['chr']]
                
                start = info['start'] - 1
                end = info['end']
                full_seq = chr_seq[start:end]
                
                rel_5ltr_start = info['5ltr_start'] - start - 1
                rel_5ltr_end = info['5ltr_end'] - start
                rel_3ltr_start = info['3ltr_start'] - start - 1
                rel_3ltr_end = info['3ltr_end'] - start
                
                # Use enhanced boundary search (including real BLAST)
                boundary_result = self.find_boundary(
                    full_seq, rel_5ltr_start, rel_5ltr_end, 
                    rel_3ltr_start, rel_3ltr_end,
                    chr_seq, start  # pass genome sequence and start position
                )
                new_boundaries, boundary_tsd_info, module2_status = boundary_result
                
                if module2_status != "pass_module2":
                    status = "false"
                    reason = module2_status
                    tsd_info = None
                else:
                    # Use enhanced status check (including real BLAST)
                    status_result = self.check_ltr_status(
                        full_seq, *new_boundaries,
                        chr_seq, start  # pass genome sequence and start position
                    )
                    status, reason, status_tsd_info = status_result
                    
                    tsd_info = status_tsd_info if status_tsd_info else boundary_tsd_info

                # Convert back to absolute genomic coordinates
                abs_boundaries = (
                    new_boundaries[0] + start + 1,
                    new_boundaries[1] + start,
                    new_boundaries[2] + start + 1,
                    new_boundaries[3] + start
                )
                
                if tsd_info:
                    tsd_size = str(tsd_info['size'])
                    tsd_upstream = tsd_info['upstream_tsd']
                    tsd_downstream = tsd_info['downstream_tsd']
                    tsd_mismatch = str(tsd_info['mismatch'])
                else:
                    tsd_size = "NA"
                    tsd_upstream = "NA"
                    tsd_downstream = "NA"
                    tsd_mismatch = "NA"
                
                results_out.write(f"{info['id']}\t{info['chr']}\t"
                              f"{info['5ltr_start']}-{info['5ltr_end']}\t"
                              f"{abs_boundaries[0]}-{abs_boundaries[1]}\t"
                              f"{info['3ltr_start']}-{info['3ltr_end']}\t"
                              f"{abs_boundaries[2]}-{abs_boundaries[3]}\t"
                              f"{status}\t{reason}\t"
                              f"{tsd_size}\t{tsd_upstream}\t{tsd_downstream}\t{tsd_mismatch}\n")
                
                if status == "pass":
                    coords = f"{abs_boundaries[0]}:{abs_boundaries[3]}"
                    tsd_header = f"tsd_size={tsd_size} tsd_seq={tsd_upstream}|{tsd_downstream}" if tsd_info else "tsd=NA"
                    fasta_out.write(f">{info['id']} chr={info['chr']} "
                                f"coords={coords} "
                                f"5ltr={abs_boundaries[0]}:{abs_boundaries[1]} "
                                f"3ltr={abs_boundaries[2]}:{abs_boundaries[3]} "
                                f"status={status} {tsd_header}\n")
                    fasta_out.write(f"{full_seq}\n")

            print(f"Processing complete! Results saved to {output_dir}")
            print(f"Only LTR-RT sequences with status=pass were written to {fasta_out_path}")
            print(f"All sequence processing results have been saved to {results_out_path}")


if __name__ == "__main__":
    genome_file = "/proj/nobackup/hpc2nstor2024-021/zhychen/data/genome/rice/all.chrs.con"
    ltr_file = "/proj/nobackup/hpc2nstor2024-028/zhychen/bin/software/ltr_filter/module1/filtered.fasta"
    output_dir = "/proj/nobackup/hpc2nstor2024-028/zhychen/bin/software/ltr_filter/module2/"
    
    
    processor = Module2(
        genome_file, 
        blast_path="blastn",  
        makeblastdb_path="makeblastdb"
    )
    processor.process_ltrs(ltr_file, output_dir)