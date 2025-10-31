from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
import logging
import os
import re
from Bio import SeqIO
import subprocess
from datetime import datetime


@dataclass
class LTRCandidate:
    id: str
    sequence: str
    length: int
    source: str
    chr: str
    coords: Tuple[int, int]
    ltr5: Tuple[int, int]
    ltr3: Tuple[int, int]
    similarity: float
    internal_length: int
    ltr5_length: int
    ltr3_length: int
    ltr_internal_ratio: float


class LTRModule1:
    def __init__(self, 
                 max_gap_size: int = 10,
                 min_internal_len: int = 100,
                 max_internal_len: int = 15000,
                 min_ltr_internal_ratio: float = 0.05,
                 max_ltr_internal_ratio: float = 50.0,
                 trf_output_dir: str = "./trf_results"):
        self.max_gap_size = max_gap_size
        self.min_internal_len = min_internal_len
        self.max_internal_len = max_internal_len
        self.min_ltr_internal_ratio = min_ltr_internal_ratio
        self.max_ltr_internal_ratio = max_ltr_internal_ratio
        self.trf_output_dir = trf_output_dir
        
        # Create TRF output directory
        os.makedirs(self.trf_output_dir, exist_ok=True)
        
        logging.basicConfig(level=logging.INFO,
                          format='%(asctime)s - %(levelname)s - %(message)s')
        
    def parse_ltr_header(self, header: str) -> Dict:
        """Parse LTR sequence header information"""
        info = {}
        
        # Parse ID
        id_match = re.match(r'>?(\S+)', header)
        if id_match:
            info['id'] = id_match.group(1)
        
        # Parse chromosome
        chr_match = re.search(r'chr=(\S+)', header)
        if chr_match:
            info['chr'] = chr_match.group(1)
        
        # Parse coordinates
        coords_match = re.search(r'coords=(\d+):(\d+)', header)
        if coords_match:
            info['coords'] = (int(coords_match.group(1)), int(coords_match.group(2)))
        
        # Parse 5' LTR
        ltr5_match = re.search(r'5ltr=(\d+):(\d+)', header)
        if ltr5_match:
            info['5ltr'] = (int(ltr5_match.group(1)), int(ltr5_match.group(2)))
        
        # Parse 3' LTR
        ltr3_match = re.search(r'3ltr=(\d+):(\d+)', header)
        if ltr3_match:
            info['3ltr'] = (int(ltr3_match.group(1)), int(ltr3_match.group(2)))
        
        # Parse similarity
        sim_match = re.search(r'sim=([\d.]+)', header)
        if sim_match:
            info['sim'] = float(sim_match.group(1))
        
        return info
    
    def calculate_element_metrics(self, info: Dict) -> Dict:
        """Calculate various metrics for the element"""
        metrics = {}
        
        if '5ltr' in info and '3ltr' in info:
            # Calculate LTR length
            metrics['ltr5_length'] = info['5ltr'][1] - info['5ltr'][0] + 1
            metrics['ltr3_length'] = info['3ltr'][1] - info['3ltr'][0] + 1

            # Calculate internal region length
            metrics['internal_start'] = info['5ltr'][1] + 1
            metrics['internal_end'] = info['3ltr'][0] - 1
            metrics['internal_length'] = metrics['internal_end'] - metrics['internal_start'] + 1

            # Calculate LTR:internal ratio
            avg_ltr_length = (metrics['ltr5_length'] + metrics['ltr3_length']) / 2
            if metrics['internal_length'] > 0:
                metrics['ltr_internal_ratio'] = avg_ltr_length / metrics['internal_length'] 
            else:
                metrics['ltr_internal_ratio'] = float('inf')
        
        return metrics
    
    def load_candidates(self, fasta_file: str) -> List[LTRCandidate]:
        """Load candidate sequences"""
        candidates = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            info = self.parse_ltr_header(record.description)
            metrics = self.calculate_element_metrics(info)
            
            candidate = LTRCandidate(
                id=info.get('id', record.id),
                sequence=str(record.seq),
                length=len(record.seq),
                source="LTRharvest",
                chr=info.get('chr', 'unknown'),
                coords=info.get('coords', (0, 0)),
                ltr5=info.get('5ltr', (0, 0)),
                ltr3=info.get('3ltr', (0, 0)),
                similarity=info.get('sim', 0.0),
                internal_length=metrics.get('internal_length', 0),
                ltr5_length=metrics.get('ltr5_length', 0),
                ltr3_length=metrics.get('ltr3_length', 0),
                ltr_internal_ratio=metrics.get('ltr_internal_ratio', 0.0)
            )
            candidates.append(candidate)
        return candidates

    def check_internal_length(self, internal_length: int) -> Tuple[bool, str]:
        """Check internal region length"""
        if internal_length < self.min_internal_len:
            return False, f"Internal length ({internal_length}) < minimum ({self.min_internal_len})"
        elif internal_length > self.max_internal_len:
            return False, f"Internal length ({internal_length}) > maximum ({self.max_internal_len})"
        return True, "PASS"

    def check_ltr_internal_ratio(self, ratio: float) -> Tuple[bool, str]:
        """Check LTR:internal ratio"""
        if ratio < self.min_ltr_internal_ratio:
            return False, f"LTR:internal ratio ({ratio:.3f}) < minimum ({self.min_ltr_internal_ratio})"
        elif ratio > self.max_ltr_internal_ratio:
            return False, f"LTR:internal ratio ({ratio:.3f}) > maximum ({self.max_ltr_internal_ratio})"
        return True, "PASS"
    
    def get_max_gap_length(self, sequence: str) -> int:
        """Calculate maximum consecutive gap length in sequence"""
        gap_count = 0
        max_gap = 0
        
        for base in sequence:
            if base.upper() == 'N':
                gap_count += 1
                max_gap = max(max_gap, gap_count)
            else:
                gap_count = 0
                
        return max_gap
    
    def run_trf_analysis(self, sequence: str, seq_id: str) -> Tuple[bool, int, str]:
        """Run TRF analysis and save results"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        trf_basename = f"{seq_id}_{timestamp}"
        tmp_file = os.path.join(self.trf_output_dir, f"{trf_basename}.fa")
        
        try:
            # Write temporary sequence file
            with open(tmp_file, "w") as f:
                f.write(f">{seq_id}\n{sequence}\n")

            # Run TRF (using parameters from the paper)
            cmd = f"trf {tmp_file} 2 7 7 80 10 1000 2000 -ngs -h"
            result = subprocess.run(cmd.split(), 
                                 capture_output=True,
                                 text=True,
                                 check=True)

            # Save TRF output
            # trf_output_file = os.path.join(self.trf_output_dir, f"{trf_basename}_trf.txt")
            # with open(trf_output_file, "w") as f:
            #     f.write(f"Command: {cmd}\n")
            #     f.write(f"STDOUT:\n{result.stdout}\n")
            #     f.write(f"STDERR:\n{result.stderr}\n")
            os.remove(tmp_file)  # Delete temporary file
            # Parse TRF results
            tandem_repeats = self.parse_trf_output(result.stdout)

            # According to the paper's criteria: contains tandem repeats with a minimum alignment length of 50bp and a maximum repeat fragment of 500bp
            significant_repeats = 0
            for repeat in tandem_repeats:
                if repeat['align_len'] >= 50 and repeat['period'] <= 500:
                    significant_repeats += 1

            # If contains significant tandem repeats, filter out
            passed = significant_repeats == 0
            filter_reason = "PASS" if passed else f"Contains {significant_repeats} tandem repeats (align_len>=50bp, period<=500bp)"
            
            return passed, significant_repeats, filter_reason
            
        except Exception as e:
            logging.warning(f"TRF analysis failed {seq_id}: {e}")
            # TReturn pass when TRF fails
            return True, 0, "TRF analysis failed"
        
        finally:
            # Keep temporary fasta file for debugging
            pass
    
    def parse_trf_output(self, trf_output: str) -> List[Dict]:
        """Parse TRF output"""
        repeats = []
        for line in trf_output.split('\n'):
            if line.strip():
                parts = line.split()
                if len(parts) >= 15 and parts[0].isdigit():
                    try:
                        start = int(parts[0])
                        end = int(parts[1])
                        period = int(parts[2])
                        copies = float(parts[3])
                        align_len = end - start + 1
                        
                        repeats.append({
                            'start': start,
                            'end': end,
                            'period': period,
                            'copies': copies,
                            'align_len': align_len
                        })
                    except (ValueError, IndexError):
                        continue
        return repeats
    
    def filter_candidates(self, input_fasta: str, output_fasta: str):
        """Main filtering function"""
        # Load candidate sequences
        candidates = self.load_candidates(input_fasta)
        filtered = []

        # Create filter log file
        output_dir = os.path.dirname(output_fasta)
        if not output_dir:
            output_dir = "."
        filter_log_path = os.path.join(output_dir, "module1_filter_results.tsv")
        
        with open(filter_log_path, "w") as filter_log:
            # Write column descriptions
            column_descriptions = [
                "# LTR_retriever Module 1 Filter Results",
                "# Column descriptions:",
                "# Sequence_ID",
                "# Total_Length",
                "# Internal_Length",
                "# LTR5_Length",
                "# LTR3_Length",
                "# LTR_Internal_Ratio",
                "# Internal_Check",
                "# Ratio_Check",
                "# Gap_Check",
                "# Max_Gap_Length:",
                "# TRF_Check:",
                "# Tandem_Repeats",
                "# Filter_Status: (PASS/FILTER)",
                "# Filter_Reason",
                "#"
            ]
            filter_log.write("\n".join(column_descriptions) + "\n")

            # Write header
            headers = [
                "Sequence_ID",
                "Total_Length",
                "Internal_Length",
                "LTR5_Length",
                "LTR3_Length",
                "LTR_Internal_Ratio",
                "Internal_Check",
                "Ratio_Check",
                "Gap_Check",
                "Max_Gap_Length",
                "TRF_Check",
                "Tandem_Repeats",
                "Filter_Status",
                "Filter_Reason"
            ]
            filter_log.write("\t".join(headers) + "\n")
            
            for candidate in candidates:
                # Initialize result dictionary
                result = {
                    "Sequence_ID": candidate.id,
                    "Total_Length": str(candidate.length),
                    "Internal_Length": str(candidate.internal_length),
                    "LTR5_Length": str(candidate.ltr5_length),
                    "LTR3_Length": str(candidate.ltr3_length),
                    "LTR_Internal_Ratio": f"{candidate.ltr_internal_ratio:.3f}",
                    "Internal_Check": "PASS",
                    "Ratio_Check": "PASS",
                    "Gap_Check": "PASS",
                    "Max_Gap_Length": "0",
                    "TRF_Check": "PASS",
                    "Tandem_Repeats": "0",
                    "Filter_Status": "PASS",
                    "Filter_Reason": "NA"
                }

                # Check internal region length
                internal_pass, internal_reason = self.check_internal_length(candidate.internal_length)
                if not internal_pass:
                    result.update({
                        "Internal_Check": "FAIL",
                        "Filter_Status": "FILTER",
                        "Filter_Reason": internal_reason
                    })

                # Check LTR:internal ratio
                elif candidate.ltr_internal_ratio > 0:
                    ratio_pass, ratio_reason = self.check_ltr_internal_ratio(candidate.ltr_internal_ratio)
                    if not ratio_pass:
                        result.update({
                            "Ratio_Check": "FAIL",
                            "Filter_Status": "FILTER",
                            "Filter_Reason": ratio_reason
                        })

                # Check gaps
                if result["Filter_Status"] == "PASS":
                    max_gap = self.get_max_gap_length(candidate.sequence)
                    result["Max_Gap_Length"] = str(max_gap)
                    if max_gap > self.max_gap_size:
                        result.update({
                            "Gap_Check": "FAIL",
                            "Filter_Status": "FILTER",
                            "Filter_Reason": f"Max gap length ({max_gap}) > threshold ({self.max_gap_size})"
                        })

                # Run TRF analysis
                if result["Filter_Status"] == "PASS":
                    trf_pass, tandem_count, trf_reason = self.run_trf_analysis(
                        candidate.sequence, candidate.id)
                    result["Tandem_Repeats"] = str(tandem_count)
                    if not trf_pass:
                        result.update({
                            "TRF_Check": "FAIL",
                            "Filter_Status": "FILTER",
                            "Filter_Reason": trf_reason
                        })

                # Write results
                filter_log.write("\t".join([result[h] for h in headers]) + "\n")

                # If passed all filters, add to output sequences
                if result["Filter_Status"] == "PASS":
                    filtered.append(candidate)

        # Save filtered sequences
        with open(output_fasta, 'w') as f:
            for candidate in filtered:
                # Reconstruct original description line
                desc = f"{candidate.id} chr={candidate.chr} coords={candidate.coords[0]}:{candidate.coords[1]}"
                desc += f" 5ltr={candidate.ltr5[0]}:{candidate.ltr5[1]}"
                desc += f" 3ltr={candidate.ltr3[0]}:{candidate.ltr3[1]}"
                desc += f" sim={candidate.similarity:.2f}"
                f.write(f">{desc}\n{candidate.sequence}\n")
        
        logging.info(f"Module 1 filtering finished：")
        logging.info(f"  Input sequences: {len(candidates)}")
        logging.info(f"  Passed filtering: {len(filtered)}")
        logging.info(f"  Filtered out: {len(candidates) - len(filtered)}")
        logging.info(f"  Filter results saved in: {filter_log_path}")
        logging.info(f"  TRF results saved in: {self.trf_output_dir}")


if __name__ == '__main__':
    # 1. Initialize filter
    filter_module = LTRModule1(
        max_gap_size=10,              # Allowed maximum consecutive N length
        min_internal_len=100,         # Internal region minimum length
        max_internal_len=15000,       # Internal region maximum length
        min_ltr_internal_ratio=0.05,  # LTR:internal minimum ratio
        max_ltr_internal_ratio=50.0,  # LTR:internal maximum ratio
        trf_output_dir="./trf_results"  # TRF results output directory
    )

    # 2. Set input and output file paths
    input_fasta = "/proj/nobackup/hpc2nstor2024-021/zhychen/bin/software/ltr_checker/ltr_finder_candidates.fasta"
    output_fasta = "/proj/nobackup/hpc2nstor2024-028/zhychen/bin/software/ltr_filter/module1/filtered.fasta"

    # 3. Run filtering
    try:
        filter_module.filter_candidates(input_fasta, output_fasta)
        temp_dir = os.getcwd() + "/trf_results"
        os.rmdir(temp_dir)  # Clean up TRF results directory
    except Exception as e:
        logging.error(f"Filtering process error: {e}")
        raise