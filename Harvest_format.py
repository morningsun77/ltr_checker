import logging
from dataclasses import dataclass
from typing import List, Dict
from Bio import SeqIO

@dataclass
class LTRElement:
    """Store information for a LTR-RT element"""
    element_id: str
    chr_id: str
    element_start: int
    element_end: int
    element_length: int
    ltr5_start: int
    ltr5_end: int
    ltr5_length: int
    ltr3_start: int
    ltr3_end: int
    ltr3_length: int
    similarity: float
    sequence: str = ""  # Store the sequence

class LTRHarvestProcessor:
    def __init__(self):
        """Initialize the processor"""
        self.logger = logging.getLogger(__name__)
        
    def create_chromosome_mapping(self, genome_file: str) -> Dict[str, str]:
        """
        Create a mapping from index to chromosome ID
        Note: The index starts from 0 and the maximum value is the number of chromosomes - 1
        """
        chr_index_mapping = {}
        try:
            # Only read chromosome IDs without loading sequences to improve efficiency
            chr_ids = [record.id for record in SeqIO.parse(genome_file, "fasta")]
            for i, chr_id in enumerate(chr_ids):
                chr_index_mapping[str(i)] = chr_id
            self.logger.info(f"Loaded {len(chr_index_mapping)} chromosome IDs from the genome file")
        except Exception as e:
            self.logger.error(f"Error reading genome file: {e}")
        return chr_index_mapping
        
    def parse_harvest_output(self, harvest_file: str, chr_index_mapping: Dict[str, str]) -> List[LTRElement]:
        """
        Parse the LTRharvest output file
        
        Args:
            harvest_file: LTRharvest output file path
            chr_index_mapping: Mapping from index to chromosome ID

        Returns:
            List of LTR-RT elements
        """
        elements = []
        element_count = 0
        
        try:
            with open(harvest_file) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                        
                    fields = line.strip().split()
                    # Ensure at least 11 columns of data
                    if len(fields) >= 11:
                        element_count += 1

                        # Handle two cases for chromosome ID
                        chr_id = fields[10]
                        if chr_id.isdigit():
                            # Case 1: Numeric index, map to chromosome ID
                            if chr_id in chr_index_mapping:
                                chr_id = chr_index_mapping[chr_id]
                            else:
                                chr_id = f"Chromosome_{chr_id}"
                        # Case 2: Already a chromosome ID, use it directly
                        
                        element = LTRElement(
                            element_id=f"LTR_{element_count}",
                            chr_id=chr_id,
                            element_start=int(fields[0]),
                            element_end=int(fields[1]),
                            element_length=int(fields[2]),
                            ltr5_start=int(fields[3]),
                            ltr5_end=int(fields[4]),
                            ltr5_length=int(fields[5]),
                            ltr3_start=int(fields[6]),
                            ltr3_end=int(fields[7]),
                            ltr3_length=int(fields[8]),
                            similarity=float(fields[9])
                        )
                        elements.append(element)
        except Exception as e:
            self.logger.error(f"Error parsing input file: {e}")

        self.logger.info(f"Parsed {len(elements)} LTR-RT elements")
        return elements

    def extract_sequences(self, elements: List[LTRElement], genome_file: str) -> List[LTRElement]:
        """Extract sequences from the genome file, optimizing memory usage"""
        elements_with_seq = []
        failed_elements = 0

        # Group elements by chromosome to reduce file read times
        chr_elements = {}
        for element in elements:
            if element.chr_id not in chr_elements:
                chr_elements[element.chr_id] = []
            chr_elements[element.chr_id].append(element)

        # Process one chromosome at a time to reduce memory usage
        for record in SeqIO.parse(genome_file, "fasta"):
            chr_id = record.id
            if chr_id in chr_elements:
                self.logger.info(f"Processing {len(chr_elements[chr_id])} elements on chromosome {chr_id}")
                chrom_seq = str(record.seq)
                chrom_len = len(chrom_seq)
                
                for element in chr_elements[chr_id]:
                    try:
                        # Ensure index is within valid range
                        if 0 <= element.element_start < chrom_len and element.element_end <= chrom_len:
                            element.sequence = chrom_seq[element.element_start:element.element_end]
                            elements_with_seq.append(element)
                        else:
                            self.logger.warning(f"Element {element.element_id} coordinates out of range: {element.element_start}-{element.element_end}, chromosome length: {chrom_len}")
                            failed_elements += 1
                    except Exception as e:
                        self.logger.error(f"Error extracting sequence for element {element.element_id}: {e}")
                        failed_elements += 1

                # Release memory
                del chr_elements[chr_id]

        # Process elements with missing chromosomes
        for chr_id, elements_list in chr_elements.items():
            self.logger.warning(f"Chromosome {chr_id} not found, skipping {len(elements_list)} elements")
            failed_elements += len(elements_list)

        self.logger.info(f"Successfully extracted sequences for {len(elements_with_seq)} elements, failed for {failed_elements}")
        return elements_with_seq

    def write_candidates_fasta(self, elements: List[LTRElement], output_file: str):
        """Generate FASTA file for candidate sequences"""
        try:
            with open(output_file, 'w') as f:
                for element in elements:
                    if not element.sequence:
                        continue

                    # Construct header with all necessary information
                    header = (f">{element.element_id} "
                             f"chr={element.chr_id} "
                             f"coords={element.element_start}:{element.element_end} "
                             f"5ltr={element.ltr5_start}:{element.ltr5_end} "
                             f"3ltr={element.ltr3_start}:{element.ltr3_end} "
                             f"sim={element.similarity:.2f}")

                    # Write sequence
                    f.write(f"{header}\n{element.sequence}\n")

            self.logger.info(f"Successfully wrote {len(elements)} candidate sequences to {output_file}")
        except Exception as e:
            self.logger.error(f"Error writing output file: {e}")

    def process_ltrharvest_results(self, harvest_file: str, genome_file: str, output_file: str, min_similarity: float = 0.0):
        """
        Complete workflow for processing LTRharvest results

        Args:
            harvest_file: LTRharvest output file path
            genome_file: Genome FASTA file path
            output_file: Output FASTA file path
            min_similarity: Minimum similarity filtering threshold (0-100)
        """
        # Step 1: Create chromosome mapping
        chr_index_mapping = self.create_chromosome_mapping(genome_file)

        # Step 2: Parse LTRharvest output
        elements = self.parse_harvest_output(harvest_file, chr_index_mapping)

        # Step 3: Filter low similarity elements (if needed)
        if min_similarity > 0:
            filtered_elements = [e for e in elements if e.similarity >= min_similarity]
            # self.logger.info(f"通过相似度过滤（阈值{min_similarity}），保留{len(filtered_elements)}/{len(elements)}个元件")
            self.logger.info(f"Similarity filtering (threshold {min_similarity}): retained {len(filtered_elements)}/{len(elements)} elements")
            elements = filtered_elements

        # Step 4: Extract sequences
        elements_with_seq = self.extract_sequences(elements, genome_file)

        # Step 5: Write to FASTA file
        self.write_candidates_fasta(elements_with_seq, output_file)
        
        return elements_with_seq
    
if __name__ == "__main__":
    processor = LTRHarvestProcessor()
    harvest = '/proj/nobackup/hpc2nstor2024-021/zhychen/rescarch/checker_res/Triticum_aestivum/harvest/genome.fa.harvest.scn'
    genome = '/proj/nobackup/hpc2nstor2024-021/zhychen/data/genome/Triticum_aestivum/Triticum_aestivum.genome.fa'
    output = '/proj/nobackup/hpc2nstor2024-021/zhychen/rescarch/checker_res/Triticum_aestivum/harvest/ltr_candidates.fasta'
    similarity = 85
        
    processor.process_ltrharvest_results(harvest, genome, output, similarity)