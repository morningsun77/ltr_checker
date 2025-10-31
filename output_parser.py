import os 
from Harvest_format import *
import re

def parse_ltr_finder_output(file_path):
    """Parse LTR_FINDER output file"""
    results = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for line in lines[1:]:  # Skip header
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 15:
                    # Extract basic LTR information
                    seq_id = parts[0]
                    location = parts[1]  # Format: start-end
                    ltr_len = parts[2]   # Format: left_len/right_len
                    element_len = parts[3]
                    similarity = parts[14] if len(parts) > 14 else "---"
                    
                    # Parse location to get retrotransposon positions
                    if '-' in location:
                        positions = location.split('-')
                        ret_start = positions[0]
                        ret_end = positions[1]
                    else:
                        ret_start = "---"
                        ret_end = "---"
                    
                    # Parse LTR lengths
                    if ',' in ltr_len:
                        ltr_lengths = ltr_len.split(',')
                        left_ltr_len = ltr_lengths[0]
                        right_ltr_len = ltr_lengths[1]
                    else:
                        left_ltr_len = ltr_len
                        right_ltr_len = ltr_len
                    
                    # Calculate LTR positions based on the standard LTR-RT structure
                    # Left LTR is at the beginning of the element
                    if ret_start != "---" and left_ltr_len != "---":
                        left_ltr_start = ret_start
                        left_ltr_end = str(int(ret_start) + int(left_ltr_len) - 1)
                    else:
                        left_ltr_start = "---"
                        left_ltr_end = "---"
                    
                    # Right LTR is at the end of the element
                    if ret_end != "---" and right_ltr_len != "---":
                        right_ltr_end = ret_end
                        right_ltr_start = str(int(ret_end) - int(right_ltr_len) + 1)
                    else:
                        right_ltr_start = "---"
                        right_ltr_end = "---"
                    
                    results.append({
                        'seq_id': seq_id,
                        'ret_start': ret_start,
                        'ret_end': ret_end,
                        'ret_len': element_len,
                        'left_ltr_start': left_ltr_start,
                        'left_ltr_end': left_ltr_end,
                        'left_ltr_len': left_ltr_len,
                        'right_ltr_start': right_ltr_start,
                        'right_ltr_end': right_ltr_end,
                        'right_ltr_len': right_ltr_len,
                        'similarity': similarity,
                        'raw_line': line.strip()
                    })
    return results

def parse_ltr_harvest_output(file_path):
    """Parse LTR_HARVEST output file"""
    results = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) >= 11:
                results.append({
                    'seq_id': parts[10],
                    'ret_start': parts[0],
                    'ret_end': parts[1],
                    'ret_len': parts[2],
                    'left_ltr_start': parts[3],
                    'left_ltr_end': parts[4],
                    'left_ltr_len': parts[5],
                    'right_ltr_start': parts[6],
                    'right_ltr_end': parts[7],
                    'right_ltr_len': parts[8],
                    'similarity': parts[9],
                    'raw_line': line.strip()
                })
    return results

def parse_ltr_detector_output(file_path):
    """Parse LtrDetector output file"""
    results = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        # Skip the two header lines
        for line in lines[2:]:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 7:
                    results.append({
                        'seq_id': parts[0],
                        'ret_start': parts[1],
                        'ret_end': parts[2],
                        'left_ltr_start': parts[3],
                        'left_ltr_end': parts[4],
                        'right_ltr_start': parts[5],
                        'right_ltr_end': parts[6],
                        'similarity': parts[7] if len(parts) > 7 else "---",
                        'raw_line': line.strip()
                    })
    return results

def convert_to_ltr_harvest_format(data, source):
    """Convert data from any source to LTR_HARVEST format"""
    lines = []
    
    # Add LTR_HARVEST headers
    lines.append("# Unified output in LTR_HARVEST format")
    lines.append("# predictions are reported in the following way")
    lines.append("# s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-id")
    lines.append("")
    
    for item in data:
        if source == "ltr_harvest":
            # Already in correct format
            lines.append(item['raw_line'])
        elif source == "ltr_finder":
            # Convert from LTR_FINDER format using parsed actual values
            ret_start = item['ret_start']
            ret_end = item['ret_end']
            ret_len = item['ret_len']
            
            # Use actual LTR positions and lengths from parsed data
            left_ltr_start = item['left_ltr_start']
            left_ltr_end = item['left_ltr_end']
            left_ltr_len = item['left_ltr_len']
            
            right_ltr_start = item['right_ltr_start']
            right_ltr_end = item['right_ltr_end']
            right_ltr_len = item['right_ltr_len']
            
            similarity = item['similarity']
            seq_id = item['seq_id']
            
            line = f"{ret_start} {ret_end} {ret_len} {left_ltr_start} {left_ltr_end} {left_ltr_len} {right_ltr_start} {right_ltr_end} {right_ltr_len} {similarity} {seq_id}"
            lines.append(line)
            
        elif source == "ltr_detector":
            # Convert from LtrDetector format
            ret_start = item['ret_start']
            ret_end = item['ret_end']
            ret_len = str(int(ret_end) - int(ret_start) + 1) if ret_start != "---" and ret_end != "---" else "---"
            
            left_ltr_start = item['left_ltr_start']
            left_ltr_end = item['left_ltr_end']
            left_ltr_len = str(int(left_ltr_end) - int(left_ltr_start) + 1) if left_ltr_start != "---" and left_ltr_end != "---" else "---"
            
            right_ltr_start = item['right_ltr_start']
            right_ltr_end = item['right_ltr_end']
            right_ltr_len = str(int(right_ltr_end) - int(right_ltr_start) + 1) if right_ltr_start != "---" and right_ltr_end != "---" else "---"
            
            similarity = item['similarity']
            seq_id = item['seq_id']
            
            line = f"{ret_start} {ret_end} {ret_len} {left_ltr_start} {left_ltr_end} {left_ltr_len} {right_ltr_start} {right_ltr_end} {right_ltr_len} {similarity} {seq_id}"
            lines.append(line)
    
    return lines

def convert_to_ltr_finder_format(data, source):
    """Convert data from any source to LTR_FINDER format"""
    lines = []
    
    # Add LTR_FINDER header
    lines.append('SeqID\tLocation\tLTR len\tInserted element len\tTSR\tPBS\tPPT\tRT\tIN (core)\tIN (c-term)\tRH\tStrand\tScore\tSharpness\tSimilarity')
    
    for item in data:
        if source == "ltr_finder":
            # Already in correct format
            lines.append(item['raw_line'])
        else:
            # Convert from other formats
            seq_id = item['seq_id']
            
            # Create location string
            if 'ret_start' in item and 'ret_end' in item:
                location = f"{item['ret_start']}-{item['ret_end']}"
            else:
                location = "---"
            
            # Get LTR length
            if source == "ltr_harvest":
                ltr_len = f"{item['left_ltr_len']}/{item['right_ltr_len']}"
            elif source == "ltr_detector":
                left_len = str(int(item['left_ltr_end']) - int(item['left_ltr_start']) + 1) if item.get('left_ltr_start') != "---" and item.get('left_ltr_end') != "---" else "---"
                right_len = str(int(item['right_ltr_end']) - int(item['right_ltr_start']) + 1) if item.get('right_ltr_start') != "---" and item.get('right_ltr_end') != "---" else "---"
                ltr_len = f"{left_len}/{right_len}"
            else:
                ltr_len = "---"
            
            # Get inserted element length
            if 'ret_len' in item:
                element_len = item['ret_len']
            elif 'ret_start' in item and 'ret_end' in item and item['ret_start'] != "---" and item['ret_end'] != "---":
                element_len = str(int(item['ret_end']) - int(item['ret_start']) + 1)
            else:
                element_len = "---"
            
            # Get similarity
            similarity = item.get('similarity', '---')
            
            # Fill in placeholders for other fields
            tsr = "---"
            pbs = "---"
            ppt = "---"
            rt = "---"
            in_core = "---"
            in_cterm = "---"
            rh = "---"
            strand = "+"
            score = "---"
            sharpness = "---"
            
            line = f"{seq_id}\t{location}\t{ltr_len}\t{element_len}\t{tsr}\t{pbs}\t{ppt}\t{rt}\t{in_core}\t{in_cterm}\t{rh}\t{strand}\t{score}\t{sharpness}\t{similarity}"
            lines.append(line)
    
    return lines

def convert_to_ltr_detector_format(data, source):
    """Convert data from any source to LtrDetector format"""
    lines = []
    
    # Add LtrDetector headers (two lines)
    lines.append("SeqID\tRetrotransposon\tLeft_LTR\tRight_LTR\t\t\tLeft_TSD\tRight_TSD\tPolypurine Tract\t\t\tTG\tCA")
    lines.append("Start\tEnd\tStart\tEnd\tStart\tEnd\tID\tStart\tEnd\tStart\tEnd\tStart\tEnd\tStrand\tPurine%\tStart\tEnd")
    
    for item in data:
        if source == "ltr_detector":
            # Already in correct format
            lines.append(item['raw_line'])
        else:
            # Convert from other formats
            seq_id = item['seq_id']
            
            # Get retrotransposon positions
            ret_start = item.get('ret_start', '---')
            ret_end = item.get('ret_end', '---')
            
            # Get LTR positions
            if source == "ltr_harvest":
                left_ltr_start = item.get('left_ltr_start', '---')
                left_ltr_end = item.get('left_ltr_end', '---')
                right_ltr_start = item.get('right_ltr_start', '---')
                right_ltr_end = item.get('right_ltr_end', '---')
            elif source == "ltr_finder":
                # Estimate LTR positions for LTR_FINDER
                if ret_start != "---" and ret_end != "---":
                    left_ltr_start = ret_start
                    left_ltr_end = str(int(ret_start) + 500)  # Estimate
                    right_ltr_start = str(int(ret_end) - 500)  # Estimate
                    right_ltr_end = ret_end
                else:
                    left_ltr_start = left_ltr_end = right_ltr_start = right_ltr_end = "---"
            else:
                left_ltr_start = left_ltr_end = right_ltr_start = right_ltr_end = "---"
            
            # Get similarity/ID
            similarity = item.get('similarity', '---')
            
            # Placeholder values for other fields
            left_tsd_start = left_tsd_end = "---"
            right_tsd_start = right_tsd_end = "---"
            ppt_start = ppt_end = "---"
            strand = "+"
            purine_perc = "---"
            tg_start = tg_end = "---"
            
            line = f"{seq_id}\t{ret_start}\t{ret_end}\t{left_ltr_start}\t{left_ltr_end}\t{right_ltr_start}\t{right_ltr_end}\t{similarity}\t{left_tsd_start}\t{left_tsd_end}\t{right_tsd_start}\t{right_tsd_end}\t{ppt_start}\t{ppt_end}\t{strand}\t{purine_perc}\t{tg_start}\t{tg_end}"
            lines.append(line)
    
    return lines

def unify_outputs(outputDir, output_format):
    """
    Unify outputs from all three software tools into a single format
    
    Parameters:
        outputDir: output directory containing the result files
        output_format: target format (ltr_harvest, ltr_finder, or ltrdetector)
    """
    print("\n" + "="*50)
    print(f"Unifying outputs to {output_format.upper()} format...")
    print("="*50)
    
    all_data = []
    
    # Parse LTR_FINDER output if it exists
    if os.path.exists(outputDir + '/ltr_finder.out'):
        print("Parsing LTR_FINDER output...")
        finder_data = parse_ltr_finder_output(outputDir + '/ltr_finder.out')
        for item in finder_data:
            item['source'] = 'ltr_finder'
        all_data.extend(finder_data)
    
    # Parse LTR_HARVEST output if it exists
    if os.path.exists(outputDir + '/ltr_harvest.out'):
        print("Parsing LTR_HARVEST output...")
        harvest_data = parse_ltr_harvest_output(outputDir + '/ltr_harvest.out')
        for item in harvest_data:
            item['source'] = 'ltr_harvest'
        all_data.extend(harvest_data)
    
    # Parse LtrDetector output if it exists
    if os.path.exists(outputDir + '/ltr_detector.out'):
        print("Parsing LtrDetector output...")
        detector_data = parse_ltr_detector_output(outputDir + '/ltr_detector.out')
        for item in detector_data:
            item['source'] = 'ltr_detector'
        all_data.extend(detector_data)
    
    # Group data by source
    data_by_source = {}
    for item in all_data:
        source = item['source']
        if source not in data_by_source:
            data_by_source[source] = []
        data_by_source[source].append(item)
    
    # Convert all data to the target format
    unified_lines = []
    
    if output_format == "ltr_harvest":
        # Add headers for LTR_HARVEST format
        unified_lines.append("# Unified output from multiple LTR annotation tools")
        unified_lines.append("# Format: LTR_HARVEST")
        unified_lines.append("# predictions are reported in the following way")
        unified_lines.append("# s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-id")
        unified_lines.append("")
        
        for source, data in data_by_source.items():
            unified_lines.append(f"# Results from {source.upper()}")
            converted = convert_to_ltr_harvest_format(data, source)
            unified_lines.extend(converted[4:])  # Skip headers from conversion
            unified_lines.append("")
            
    elif output_format == "ltr_finder":
        # Add header for LTR_FINDER format
        unified_lines.append('SeqID\tLocation\tLTR len\tInserted element len\tTSR\tPBS\tPPT\tRT\tIN (core)\tIN (c-term)\tRH\tStrand\tScore\tSharpness\tSimilarity\tSource')
        
        for source, data in data_by_source.items():
            converted = convert_to_ltr_finder_format(data, source)
            for line in converted[1:]:  # Skip header
                unified_lines.append(line + f"\t{source.upper()}")
                
    elif output_format == "ltrdetector":
        # Add headers for LtrDetector format
        unified_lines.append("SeqID\tRetrotransposon\tLeft_LTR\tRight_LTR\t\t\tLeft_TSD\tRight_TSD\tPolypurine Tract\t\t\tTG\tCA\tSource")
        unified_lines.append("Start\tEnd\tStart\tEnd\tStart\tEnd\tID\tStart\tEnd\tStart\tEnd\tStart\tEnd\tStrand\tPurine%\tStart\tEnd\t")
        
        for source, data in data_by_source.items():
            converted = convert_to_ltr_detector_format(data, source)
            for line in converted[2:]:  # Skip headers
                unified_lines.append(line + f"\t{source.upper()}")
    
    # Write unified output
    output_file = outputDir + '/unified_output.txt'
    with open(output_file, 'w') as f:
        for line in unified_lines:
            f.write(line + '\n')

def generate_ltr_fasta(outputDir, genome_file):
    """
    Generate LTR-RT FASTA file from any output format using Harvest.py functions
    
    Parameters:
        outputDir: output directory containing the result files
        genome_file: original genome file path
    """
    print("\n" + "="*50)
    print("Generating LTR-RT FASTA file...")
    print("="*50)
    
    # Set up logging
    logging.basicConfig(level=logging.INFO)
    
    # Initialize the processor
    processor = LTRHarvestProcessor()
    
    # Temporary file for unified harvest format output
    temp_harvest_file = outputDir + '/temp_unified_harvest.txt'
    
    # Collect all results and convert to harvest format
    all_elements = []
    element_count = 0
    
    # Parse and convert each output file to harvest format
    if os.path.exists(outputDir + '/ltr_finder.out'):
        print("Processing LTR_FINDER output...")
        finder_data = parse_ltr_finder_output(outputDir + '/ltr_finder.out')
        for item in finder_data:
            item['source'] = 'ltr_finder'
        # Convert to harvest format
        converted = convert_to_ltr_harvest_format(finder_data, 'ltr_finder')
        all_elements.extend(converted[4:])  # Skip headers
    
    if os.path.exists(outputDir + '/ltr_harvest.out'):
        print("Processing LTR_HARVEST output...")
        harvest_data = parse_ltr_harvest_output(outputDir + '/ltr_harvest.out')
        for item in harvest_data:
            all_elements.append(item['raw_line'])
    
    if os.path.exists(outputDir + '/ltr_detector.out'):
        print("Processing LtrDetector output...")
        detector_data = parse_ltr_detector_output(outputDir + '/ltr_detector.out')
        for item in detector_data:
            item['source'] = 'ltr_detector'
        # Convert to harvest format
        converted = convert_to_ltr_harvest_format(detector_data, 'ltr_detector')
        all_elements.extend(converted[4:])  # Skip headers
    
    # Also check for unified output file
    if os.path.exists(outputDir + '/unified_output.txt'):
        print("Processing unified output...")
        with open(outputDir + '/unified_output.txt', 'r') as f:
            lines = f.readlines()
            for line in lines:
                # Skip comment and header lines
                if not line.startswith('#') and line.strip() and not 'predictions are reported' in line and not 's(ret)' in line:
                    # Check if it's already in harvest format
                    parts = line.strip().split()
                    if len(parts) >= 11 and parts[0].isdigit():
                        all_elements.append(line.strip())
    
    # Write temporary harvest format file
    with open(temp_harvest_file, 'w') as f:
        f.write("# Temporary unified harvest format for FASTA generation\n")
        f.write("# s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-id\n")
        for element in all_elements:
            if element and not element.startswith('#'):
                f.write(element + '\n')
    
    # Generate FASTA file using Harvest.py processor
    output_fasta = outputDir + '/ltr_candidates.fasta'
    
    try:
        # Process the harvest format file to generate FASTA
        elements_with_seq = processor.process_ltrharvest_results(
            temp_harvest_file, 
            genome_file, 
            output_fasta,
            min_similarity=0.0  # No filtering, keep all elements
        )
        
        print(f"\nSuccessfully generated LTR-RT FASTA file: ltr_candidates.fasta")
        print(f"Total LTR-RT sequences: {len(elements_with_seq)}")
        
    except Exception as e:
        print(f"Error generating FASTA file: {e}")
        print("Attempting alternative method...")
        
        # Alternative method: directly parse and extract sequences
        elements = []
        with open(temp_harvest_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                fields = line.strip().split()
                if len(fields) >= 11:
                    element_count += 1
                    chr_id = fields[10]
                    
                    # Handle chromosome ID that might be in format Chr1:start-end
                    original_chr_id = chr_id
                    offset = 0
                    if ':' in chr_id:
                        chr_parts = chr_id.split(':')
                        chr_id = chr_parts[0]
                        if '-' in chr_parts[1]:
                            offset = int(chr_parts[1].split('-')[0])
                    
                    element = LTRElement(
                        element_id=f"LTR_{element_count}",
                        chr_id=chr_id,
                        element_start=int(fields[0]) + offset,
                        element_end=int(fields[1]) + offset,
                        element_length=int(fields[2]),
                        ltr5_start=int(fields[3]) + offset,
                        ltr5_end=int(fields[4]) + offset,
                        ltr5_length=int(fields[5]),
                        ltr3_start=int(fields[6]) + offset,
                        ltr3_end=int(fields[7]) + offset,
                        ltr3_length=int(fields[8]),
                        similarity=float(fields[9])
                    )
                    elements.append(element)
        
        # Extract sequences
        elements_with_seq = processor.extract_sequences(elements, genome_file)
        
        # Write FASTA
        processor.write_candidates_fasta(elements_with_seq, output_fasta)
        
        print(f"\nSuccessfully generated LTR-RT FASTA file using alternative method: ltr_candidates.fasta")
        print(f"Total LTR-RT sequences: {len(elements_with_seq)}")
    
    finally:
        # Clean up temporary file
        if os.path.exists(temp_harvest_file):
            os.remove(temp_harvest_file)
    
    print("="*50)

def fasta_to_gff(fasta_file, gff_file):
    """
    Parse FASTA file with LTR information and generate GFF3 format file.
    
    Args:
        fasta_file: Input FASTA file path
        gff_file: Output GFF file path
    
    Expected FASTA header format:
        >LTR_2 chr=Chr1 coords=465118:477981 5ltr=465118:466192 3ltr=476907:477981 sim=1.00
    """
    
    with open(fasta_file, 'r') as f_in, open(gff_file, 'w') as f_out:
        # Write GFF3 header
        f_out.write("##gff-version 3\n")
        
        for line in f_in:
            line = line.strip()
            
            # Process only header lines (starting with >)
            if not line.startswith('>'):
                continue
            
            # Remove leading '>' character
            line = line.lstrip('>')
            
            # Extract fields using regular expression
            pattern = r'(\S+)\s+chr=(\S+)\s+coords=(\d+):(\d+)\s+5ltr=(\d+):(\d+)\s+3ltr=(\d+):(\d+)\s+sim=([\d.]+)'
            match = re.match(pattern, line)
            
            if not match:
                print(f"Warning: Cannot parse line: {line}")
                continue
            
            # Parse matched groups
            ltr_id, chrom, start, end, ltr5_start, ltr5_end, ltr3_start, ltr3_end, similarity = match.groups()
            
            # Write LTR retrotransposon main record
            f_out.write(f"{chrom}\t"
                       f"LTR_retrotransposon\t"
                       f"{start}\t"
                       f"{end}\t"
                       f".\t"
                       f"+\t"
                       f".\t"
                       f"ID={ltr_id};Name={ltr_id};ltr_similarity={similarity}\n")
            
            # Write 5'LTR record
            f_out.write(f"{chrom}\t"
                       f"five_prime_LTR\t"
                       f"{ltr5_start}\t"
                       f"{ltr5_end}\t"
                       f".\t"
                       f"+\t"
                       f".\t"
                       f"ID={ltr_id}_5LTR;Parent={ltr_id}\n")
            
            # Write 3'LTR record
            f_out.write(f"{chrom}\t"
                       f"three_prime_LTR\t"
                       f"{ltr3_start}\t"
                       f"{ltr3_end}\t"
                       f".\t"
                       f"+\t"
                       f".\t"
                       f"ID={ltr_id}_3LTR;Parent={ltr_id};ltr_similarity={similarity}\n")
