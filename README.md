# LTR_Checker

**LTR-checker: deep-learning (DL) guided ensemble structural identification of LTR retrotransposon in plant genomes**

## Overview

LTR_Checker is a sophisticated bioinformatics pipeline designed for the detection, annotation, and analysis of Long Terminal Repeat (LTR) retrotransposons in genomic sequences. The tool integrates machine learning-based detection with multiple established LTR detection software and employs a rigorous six-module filtering and validation system to ensure high-quality results.

## Features

- **Deep Learning Detection**: Uses a PyTorch-based convolutional neural network (CNN) model for initial LTR-RT detection
- **Multi-Software Integration**: Supports LTR_FINDER, LTR_HARVEST, and LtrDetector for comprehensive annotation
- **Six-Module Quality Control Pipeline**: Comprehensive filtering and validation system
- **Protein Domain Classification**: Automatic classification into Copia, Gypsy, or Unknown superfamilies
- **Sequence Similarity Filtering**: Removes redundant and highly similar elements
- **Final Validation**: Multi-step cleanup to eliminate false positives

## Dependencies

### Required Software

**Python 3.9+** with the following packages:
- PyTorch 1.10.1+ (with CUDA support optional)
- Biopython 1.85+
- NumPy 1.24.3+
- Pandas 2.2.3+
- Matplotlib 3.9.4+
- Seaborn 0.13.2+
- Scikit-learn
- tqdm 4.67.1+

### External Tools

- **LTR_FINDER**: For LTR retrotransposon detection
- **LTR_HARVEST** (part of GenomeTools): For alternative LTR detection
- **LtrDetector**: For additional LTR detection capabilities
- **BLAST+**: For sequence alignment and comparison
- **HMMER**: For protein domain annotation (hmmsearch)
- **TRF (Tandem Repeats Finder)**: For tandem repeat detection

### Database Files (Included)

- `REXdb_protein_database_viridiplantae_v3.0.hmm`: HMM profiles for plant LTR-RT protein domains
- `TEfam.hmm`: Transposable element family profiles
- `Tpases020812DNA.txt` and `Tpases020812LINE.txt`: Transposase reference sequences

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/yourusername/ltr_checker.git
cd ltr_checker
```

### 2. Install Python Dependencies
```bash
pip install -r requirements.txt
```

### 3. Install External Tools

Ensure that all external tools (LTR_FINDER, GenomeTools, LtrDetector, BLAST+, HMMER, TRF) are installed and accessible in your system PATH.

### 4. Prepare Database Files

The required HMM profiles and reference sequences should be placed in the `dataset/` directory.

## Usage

### Basic Command
```bash
python ltr_checker.py --genome input.fasta --output output_dir --software ltr_finder --threads 8
```

### Command-Line Options
```
Required Arguments:
  --genome FASTA          Input genome file in FASTA format
  --output DIR            Output directory for results

Optional Arguments:
  --software {ltr_finder,ltr_harvest,ltrdetector}
                          LTR detection software to use (default: ltr_finder)
  --threads INT           Number of threads for parallel processing (default: 1)
  --model PATH            Path to pre-trained CNN model (default: model/ltr_cnn.pth)
  --min-length INT        Minimum LTR length (default: 100)
  --max-length INT        Maximum LTR length (default: 6000)
  --min-distance INT      Minimum distance between LTRs (default: 1000)
  --max-distance INT      Maximum distance between LTRs (default: 15000)
  --similarity FLOAT      Minimum LTR similarity threshold (default: 0.8)
  --tsd-length INT        Target site duplication length (default: 5)
  --skip-modules LIST     Comma-separated list of modules to skip (e.g., "1,3,5")
  --keep-intermediate     Keep intermediate files for debugging
```

### Example Workflows

#### 1. Basic Detection with LTR_FINDER
```bash
python ltr_checker.py \
    --genome genome.fasta \
    --output ./results \
    --software ltr_finder \
    --threads 16
```

#### 2. Multi-Software Consensus
```bash
# Run with each software
python ltr_checker.py --genome genome.fasta --output ltr_finder.out --software ltr_finder
python ltr_checker.py --genome genome.fasta --output ltr_harvest.out --software ltr_harvest
python ltr_checker.py --genome genome.fasta --output ltrdetector.out --software ltrdetector

# Merge results (implementation-specific)
```

#### 3. Custom Parameters
```bash
python ltr_checker.py \
    --genome large_genome.fasta \
    --output results/custom \
    --software ltr_harvest \
    --threads 32 \
    --min-length 150 \
    --max-length 5000 \
    --similarity 0.85 \
    --keep-intermediate
```

## Pipeline Workflow

### Stage 1: Deep Learning Detection

The CNN model scans the genome sequence using sliding windows to identify potential LTR-RT candidates.

### Stage 2: Software Annotation

Selected software (LTR_FINDER/LTR_HARVEST/LtrDetector) performs detailed structural annotation of candidates.

### Stage 3: Six-Module Quality Control

**Module 1: Initial Quality Filtering**
- Validates basic structural requirements
- Checks LTR length and distance constraints
- Filters by LTR similarity threshold

**Module 2: LTR Boundary Refinement**
- Refines LTR boundaries using sequence alignment
- Adjusts coordinates for optimal accuracy

**Module 3: Nested Insertion Detection**
- Identifies and filters nested LTR-RTs
- Prevents false positives from nested elements

**Module 4: Protein Domain Annotation**
- Performs HMMER searches against REXdb and TEfam
- Classifies elements into Copia, Gypsy, or Unknown superfamilies
- Filters elements lacking characteristic protein domains

**Module 5: Sequence Similarity Filtering**
- Uses BLAST to identify redundant elements
- Removes highly similar sequences
- Reduces dataset to unique representatives

**Module 6: Final Cleanup and Validation**
- Eliminates false positives from DNA transposons and LINEs
- Validates TSD (Target Site Duplication) presence
- Filters tandem repeats using TRF
- Performs final sequence quality checks

## Output Files
```
output_dir/
├── final_ltr_library.fasta          # High-confidence LTR-RT sequences
├── final_ltr_library.gff3           # Genomic coordinates and annotations
├── classification_summary.txt       # Superfamily classification statistics
├── protein_domains.txt              # Detailed protein domain annotations
├── quality_metrics.txt              # Quality control statistics
├── intermediate/                    # Intermediate files (if --keep-intermediate)
│   ├── module1_output/
│   ├── module2_output/
│   └── ...
└── logs/                            # Processing logs
    └── ltr_checker.log
```

## Performance Considerations

- **Memory Usage**: Proportional to genome size and number of threads
- **Runtime**: Depends on genome size, number of candidates, and selected software
- **GPU Support**: Optional CUDA support for CNN inference
- **Disk Space**: Requires temporary space for intermediate files

## Support

For questions, issues, or feature requests, please contact the author or refer to the source code documentation.

## License

Please refer to the license information provided with the software distribution.                                                                                            │ │
