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

**Python 3.9+** with the following main packages:
- PyTorch 1.10.1+ (with CUDA support optional)
- Biopython 1.85+
- NumPy 1.24.3+
- Pandas 2.2.3+
- Matplotlib 3.9.4+
- Seaborn 0.13.2+
- Scikit-learn
- tqdm 4.67.1+

### External Tools

- **TRF (Tandem Repeats Finder) 4.10.0** - Tandem repeat detection
- **BLAST+ (RMBlast variant) 2.14.1+** - Sequence alignment and comparison
- **HMMER 3.4** - Protein domain annotation (hmmsearch)
- **EMBOSS 6.6.0.0** - Sequence analysis utilities
- **RepeatMasker 4.1.8** - Repeat sequence masking
- **LTR_FINDER 1.07** - LTR retrotransposon detection
- **GenomeTools 1.6.6** (includes LTR_HARVEST) - Alternative LTR retrotransposon detection
- **LtrDetector** - Additional LTR retrotransposon detection capabilities
- **MUSCLE 5.1.0**

### Database Files (Included)

- `REXdb_protein_database_viridiplantae_v3.0.hmm` - HMM profiles for plant LTR-RT protein domains
- `TEfam.hmm` - Transposable element family profiles
- `Tpases020812DNA.txt` and `Tpases020812LINE.txt` - Transposase reference sequences

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/yourusername/ltr_checker.git
cd ltr_checker
```

### 2. Install Python Dependencies
```bash
conda create -n ltr_checker python=3.9
conda activate ltr_checker
pip install -r requirements.txt
```

### 3. Install External Tools

Ensure all external tools listed above are installed and accessible in your system PATH.

### 4. Prepare Database Files

Place the required HMM profiles and reference sequences in the `dataset/` directory.

## Usage

### Basic Command
```bash
python ltr_checker.py --genome input.fasta --output output_dir
```

### Command-Line Options
```
Required Arguments:
  --genome FASTA          Input genome file in FASTA format
  --output DIR            Output directory for results

Optional Arguments:
  -h, --help            show this help message and exit
  --genome GENOME       genome file
  --threads THREADS     number of threads, default=10
  --output OUTPUT       output path
  --stride STRIDE       stride, default=10000
  --max MAX             max separation distance of two LTR,default=15000
  --min MIN             min separation distance of two LTR, default=1000
  --max_ltr MAX_LTR     max length of LTR, default=7000
  --min_ltr MIN_LTR     min length of LTR, default=100
  --tgca TGCA           whether need TGCA, default=no
  --tsd TSD             whether need TSD, default=no
  --model MODEL         path of model
  --split SPLIT         chromosome segmentation number, default=100
  --device DEVICE       cpu or cuda, default=cpu
  --software SOFTWARE   software to use: ltr_finder, ltr_harvest, ltrdetector, or all, default=ltr_finder
  --identity IDENTITY   minimum identity between 5' and 3' LTRs for LtrDetector, default=85
  --output_format OUTPUT_FORMAT, ltr_harvest, ltr_finder, ltrdetector and all.
                        Unified output format when using --software all: ltr_harvest, ltr_finder, or ltrdetector, default=ltr_harvest
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
