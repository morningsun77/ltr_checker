# LTR_checker

**LTR_checker: deep-learning (DL) guided ensemble structural identification of LTR retrotransposon in plant genomes**

## Overview

LTR_Checker is a sophisticated bioinformatics pipeline designed for the detection, annotation, and analysis of Long Terminal Repeat retrotransposons (LTR-RTs) in genomic sequences. The tool integrates CNN detection with multiple established LTR-RT identification methods and employs a rigorous six-module filtering and validation procedure to ensure high-quality results.

## Dependencies

**Python 3.9+** with the following main packages:
- PyTorch 1.11.0 (with CUDA support optional)
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

### 3. Install PyTorch
```bash
pip install torch==1.11.0+cu113 torchvision==0.12.0+cu113 torchaudio==0.11.0 --extra-index-url https://download.pytorch.org/whl/cu113
```

### 4. Install External Tools

Ensure all external tools listed above are installed and accessible in your system path.

### 5. Prepare Database Files

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
  -h, --help            show help message and exit
  --threads THREADS     number of threads, default=10
  --stride STRIDE       stride size for sliding windows in LTR-RT identification, in bp, default=10000
  --max MAX             max separation distance of two LTR,default=15000
  --min MIN             min separation distance of two LTR, default=1000
  --max_ltr MAX_LTR     max length of LTR, default=7000
  --min_ltr MIN_LTR     min length of LTR, default=100
  --tgca TGCA           whether require the presence of TGCA, default=no
  --tsd TSD             whether require the presence of TSD, default=no
  --model MODEL         path to the CNN model for LTR-RT identification
  --split SPLIT         how many chromosme segments you want to split into, default=2
  --device DEVICE       cpu or cuda, default=cpu
  --method MRTHOD       method to use: ltr_finder, ltr_harvest, ltrdetector, or all, default=ltr_finder
  --identity IDENTITY   minimum identity between 5' and 3' LTRs for LtrDetector, default=85
  --filter FILTER       whether to run filtering modules (Module 1-6), options: yes/no, default=yes
  --output_format OUTPUT_FORMAT, ltr_harvest, ltr_finder, ltrdetector and all.
                        Unified output format when using --method all: ltr_harvest, ltr_finder, or ltrdetector, default=ltr_harvest
```

### Example

#### 1. Single or Ensemble identification
```bash
# Run with each method, the single identification
python ltr_checker.py --genome test/test.fa --output output_dir --method ltr_finder
python ltr_checker.py --genome test/test.fa --output output_dir --method ltr_harvest
python ltr_checker.py --genome test/test.fa --output output_dir --method ltrdetector

# Run with all methods, the ensemble identification
python ltr_checker.py --genome test/test.fa --output output_dir --method all
```
**Note**: A single method tool can be selected, or all three methods tools can be combined together. Using `--method all` provides more comprehensive results by integrating predictions from LTR_FINDER, LTR_HARVEST, and LtrDetector, which can improve detection accuracy and coverage.

#### 2. Skip Quality Control Filtering
```bash
python ltr_checker.py \
    --genome test/test.fa \
    --output ./results \
    --method all \
    --filter no \
```

**Note**: The `--filter` parameter controls whether to apply the six-module quality control pipeline (default: `yes`). Setting `--filter no` skips all filtering modules and directly outputs results from the selected method (ltr_finder.out, ltr_harvest.out, and/or ltr_detector.out). This option is faster but may include more false positives.

## Pipeline Workflow

### Stage 1: Deep Learning Detection

The CNN model scans the genome sequence using sliding windows to identify potential LTR-RT candidates.

### Stage 2: Method Annotation

Selected method (LTR_FINDER/LTR_HARVEST/LtrDetector) performs detailed structural annotation of candidates.

### Stage 3: Six-Module Quality Control

**Module 1: Removal of gapped sequences, tandem repeats, and length validation**
**Module 2: boundary refinement and target site duplication (TSD) verification**
**Module 3: calculation of LTR divergence and filtering**
**Module 4: Protein domain identification and classification**
**Module 5: Filtering candidate with significant non-LTR sequence**
**Module 6: Validation of nested insertion**


## Output Files
```
output_dir/
├── ltr_candidates.fasta             # Initial LTR-RT candidates from CNN detection
├── ltr_finder.out                   # LTR_FINDER results (if --method ltr_finder or all)
├── ltr_harvest.out                  # LTR_HARVEST results (if --method ltr_harvest or all)
├── ltr_detector.out                 # LtrDetector results (if --method ltrdetector or all)
├── unified_output.txt               # All methods results (if --method all)
├── final_ltr_rt.fasta               # High-confidence LTR-RT sequences (final output)
├── final_ltr_rt.gff3                # Genomic coordinates and annotations
├── module1/                         
├── module2/                        
├── module3/                         
├── module4/                         
├── module5/                         
├── module6/                         
```

## Performance Considerations

- **Memory Usage**: Proportional to genome size and number of threads
- **Runtime**: Depends on genome size, number of candidates, and selected method
- **GPU Support**: Optional CUDA support for CNN inference
- **Disk Space**: Requires temporary space for intermediate files

## Support

For questions, issues, or feature requests, please contact Zhaoyang Chen(zhao-yamg.chen@umu.se) or Jianfeng Mao(jianfeng.mao@umu.se), or refer to the source code documentation.
