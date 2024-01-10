# LTR-checker: deep-learning guided structural identification of LTR retrotransposon providing decreased computation demand and high flexibility
1. **Install**

conda create -n ltr_checker python=3.8  
  
conda activate ltr_checker  
  
pip3 install -r requirements.txt  

conda install ltr_finder -c bioconda  

2. **Workflow**

![image](https://github.com/morningsun77/ltr_checker/blob/main/Image/workflow%20of%20ltr_checker.png)
  
3. **Usage**

```
python3 ltr_finder.py --genome <genome_file> --threads <num_threads> --output <output_path> --stride <stride_value> --max <max_distance> --min <min_distance> --max_ltr <max_ltr_length> --min_ltr <min_ltr_length> --tgca <tgca_option> --tsd <tsd_option> --model <model_path> --split <num_segments>
```

4. **Options**
   
**·** --genome: Path to the input genome file (required).

**·** --threads: Number of threads to use (default=10).

**·** --output: Path to the output directory (required).

**·** --stride: Stride value for processing (default=10000).

**·** --max: Maximum separation distance of two LTR elements (default=15000).

**·** --min: Minimum separation distance of two LTR elements (default=1000).

**·** --max_ltr: Maximum length of LTR elements (default=7000).

**·** --min_ltr: Minimum length of LTR elements (default=100).

**·** --tgca: Specify whether TGCA is needed (default=no).

· --tsd: Specify whether TSD is needed (default=no).

· --model: Path to the machine learning model file (required).

· --split: Number of chromosome segmentation (default=100).


