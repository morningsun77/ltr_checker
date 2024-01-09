# ltr_checker
# A de-novo tool to annotate LTR-RT
1.Install  
conda create -n ltr_checker python=3.8  
  
conda activate ltr_checker  
  
pip3 install -r requirements.txt  

2.Workflow

![image](https://github.com/morningsun77/ltr_checker/blob/main/Image/workflow%20of%20ltr_checker.png)
  
3.Usage

```
python3 ltr_finder.py --genome <genome_file> --threads <num_threads> --output <output_path> --stride <stride_value> --max <max_distance> --min <min_distance> --max_ltr <max_ltr_length> --min_ltr <min_ltr_length> --tgca <tgca_option> --tsd <tsd_option> --model <model_path> --split <num_segments>
```

