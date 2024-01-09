# ltr_checker
# A de-novo tool to annotate LTR-RT
1.Install  
conda create -n ltr_checker python=3.8  
  
conda activate ltr_checker  
  
pip3 install -r requirements.txt  

2.Workflow

![image](https://github.com/morningsun77/ltr_checker/blob/main/Image/workflow%20of%20ltr_checker.png)
  
3.Usage

·python ltr_finder.py --genome example_genome.fasta --threads 8 --output output_directory --stride 5000 --max 20000 --min 5000 --max_ltr 8000 --min_ltr 200 --tgca yes --tsd no --model model_file.pkl --split 50·
