#### We provide an example dataset to run the whole pipeline. 

Please download the example data <a href="https://zenodo.org/record/5117764">from Zenodo</a>.
```
AquilaSV_exampledata
|-selected.bam (hg19)
|-selected.bam.bai
|
|-test_freebayes.vcf (hg19)
|
|-genome_hg19.fa         
```

Run the whole pipeline:
```
python AquilaSV/bin/AquilaSV_step1.py --bam_file selected.bam --vcf_file test_freebayes.vcf --chr_num 3 --out_dir test_sv

python AquilaSV/bin/AquilaSV_step2.py --out_dir test_sv --chr_num 3 --reference genome_hg19.fa

python AquilaSVbin/AquilaSV_step3.py  --assembly_dir test_sv  --ref_file genome_hg19.fa  --chr_num 3 
```
