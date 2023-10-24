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

Run the whole pipeline (You need to provide the OER folder generated in step 0):
```
python3 RegionIndel/bin/RegionIndel_step1.py  --bam_file selected.bam --vcf_file test_freebayes.vcf --chr_num 3 --out_dir test_sv --OER_dir ./OER/part2


python3 RegionIndel/bin/RegionIndel_step2.py --out_dir test_sv --chr_num 3 


python3 RegionIndel/bin/RegionIndel_step3.py  --assembly_dir test_sv  -o_dir test_sv --ref_file genome_hg19.fa  --chr_num 3 


python3 RegionIndel/bin/remove_redundancy.py   \
-i ./test_sv/RegionIndel_Step3_Results/RegionIndel_Contig_final_sorted_sv.vcf  \
-o ./test_sv/Remove_redundancy/
```
