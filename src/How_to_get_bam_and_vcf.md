### To generate BAM file for the target region:

#### step1. generate the bed file "selected.bed". 
The bed file contains the chromosome name, start position and end position of your target region. eg. `chr1 123456 223456`. 

#### step2. generate the "selected.bam" by a WGS bam file (or a bam file from one specific chromosome) and the bed file. 

`samtools -b -L selected.bed wgs.bam > selected.bam`

`samtools index selected.bam`

:octocat: To generate the original WGS bam file for 10X or stLFR linked-reads, you can run LongRanger/EMA/BWA. 
Check <a href="https://github.com/maiziex/Aquila_stLFR/blob/master/src/How_to_get_bam_and_vcf.md">here</a> for more details. 


### To generate VCF file through FreeBayes:
`freebayes -f genome_hg19.fa wgs.bam > test_freebayes.vcf`

or 

`freebayes -f genome_hg19.fa selected.bam > test_freebayes.vcf`

Once you generate the `selected.bam` and `test_freebayes.vcf` in your working directory, you can use use it to run AquilaSV pipeline.
