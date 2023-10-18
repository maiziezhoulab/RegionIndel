
#-----------------Parameter------------

bed_file="test.bed"  # Specify your BED file
wgsbam=possorted_bam.bam #your whole genome bam file
wgsvcf=possorted_bam_freebayes.vcf #your whole genome VCF file
oerdir=./OER/part2   # the OER file generated in step 0
reference=genome_hg19.fa # reference file
RegionInde_dir=./RegionIndel/

#-----------------Pipeline------------

# Read the BED file line by line and extract chr, start, and end

while IFS=$'\t' read -r chr start end rest; do

	echo "Processing region--------Chromosome: $chr, Start: $start, End: $end"

	# specify output folder
	outdir=./out_${chr}_${start}_${end}/
	# create output folder
	mkdir -p $outdir

	# extract region bam file
	samtools view $wgsbam -Sb ${chr}:${start}-${end} > $outdir/selected.bam
	# create index for region bam file
	samtools index $outdir/selected.bam


	# step1 phasing
	python3 ${RegionInde_dir}/bin/RegionIndel_step1.py  \
	--bam_file $outdir/selected.bam \
	--vcf_file ${wgsvcf} \
	--chr_num $chr \
	--out_dir $outdir \
	--fd $oerdir


	# step2 assmebly
	python3 ${RegionInde_dir}/bin/RegionIndel_step2.py \
	--out_dir $outdir/ \
	--chr_num $chr \
	--reference $reference


	# step3 variant call
	python3 ${RegionInde_dir}/bin/RegionIndel_step3.py  \
	--assembly_dir $outdir/  \
	--ref_file $reference  \
	--chr_num $chr

	# step4 remove redundancy
	python3 ${RegionInde_dir}/bin/remove_redundancy.py.py   \
	-i $outdir/RegionIndel_Step3_Results/RegionIndel_Contig_final_sorted_sv.vcf  \
	-o $outdir/Remove_redundancy/


done < "$bed_file"

# collect and merge region VCF
cat $outdir/Remove_redundancy/RegionIndel_no_redundancy.vcf|grep '#' > RegionIndel_merged.vcf
cat out_chr*_*_*/Remove_redundancy/RegionIndel_no_redundancy.vcf|grep -v '#' >> RegionIndel_merged.vcf

# remove redundancy for region VCF
python3 ${RegionInde_dir}/bin/remove_redundancy.py.py   \
-i  RegionIndel_merged.vcf \
-o ./Remove_redundancy/













