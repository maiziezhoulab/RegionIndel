import pdb
#pdb.set_trace()
from collections import defaultdict
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import os
import sys

parser = ArgumentParser(description="Variants Calling:")
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=22)
parser.add_argument('--out_dir','-o_dir', help="Directory to store outputs", default='results_/')

args = parser.parse_args()

import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")


def Extract_SNV_info(hap_file,output_file,chr_num):
    f = open(hap_file,"r")
    contig_dict = defaultdict(list)
    curr = 0
    count_2 = 0
    for line in f:

        curr += 1
        data = line.rsplit()
        if data[0] == "V":
            mapq = int(data[5])
            repeat = int(data[4])
            cur_chr_num = data[1]
            if cur_chr_num == "chrX":
                cur_chr_num = "chr23"
            if cur_chr_num == "chr" + str(chr_num): 
                cur_chr_num = data[1]
                ref_start = int(data[2])
                ref_end = int(data[3])
                ref = data[6]
                alt = data[7]
                contig_num = data[8]
                contig_start = int(data[9])
                contig_end = int(data[10])
                strand_dir = data[-1]
                contig_dict[(cur_chr_num,ref_start)].append([cur_chr_num,ref_start, ref_end, ref, alt, contig_num,contig_start,contig_end,strand_dir])

    contig_dict_filter = defaultdict(list)
    count = 0
    for key, val in contig_dict.items():

        contig_dict_filter[key] = val[0]
    
    SNV_dict = defaultdict(list)
    for key, val in contig_dict_filter.items():
        ref = val[3] 
        alt = val[4]
        if len(ref) == 1 and len(alt) == 1:
            SNV_dict[key] = val
        
    pickle.dump(SNV_dict, open(output_file,"wb"))

    return SNV_dict


def Extract_ref_info(hap_file,output_file,chr_num):
    f = open(hap_file,"r")
    contig_dict = defaultdict(list)
    curr = 0
    count_2 = 0
    ref_dict = defaultdict(int)
    for line in f:

        curr += 1
        data = line.rsplit()
        if data[0] == "R":
            start_ = int(data[2])
            end_ = int(data[3])
            for _step in range(start_,end_+1):
                ref_dict[_step] = 1

    pickle.dump(ref_dict, open(output_file,"wb"))

    return ref_dict


def compare_two_haploid_SV(SV_dict_contig_1_file, SV_dict_contig_2_file,ref_dict_contig_1_file,ref_dict_contig_2_file,chr_num,out_dir):
    SV_dict_contig_1 = pickle.load(open(SV_dict_contig_1_file,"rb"))
    SV_dict_contig_2 = pickle.load(open(SV_dict_contig_2_file,"rb"))
    ref_dict_contig_1 = pickle.load(open(ref_dict_contig_1_file,"rb"))
    ref_dict_contig_2 = pickle.load(open(ref_dict_contig_2_file,"rb"))
    count = 0
    count_same = 0
    count_same_wrong = 0
    count_diff = 0
    snv_homo = defaultdict(list)
    snv_hetero = defaultdict(list)
    snv_total = defaultdict(list)
    count_hetero_1 = 0
    count_hetero_2 = 0

    all_candidate_SNV_pos = set(list(SV_dict_contig_1.keys())+list(SV_dict_contig_2.keys()))

    for key in all_candidate_SNV_pos:
        if (key in SV_dict_contig_1) and (key in SV_dict_contig_2):
            ### homo SNV
            val = SV_dict_contig_1[key]
            val_2 = SV_dict_contig_2[key]
            chr_num_1 = val[0]
            chr_num_2 = val_2[0]
            start_1 = val[1]
            end_1 = val[2]
            start_2 = val_2[1]
            end_2 = val_2[2]
            ref_1 = val[3]
            ref_2 = val_2[3]
            alt_1 = val[4]
            alt_2 = val_2[4]

            contig_num_1 = val[-4]
            contig_start_1 = val[-3]
            contig_end_1 = val[-2]
            contig_strand_dir_1 = val[-1]
            contig_num_2 = val_2[-4]
            contig_start_2 = val_2[-3]
            contig_end_2 = val_2[-2]
            contig_strand_dir_2 = val_2[-1]

            snv_homo[key] = [chr_num_1,start_1,end_1,ref_1,alt_1,contig_num_1,contig_start_1,contig_end_1,contig_strand_dir_1,contig_num_2,contig_start_2,contig_end_2,contig_strand_dir_2]  
            snv_total[key] = [chr_num_1,start_1,end_1,ref_1,alt_1,contig_num_1,contig_start_1,contig_end_1,contig_strand_dir_1,contig_num_2,contig_start_2,contig_end_2,contig_strand_dir_2]  

        elif (key in SV_dict_contig_1) and (key[1] in ref_dict_contig_2):

            ### heter SNV (1/0)
            val = SV_dict_contig_1[key]
            chr_num_1 = val[0]
            start_1 = val[1]
            end_1 = val[2]
            ref_1 = val[3]
            alt_1 = val[4]
            contig_num_1 = val[-4]
            contig_start_1 = val[-3]
            contig_end_1 = val[-2]
            contig_strand_dir_1 = val[-1]


            snv_hetero[key] = [val[0],start_1,end_1,ref_1,alt_1,contig_num_1,contig_start_1,contig_end_1,contig_strand_dir_1,'0_PS0:0_hp0',0,0,0]
            snv_total[key] = [val[0],start_1,end_1,ref_1,alt_1,contig_num_1,contig_start_1,contig_end_1,contig_strand_dir_1,'0_PS0:0_hp0',0,0,0]

        elif (key in SV_dict_contig_2) and (key[1] in ref_dict_contig_1):

            ### heter SNV (0/1)
            val_2 = SV_dict_contig_2[key]
            start_2 = val_2[1]
            end_2 = val_2[2]
            ref_2 = val_2[3]
            alt_2 = val_2[4]
            contig_num_2 = val_2[-4]
            contig_start_2 = val_2[-3]
            contig_end_2 = val_2[-2]
            contig_strand_dir_2 = val_2[-1]
            count_hetero_2 += 1
            snv_hetero[key] = [val_2[0],start_2,end_2,ref_2,alt_2,'0_PS0:0_hp0',0,0,0,contig_num_2,contig_start_2,contig_end_2,contig_strand_dir_2]
            snv_total[key] = [val_2[0],start_2,end_2,ref_2,alt_2,'0_PS0:0_hp0',0,0,0,contig_num_2,contig_start_2,contig_end_2,contig_strand_dir_2]

    pickle.dump(snv_homo,open(out_dir + "snv_homo_chr" + str(chr_num) + ".p","wb"))
    pickle.dump(snv_hetero,open(out_dir + "snv_hetero_chr" + str(chr_num) + ".p","wb"))
    pickle.dump(snv_total,open(out_dir + "snv_total_chr" + str(chr_num) + ".p","wb"))
    return (snv_homo,snv_hetero)




def extract_SV(lib_file,lib_prefix,indel_tag,h_tag):
    lib_dict = pickle.load(open(lib_file,"rb"))
    fw_lib_contig1 = open(lib_prefix + "_contig1.bed","w") 
    fw_lib_contig2 = open(lib_prefix + "_contig2.bed","w") 
    fw_lib_ref = open(lib_prefix + "_ref.bed","w") 
    count = 1
    for key, val in lib_dict.items():
        fw_lib_contig1.writelines(str(val[3]) + "\t" + str(val[4]) + "\t" + str(val[5]) + "\t" + indel_tag + "\t" + h_tag+ "\t" + str(count) + "\n")
        fw_lib_contig2.writelines(str(val[6]) + "\t" + str(val[7]) + "\t" + str(val[8]) + "\t" + indel_tag + "\t" + h_tag + "\t" + str(count)+ "\n")
        fw_lib_ref.writelines("chr" + str(val[0]) + "\t" + str(val[1]) + "\t" + str(val[2]) + "\t" + indel_tag + "\t" + h_tag + "\t" + str(count) + "\n")
        count += 1


def Write_vcf_for_SV(hetero_file,homo_file,fw,count_id):
    homo = pickle.load(open(homo_file,"rb"))
    hetero = pickle.load(open(hetero_file,"rb"))
    for key, val in homo.items():
        chr_num = str(val[0])
        start_ = val[1]
        end_ = val[2]
        ref = val[3]
        alt = val[4]
        GT = "1/1"
        contig_num_1 = val[5]
        contig_start_1 = val[6]
        contig_end_1 = val[7]
        contig_strand_dir_1 = val[8]
        contig_num_2 = val[9]
        contig_start_2 = val[10]
        contig_end_2 = val[11]
        contig_strand_dir_2 = val[12]
        new_contig_num_1 = contig_num_1.replace(":","_")
        new_contig_num_2 = contig_num_2.replace(":","_")
        contig_info = str(new_contig_num_1) + "_" + str(contig_start_1) + "_" + str(contig_end_1) + "_" + contig_strand_dir_1 + "_" +  str(new_contig_num_2) + "_" + str(contig_start_2) + "_" + str(contig_end_2) + "_" + contig_strand_dir_2 

        fw.writelines(chr_num + "\t" + str(start_) + "\t" + "event" + str(count_id) + "\t" + ref + "\t" + alt + "\t" + "." + "\t" + "PASS" + "\t" + "SVTYPE=SNP" + "\t" + "GT:Contig" + "\t"  + "1/1" + ":" + contig_info +  "\n")
        count_id += 1

    for key, val in hetero.items():
        chr_num = str(val[0])
        start_ = val[1]
        end_ = val[2]
        ref = val[3]
        alt = val[4]
        GT = "0/1"
        contig_num_1 = val[5]
        contig_start_1 = val[6]
        contig_end_1 = val[7]
        contig_strand_dir_1 = val[8]
        contig_num_2 = val[9]
        contig_start_2 = val[10]
        contig_end_2 = val[11]
        contig_strand_dir_2 = val[12]
        new_contig_num_1 = contig_num_1.replace(":","_")
        new_contig_num_2 = contig_num_2.replace(":","_")
        contig_info = str(new_contig_num_1) + "_" + str(contig_start_1) + "_" + str(contig_end_1) + "_" + str(contig_strand_dir_1) + "_" +  str(new_contig_num_2) + "_" + str(contig_start_2) + "_" + str(contig_end_2) + "_" + str(contig_strand_dir_2) 

        fw.writelines(chr_num + "\t" + str(start_) + "\t" + "event" + str(count_id) + "\t" + ref + "\t" + alt + "\t" + "." + "\t" + "PASS" + "\t" + "SVTYPE=SNP" + "\t" + "GT:Contig" + "\t"  + "0/1"  + ":"  + contig_info + "\n")
        count_id += 1

    return (fw,count_id) 


if __name__ == "__main__":
    out_dir = args.out_dir
    if os.path.exists(out_dir):
        logger.info("using existing output folder: " + out_dir)
    else:
        os.makedirs(out_dir)
    chr_start = args.chr_start
    chr_end = args.chr_end

    count_total_1 = 0
    count_total_2 = 0
    logger.info("---------SNV candidates count")
    for chr_num in range(chr_start,chr_end + 1):
        SNV_dict_contig_1 = Extract_SNV_info(out_dir + "RegionIndel_Contig_chr" + str(chr_num) + "_hp1.var.txt",out_dir + "SNV_dict_contig_1_chr" + str(chr_num) + ".p",chr_num)
        SNV_dict_contig_2 = Extract_SNV_info(out_dir + "RegionIndel_Contig_chr" + str(chr_num) + "_hp2.var.txt",out_dir + "SNV_dict_contig_2_chr" + str(chr_num) + ".p",chr_num)
        logger.info("chr" + str(chr_num) + ":")
        logger.info("haplotype 1: "+str(len(SNV_dict_contig_1)))
        logger.info("haplotype 2: "+str(len(SNV_dict_contig_2)))
        count_total_1 += len(SNV_dict_contig_1)
        count_total_2 += len(SNV_dict_contig_2)

    logger.info("===Total:")
    logger.info("haplotype 1: %d"%count_total_1)
    logger.info("haplotype 2: %d"%count_total_2)

    for chr_num in range(chr_start,chr_end + 1):
        ref_dict_contig_1 = Extract_ref_info(out_dir + "RegionIndel_Contig_chr" + str(chr_num) + "_hp1.var.txt",out_dir + "ref_dict_contig_1_chr" + str(chr_num) + ".p",chr_num)
        ref_dict_contig_2 = Extract_ref_info(out_dir + "RegionIndel_Contig_chr" + str(chr_num) + "_hp2.var.txt",out_dir + "ref_dict_contig_2_chr" + str(chr_num) + ".p",chr_num)

    count_hetero_total = 0
    count_homo_total = 0
    # count_hetero_1_total = 0
    # count_hetero_2_total = 0
    # count_hetero_final = 0
    # count_same_total = 0
    # count_same_wrong_total = 0
    for chr_num in range(chr_start,chr_end + 1):
        #snv_homo,snv_hetero,count_same,count_hetero,count_same_wrong,count_hetero_1,count_hetero_2 = compare_two_haploid_SV(out_dir + "SNV_dict_contig_1_chr"+str(chr_num) + ".p", out_dir + "SNV_dict_contig_2_chr" + str(chr_num) + ".p", out_dir + "ref_dict_contig_1_chr" + str(chr_num) + ".p", out_dir + "ref_dict_contig_2_chr" + str(chr_num) + ".p", chr_num,out_dir)
        snv_homo,snv_hetero = compare_two_haploid_SV(out_dir + "SNV_dict_contig_1_chr"+str(chr_num) + ".p", out_dir + "SNV_dict_contig_2_chr" + str(chr_num) + ".p", out_dir + "ref_dict_contig_1_chr" + str(chr_num) + ".p", out_dir + "ref_dict_contig_2_chr" + str(chr_num) + ".p", chr_num,out_dir)

        count_hetero_total += len(snv_hetero)
        count_homo_total += len(snv_homo)

    logger.info("---------After genotyping")
    logger.info("total_hetero:%d"%count_hetero_total)
    logger.info("total_homo:%d"%count_homo_total)
    logger.info("total_SNV:%d"%(count_hetero_total+count_homo_total))


    vcf_output = out_dir + "RegionIndel_SNPs_chr" + str(chr_num) + ".vcf"
    count_id = 1
    fw = open(vcf_output,"w")
    #fw.writelines("#CHROM" + "\t" + "POS" + "\t" + "ID" + "\t" + "REF" + "\t" + "ALT" + "\t" + "QUAL" + "\t" + "FILTER" + "\t" + "INFO" + "\t" + "FORMAT" + "\t"  + "SAMPLE" + "\n")

    for chr_num in range(chr_start,chr_end + 1):
        fw,count_id = Write_vcf_for_SV(out_dir + "snv_hetero_chr" + str(chr_num) + ".p",out_dir + "snv_homo_chr" + str(chr_num) + ".p", fw,count_id)















