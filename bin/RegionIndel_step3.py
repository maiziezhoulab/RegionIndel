import pdb
#pdb.set_trace()
import os
from argparse import ArgumentParser
from subprocess import Popen, PIPE
from multiprocessing import Pool,cpu_count,active_children,Manager
import time
import sys
import json
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 
config_dc = json.load(open(code_path+'/config.txt'))
paftools_path = config_dc['paftools']
k8_path = config_dc['k8']

parser = ArgumentParser(description="Run depth all:")
parser.add_argument('--chr_num','-start',type=int,help="chromosome number for target variant or region", required=True)
parser.add_argument('--var_size','-v',type=int,help="variant size, cut off size for indel and SV, default = 1", default=1)
parser.add_argument('--num_of_threads','-t',type=int,help="number of threads, default = 1", default=1)
parser.add_argument('--assembly_dir','-i_dir', help="Required parameter, folder to store RegionIndel assembly results at RegionIndel assembly steps",required=True)
parser.add_argument('--out_dir','-o_dir', help="Directory to store outputs, default = ./RegionIndel_Step3_Results", default='RegionIndel_Step3_Results/')
parser.add_argument('--ref_file','-r', help="Required parameter, human reference fasta file",required=True)
parser.add_argument('--delete_temp_file','-d', action='store_true')

global num_of_threads
args = parser.parse_args()
#all_regions_flag = args.all_regions_flag
all_regions_flag = 1
#clean_flag = args.clean_flag
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 
__author__ = "Maizie&Can@Vandy"

import logging
## set logger
logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
            level=logging.INFO,
                datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")


def Split_supercontig_by_haplotype(fasta_file,xin):
    f = open(fasta_file,"r")
    output_prefix = fasta_file.split(".fasta")[0]
    hp1_file = output_prefix + "_hp1.fasta"
    hp2_file = output_prefix + "_hp2.fasta"
    fw_hp1 = open(hp1_file,"w")
    fw_hp2 = open(hp2_file,"w")
    count = 0
    for line in f:
        data = line.rsplit()
        if count%2 == 0:
            hp_flag = data[0].split("_")[-1]
            if hp_flag == "hp1":
                fw_hp1.writelines(line)
            elif hp_flag == "hp2":
                fw_hp2.writelines(line)
        elif count%2 == 1:
            if hp_flag == "hp1":
                fw_hp1.writelines(line)
            elif hp_flag == "hp2":
                fw_hp2.writelines(line)
        count += 1
    f.close()
    fw_hp1.close()
    fw_hp2.close()
    #print("done")


def Split_haplotype(chr_start,chr_end,in_dir):
    pool = Pool(processes=chr_end - chr_start +1)
    for chr_num in range(chr_start,chr_end + 1):
        input_file = in_dir + "RegionIndel_Contig_chr" + str(chr_num) + ".fasta"
        pool.apply_async(Split_supercontig_by_haplotype,(input_file,"xin"))
    pool.close()
    while len(active_children()) > 1:
        time.sleep(0.5)
    pool.join()
    logger.info("finish splitting haplotype")


def get_paf(ref_file,hap_file,hap_paf,xin):
    try:
        use_cmd =  "minimap2 -cx asm5 -t8 --cs " + ref_file + " " + hap_file  + " > " + hap_paf
    except:
        use_cmd =  code_path + "minimap2/" + "minimap2 -cx asm5 -t22 --cs " + ref_file + " " + hap_file  + " > " + hap_paf
    logger.info(use_cmd)
    Popen(use_cmd,shell=True).wait()


def sort_paf(hap_paf,hap_paf_sorted,xin):
    use_cmd = "sort -k6,6 -k8,8n " + hap_paf + " > "  + hap_paf_sorted  
    Popen(use_cmd,shell=True).wait()


def get_var(hap_paf_sorted,hap_var_txt,xin):
    #code_path1 = "/data/maiziezhou_lab/CanLuo/Software/miniconda3/lib/python3.8/site-packages/bin/"
    #os.chmod(k8_path +"/k8-Linux",0o777)
    #change_mode_cmd = "chmod 777 "+code_path1
    #Popen(change_mode_cmd,shell=True).wait()
    use_cmd = k8_path + "/k8-Linux "  + paftools_path  +  "/paftools.js " + " call -l 1 -L 1 -q 20 " +  hap_paf_sorted + " > " +  hap_var_txt 
    Popen(use_cmd,shell=True).wait()
    

def assembly_based_variants_call_paf(chr_start,chr_end,ref_file,in_dir,out_dir):
    for chr_num in range(chr_start,chr_end + 1):
        hap1_file = in_dir + "RegionIndel_Contig_chr" + str(chr_num) + "_hp1.fasta"
        hap2_file = in_dir + "RegionIndel_Contig_chr" + str(chr_num) + "_hp2.fasta"
        hap1_paf = out_dir + "RegionIndel_Contig_chr" + str(chr_num) + "_hp1.paf"
        hap2_paf = out_dir + "RegionIndel_Contig_chr" + str(chr_num) + "_hp2.paf"
        logger.info("start getting paf...")
        pool = Pool(processes=2)
        pool.apply_async(get_paf,(ref_file,hap1_file,hap1_paf,"xin"))
        pool.apply_async(get_paf,(ref_file,hap2_file,hap2_paf,"xin"))
        pool.close()
        while len(active_children()) > 1:
            time.sleep(0.5)
        pool.join()

    logger.info("finish extracting paf files")


def assembly_based_variants_call_sort(chr_start,chr_end,out_dir):
    for chr_num in range(chr_start,chr_end + 1):
        hap1_paf = out_dir + "RegionIndel_Contig_chr" + str(chr_num) + "_hp1.paf"
        hap2_paf = out_dir + "RegionIndel_Contig_chr" + str(chr_num) + "_hp2.paf"
        hap1_paf_sorted = out_dir + "RegionIndel_Contig_chr" + str(chr_num) + "_hp1.paf.sorted" 
        hap2_paf_sorted = out_dir + "RegionIndel_Contig_chr" + str(chr_num) + "_hp2.paf.sorted" 
        pool = Pool(processes=2)
        pool.apply_async(sort_paf,(hap1_paf,hap1_paf_sorted,"xin"))
        pool.apply_async(sort_paf,(hap2_paf,hap2_paf_sorted,"xin"))
        pool.close()
        while len(active_children()) > 1:
            time.sleep(0.5)
        pool.join()

    logger.info("finish sorting paf files")
 

def assembly_based_variants_call(chr_start,chr_end,out_dir):
    for chr_num in range(chr_start,chr_end + 1):
        hap1_paf_sorted = out_dir + "RegionIndel_Contig_chr" + str(chr_num) + "_hp1.paf.sorted"
        hap2_paf_sorted = out_dir + "RegionIndel_Contig_chr" + str(chr_num) + "_hp2.paf.sorted"
        hap1_var_txt = out_dir + "RegionIndel_Contig_chr" + str(chr_num) + "_hp1.var.txt"
        hap2_var_txt = out_dir + "RegionIndel_Contig_chr" + str(chr_num) + "_hp2.var.txt"
        pool = Pool(processes=2)
        pool.apply_async(get_var,(hap1_paf_sorted,hap1_var_txt,"xin"))
        pool.apply_async(get_var,(hap2_paf_sorted,hap2_var_txt,"xin"))
        pool.close()
        while len(active_children()) > 1:
            time.sleep(0.5)
        pool.join()

    logger.info("finish calling variants from paf files")


def Run_del_or_ins(all_cmd,xin):
    Popen(all_cmd,shell=True).wait()


def Call_SNV_info_from_contigs(chr_start,chr_end,out_dir,num_of_threads):
    total_num = chr_end - chr_start +1
    pool = Pool(num_of_threads)
    count = 1
    for chr_num in range(chr_start,chr_end + 1):
        all_cmd = "python3 " + code_path + "Extract_SNV_info_from_contigs_forcontiginfo_forall.py  --chr_start " + str(chr_num) + " --chr_end " + str(chr_num) + " --out_dir " + out_dir 
        count += 1
        pool.apply_async(Run_del_or_ins,(all_cmd,"xin"))
        if (count - 1)%num_of_threads == 0 or (count - 1) == total_num:
            pool.close()
            while len(active_children()) > 1:
                time.sleep(0.5)
            pool.join()

            if (count - 1) == total_num:
                print("finished chr" + str(chr_num))
            else:
                pool = Pool(num_of_threads)
    curr_file = out_dir + "RegionIndel_SNPs.vcf"
    exists = os.path.exists(curr_file)
    if exists:
        Popen("rm -rf " + curr_file,shell=True).wait()
    for chr_num in range(chr_start,chr_end + 1):
        one_vcf = out_dir + "RegionIndel_SNPs_chr" + str(chr_num) + ".vcf"
        cat_cmd = "cat " + one_vcf  + " >> " + out_dir + "RegionIndel_SNPs.vcf"
        Popen(cat_cmd,shell=True).wait()
    print("finish SNV calling")


def Call_SV_del_from_contigs(chr_start,chr_end,out_dir,num_of_threads,v_size):
    total_num = chr_end - chr_start +1
    pool = Pool(num_of_threads)
    count = 1
    for chr_num in range(chr_start,chr_end + 1):
        if all_regions_flag == 1:
            all_cmd = "python3 " + code_path + "Extract_DEL_allregions.py --chr_start " + str(chr_num) + " --chr_end " + str(chr_num) + " --out_dir " + out_dir + " --var_size " + str(v_size)
       # else:
         #   all_cmd = "python3 " + code_path + "Extract_SV_info_from_contigs_use_overlap_for_del_forcontiginfo.py --chr_start " + str(chr_num) + " --chr_end " + str(chr_num) + " --out_dir " + out_dir + " --var_size " + str(v_size)

        count += 1
        pool.apply_async(Run_del_or_ins,(all_cmd,"xin"))
        if (count - 1)%num_of_threads == 0 or (count - 1) == total_num:
            pool.close()
            while len(active_children()) > 1:
                time.sleep(0.5)
            pool.join()

            if (count - 1) == total_num:
                print("finished chr" + str(chr_num))
            else:
                pool = Pool(num_of_threads)
    curr_file = out_dir + "RegionIndel_DEL.vcf"
    exists = os.path.exists(curr_file)
    if exists:
        Popen("rm -rf " + curr_file,shell=True).wait()
    for chr_num in range(chr_start,chr_end + 1):
        one_vcf = out_dir + "RegionIndel_DEL_chr" + str(chr_num) + ".vcf"
        cat_cmd = "cat " + one_vcf  + " >> " + out_dir + "RegionIndel_DEL.vcf"
        Popen(cat_cmd,shell=True).wait()

    print("finish DEL calling")


def Call_SV_ins_from_contigs(chr_start,chr_end,out_dir,num_of_threads,v_size):
    total_num = chr_end - chr_start +1
    pool = Pool(num_of_threads)
    count = 1
    for chr_num in range(chr_start,chr_end + 1):
        if all_regions_flag == 1:
            all_cmd = "python3 " + code_path + "Extract_INS_allregions.py --chr_start " + str(chr_num) + " --chr_end " + str(chr_num) + " --out_dir " + out_dir + " --var_size " + str(v_size)
        #else:
            
            #all_cmd = "python3 " + code_path + "Extract_SV_info_from_contigs_use_shift_for_ins2_forcontiginfo.py --chr_start " + str(chr_num) + " --chr_end " + str(chr_num) + " --out_dir " + out_dir + " --var_size " + str(v_size)

        count += 1
        #Run_del_or_ins(all_cmd,"xin")
        pool.apply_async(Run_del_or_ins,(all_cmd,"xin"))
        if (count - 1)%num_of_threads == 0 or (count - 1) == total_num:
            pool.close()
            while len(active_children()) > 1:
                time.sleep(0.5)
            pool.join()
            if (count - 1) == total_num:
                print("finished chr" + str(chr_num))
            else:
                pool = Pool(num_of_threads)
    curr_file = out_dir + "RegionIndel_INS.vcf"
    exists = os.path.exists(curr_file)
    if exists:
        Popen("rm -rf " + curr_file,shell=True).wait()
    for chr_num in range(chr_start,chr_end + 1):
        one_vcf = out_dir + "RegionIndel_INS_chr" + str(chr_num) + ".vcf"
        cat_cmd = "cat " + one_vcf  + " >> " + out_dir + "RegionIndel_INS.vcf"
        Popen(cat_cmd,shell=True).wait()

    print("finish INS calling")



def Merge_all_variants(out_dir):
    snp_vcf = out_dir + "RegionIndel_SNPs.vcf"
    del_vcf = out_dir + "RegionIndel_DEL.vcf"
    ins_vcf = out_dir + "RegionIndel_INS.vcf"
    final_vcf = out_dir + "RegionIndel_final.vcf"
    final_sort_vcf = out_dir + "RegionIndel_final_sorted.vcf"
    use_cmd = "cat " + snp_vcf + " "  + del_vcf + " " + ins_vcf + " > " + final_vcf 
    sort_cmd = "cat " + final_vcf + " | awk '$1 ~ /^#/ {print $0;next} {print $0 | \"LC_ALL=C sort -k1,1 -k2,2n\"}' > " + final_sort_vcf
    Popen(use_cmd,shell=True).wait()
    Popen(sort_cmd,shell=True).wait()
    print("finish merging all variants")


def Clean_data(out_dir):
    rm_cmd = "rm " + out_dir + "/*chr*"
    Popen(rm_cmd,shell=True).wait()

def delete_files(del_dir):
    del_cmd = "rm -r "+del_dir
    Popen(del_cmd,shell=True).wait()
    return()


def main():
    if len(sys.argv) == 1:
        Popen("python3 " + "RegionIndel_step3.py -h",shell=True).wait()
    else:
        out_dir = args.out_dir + "/RegionIndel_Step3_Result/"
        if os.path.exists(out_dir):
            print("using existing output folder: " + out_dir)
        else:
            os.makedirs(out_dir)
        chr_start = args.chr_num
        chr_end = args.chr_num
        v_size = args.var_size
        num_of_threads = args.num_of_threads
        in_dir = args.assembly_dir + "/" + "Assembly_Contigs_files/"  
        ref_file = args.ref_file
        deletion_mode=args.delete_temp_file


        Split_haplotype(chr_start,chr_end,in_dir)
        assembly_based_variants_call_paf(chr_start,chr_end,ref_file,in_dir,out_dir)
        assembly_based_variants_call_sort(chr_start,chr_end,out_dir)
        assembly_based_variants_call(chr_start,chr_end,out_dir)
        Call_SNV_info_from_contigs(chr_start,chr_end,out_dir,num_of_threads)
        Call_SV_del_from_contigs(chr_start,chr_end,out_dir,num_of_threads,v_size)
        Call_SV_ins_from_contigs(chr_start,chr_end,out_dir,num_of_threads,v_size)
        Merge_all_variants(out_dir)



        #prepare files
        cmd = "bgzip -c %sRegionIndel_final_sorted.vcf >  %sRegionIndel_final_sorted.vcf.gz;tabix -fp vcf %sRegionIndel_final_sorted.vcf.gz"%(out_dir,out_dir,out_dir)
        Popen(cmd,shell=True).wait()
        #reformat
        cmd = "python3  %sReformat.py -r  %s -i %sRegionIndel_final_sorted.vcf -o %sRegionIndel_Contig_final.vcf --add_header 19 --base_norm --gz_tbi"%(code_path,ref_file,out_dir,out_dir)
        Popen(cmd,shell=True).wait()



        if deletion_mode:
            files = [out_dir+'/'+file_  for file_ in os.listdir(out_dir) if "Contig_final" not in file_]
            for path in files:
                os.system("rm "+path)







if __name__ == "__main__":
    main()
