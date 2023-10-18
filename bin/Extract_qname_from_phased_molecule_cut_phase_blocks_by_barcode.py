
from collections import defaultdict
import pickle
import os
from argparse import ArgumentParser

import pysam
from subprocess import Popen
from retrieve_bc import *

from tqdm import tqdm


def reverse_complement(seq):
    tab = str.maketrans("ACTGN", "TGACN")
    return seq.translate(tab)[::-1]


def bam2fq_bad_bc(bamfile, qname_pb_hp_dict , output_dir):

    samfile = pysam.AlignmentFile(bamfile)
    read1_dc = {}
    read2_dc = {}

    for read in samfile.fetch():
        if read.is_secondary is False:
            bc = -1
            for tup in read.get_tags():
                if tup[0] == 'BX':
                    bc = tup[1]
                    break 

            if bc == '0_0_0':
                if read.is_read1:
                    if read.is_reverse:
                        seq = reverse_complement(read.seq)
                        qual = read.qual[::-1]
                    else:
                        seq = read.seq 
                        qual = read.qual
                    read1_dc[read.qname] = [seq, qual]
                elif read.is_read2:
                    if read.is_reverse:
                        seq = reverse_complement(read.seq)
                        qual = read.qual[::-1]
                    else:
                        seq = read.seq 
                        qual = read.qual
                    read2_dc[read.qname] = [seq, qual]

    for qname in read1_dc:
        if qname in read2_dc:
            seq1,qual1 = read1_dc[qname]
            seq2,qual2 = read2_dc[qname]
            if len(seq1)==len(seq2):
                lines = '@'+qname+'\tBX:Z:0_0_0\n'+\
                seq1+'\n'+\
                '+\n'+\
                qual1+'\n'+\
                '@'+qname+'\tBX:Z:0_0_0\n'+\
                seq2+'\n'+\
                '+\n'+\
                qual2+'\n'

                for hp_list in qname_pb_hp_dict[qname]:
                    for key in hp_list:
                        _PS_flag = key[0]
                        _HP_flag = key[1]
                        outfile = output_dir +"/fastq_by_" + str(_PS_flag)  + "_" + _HP_flag + ".fastq"
                        with open(outfile,'a') as f:
                            f.write(lines)

    return  


def extract_bad_reads(bam_file, out_dir, num_threads, qname_pb_hp_dict ):
    bam_sorted_file = out_dir + "/bam_sort_by_qname.bam"


    sort_bam_cmd = "samtools sort -@ " + str(num_threads) + " -n " + bam_file + " -o " + bam_sorted_file

    Popen(sort_bam_cmd,shell=True).wait()

    bam2fq_bad_bc(bam_file,  qname_pb_hp_dict, out_dir )


def Extract_mole_num_from_phased_mole(phased_file,chr_num):
    phase_dict = defaultdict(list)
    f = open(phased_file,"r")
    for line in f:
        data = line.rsplit()
        PS_flag_info = data[-2]
        PS_flag = int(PS_flag_info.split(":")[1])
        HP_flag = data[-1]
        mole_num = int(data[6])
        phase_dict[(PS_flag,HP_flag)].append(mole_num)

    return phase_dict

def load_qname_bc_dict(bamfile):

    samfile = pysam.AlignmentFile(bamfile)
    dc = {}
    for read in samfile.fetch():
        bc = -1
        for tup in read.get_tags():
            if tup[0] == 'BX':
                bc = tup[1]
                break 

        if (bc!=-1) :
            dc[read.qname] = bc 

    return dc 

def combine_dict(phased_dict_mole_num,mole_qname_dict,qname_bc_dict):
    bc_pb_hp_dict = defaultdict(list)
    qname_pb_hp_dict = defaultdict(list)
    for key, mole_num_list in phased_dict_mole_num.items():
        for mole_num in mole_num_list:
            qname_list = mole_qname_dict[mole_num]
            for qname in qname_list:
                qname_pb_hp_dict[qname].append(key)
                bc = qname_bc_dict[qname]
                bc_pb_hp_dict[bc].append(key)
    return bc_pb_hp_dict, qname_pb_hp_dict





def write_fastqs(bc_pb_hp_dict, output_dir):

    for bc,hp_list in tqdm(bc_pb_hp_dict.items(), desc = 'assign barcode'):
        if bc!='0_0_0':
            for key in hp_list:
                _PS_flag = key[0]
                _HP_flag = key[1]
                outfile = output_dir +"/fastq_by_" + str(_PS_flag)  + "_" + _HP_flag + ".fastq"
                extract_one_barcode(bc, outfile, write_mode = 'a')
                


def Extract_start(output_dir,phased_h5_file,mole_qname_dict_file,bamfile, chr_num, num_threads):
    phased_dict_mole_num = Extract_mole_num_from_phased_mole(phased_h5_file,chr_num)
    mole_qname_dict = pickle.load(open(mole_qname_dict_file,"rb"))
    qname_bc_dict = load_qname_bc_dict(bamfile)
    bc_pb_hp_dict, qname_pb_hp_dict = combine_dict(phased_dict_mole_num,mole_qname_dict,qname_bc_dict)


    if os.path.exists(output_dir):
        os.system("rm -r " + output_dir)

    os.system("mkdir -p " + output_dir)

    extract_bad_reads(bam_file, output_dir, num_threads, qname_pb_hp_dict )
    write_fastqs(bc_pb_hp_dict, output_dir)







work_dir = "/data/maiziezhou_lab/CanLuo/AquilaSV_Project/stLFR_no_OER/stLFR_chr22_add_OER_result/chr22_28485965_HG3_PB_HySA_29009"
chr_num = 22
num_threads = 22
sample_name = 'target'
phased_file_dir = work_dir + "/results_phased_probmodel/"
phased_h5_file = phased_file_dir + "chr" + str(chr_num) + ".phased_final"
h5_folder = work_dir+'/H5_for_molecules/'
mole_qname_dict_file = h5_folder + sample_name + "_chr" + str(chr_num) +  "_qname.p"
bamfile = work_dir+'/selected.bam'


output_dir = work_dir + "/fastqs/"

print(output_dir)




Extract_start(output_dir,phased_h5_file,mole_qname_dict_file,bamfile, chr_num, num_threads)











