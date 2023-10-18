#!/usr/bin/env python
import pdb
#pdb.set_trace()
from subprocess import Popen
from argparse import ArgumentParser
import os
import sys
import pickle
import os.path
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 
__author__ = "Maizie&Can@Vandy"
parser = ArgumentParser(description="Author: maiziezhoulab@gmail.com\n",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--chr_num','-chr',type=int,help="chromosome number for target variant or region", required=True)
parser.add_argument('--out_dir','-o', help="Directory to store assembly results, default = ./AquilaSV_results",default="./AquilaSV_results")
#parser.add_argument('--reference','-ref', help="reference fasta file, run ./install to download it",required=True)
parser.add_argument('--num_threads','-t',type=int,help="number of threads, default = 10, this correponds to number of small files get assembled simulateoulsy", default=10)
parser.add_argument('--num_threads_spades','-t_spades',type=int,help="number of threads for spades, default = 5", default=5)
parser.add_argument('--block_len_use','-bl',type=int,help="phase block len threshold, default = 100000",default=100000)

args = parser.parse_args()

def read_ref(fasta_file,chr_num,out_dir):
    f = open(fasta_file,"r")
    count = 0
    ref_seq = ""
    for line in f:
        if count > 0:
            data = line.rsplit()
            ref_seq += data[0]
        count += 1
    print("total_len for chr" + str(chr_num))
    print(len(ref_seq))
    pickle.dump(ref_seq, open(out_dir + "ref_seq_chr" + str(chr_num) +  ".p","wb"))


def extract_ref_chr(ref_file,chr_num,out_dir):
    fw = open(out_dir + "genome_ref_chr" + str(chr_num) + ".fasta","w")
    f = open(ref_file,"r")
    flag = 0
    total_len = 0
    for line in f:
        data = line.rsplit()
        if data[0] == ">chr" + str(chr_num):
            fw.writelines(">" + str(chr_num) + "\n")
            flag = 1
        elif data[0][0] == ">" and flag == 1:
            break
        else:
            if flag == 1:
                total_len += len(data[0])
                fw.writelines(data[0] + "\n")
    print("chr" + str(chr_num) + ":")
    print(total_len)


def local_assembly_for_small_chunks(chr_start,chr_end,num_threads,num_threads_spades,Local_Assembly_dir,Assembly_Contigs_dir,assembler):
    if assembler=="spades":
        use_cmd = "python3 " + code_path + "Run_spades_final_MT_2_all_noec_deltemp.py" + " --num_threads " + str(num_threads) + " --num_threads_spades " + str(num_threads_spades) + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end) + " --out_dir " + Local_Assembly_dir + " --minicontig_dir " + Assembly_Contigs_dir
        Popen(use_cmd,shell=True).wait()
    elif assembler=='velvet':
        use_cmd = "python " + code_path + "Run_velvet_final_MT_2_all_noec_deltemp.py" + " --num_threads " + str(num_threads) + " --num_threads_spades " + str(num_threads_spades) + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end) + " --out_dir " + Local_Assembly_dir + " --minicontig_dir " + Assembly_Contigs_dir
        Popen(use_cmd,shell=True).wait()



def Complete_contiguity(chr_start,chr_end,Assembly_Contigs_dir,phase_blocks_cut_highconf_dir,cut_threshold,ref_dir,ref_idx_file):
    use_cmd = "python " + code_path + "Concatenate_contigs_from_microcontigs.py" + " --chr_start "  + str(chr_start) + " --chr_end " + str(chr_end) + " --out_dir " + Assembly_Contigs_dir + " --phase_cut_folder " + phase_blocks_cut_highconf_dir + " --cut_threshold " + str(cut_threshold) + " --ref_file " + ref_idx_file + " --ref_dir " + ref_dir
    Popen(use_cmd,shell=True).wait()
  

def fastq_to_bam(fasta_path,bam_path,ref_path):
    file_name = bam_path.split('/')[-1].split('.')[0]
    file_dir = os.dirname(bam_path)
    cmd1 = "minimap2 -a %s %s %s/%s.sam;"%(ref_path,fasta_path,file_dir,file_name)
    cmd2 = "samtools view -Sb %s/%s.sam | samtools sort > %s;"%(file_dir,file_name,bam_path)
    cmd3 = "samtools index %s;"%(bam_path)
    cmd4 = "rm %s/%s.sam;"%(file_dir,file_name)
    Popen(cmd1+cmd2+cmd3+cmd4,shell=True).wait()

#if __name__ == "__main__":
def main():
    if len(sys.argv) == 1:
        Popen("python " + "AquilaSV_step2.py -h",shell=True).wait()
    else:
        chr_start = args.chr_num
        chr_end = args.chr_num
        #ref_file = args.reference
        cut_threshold = args.block_len_use
        num_threads = int(args.num_threads)
        num_threads_spades = int(args.num_threads_spades)
        #ref_dir = args.out_dir + "/ref_dir/"
        #if os.path.exists(ref_dir):
        #    print("using existing output folder: " + ref_dir)
        #else:
        #    os.makedirs(ref_dir)
        #if ~os.path.exists(ref_dir + "ref.mmi"):
        #    try:
        #        mk_ref_idx = "minimap2 -d " + ref_dir + "ref.mmi "  + ref_file
        #    except:
        #        mk_ref_idx = code_path + "minimap2/" + "minimap2 -d " + ref_dir + "ref.mmi "  + ref_file

        #    Popen(mk_ref_idx,shell=True).wait()
        #for chr_num in range(chr_start,chr_end + 1):
        #    extract_ref_chr(ref_file,chr_num,ref_dir)
        #extract_ref_chr(ref_file,"X",ref_dir)
        #for chr_num in range(chr_start,chr_end + 1):
        #    read_ref(ref_dir + "genome_ref_chr" + str(chr_num) + ".fasta",chr_num,ref_dir)
        #read_ref(ref_dir + "genome_ref_chrX.fasta",23,ref_dir)

        #ref_idx_file = ref_dir + "ref.mmi"
        Assembly_Contigs_dir = args.out_dir + "/Assembly_Contigs_files/"
        phase_blocks_cut_highconf_dir = args.out_dir + "/phase_blocks_cut_highconf/"
        Local_Assembly_dir = args.out_dir + "/Local_Assembly_by_chunks/"
        #assembler=args.assembler
        assembler='spades'


        local_assembly_for_small_chunks(chr_start,chr_end,num_threads,num_threads_spades,Local_Assembly_dir,Assembly_Contigs_dir,assembler)

       # cmd1 = "minimap2 -a  %s %sAquila_cutPBHC_minicontig_chr%d.fasta > %sAquila_cutPBHC_minicontig_chr%d.sam;"%(ref_idx_file,Assembly_Contigs_dir,chr_start,Assembly_Contigs_dir,chr_start)
       # cmd2 = "samtools view -Sb %sAquila_cutPBHC_minicontig_chr%d.sam | samtools sort > %sAquila_Contig_chr%d.bam;"%(Assembly_Contigs_dir,chr_start,Assembly_Contigs_dir,chr_start)
       # cmd3 = "samtools index %sAquila_Contig_chr%d.bam;rm %sAquila_cutPBHC_minicontig_chr%d.sam;"%(Assembly_Contigs_dir,chr_start,Assembly_Contigs_dir,chr_start)
       # with open("test.sh",'w') as f:
       #     f.write(cmd1+cmd2+cmd3)
       # Popen(cmd1+cmd2+cmd3,shell=True).wait()
       # del_cmd = "rm -r -v %schr*_files_cutPBHC/*_spades_assembly"%(Local_Assembly_dir)
       # Popen(del_cmd,shell=True).wait()
       # ###Complete_contiguity(chr_start,chr_end,Assembly_Contigs_dir,phase_blocks_cut_highconf_dir,cut_threshold,ref_dir,ref_idx_file)

if __name__ =="__main__":
    main()
    
