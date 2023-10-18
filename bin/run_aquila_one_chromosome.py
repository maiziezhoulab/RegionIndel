'''

config file


{
"spades":"/data/maiziezhou_lab/CanLuo/Software/AquilaSV/bin/SPAdes-3.13.0-Linux/",
"k8":"/data/maiziezhou_lab/CanLuo/Software/AquilaSV/bin/k8-0.2.4/",
"paftools":"/data/maiziezhou_lab/CanLuo/Software/AquilaSV/bin/paftools/",
"header":"/data/maiziezhou_lab/CanLuo/Software/AquilaSV/bin/header/",
"feedback":"/data/maiziezhou_lab/CanLuo/Link_Reads_Project/oer_test3/feedback_v3/for_download/",
"bamfile":"/data/maiziezhou_lab/Datasets/L5_NA24385_10x/hg19_longranger_align/by_chr/possorted_chr*.bam",
"vcf":"/data/maiziezhou_lab/Datasets/L5_NA24385_10x/Freebayes_results_hg19/L5_hg19_ref_freebayes.vcf",
"reference":"/data/maiziezhou_lab/CanLuo/long_reads_project/DipPAV2/hg19_ref_by_chr/hg19_chr*.fa",
"bench":"/data/maiziezhou_lab/CanLuo/Link_Reads_Project/bench.vcf"
}



'''


from argparse import ArgumentParser
import os 
from subprocess import Popen
from joblib import Parallel, delayed
from tqdm import tqdm
import json
code_path = os.path.dirname(os.path.abspath( __file__ ))+'/'

parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--chr_num','-chr',type = int)
parser.add_argument('--output_dir','-o')
parser.add_argument('--flanking','-fl',type = int, default = 25000)
parser.add_argument('--config_path','-conf')
parser.add_argument('--delete_temp_file','-d', action='store_true')
parser.add_argument('--n_thread','-t', type = int, default = 10)
args = parser.parse_args()
chr_num = args.chr_num
flanking = args.flanking
output_dir = args.output_dir
n_thread = args.n_thread
config_path = args.config_path

import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")

os.system("mkdir -p "+output_dir)
# config_dc = json.load(open(code_path+'/config.txt'))
config_dc = json.load(open(config_path))
with open(config_dc['bench'],'r') as f:
    s = f.readlines()
bam_file = config_dc['bamfile'].replace("*",str(chr_num))
if not os.path.exists(bam_file+'.bai'):
    logger.info("create index for "+bam_file)
    os.system("samtools index "+bam_file)
chr_name = 'chr'+str(chr_num)
lines = []
for line in s:
    chr_num_line = line.split()[0]
    if (chr_num_line==chr_name): # and ('SVTYPE=INS' in line):
        lines.append(line.split()[:3])

logger.info("chr%d: %d target regions"%(chr_num, len(lines)))


log_file = output_dir+'/log.txt'
# load already run
if os.path.exists(log_file):
    with open(log_file,'r') as f:
        finished_regions = f.read().split('\n')[:-1]
else:
    finished_regions = []
    os.system("touch "+log_file)

lines_not_finished = [ line for line in lines if ("_".join(line) not in finished_regions) ]

logger.info("not finished tasks:"+str(lines_not_finished))
def sv_one_region(line):
    sv_dir = output_dir+'/'+"_".join(line)+'/'
    pos = int(line[1])
    logger.info("Region : %s"%(' '.join(line) ))
    cmd = "python3 "+code_path+"//run_aquila_one_region.py -chr %d -pos %d -o %s -fl %d -conf %s"%(
        chr_num,
        pos,
        sv_dir,
        flanking,
        config_path)
    logger.info(cmd)
    Popen(cmd,shell = True).wait()
    os.system("echo %s >> %s"%("_".join(line),log_file))
    return 

results = Parallel(n_jobs=n_thread)(delayed(sv_one_region)(line) for line in tqdm(lines_not_finished, desc = "variant call progress"))












