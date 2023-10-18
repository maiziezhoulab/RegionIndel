from subprocess import Popen
from argparse import ArgumentParser
import os
import json
code_path = os.path.dirname(os.path.abspath( __file__ ))+'/'


parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--chr_num','-chr',type = int)
parser.add_argument('--position','-pos',type = int)
parser.add_argument('--output_dir','-o')
parser.add_argument('--config_path','-conf')
parser.add_argument('--flanking','-fl',type = int, default = 25000)
parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
chr_num = args.chr_num
pos = args.position
output_dir = args.output_dir
flanking = args.flanking
config_path = args.config_path

# config_dc = json.load(open(code_path+'/config.txt'))
config_dc = json.load(open(config_path))


import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")

# chr_num = 21
# pos = 34144186
# flanking = 50000

os.system("mkdir -p "+output_dir)
start = pos-flanking
end = pos+ flanking


bam_file = config_dc['bamfile'].replace("*",str(chr_num))
vcf_file = config_dc['vcf']
#Aquila_dir = "/data/maiziezhou_lab/CanLuo/Software/AquilaSV/bin_v1/"
Aquila_dir = code_path
selected_bam_file = output_dir + '/selected.bam'
reference_file = config_dc['reference'].replace('*',str(chr_num))

cmd = "samtools view -b %s chr%d:%d-%d > %s; samtools index %s"%(
	bam_file,
	chr_num,
	start,
	end,
	selected_bam_file,
	selected_bam_file)
Popen(cmd, shell = True).wait()

cmd = "python3 "+Aquila_dir+"/AquilaSV_step1.py \
--bam_file %s \
--vcf_file %s \
--chr_num %d \
--out_dir %s -conf %s"%(
	selected_bam_file,
	vcf_file,
	chr_num,
	output_dir,config_path)
Popen(cmd, shell = True).wait()



cmd = "python3 "+Aquila_dir+"/AquilaSV_step2.py --out_dir %s --chr_num %d "%(
	output_dir,
	chr_num,)
Popen(cmd, shell = True).wait()

cmd = "python3 "+Aquila_dir+"/AquilaSV_step3.py  --assembly_dir %s  --ref_file %s --chr_num %d -o %s -d"%(
	output_dir,
	reference_file,
	chr_num, output_dir)
Popen(cmd, shell = True).wait()


