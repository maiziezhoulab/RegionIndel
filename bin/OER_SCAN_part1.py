from argparse import ArgumentParser
import pysam
from collections import Counter
import os
from subprocess import Popen
from tqdm import tqdm
from joblib import Parallel, delayed
import pickle
parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--input_path','-i')
parser.add_argument('--output_dir','-o')
parser.add_argument('--n_thread','-t', type = int, default = 23)
parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
input_path = args.input_path
output_dir = args.output_dir
n_thread = args.n_thread

import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
level=logging.INFO,
datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")

os.system("mkdir -p " + output_dir)

samfile = pysam.AlignmentFile(input_path)
dc_OER = {}
cnt = 0
retrieve_list = ["chr"+str(i) for i in range(1,23)]
retrieve_list.append("hs37d5")
print(retrieve_list)
dc_OER_seq = {}
for read in samfile:
	cnt +=1
	chr_name = read.reference_name
	mate_ref = read.next_reference_name
	if  mate_ref != chr_name:
		if mate_ref in retrieve_list:
			if read.is_read1:
				read_pair = 1 
			else:
				read_pair = 2
			seq = read.get_forward_sequence()
			# qual = read.qual
			# if read.is_reverse:
			# 	qual = qual[::-1]
			qname = read.qname
			mate_pos = read.next_reference_start
			#if qname in dc_OER:
			#	old_mate_ref = dc_OER[qname]['mate_ref']
			#	if old_mate_ref != mate_ref:
			#		logger.info(qname+old_mate_ref+mate_ref)
			dc_OER[qname] = {'cur_ref': chr_name, 'cur_pos': read.pos,
			'mate_ref':mate_ref,'mate_pos':mate_pos,'mate_pair':3-read_pair}
			dc_OER_seq[qname] = (seq, read_pair)

logger.info("# total reads: %d"%cnt)
logger.info("# total OER: %d"%(len(dc_OER)))

## collect all other chr

b = [dc['mate_ref']  for name,dc in dc_OER.items()]

print(Counter(b))


dc_by_chr = {}
for qname,dc in tqdm(dc_OER.items()):
	mate_ref = dc['mate_ref']
	mate_pos = dc['mate_pos']
	# mate_pair = dc['mate_pair']
	if mate_ref in dc_by_chr:
		dc_by_chr[mate_ref].append((mate_pos,qname))
	else:
		dc_by_chr[mate_ref]=[(mate_pos,qname)]


for mate_ref in retrieve_list:
	if mate_ref in dc_by_chr:
		logger.info("write "+mate_ref)
		qlist = dc_by_chr[mate_ref]
		dc = {}
		logger.info("qnames: %d"%(len(qlist)))
		for mate_pos,qname in qlist:
			if mate_pos in dc:
				dc[mate_pos].append(qname)
			else:
				dc[mate_pos] = [qname]
		logger.info("poss: %d"%(len(dc)))
		out_path = output_dir+'/OER_%s_mate%s.p'%(chr_name,mate_ref)
		pickle.dump(dc, open(out_path,'wb'))




















