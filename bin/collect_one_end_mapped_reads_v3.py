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
parser.add_argument('--OER_dir','-oer')
parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
input_path = args.input_path
output_dir = args.output_dir
n_thread = args.n_thread
OER_dir = args.OER_dir

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
			qual = read.qual
			if read.is_reverse:
				qual = qual[::-1]
			qname = read.qname
			mate_pos = read.next_reference_start
			dc_OER[qname] = {'cur_ref': chr_name, 'cur_pos': read.pos,
			'mate_ref':mate_ref,'mate_pos':mate_pos}
			dc_OER_seq[qname] = (seq, read_pair)

logger.info("# total reads: %d"%cnt)
logger.info("# total OER: %d"%(len(dc_OER)))

## collect all other chr

b = [dc['mate_ref']  for name,dc in dc_OER.items()]

print(Counter(b))
dc_chr = {}




def search_one_chr(mate_ref,qnames,chr_name,OER_dir):
	dc_mate_one_chr = {}
	logger.info("search "+ mate_ref)
	qnames = set(qnames)
	pickle_path = OER_dir + "/feedback_from_%s_to_%s.p"%(mate_ref, chr_name)
	dc_mate = pickle.load(open(pickle_path,'rb'))
	for qname in qnames:
		if qname in dc_mate:
			dc_mate_one_chr[qname] = dc_mate[qname]
		else:
			print(qname)
	return dc_mate_one_chr



dc_by_chr = {}
for qname,dc in tqdm(dc_OER.items()):
	mate_ref = dc['mate_ref']
	if mate_ref in dc_by_chr:
		dc_by_chr[mate_ref].append(qname)
	else:
		dc_by_chr[mate_ref]=[qname]




sequences = Parallel(n_jobs=n_thread)(delayed(search_one_chr )(mate_ref, dc_by_chr[mate_ref],chr_name,OER_dir) for mate_ref in tqdm(dc_by_chr))



dc_mate_seq = {}
for dc in sequences:
	dc_mate_seq.update(dc)

logger.info("raw OER: %d"%(len(dc_OER)))
logger.info("retrieve mates: %d"%len(dc_mate_seq))

out_path = output_dir+'/OER.fa'
with open(out_path,'w') as f:
	for qname,seq2 in dc_mate_seq.items():
		seq1, read_pair1 = dc_OER_seq[qname]
		
		if read_pair1 == 1:
			f.write(">"+qname+'\n'+seq1+'\n')
			f.write(">"+qname+'\n'+seq2+'\n')
		else:
			f.write(">"+qname+'\n'+seq2+'\n')
			f.write(">"+qname+'\n'+seq1+'\n')
















