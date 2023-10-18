from argparse import ArgumentParser
import pysam
from collections import Counter
import os
from subprocess import Popen
from tqdm import tqdm
from joblib import Parallel, delayed
import pickle

parser = ArgumentParser(description="",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--chr_name','-chrn')
parser.add_argument('--input_path','-i')
parser.add_argument('--output_dir','-o')
parser.add_argument('--rq_dir','-rq', help = "part1 dir")
parser.add_argument('--n_thread','-t', type = int, default = 20)
parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
chr_name = args.chr_name
input_path = args.input_path
rq_dir = args.rq_dir
output_dir = args.output_dir
n_thread = args.n_thread

import logging
## set logger
logging.basicConfig(
format='%(asctime)s %(levelname)-8s %(message)s',
level=logging.INFO,
datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(" ")


## collect all request card
dc_all = {}
dc_all_name = {}
#rq_dir = "/data/maiziezhou_lab/CanLuo/Link_Reads_Project/oer_test3/"
card_list = []
request_chr_name_list = []
for i in range(1,23):
	request_chr_name = 'chr'+str(i)
	logger.info("load request from "+request_chr_name+'...')
	card_path = rq_dir + "/OER_%s_mate%s.p"%(request_chr_name,chr_name)
	if not os.path.exists(card_path):
		logger.info(card_path + ' does not exist')
		continue
	else:
		dc = pickle.load(open(card_path,'rb'))
		card_list.append(dc)
		request_chr_name_list.append(request_chr_name)
		for pos,qnames in tqdm(dc.items()):
			# qname_chr_dc = {}
			# for q in qnames:
			# 	qname_chr_dc[q] = request_chr_name

			if pos in dc_all:
				dc_all[pos][request_chr_name] = qnames
				dc_all_name[pos].extend(qnames)
			else:
				dc_all[pos] = {request_chr_name:qnames}
				dc_all_name[pos] = qnames

for pos in dc_all_name:
	dc_all_name[pos] = list(set(dc_all_name[pos]))


for pos in dc_all:
	dc1 = {}
	dc = dc_all[pos]
	for chr_name,qnames in dc.items():
		for q in qnames:
			if q in dc1:
				dc1[q].append(chr_name)
			else:
				dc1[q] = [chr_name]
	dc_all[pos] = dc1

pickle.dump(dc_all_name, open(rq_dir+'/request_to_%s.p'%chr_name,'wb'))
pickle.dump(dc_all, open(rq_dir+'/request_to_%s_with_from_info.p'%chr_name,'wb'))
logger.info("num of pos: %d"%(len(dc_all)))
pos_list = sorted(list(dc_all.keys()))
##### loop through samfile





os.system("mkdir -p " + output_dir)

samfile = pysam.AlignmentFile(input_path)
dc_OER = {}
cnt = 0
mate_list = set(["chr"+str(i) for i in range(1,23)])

for mate in mate_list:
	dc_OER[mate] = {}
print(mate_list)
pos_pointer = 0
dc_mate = {}
logger.info("start pos: %d end pos: %d"%(pos_list[0],pos_list[-1]))
for read in samfile.fetch(args.chr_name,pos_list[0],pos_list[-1]+1):
	cnt +=1
	pos = read.pos 

	if cnt%10000000==0:
		logger.info(str(cnt))

	comp_pos = pos_list[pos_pointer]
	# if pos<comp_pos:
	# 	continue
	if pos==comp_pos:
		q = read.qname
		qnames = dc_all_name[pos]
		if q in qnames:
			seq = read.get_forward_sequence()
			if pos in dc_mate:
				if q not in dc_mate[pos]:
					dc_mate[pos][q] = seq
			else:
				dc_mate[pos] = { q:seq}
	elif (pos_pointer+1)<= len(pos_list)-1:
		if pos == pos_list[pos_pointer+1]:
			q = read.qname
			qnames = dc_all_name[pos]
			if q in qnames:
				seq = read.get_forward_sequence()
				if pos in dc_mate:
					if q not in dc_mate[pos]:
						dc_mate[pos][q] = seq
					# else:
					# 	old_seq = dc_mate[pos][q]
					# 	new_seq = seq
					# 	if len(new_seq)>len(old_seq):
					# 		dc_mate[pos][q] = new_seq
				else:
					dc_mate[pos] = {q:seq}
			pos_pointer+=1



logger.info("# find pos: %d"%len(dc_mate))
pickle.dump(dc_mate,open(output_dir+'/feedback_from_%s.p'%(args.chr_name),'wb'))


dc_by_chr = {}

for pos in tqdm(dc_mate,desc="split result back to chromosomes"):
	dc = dc_mate[pos]
	find_qnames = list(dc.keys())
	requested_qnames = dc_all_name[pos]
	if len(find_qnames)<len(requested_qnames):
		logger.info(set(requested_qnames)-set(find_qnames))
	for q in find_qnames:
		seq = dc[q]
		chrs = dc_all[pos][q]
		for chr_name in chrs:
			if chr_name in dc_by_chr:
				dc_by_chr[chr_name][q] = seq 
			else:
				dc_by_chr[chr_name] = {q:seq}

#### 
logger.info("write feedback by chromsome")

for chr_name in tqdm(dc_by_chr):
	dc = dc_by_chr[chr_name]
	out_path = output_dir+'/feedback_from_%s_to_%s.p'%(args.chr_name,chr_name)
	pickle.dump(dc, open(out_path,'wb'))
