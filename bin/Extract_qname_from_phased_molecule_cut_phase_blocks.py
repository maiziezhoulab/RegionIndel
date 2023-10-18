#### modify for no cut
import pdb
#pdb.set_trace()
from collections import defaultdict
import pickle
import os
from argparse import ArgumentParser

def save_pickle_file(dict1,fname):
    for value in dict1:
        dict1[value] = dict(dict1[value])
    my_dict = dict(dict1)
    with open(fname, "wb") as f:
        pickle.dump(my_dict,f) 


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

    #pickle.dump(phase_dict,open("phase_dict_chr" + str(chr_num) + ".p","wb"))
    print("done~")
    return phase_dict


# def Extract_qname_old(phased_dict_mole_num,mole_qname_dict,barcoded_fastq_file,chr_num,output_dir):
#     phased_dict_qname_2 = defaultdict(lambda: defaultdict(int))
#     qname_total_dict = defaultdict(int)
#     curr = 0
#     for key, mole_num_list in phased_dict_mole_num.items():
#         curr += 1
#         for mole_num in mole_num_list:
#             qname_list = mole_qname_dict[mole_num]
#             for qname in qname_list:
#                 phased_dict_qname_2[qname][key] = 1
#     #save_pickle_file(phased_dict_qname_2,"phased_dict_qname_2.p")
    
#     f = open(barcoded_fastq_file,"r")
#     count = 0
#     flag = 0
#     for line in f:
#         data = line.rsplit()
#         if count%8 == 0:
#             qname_curr = data[0][1:]
#             if phased_dict_qname_2[qname_curr] != {}:
#                 flag = 1
#                 PS_flag_info = phased_dict_qname_2[qname_curr]
#                 for _PS_HP, _val in PS_flag_info.items():
#                     _PS_flag = _PS_HP[0]
#                     _HP_flag = _PS_HP[1]
#                 file_curr = output_dir +"fastq_by_" + str(_PS_flag)  + "_" + _HP_flag + ".fastq"
#                 if os.path.isfile(file_curr):
#                     fw_curr = open(file_curr,"a")
#                     fw_curr.write(line)
#                 else:
#                     fw_curr = open(file_curr,"w")
#                     fw_curr.write(line)
#                 del phased_dict_qname_2[qname_curr]

#             else:
#                 del phased_dict_qname_2[qname_curr]

#         elif count%8 == 7:
#             if flag == 1:
#                 flag = 0
#                 fw_curr.write(line)
#                 fw_curr.close()

#         else:
#             if flag == 1:
#                 fw_curr.write(line)

#         count += 1

#     print("finished extracting...")

def Extract_qname(phased_dict_mole_num,mole_qname_dict,barcoded_fastq_file,chr_num,output_dir):
    phased_dict_qname_2 = defaultdict(lambda: defaultdict(int))
    qname_total_dict = defaultdict(int)
    curr = 0
    for key, mole_num_list in phased_dict_mole_num.items():
        curr += 1
        for mole_num in mole_num_list:
            qname_list = mole_qname_dict[mole_num]
            for qname in qname_list:
                phased_dict_qname_2[qname][key] = 1

    with open(barcoded_fastq_file,"r") as f:
        fastq_lines = f.readlines()

    block_list = [''.join(fastq_lines[i:i+4]) for i in range(0,len(fastq_lines),4)]
    name_list = [fastq_lines[i][1:-1] for i in range(0,len(fastq_lines),4)]
    print(name_list[:10])
    print(len(name_list))
    for i in range(len(name_list)):
        qname_curr = name_list[i]
        line = block_list [i]
        if phased_dict_qname_2[qname_curr] != {}:
            flag = 1
            PS_flag_info = phased_dict_qname_2[qname_curr]
            for _PS_HP, _val in PS_flag_info.items():
                _PS_flag = _PS_HP[0]
                _HP_flag = _PS_HP[1]
            file_curr = output_dir +"fastq_by_" + str(_PS_flag)  + "_" + _HP_flag + ".fastq"
            if os.path.isfile(file_curr):
                fw_curr = open(file_curr,"a")
                fw_curr.write(line)
            else:
                fw_curr = open(file_curr,"w")
                fw_curr.write(line)
            fw_curr.close()
            # del phased_dict_qname_2[qname_curr]
    print(output_dir)
    fastq_list = [output_dir+path for path in os.listdir(output_dir) if '.fastq' in path]
    
    #print('hey',fastq_list)
    #for fastq_path in fastq_list:
    #    print('hahaha')
    
    # for fastq_path in fastq_list:
    #     print(fastq_path)
    #     with open(fastq_path,'r') as f:
    #         s = f.readlines()
    #     print(s[:10])
    #     name_list = [s[i] for i in range(0,len(s),4)]
    #     block_list = [''.join(s[i:i+4]) for i in range(0,len(s),4)]
    #     idx = np.argsort(name_list)
    #     ct = np.unique(name_list,return_counts = True)
    #     unpaired_name = set(ct[0][ct[1]!=2])
    #     print('up:',unpaired_name)
    #     block_list = np.array(block_list)[idx]
    #     name_list = np.array(name_list)[idx]
    #     block_list = [ block_list[i] for i in range(len(block_list)) if name_list[i] not in unpaired_name]
    #     print(len(block_list))
    #     with open(fastq_path,'w') as f:
    #         f.writelines(block_list)
        
    
    print("finished extracting...")

def Extract_start(output_dir,chr_num,phased_h5_file,mole_qname_dict_file,chr_fastq,xin):
    phased_dict_mole_num = Extract_mole_num_from_phased_mole(phased_h5_file,chr_num)
    mole_qname_dict = pickle.load(open(mole_qname_dict_file,"rb"))
    Extract_qname(phased_dict_mole_num,mole_qname_dict,chr_fastq,chr_num,output_dir)
