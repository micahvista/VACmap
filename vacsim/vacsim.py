import numpy as np
from Bio.Seq import Seq
from heapdict import heapdict
import mappy as mp
import pandas as pd
import sys
import pysam

def chioce_contig(contig_list, contigprob_list):
    return np.random.choice(contig_list, p = contigprob_list)
def prepare_usable_loc_info(contig2seq):
    contig_list = []
    contigprob_list = []
    contig2usable_interval = dict()
    allsize = 0
    prestart = 0
    for contig in contig2seq:
        hd = heapdict()
        iloc = -1

        Nsize = 0
        good = False
        print(contig)
        for c in contig2seq[contig]:
            iloc += 1
            if(c not in 'ATGC'):
                Nsize += 1
                if(good == True):
                    print(prestart, iloc)
                    if(iloc - prestart > 200):
                        hd[(prestart, iloc)] = 1/(iloc - prestart)
                    good = False

            else:
                if(good == False):
                    prestart = iloc
                    good = True

        if(good == True):
            print(prestart, len(contig2seq[contig]))
            hd[(prestart, len(contig2seq[contig]))] = 1/(len(contig2seq[contig]) - prestart)
        contig_list.append(contig)
        contigprob_list.append(len(contig2seq[contig])-Nsize)
        allsize += contigprob_list[-1]
        contig2usable_interval[contig] = hd
    contigprob_list = np.array(contigprob_list)/allsize   
    return contig2usable_interval, contig_list, contigprob_list
def repalce(from_seq, from_start, from_end, reverse, random):#return replaced string, include smaller, larger are exclude
    acid = ['A', 'T', 'G', 'C']
    if(from_seq == ''):#remove this subseq
        return ''
    if(random == True):
        return ''.join([np.random.choice(acid) for i in range(from_end - from_start)])
    if(reverse == 0):
        return from_seq[from_start: from_end]
    else:
        return str(Seq(from_seq[from_start: from_end]).reverse_complement())
def exchange(seq_a, start_a, end_a, reverse_a, seq_b, start_b, end_b, reverse_b):
    if(reverse_a == 0):
        a_pos_seq = seq_b[start_b: end_b]
    else:
        a_pos_seq = str(Seq(seq_b[start_b: end_b]).reverse_complement())
    if(reverse_b == 0):
        b_pos_seq = seq_a[start_a: end_a]
    else:
        b_pos_seq = str(Seq(seq_a[start_a: end_a]).reverse_complement())
    return a_pos_seq, b_pos_seq
def insert(from_seq, from_start, from_end, reverse, times):
    acid = ['A', 'T', 'G', 'C']
    if(from_seq == ''):#remove this subseq
        return ''.join([np.random.choice(acid) for i in range(from_end - from_start)])

    if(reverse == 0):
        return from_seq[from_start: from_end]*times
    else:
        return str(Seq(from_seq[from_start: from_end]).reverse_complement())*times     
        
def decode_sim_sv_info(sim_sv_info):
    contig_1, contig_2 = 'unkown', 'unkown'
    match = 200
    svstart_1, svstart_2 = match, match
    times = int(sim_sv_info.split(',')[-1])
    sim_sv_info = sim_sv_info.split(',')[:-1]
    sv_list = []
    for i in range(times):
        for op in sim_sv_info:
            svtype = op.split(':')[0]
            svlen = np.random.randint(int(op.split(':')[1]), int(op.split(':')[2]))
            if(svtype == 'DEL'):#DEL:contig:st:end
                sv_list.append([svtype, contig_1, svstart_1, svstart_1 + svlen])
                svstart_1 += svlen
            elif(svtype == 'INS'):#INS:contig:st:size
                sv_list.append([svtype, contig_1, svstart_1, svlen])
            elif(svtype == 'DUP'):#DUP:contig:st:end:reverse:times
                reverse = int(op.split(':')[3])
                dup_times = int(op.split(':')[4])
                sv_list.append([svtype, contig_1, svstart_1, svstart_1 + svlen, reverse, dup_times])
                svstart_1 += svlen
            elif(svtype == 'INV'):#INV:contig:st:end
                sv_list.append([svtype, contig_1, svstart_1, svstart_1 + svlen])
                svstart_1 += svlen
            elif('TRA'):#TRA:contig:st:end:contig:st:end:reverse
                reverse = int(op.split(':')[3])
                sv_list.append([svtype, contig_1, svstart_1, svstart_1 + svlen, contig_2, svstart_2, svstart_2 + svlen, reverse])
                svstart_1 += svlen
                svstart_2 += svlen + match
            svstart_1 += match
    return svstart_1, svstart_2, sv_list
def add_SV(contig_1, contig_2, refstart_1, refstart_2, oneinfo, contig2seq, unique_id):
    if(oneinfo[0] == 'DEL'):
        return [[contig_1, refstart_1+oneinfo[2], refstart_1+oneinfo[3], '', 'DEL', unique_id, ['DEL', contig_1, refstart_1+oneinfo[2], refstart_1+oneinfo[3]]]]
    
    if(oneinfo[0] == 'INS'):
        alt_seq = insert('', refstart_1+oneinfo[2], refstart_1+oneinfo[2] + oneinfo[3], reverse = 0, times = 1)
        return [[contig_1, refstart_1+oneinfo[2], refstart_1+oneinfo[2], alt_seq, 'INS', unique_id, ['INS', contig_1, refstart_1+oneinfo[2], oneinfo[3]]]]
    
    if(oneinfo[0] == 'INV'):
        alt_seq = repalce(contig2seq[contig_1], refstart_1+oneinfo[2], refstart_1+oneinfo[3], reverse = 1, random = False)
        return [[contig_1, refstart_1+oneinfo[2], refstart_1+oneinfo[3], alt_seq, 'INV', unique_id, ['INV', contig_1, refstart_1+oneinfo[2],  refstart_1+oneinfo[3]]]]
    
    if(oneinfo[0] == 'DUP'):
        alt_seq = insert(contig2seq[contig_1], refstart_1+oneinfo[2], refstart_1+oneinfo[3], reverse = oneinfo[4], times = oneinfo[5])
        return [[contig_1, refstart_1+oneinfo[3], refstart_1+oneinfo[3], alt_seq, 'DUP', unique_id, ['DUP', contig_1, refstart_1+oneinfo[2], refstart_1+oneinfo[3], oneinfo[-2], oneinfo[-1]]]]
    
    if(oneinfo[0] == 'TRA'):
        a_pos_seq, b_pos_seq = exchange(contig2seq[contig_1], refstart_1+oneinfo[2], refstart_1+oneinfo[3], oneinfo[7], contig2seq[contig_2], refstart_2+oneinfo[5], refstart_2+oneinfo[6], oneinfo[7])
        tmpinfo = ['TRA', contig_1, refstart_1+oneinfo[2], refstart_1+oneinfo[3], contig_2, refstart_2+oneinfo[5], refstart_2+oneinfo[6], oneinfo[7]]
        return [[contig_1, refstart_1+oneinfo[2], refstart_1+oneinfo[3], a_pos_seq, 'TRA', unique_id, tmpinfo], [contig_2, refstart_2+oneinfo[5], refstart_2+oneinfo[6], b_pos_seq, 'TRA', unique_id, tmpinfo]]
def chioce_contig(contig_list, contigprob_list):
    return np.random.choice(contig_list, p = contigprob_list)
def get_start_loc(svspan_1, svspan_2, contig_list, contigprob_list, contig2usable_interval):
    edge_size = 200
    pass_flag = False
    count = 0
    while(True):
        count += 1
        if(count > 50):
            return pass_flag, '-1', -1, '-1', -1
        while(True):
            contig_1 = chioce_contig(contig_list, contigprob_list)
            contig_2 = chioce_contig(contig_list, contigprob_list)
            if(contig_1 != contig_2):
                break

        item = contig2usable_interval[contig_1].popitem()
        if(((item[0][1] - item[0][0]) - 2*edge_size) < svspan_1):
            contig2usable_interval[contig_1][item[0]] = item[1]
            continue
        if(item[0][0] + edge_size >= item[0][1] - svspan_1 - edge_size):
            continue
        refstart_1 = np.random.randint(item[0][0] + edge_size, item[0][1] - svspan_1 - edge_size)
        item_1 = item

        item = contig2usable_interval[contig_2].popitem()
        if(((item[0][1] - item[0][0]) - 100) < svspan_2):
            contig2usable_interval[contig_2][item[0]] = item[1]
            continue
        if(item[0][0] + edge_size >= item[0][1] - svspan_2 - edge_size):
            continue
        refstart_2 = np.random.randint(item[0][0] + edge_size, item[0][1] - svspan_2 - edge_size)
        item_2 = item
        pass_flag = True
        
        contig2usable_interval[contig_1][(item_1[0][0], refstart_1)] = 1 / (refstart_1 - item_1[0][0])
        contig2usable_interval[contig_1][(refstart_1 + svspan_1, item_1[0][1])] = 1 / (item_1[0][1] - refstart_1 - svspan_1)
        
        contig2usable_interval[contig_2][(item_2[0][0], refstart_2)] = 1 / (refstart_2 - item_2[0][0])
        contig2usable_interval[contig_2][(refstart_2 + svspan_2, item_2[0][1])] = 1 / (item_2[0][1] - refstart_2 - svspan_2)
        return pass_flag, contig_1, refstart_1, contig_2, refstart_2
        
def create_sim_sv_info_list(simulate_info_text):
    sim_sv_info_list = []
    for line in simulate_info_text.split('\n'):
        times = int(line.split('\t')[1])
        svstyle = line.split('\t')[0]
        for i in range(times):
            sim_sv_info_list.append(svstyle)
    return sim_sv_info_list
def create_sim_sv_info_list(simulate_info_text):
    sim_sv_info_list = []
    for line in simulate_info_text.split('\n'):
        if(line == ''):
            continue
        times = int(line.split('\t')[1])
        svstyle = line.split('\t')[0]
        for i in range(times):
            sim_sv_info_list.append(svstyle)
    return sim_sv_info_list
def getfirst_python(x):
    return x[0]
def getsecond_python(x):
    return x[1]
def random_create_sim_sv_info_list(eventlist = ['DEL', 'INS', 'INV', 'DUP'], sizerange = [100, 1000], eventmaxcount = 20, svcount = 20000):
    sim_sv_info_list = []
    lenstring = ':'+str(sizerange[0])+':'+str(sizerange[1])
    for line in range(svcount):
        eventcount = np.random.randint(1, max(eventmaxcount, 2))
        simulatedevent = 0
        svstyle = ''
        while(simulatedevent < eventcount):
            simevent = np.random.choice(eventlist)
            if(simevent in ('DEL', 'INS', 'INV')):
                svstyle += simevent + lenstring + ','
                simulatedevent += 1
            elif(simevent == 'DUP'):
                repeat = 1
                #repeat = np.random.randint(1, eventcount - simulatedevent + 1)
                simulatedevent += repeat
                svstyle += simevent+lenstring+':0:'+str(repeat)+','
            elif(simevent == 'TRA'):
                simulatedevent += 1
                reverse = str(np.random.randint(0, 1))
                svstyle += simevent + lenstring +':'+ reverse  + ','
                
        sim_sv_info_list.append(svstyle+'1')
    return sim_sv_info_list
def count_event(decoded_info):
    eventcount = 0
    for line in decoded_info[-1]:
        if(line[0] == 'DUP'):
            eventcount += line[-1]
        else:
            eventcount += 1
    return eventcount

def decode_sim_sv_info(sim_sv_info):
    contig_1, contig_2 = 'unkown', 'unkown'
    match = 0
    svstart_1, svstart_2 = match, match
    times = int(sim_sv_info.split(',')[-1])
    sim_sv_info = sim_sv_info.split(',')[:-1]
    sv_list = []
    preop = ''
    #print(sim_sv_info)
    for i in range(times):
        for op in sim_sv_info:
            svtype = op.split(':')[0]
            svlen = np.random.randint(int(op.split(':')[1]), int(op.split(':')[2]))
            if(svtype == 'DEL'):#DEL:contig:st:end
                sv_list.append([svtype, contig_1, svstart_1, svstart_1 + svlen])
                svstart_1 += svlen
            elif(svtype == 'INS'):#INS:contig:st:size
                sv_list.append([svtype, contig_1, svstart_1, svlen])
            elif(svtype == 'DUP'):#DUP:contig:st:end:reverse:times
                reverse = int(op.split(':')[3])
                dup_times = int(op.split(':')[4])
                sv_list.append([svtype, contig_1, svstart_1, svstart_1 + svlen, reverse, dup_times])
                svstart_1 += svlen
            elif(svtype == 'INV'):#INV:contig:st:end
                sv_list.append([svtype, contig_1, svstart_1, svstart_1 + svlen])
                svstart_1 += svlen
            elif(svtype == 'TRA'):#TRA:contig:st:end:contig:st:end:reverse
                reverse = int(op.split(':')[3])
                sv_list.append([svtype, contig_1, svstart_1, svstart_1 + svlen, contig_2, svstart_2, svstart_2 + svlen, reverse])
                svstart_1 += svlen
                svstart_2 += svlen + match
            elif(svtype == 'NML'):
                svstart_1 += svlen
                if(preop == 'TRA'):
                    svstart_2 += svlen
            preop = svtype    
            svstart_1 += match
    return svstart_1, svstart_2, sv_list
def get23_python(x):
    return x[1], x[2]
import pandas as pd
def tovcf(svlist, contig2length, outputpath):
    top = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">\n"""
    body = ''
    for contig in contig2length:
        body += "##contig=<ID="+contig+",length="+str(int(contig2length[contig]))+">\n"
    tail = """##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV:DEL=Deletion, INS=Insertion">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=RE,Number=1,Type=Integer,Description="Number of read support this record">
##INFO=<ID=bp,Number=1,Type=Integer,Description="Breakpoint">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t.\n"""

    myvcf = top+body+tail
    genomapper = {0:'0/0', 1:'0/1', 2:'1/1'}
    

    re = 999    
    for rec in pd.DataFrame(svlist).drop_duplicates().sort_values([0, 1]).values:


        contig = rec[0]
        start = int(rec[1])
        svlen = int(rec[2])
        bp = str(rec[3][0])
        uid = str(rec[3][1]) 
        code = int(rec[-1])
        geno = '.'
        
        if(code == 2):
            recinfo = 'SVLEN=' + str(svlen)+';SVTYPE=' + 'DEL'+';END='+str(start+svlen)+';bp='+bp+'\tGT\t'+geno+'\n'
            myvcf += (contig +'\t'+ str(start)+'\t'+ uid +'\t'+ '.'+'\t'+ '.'+'\t'+ str(re)+'\t'+ 'PASS'+'\t'+recinfo)
        elif(code == 1):
            recinfo = 'SVLEN=' + str(svlen)+';SVTYPE=' + 'INS'+';END='+str(start+1)+';bp='+bp+'\tGT\t'+geno+'\n'
            myvcf += (contig +'\t'+ str(start)+'\t'+ uid +'\t'+ '.'+'\t'+ '.'+'\t'+ str(re)+'\t'+ 'PASS'+'\t'+recinfo)
        elif(code == 3):
            recinfo = 'SVLEN=' + str(svlen)+';SVTYPE=' + 'INV'+';END='+str(start+svlen)+';bp='+bp+'\tGT\t'+geno+'\n'
            myvcf += (contig +'\t'+ str(start)+'\t'+ uid +'\t'+ '.'+'\t'+ '.'+'\t'+ str(re)+'\t'+ 'PASS'+'\t'+recinfo)
        elif(code == 4):
            recinfo = 'SVLEN=' + str(svlen)+';SVTYPE=' + 'DUP'+';END='+str(start+svlen)+';bp='+bp+'\tGT\t'+geno+'\n'
            myvcf += (contig +'\t'+ str(start)+'\t'+ uid +'\t'+ '.'+'\t'+ '.'+'\t'+ str(re)+'\t'+ 'PASS'+'\t'+recinfo)
        elif(code == 5):
            recinfo = 'SVLEN=' + str(999)+';SVTYPE=' + 'BND'+';END='+str(svlen)+';CHR2='+str(rec[4])+';bp='+bp+'\tGT\t'+geno+'\n'
            myvcf += (contig +'\t'+ str(start)+'\t'+ uid +'\t'+ '.'+'\t'+ '.'+'\t'+ str(re)+'\t'+ 'PASS'+'\t'+recinfo)


    with open(outputpath, "w") as f:
        f.write(myvcf)
        
def random_create_sim_sv_info_list(eventlist = ['DEL:100:200,NML:100:200', 'INS:100:200', 'INV:100:200', 'DUP:200:300', 'TRA:200:300,NML:100:200'], eventcountrange = [1, 20], svcount = 20000):
    sim_sv_info_list = [] 
    allowset = set(('INS', 'DUP'))
    for line in range(svcount):
        eventcount = np.random.randint(max(eventcountrange[0], 1), max(eventcountrange[1], 1))
        simulatedevent = 0
        svstyle = ''
        simevent = ''
        simSVevent = ''
        while(simulatedevent < eventcount):
            count = 0
            while(True):
                count += 1
                if(count > 100):
                    print('BAD ', eventlist)
                    return 
                simeventlist = np.random.choice(eventlist)
                if(simeventlist.split(',')[-1] == ''):
                    simeventlist.pop(-1)
                if(simeventlist.split(',')[-1].split(':')[0] == simevent and simevent in ('DEL', 'INS')):
                    continue
                if(simSVevent in ('DEL', 'INS', 'DUP')):
                    if(simevent == 'NML'):
                        if(set((simeventlist.split(',')[0].split(':')[0], simSVevent)) != allowset and simeventlist.split(',')[0].split(':')[0] in ('DEL', 'INS', 'DUP')):
                            if(simeventlist.split(',')[0].split(':')[0] != simSVevent):
                                continue
                    elif(simeventlist.split(',')[0].split(':')[0] == 'NML'):
                        if(len(simeventlist.split(',')) > 1 and simeventlist.split(',')[1] != ''):
                            if(set((simeventlist.split(',')[0].split(':')[0], simSVevent)) != allowset and simeventlist.split(',')[1].split(':')[0] in ('DEL', 'INS', 'DUP')):
                                if(simeventlist.split(',')[1].split(':')[0] != simSVevent):
                                    continue 
                break
            for simeventinfo in simeventlist.split(','):
                if(simeventinfo == ''):
                    continue
                simevent, sizestart, sizeend = tuple(simeventinfo.split(':'))
                lenstring = ':'+sizestart+':'+sizeend
                if(simevent in ('DEL', 'INS', 'INV', 'NML')):
                    svstyle += simevent + lenstring + ','
                elif(simevent == 'DUP'):
                    repeat = 1
                    #repeat = np.random.randint(1, eventcount - simulatedevent + 1)
                    simulatedevent += repeat
                    svstyle += simevent+lenstring+':0:'+str(repeat)+','
                elif(simevent == 'TRA'):
                    simulatedevent += 1
                    reverse = str(np.random.randint(0, 2))
                    svstyle += simevent + lenstring +':'+ reverse  + ','
                if(simevent != 'NML'):
                    simSVevent = simevent
                    simulatedevent += 1
        sim_sv_info_list.append(svstyle+'1')
    return sim_sv_info_list
def decode_parameterfile(s):
    simulate_info_text = ''
    decoded_sim_sv_info_list = []
    for line in s.split('\n'):
        line = (''.join([c for c in line if(c != ' ')]))
        if(line == ''):
            continue
        if(line.split('{')[0] not in ('Specified', 'Random')):
            continue
        if(line.split('{')[0] == 'Specified'):
            line = line[10:-1]
            number = (line.split(';')[-1].split('=')[1])
            simulate_info_text += (line.split(';')[0]+',1\t'+number+'\n')
        else:
            line = line[7:-1]
            parameterdict = dict()
            for item in line.split(';'):
                parameterdict[item.split('=')[0]] = eval(item.split('=')[1])
            
            sim_sv_info_list = random_create_sim_sv_info_list(eventlist = parameterdict['eventset'], eventcountrange = parameterdict['eventcount'], svcount = parameterdict['number'])
            for sim_sv_info in sim_sv_info_list:
                decoded_sim_sv_info_list.append(decode_sim_sv_info(sim_sv_info))
            
    #print(simulate_info_text)
    for sim_sv_info in create_sim_sv_info_list(simulate_info_text):
        decoded_sim_sv_info_list.append(decode_sim_sv_info(sim_sv_info))

    decoded_sim_sv_info_list.sort(key = getfirst_python)
    decoded_sim_sv_info_list = decoded_sim_sv_info_list[::-1]
    return decoded_sim_sv_info_list
def insertSVandoutput(parameterfilepath, inputgenomepath, altedgenomepath, outputvcfpath, heterozygous_ratio, mode):
    with open(parameterfilepath, 'r') as file:
        parameter_file = file.read()
    decoded_sim_sv_info_list = decode_parameterfile(parameter_file)
    print('Loading genome to memory')
    sim_contig2seq = dict()
    for rec in mp.fastx_read(inputgenomepath):
        #if(rec[0] not in ('chr21', 'chr22')):
            #continue
            #break
        sim_contig2seq[rec[0]] = rec[1].upper()
    print('Loading genome to memory completed ', len(sim_contig2seq))
    print('Checking N in genome')
    contig2usable_interval, contig_list, contigprob_list = prepare_usable_loc_info(sim_contig2seq)
    print('Checking N in genome completed')
    unique_id = 1
    contig2op_list = dict()
    print('Trying to insert variations')
    hetdict = dict()
    for decoded_info in decoded_sim_sv_info_list: 
        svend_1, svend_2, onesvinfo = decoded_info
        pass_flag, contig_1, refstart_1, contig_2, refstart_2 = get_start_loc(svend_1, svend_2, contig_list, contigprob_list, contig2usable_interval)
        if(pass_flag != True):
            print('Failed to simluate, no space availible')
            return
        eventcount = count_event(decoded_info)
        for oneinfo in onesvinfo:
            #print(oneinfo)
            if(np.random.randint(100) > (heterozygous_ratio*100)):
                hetdict[unique_id] = False
            else:
                hetdict[unique_id] = True
            for oneop in add_SV(contig_1, contig_2, refstart_1, refstart_2, oneinfo, sim_contig2seq, unique_id):
                oneop[-1].append(eventcount)
                oneop[-1].append(unique_id)
                try:
                    contig2op_list[oneop[0]].append(oneop)
                except:
                    contig2op_list[oneop[0]] = [oneop]

        unique_id += 1

    vcflist = []
    modifed_contig_dict = dict()
    for contig in contig2op_list:
        contig2op_list[contig].sort(key = getsecond_python)
        l = []
        l1 = []
        prestart = 0
        for item in contig2op_list[contig]:
        
            l.append(sim_contig2seq[contig][prestart: item[1]])
            l1.append(sim_contig2seq[contig][prestart: item[1]])
        
            
            l.append(item[3])
            unique_id = item[-1][-1]
            if(hetdict[unique_id] == False):
                l1.append(item[3])
            else:
                l1.append(sim_contig2seq[contig][item[1]: item[2]])
            
        
            prestart = item[2]
            vcflist.append(item[-1])
        l.append(sim_contig2seq[contig][prestart: ])
        l1.append(sim_contig2seq[contig][prestart: ])
        if(heterozygous_ratio > 0):
            modifed_contig_dict[contig+'_hap1'] = ''.join(l)
            modifed_contig_dict[contig+'_hap2'] = ''.join(l1)
        else:
            modifed_contig_dict[contig] = ''.join(l)
    print('Insert variations completed')
    print('writing alted genome to '+ altedgenomepath, len(modifed_contig_dict))
    ofile = open(altedgenomepath, "w")
    for contig in modifed_contig_dict:
        ofile.write(">" + contig + "\n" +modifed_contig_dict[contig] + "\n")
    ofile.close()

    vcflist.sort(key = get23_python)
    svtype2code = {'INS': 1, 'DEL': 2, 'INV': 3, 'DUP': 4, 'TRA': 5, 'BND': 5}
    tovcflist = []
    for line in vcflist:
        if(line[0] in ('DEL', 'INV', 'DUP')):
            tovcflist.append([line[1], line[2], line[3] - line[2], tuple(line[-2:]),-1, svtype2code[line[0]]])
        elif(line[0] == 'INS'):
            tovcflist.append([line[1], line[2], line[3], tuple(line[-2:]), -1,1])
        elif(line[0] == 'TRA'):
            tovcflist.append([line[1], line[2], line[5], tuple(line[-2:]), line[4], 5])
            tovcflist.append([line[1], line[3], line[6], tuple(line[-2:]), line[4], 5])

    contig2length = {}
    for contig in sim_contig2seq:
        contig2length[contig] = len(sim_contig2seq[contig])
    print('writing vcf file to '+ outputvcfpath)
    tovcf(tovcflist, contig2length, outputvcfpath)
    
    if(mode == 'reference'):
        print('Reference mode')
        bias = 0
        header = pysam.VariantFile(outputvcfpath, mode = 'r').header.copy()
        l = []
        fail = False
        for rec in pysam.VariantFile(outputvcfpath, mode = 'r'):
            #print(rec.contig, rec.start, rec.stop, rec.info['SVTYPE'], abs(rec.info['SVLEN'][0]))
            if(heterozygous_ratio > 0):
                rec.contig = rec.contig + '_hap1'
            rec.start = rec.start + bias


            if(rec.info['SVTYPE'] == 'INS'):
                rec.info['SVTYPE'] = 'DEL'
                rec.stop = rec.start + abs(rec.info['SVLEN'][0])
                bias += abs(rec.info['SVLEN'][0])

            elif(rec.info['SVTYPE'] == 'DEL'):
                rec.info['SVTYPE'] = 'INS'
                rec.stop = rec.start + 1 
                bias -= abs(rec.info['SVLEN'][0])
            else:
                print('Fail in reference mode: unsupported SVTYPE', rec.info['SVTYPE'])
                fail = True
                break
            l.append(rec)
        if(fail == True):
            return 1
        f = pysam.VariantFile(outputvcfpath, mode = 'w', header = header)
        for rec in l:
            f.write(rec.copy())
        f.close()
    print('All completed')
    return 0


pstring = sys.argv
forvaule = False
pdict = {}
op = ''
for item in pstring[1:]:
    if(forvaule == True):
        forvaule = False
        if(op in ('parameterfilepath', 'inputgenomepath', 'altedgenomepath', 'outputvcfpath', 'heterozygous_ratio', 'mode')):
            pdict[op] = item
        
    if('-' == item[:1]):
        forvaule = True
        op = item[1:]
hit = False
try:
    parameterfilepath = pdict['parameterfilepath']
    inputgenomepath = pdict['inputgenomepath']
    altedgenomepath = pdict['altedgenomepath']
    outputvcfpath = pdict['outputvcfpath']
    if('heterozygous_ratio' in pdict):
        heterozygous_ratio = float(pdict['heterozygous_ratio'])
    else:
        heterozygous_ratio = 0.8
    if('mode' in pdict):
        mode = pdict['mode']
    else:
        mode = 'read'
    hit = True
except:
    pass
    
if(hit == True):
    insertSVandoutput(parameterfilepath, inputgenomepath, altedgenomepath, outputvcfpath, heterozygous_ratio, mode)
else:
    print('Usage')
    print('python vacsim.py -parameterfilepath parameterfilepath -inputgenomepath inputgenomepath -altedgenomepath altedgenomepath -outputvcfpath outputvcfpath -heterozygous_ratio 0.8')
     
