'''   
 * @Title: VACmap 
 * @author: Ding Hongyu
 * @date: 17 6 2024
 * @version V1.0.0
'''
import edlib

import pandas as pd
import pysam
import pysam
import numpy as np
import heapq
from numba import njit, jit

from Bio import SeqIO

import time 

from numba.typed import Dict
from numba.typed import List
import matplotlib.pyplot as plt
import copy
from Bio.Seq import Seq
import re

from cigar import Cigar
from io import StringIO
import multiprocessing
import os
import vacmap_index as mp
import gzip
import array
import sys
import logging
def print_log(*a):
    a = str(a)
    logging.info(a)

@njit
def pos2contig(pos, contig2start):
    for precontig in contig2start:
        break
    for contig in contig2start:
        
        if(pos < contig2start[contig]):
            break
        precontig = contig
    return precontig
@njit
def getfirst(x):
    return x[0]
@njit
def getsecond(x):
    return x[1]
def sort_by_length(x):
    return x[4] - x[3]
def get_overlapsize(readloc_set_a, readloc_set_b):
    return len(readloc_set_a&readloc_set_b)/min(len(readloc_set_a), len(readloc_set_b))
def create_header(reffilepath):
    reader = SeqIO.parse(reffilepath, 'fasta')
    SQlist = []
    contig2iloc = {}
    iloc = -1
    for rec in reader:
        name = rec.id
        seq = str(rec.seq)
        iloc += 1
        SQlist.append({'LN': len(seq), 'SN': name})
        contig2iloc[name] = iloc
    header = {'HD': {'VN': '1.0'},
             'SQ': SQlist}
    return header, contig2iloc


def cigartocigartuple(cigar):
    cigartuple = []
    numberset = set()
    for i in range(10):
        numberset.add(str(i))
    number = ''
    for c in cigar:
        if(c in numberset):
            number += c
        else:
            cigartuple.append([c, int(number)])
            number = ''
    return cigartuple
def getreadrefend(cigarstring, readloc, refloc, reportsize = 50):
    index_object_result = cigartocigartuple(cigarstring)
    dellist = []
    inslist = []
    for item in index_object_result:
        if(item[0] in ['D', 'M', 'X', '=']):
            if(item[0] == 'D' and item[1] > reportsize):
                dellist.append([readloc, refloc, item[1]])
            refloc += item[1]
        if(item[0] in ['I', 'S', 'M', 'X', '=']):
            if(item[0] == 'I' and item[1] > reportsize):
                inslist.append([readloc, refloc, item[1]])
            readloc += item[1]

    print_log(readloc, refloc)
    return dellist, inslist
@njit
def get_length(x):
    return len(x)
@njit
def group_anchor(one_mapinfo, bias = 2000):
    one_mapinfo.sort(key = getsecond)

    cluster_list = List([List([one_mapinfo[0]])])
    for item in one_mapinfo[1:]:
        if(abs(item[1] - cluster_list[-1][-1][1]) > bias):
            if(len(cluster_list[-1]) < 3):
                cluster_list.pop(-1)

            cluster_list.append(List([item]))
        else:
            cluster_list[-1].append(item)

    if(len(cluster_list[-1]) < 3):
        cluster_list.pop(-1)
    cluster_list.sort(key = get_length)
    return cluster_list

@njit
def get_readloc_set_bin(one_mappos, bin_size):
    readloc_set = set()
    for item in one_mappos:
        readloc_set.add(item[0]//bin_size)
    return readloc_set

@njit
def get_chaininfo_merged(S, P, take_index, one_mapinfo, detail_dict, count, min_scores = 15., kmersize = 15, maxgap = 500):#count were used
    
    unusedindex_set = set(list(range(len(S))))
    usedindex_set = set()
    chaininfo = List()

    while(True):
        if(S[take_index] > min_scores):
            #print_log(S[take_index])
            path = List()
            
            count += 1
            #print_log(len(usedindex_set), len(unusedindex_set))
            onechaininfo = List([S[take_index], one_mapinfo[take_index][0], one_mapinfo[take_index][1], count])
            usedindex_set.add(take_index)
            unusedindex_set.remove(take_index)
            path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))



            while(True):
                if((P[take_index] == 9999999) or (P[take_index] in usedindex_set)):
                    break
                pre_take_index = take_index
                take_index = P[take_index]
                usedindex_set.add(take_index)
                unusedindex_set.remove(take_index)
                path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
                
                
            onechaininfo.insert(1, one_mapinfo[take_index][0])
            onechaininfo.insert(3, one_mapinfo[take_index][1])
            #onechaininfo.insert(5, take_index)
            onechaininfo[0] -= S[take_index]
            
            chaininfo.append(onechaininfo)
            detail_dict[count] = path[::-1]
            break
            max_scores = kmersize
            if(len(unusedindex_set) == 0):
                break
            if((onechaininfo[0] < 20.) or (len(path) < 3)):
                chaininfo.pop(-1)
            for tmp_index in unusedindex_set:
                if(max_scores <= S[tmp_index]):
                    max_scores = S[tmp_index]
                    take_index = tmp_index
            #print_log(len(usedindex_set), len(unusedindex_set))
            #print_log()
            break
        else:
            break
    #return dropoverlap(chaininfo), detail_dict, count
    return chaininfo, detail_dict, count



@njit
def getClosest(val1, val2, target, iloc1, iloc2):
 
    if(target - val1 >= val2 - target):
        return val2 - target, target - val1, iloc2, iloc1
    else:
        return target - val1, val2 - target, iloc1, iloc2
@njit
def seq2hashtable_multi_test(onelookuptable, seq, start, kmersize):

    for iloc in range(0, len(seq) - kmersize + 1, 1):
        hashedkmer = hash(seq[iloc:iloc+kmersize])
        if(hashedkmer in onelookuptable):

            onelookuptable[hashedkmer].append(start + iloc)
        else:
            onelookuptable[hashedkmer] = List([start + iloc])



@njit
def test_indel_present(cigar, report_size = 10):#tested

    zero = ord('0')
    INSsyb = ord('I') - zero
    SOFTsyb = ord('S') - zero
    HARDsyb = ord('H') - zero
    PADsyb = ord('P') - zero
    DELsyb = ord('D') - zero
    SKIPsyb = ord('N') - zero
    delsyb = ord('^') - zero
    
    DEL = False
    INS = False
    
    BAD = False
    
    count = 0
    
    typed_cigar = List()
    [typed_cigar.append(ord(item) - zero) for item in cigar] 
    
    number = 0
    for item in typed_cigar:
        if(item < 10):
            number = number * 10 + item
            continue

        elif((number >= report_size) and (item in (DELsyb, INSsyb))):
            count += 1
            if(item == DELsyb):
                
                DEL = True
                if(INS == True):
                    BAD = True
                    break
            else:
                INS = True
                if(DEL == True):
                    BAD = True
                    break
        number = 0
    if(count > 1):
        BAD = True
    return BAD, DEL, INS


@njit
def get_reverse_complement(seq):
    s2n = dict()
    s2n['A'] = 'T'
    s2n['T'] = 'A'
    s2n['G'] = 'C'
    s2n['C'] = 'G'
    s2n['N'] = 'N'
    r_seq = ''
    for c in seq[::-1]:
        r_seq += s2n[c]
    return r_seq

@njit
def List_merge(List_of_List):
    merged_List = List_of_List[0][:]
    for one_List in List_of_List[1:]:
        for item in one_List:
            merged_List.append(item)
    return merged_List
@njit
def get_target_query_for_drop(a_1, a_2, testseq, kmersize, pos2contig, contig2start, contig2seq, reverse = False, rc = False, debug = False):

    if(a_1[2] == 1):
        query_st = a_1[0]
        query_en = a_2[0] + a_2[3]

    else:
        query_st = a_1[0]
        query_en = a_2[0] + a_2[3]

    query = testseq[query_st: query_en]
    testcontig = pos2contig(a_1[1], contig2start)
    if(a_1[2] == 1):
        target_st = a_1[1] - contig2start[testcontig]
        target_en = a_2[1] + a_2[3] - contig2start[testcontig]
    else:
        target_en = a_1[1] + min(kmersize, a_1[3]) - contig2start[testcontig]
        target_st = a_2[1] - a_2[3] + min(kmersize, a_2[3]) - contig2start[testcontig]
    target = contig2seq[testcontig][target_st: target_en]
    if(reverse == False and rc == False):
        return target, query, target_st+contig2start[testcontig], target_en+contig2start[testcontig], query_st, query_en
    elif(rc == True and reverse == True):
        return get_reverse_complement(target)[::-1], query[::-1], target_st+contig2start[testcontig], target_en+contig2start[testcontig], query_st, query_en
    elif(reverse == True):
        return target[::-1], query[::-1], target_st+contig2start[testcontig], target_en+contig2start[testcontig], query_st, query_en
    else:
        return get_reverse_complement(target), query, target_st+contig2start[testcontig], target_en+contig2start[testcontig], query_st, query_en

@njit
def get_refseq(contig, start, end, contig2seq):
    return contig2seq[contig][start: end]
@njit
def get_reflen(contig, contig2seq):
    return len(contig2seq[contig])
@njit
def get_contig2start(contig, contig2start):
    return contig2start[contig]

import array

def compute_NM_tag(query, target):
    return edlib.align(query = query, target = target, task = 'distance')['editDistance']
def get_bam_dict1(mapinfo, query, qual, contig2iloc, contig2seq):
    #'hhk',         ,  '1', '+', 11, 9192, 2767041, 2776138, 60
    #      0            1    2   3    4      5         6      7
    #'18_19897150_+', '18', '+', 0, 4776, 19832244, 19837393, 1]
    for line in mapinfo:
        if(len(Cigar(line[-1])) != len(query)):
            print_log()
            print_log(line[0], len(Cigar(line[-1])), len(query))
            print_log()
            return 
            break
    mq = mapinfo[-1][7]
    if(mq != 0):
        mq = 60
    rc_query = str(Seq(query).reverse_complement())
    mapinfo.sort(key = sort_by_length)
    mapinfo = mapinfo[::-1]
    iloc2nm = dict()
    tmpiloc = -1
    for item in mapinfo:
        tmpiloc += 1
        if(item[2] == '+'):
            nm = compute_NM_tag(query[item[3]: item[4]], get_refseq(item[1], item[5], item[6], contig2seq))
        else:
            nm = compute_NM_tag(rc_query[item[3]: item[4]], get_refseq(item[1], item[5], item[6], contig2seq))
        iloc2nm[tmpiloc] = nm
    if((qual != None) and (len(qual) == len(query))):
        query_qualities = fastq_q2b(qual)
        rc_query_qualities = query_qualities[::-1]
    a_list = []
    for iloc in range(len(mapinfo)):
        bam_dict = dict()
        primary = mapinfo[iloc]
        bam_dict['readid'] = primary[0]
        bam_dict['contig'] = primary[1]
        if(iloc == 0):
            base_value = 0
        else:
            base_value = 2048
        if(primary[2] == '+'):
            bam_dict['flag'] = 0 + base_value

        else:
            bam_dict['flag'] = 16 + base_value


        bam_dict['refstart'] = primary[5]


        bam_dict['cigar'] = primary[8]

        if(len(mapinfo) > 1):
            salist = []
            tmpiloc = -1
            for item in mapinfo:
                tmpiloc += 1
                if(tmpiloc == iloc):
                    continue
                nm = iloc2nm[tmpiloc]
                salist.append(''.join((item[1], ',', str(item[5]+1), ',', item[2], ',', item[8], ',', str(mq), ',', str(nm)+';')))

            bam_dict['sa'] = ''.join(salist)


        a = pysam.AlignedSegment()
        a.query_name = bam_dict['readid']

        a.flag = bam_dict['flag']
        a.reference_id = contig2iloc[bam_dict['contig']]
        a.reference_start = bam_dict['refstart']
        item = primary

        a.mapping_quality = mq
        a.cigarstring = bam_dict['cigar']
        if(item[2] == '+'):
            a.query_sequence = query
            a.template_length = len(query)
            if((qual != None) and (len(qual) == len(query))):
                a.query_qualities = query_qualities
        else:
            a.query_sequence = rc_query
            a.template_length = len(rc_query)
            if((qual != None) and (len(qual) == len(query))):
                a.query_qualities = rc_query_qualities
        if('sa' in bam_dict):
            a.tags = [('NM', iloc2nm[iloc]), ("SA", bam_dict['sa'])]
        else:
            a.tags = [('NM', iloc2nm[iloc])]
        a_list.append(a)
    return a_list

@njit
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_numpy(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    g_max_scores = 0.
    g_max_index = -1
     



    n = len(one_mapinfo)
    S = np.zeros(n, np.float64)
    P = np.zeros(n, np.int64)
    start = 0
    
    testdict = dict()
    testdict[-1] = 0.
    testdict.pop(-1)

    
    for i in range(n):
        P[i] = 9999999
        max_scores = one_mapinfo[i][3]
        pre_index = 9999999


        #while(((i - start) > 100) and ((one_mapinfo[i][0] - one_mapinfo[start][0]) > max_gap)):
            #start += 1
        j = i
        skipcount = 0
        #while(j > start):
            #j -= 1
        if(len(testdict) > 50):

            cutoff = g_max_scores - max(2.5*skipcost, 60.)
            for loc_in_one_mapinfo in testdict:
                if(S[loc_in_one_mapinfo] < cutoff and (one_mapinfo[i][0] - one_mapinfo[loc_in_one_mapinfo][0]>100)):
                    
                    testdict.pop(loc_in_one_mapinfo)
                    #else:
                        #break




               
        
        for j in testdict:
            
                    
            #(12240, 2791460820, 1, 21) (12251, 2791460711, 1, 108)
            nocost = False
            filtered = True
            if(one_mapinfo[i][0] == one_mapinfo[j][0]):
                continue

                      #      x_2 - x_1 - a_1
            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            if((readgap < 0)):
                continue
                #if((one_mapinfo[i][2] == one_mapinfo[j][2]) and ((readgap + one_mapinfo[i][3]) > 10)):
                    #pass
                #else:
                    #continue

            while(True):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    break

                if(one_mapinfo[i][2] == 1):
                              #    y_2 - y_1 - a_1
                    refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                else:
                             #     y_1 - y_2 - a_1
                    refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                
    
                if((one_mapinfo[i][2] == 1) and (refgap < 0)):
                    break
                if(one_mapinfo[i][2] == -1):
                    if(refgap > 0):
                        break
                    else:
                        refgap = abs(refgap)

                gapcost = abs(abs(readgap) - refgap)
                if(((readgap ) > maxgap) or (gapcost > maxdiff)):
                    break

                if(gapcost != 0):
                    gapcost = 0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost)
                    #gapcost = 0.01 * kmersize * gapcost + 2.0 * np.log2(gapcost)
        
        
                test_scores = S[j] + (min(readgap, 0) + one_mapinfo[i][3]) - gapcost
                filtered = False
                break
            if(filtered == True):
                """if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(magic > abs(one_mapinfo[i][1] - one_mapinfo[j][1])):
                        continue"""
                test_scores = S[j] - skipcost + one_mapinfo[i][3]# - np.log2(readgap)# - np.log2(abs(abs(refgap) - readgap)+1)
                #continue
            if(test_scores >= max_scores):
                max_scores = test_scores
                pre_index = j

       
            
                
        S[i] = max_scores
        P[i] = pre_index 
        
        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
            
 
    
            

        testdict[i] = max_scores
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    while(True):
        if((P[take_index] == 9999999)):
            break
        take_index = P[take_index]
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

    return g_max_scores, np.array(path)#, testdict


@njit
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_linear(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    g_max_scores = 0.
    g_max_index = -1
     



    n = len(one_mapinfo)
    S = np.zeros(n, np.float64)
    P = np.zeros(n, np.int64)
    start = 0
    
    testdict = dict()
    testdict[-1] = 0.
    testdict.pop(-1)

    
    for i in range(n):
        P[i] = 9999999
        max_scores = one_mapinfo[i][3]
        pre_index = 9999999


        #while(((i - start) > 100) and ((one_mapinfo[i][0] - one_mapinfo[start][0]) > max_gap)):
            #start += 1
        j = i
        skipcount = 0
        #while(j > start):
            #j -= 1
        if(len(testdict) > 50):

            cutoff = g_max_scores - max(2.5*skipcost, 60.)
            for loc_in_one_mapinfo in testdict:
                if(S[loc_in_one_mapinfo] < cutoff and (one_mapinfo[i][0] - one_mapinfo[loc_in_one_mapinfo][0]>100)):
                    
                    testdict.pop(loc_in_one_mapinfo)
                    #else:
                        #break




               
        
        for j in testdict:
            
                    
            #(12240, 2791460820, 1, 21) (12251, 2791460711, 1, 108)
            nocost = False
            filtered = True
            if(one_mapinfo[i][0] == one_mapinfo[j][0]):
                continue

                      #      x_2 - x_1 - a_1
            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            if((readgap < 0)):
                continue
                #if((one_mapinfo[i][2] == one_mapinfo[j][2]) and ((readgap + one_mapinfo[i][3]) > 10)):
                    #pass
                #else:
                    #continue

            while(True):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    break

                if(one_mapinfo[i][2] == 1):
                              #    y_2 - y_1 - a_1
                    refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                else:
                             #     y_1 - y_2 - a_1
                    refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                
    
                if((one_mapinfo[i][2] == 1) and (refgap < 0)):
                    break
                if(one_mapinfo[i][2] == -1):
                    if(refgap > 0):
                        break
                    else:
                        refgap = abs(refgap)

                gapcost = abs(abs(readgap) - refgap)
                if(((readgap ) > maxgap) or (gapcost > maxdiff)):
                    break

                if(gapcost != 0):
                    gapcost = 0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost)
                    #gapcost = 0.01 * kmersize * gapcost + 2.0 * np.log2(gapcost)
        
        
                test_scores = S[j] + (min(readgap, 0) + one_mapinfo[i][3]) - gapcost
                filtered = False
                break
            if(filtered == True):

                #test_scores = S[j] - skipcost + one_mapinfo[i][3]# - np.log2(readgap)# - np.log2(abs(abs(refgap) - readgap)+1)
                continue
            if(test_scores >= max_scores):
                max_scores = test_scores
                pre_index = j

       
            
                
        S[i] = max_scores
        P[i] = pre_index 
        
        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
            
 
    
            

        testdict[i] = max_scores
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    while(True):
        if((P[take_index] == 9999999)):
            break
        take_index = P[take_index]
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

    return g_max_scores, path#, testdict



@njit
def findClosest(arr, target):
    n = len(arr)
    if(target <= arr[0][0]):
        return arr[0][0] - target, arr[0][0] - target, 0, 0
    if(target >= arr[n - 1][0]):
        return target - arr[n - 1][0], target - arr[n - 1][0], n-1, n-1

    i = 0
    j = n
    mid = 0
    while(i < j):
        mid = (i + j) // 2

        if(arr[mid][0] == target):
            return abs(arr[mid][0] - target), abs(arr[mid][0] - target), mid, mid

        if(target < arr[mid][0]) :

            if(mid > 0 and target > arr[mid - 1][0]):
                return abs(arr[mid - 1][0] - target), abs(arr[mid][0] - target), mid-1, mid

            j = mid

        else :
            if(mid < n - 1 and target < arr[mid + 1][0]):
                return abs(arr[mid][0] - target), abs(arr[mid + 1][0] - target), mid, mid+1

            i = mid + 1

    return abs(arr[mid][0] - target), abs(arr[mid][0] - target), mid, mid


@njit
def drop_misplaced_alignment_test(alignment_list, iloc, debug):
    if((alignment_list[iloc][0][2] == alignment_list[iloc + 1][0][2]) and (alignment_list[iloc][0][2] == alignment_list[iloc + 2][0][2])):


        prealignment_size = alignment_list[iloc][-1][0] + alignment_list[iloc][-1][3] - alignment_list[iloc][0][0]
        midalignment_size = alignment_list[iloc+1][-1][0] + alignment_list[iloc+1][-1][3] - alignment_list[iloc+1][0][0]
        if(midalignment_size > 1000):
            return False
        nowalignment_size = alignment_list[iloc+2][-1][0] + alignment_list[iloc+2][-1][3] - alignment_list[iloc+2][0][0]

        mid_gap = alignment_list[iloc+2][0][0] - (alignment_list[iloc][-1][0] + alignment_list[iloc][-1][3])
        #if(midalignment_size < 500 and (min(mid_gap, midalignment_size) / max(mid_gap, midalignment_size)) < 0.7):
            #if(debug): print_log('(min(mid_gap, midalignment_size) / max(mid_gap, midalignment_size)) < 0.7')
            #alignment_list.pop(iloc+1)
            #if(debug): print_log('misplaced alignment removed ', midalignment_size)
            #return True

        preitem = alignment_list[iloc][-1]
        nowitem = alignment_list[iloc + 1][0]
        readgap = nowitem[0] - preitem[0] - preitem[3]
        if(preitem[2] == 1):
            refgap = nowitem[1] - preitem[1] - preitem[3]

        else:
            refgap = preitem[1]  - nowitem[1] - nowitem[3]

        if(abs(refgap) < 100000):
            DEL = 0
            INS = 0
            if((readgap - refgap) < -30):
                DEL += 1
            elif((readgap - refgap) > 30):
                INS += 1
            else:

                return False
            gap_1 = abs(readgap - refgap)
            preitem = alignment_list[iloc + 1][-1]
            nowitem = alignment_list[iloc + 2][0]

            readgap = nowitem[0] - preitem[0] - preitem[3]
            if(preitem[2] == 1):
                refgap = nowitem[1] - preitem[1] - preitem[3]

            else:
                refgap = preitem[1]  - nowitem[1] - nowitem[3]

            if(abs(refgap) < 100000):
                if((readgap - refgap) < -30):
                    DEL += 1
                elif((readgap - refgap) > 30):
                    INS += 1
                else:

                    return False
                gap_2 = abs(readgap - refgap)
                if(DEL == 1 and INS == 1 and (midalignment_size < 500 or (max(gap_1, gap_2)/midalignment_size) > 0.5)):

                    alignment_list.pop(iloc+1)
                    return True
    return False

@njit
def get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start):
    testcontig = pos2contig(preitem[1], contig2start)
    refbias = contig2start[testcontig]
    query_st, query_en = preitem[0], nowitem[0]
    if(preitem[2] == 1):
        target_st, target_en = preitem[1], nowitem[1]
        query = testseq[preitem[0]: nowitem[0]]
        target = contig2seq[testcontig][preitem[1] - refbias: nowitem[1] - refbias]

    else:
        target_st, target_en = nowitem[1] + min(kmersize, nowitem[3]), preitem[1] + min(kmersize, preitem[3])
        query = rc_testseq[testseq_len - nowitem[0]: testseq_len - preitem[0]]
        target = contig2seq[testcontig][nowitem[1] - refbias + min(kmersize, nowitem[3]): preitem[1] - refbias+min(kmersize, preitem[3])]
    return target, query, target_st, target_en, query_st, query_en



    
            





    
def extend_edge_test(testseq, testseq_len, alignment_list, kmersize, pos2contig, contig2start, contig2seq, san, debug):



    zdrop_value = 400  

    mismatch_value = -4

    gap_open_extend = -4
    gap_extend_extend = -4
    zdrop_value_extend = 50
    bw = 100




    max_extend_size = 20000
    if(debug == True): print_log('Extend')
    if(debug == True): print_log('len(alignment_list): ', len(alignment_list))
    onealignment_index = -1
    while(True):
        onealignment_index += 1
        if(debug == True): print_log('onealignment_index', onealignment_index)
        if(onealignment_index >= len(alignment_list)):
            break
        onealignment = alignment_list[onealignment_index]
        if(debug == True): print_log('Alignment: start end ', onealignment[0][0], onealignment[-1][0] + onealignment[-1][3])
        if(debug): print_log(onealignment[0], onealignment[-1])
        if(onealignment[0][0]>0):
            while(True):

                if(onealignment[0][0] == 0):
                    break
                pre_onealignment_index = max(onealignment_index - san, 0)
                if(onealignment_index == 0 or  onealignment_index - san < 0):
                    looksize = onealignment[0][0] - 0
                else:
                    looksize = onealignment[0][0] - (alignment_list[pre_onealignment_index][-1][0] + alignment_list[pre_onealignment_index][-1][3])

                if(looksize < 0):

                    print_log('error')
                    print_log(alignment_list[pre_onealignment_index][-1], onealignment[0])


                preitem = onealignment[0]
                nowitem = onealignment[1]
                if(debug): print_log('enlong start')
                if(debug): print_log(preitem)
                if(debug): print_log(nowitem)

                target, query, target_st, target_en, query_st, query_en = get_target_query_for_drop(preitem, nowitem, testseq, kmersize, pos2contig, contig2start, contig2seq, reverse = False, rc = False)

                testcontig = pos2contig(preitem[1], contig2start)

                if(preitem[2] == 1):
                    
                    #looksize = min(looksize, target_st - contig2start[testcontig])
                    looksize = min(looksize, target_st - get_contig2start(testcontig, contig2start))
                    
                    if(looksize > max_extend_size):
                        looksize = max_extend_size
                        skip_loop = False
                    else:
                        skip_loop = True
                    if(looksize == 0):
                        if(debug): print_log('looksize == 0, skiped')
                        break
                    query = testseq[max(query_st - looksize, 0): query_st][::-1]
                    
                    #target = contig2seq[testcontig][target_st-contig2start[testcontig] - len(query): target_st-contig2start[testcontig]][::-1]
                    target = get_refseq(testcontig, target_st-get_contig2start(testcontig, contig2start) - len(query), target_st-get_contig2start(testcontig, contig2start), contig2seq)[::-1]
                    
                    if(debug): print_log(len(target), len(query), target_st, target_en, query_st, query_en)


                    #zdroped, max_extend, cigarstring, Score = fast_globalms_align_extend(target, query, match = 2, mismatch = -4, gap_open = gap_open_extend, gap_extend = gap_extend_extend, zdrop = zdrop_value_extend, loc_only = True)
                    cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 4, gap_open_2 = 4, gap_extend_2 = 4, bw=bw, zdropvalue=zdrop_value_extend)
          
                    
                    max_extend = (t_e, q_e)
                    onealignment[0] = ((query_st - max_extend[1], target_st - max_extend[0], 1, 0))
                else:
                    
                    #looksize = min(looksize, contig2start[testcontig] + len(contig2seq[testcontig]) - target_en)
                    looksize = min(looksize, get_contig2start(testcontig, contig2start) + get_reflen(testcontig, contig2seq) - target_en)
                    
                    if(looksize > max_extend_size):
                        looksize = max_extend_size
                        skip_loop = False
                    else:
                        skip_loop = True
                    if(looksize == 0):
                        if(debug): print_log('looksize == 0, skiped')
                        break
                    query = testseq[max(query_st - looksize, 0): query_st][::-1]
                    
                    #target = get_reverse_complement(contig2seq[testcontig][target_en-contig2start[testcontig]: target_en + len(query)-contig2start[testcontig]])[::-1]
                    target = str(Seq(get_refseq(testcontig, target_en-get_contig2start(testcontig, contig2start), target_en + len(query)-get_contig2start(testcontig, contig2start), contig2seq)).reverse_complement())[::-1]
                    
                    if(debug): print_log(len(target), len(query), target_st, target_en, query_st, query_en)
                        
  
                    #zdroped, max_extend, cigarstring, Score = fast_globalms_align_extend(target, query, match = 2, mismatch = -4, gap_open = gap_open_extend, gap_extend = gap_extend_extend, zdrop = zdrop_value_extend, loc_only = True)
                    cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 4, gap_open_2 = 4, gap_extend_2 = 4, bw=bw, zdropvalue=zdrop_value_extend)

                    
                    max_extend = (t_e, q_e)
                    onealignment[0] = ((query_st - max_extend[1], target_en + max_extend[0], -1, 0))
                if(debug): print_log(max_extend)
                if(debug): print_log(onealignment[0])
                if(debug): print_log(target)
                if(debug): print_log(query)


                break
        else:
            tmpitem = alignment_list[onealignment_index][0]
            if(alignment_list[onealignment_index][0][2] == 1):
                alignment_list[onealignment_index][0] = ((tmpitem[0], tmpitem[1], 1, 0))
            else:
                alignment_list[onealignment_index][0] = ((tmpitem[0], tmpitem[1] + tmpitem[3], -1, 0))
            
        if((onealignment[-1][0]+onealignment[-1][3])<len(testseq)):
            while(True):
                if((onealignment[-1][0]+onealignment[-1][3])>=len(testseq)):
                    break
                pre_onealignment_index = min(onealignment_index + san, len(alignment_list))
                if(pre_onealignment_index == len(alignment_list)):
                    looksize = testseq_len - (alignment_list[onealignment_index][-1][0] + alignment_list[onealignment_index][-1][3])
                else:
                    looksize = alignment_list[pre_onealignment_index][0][0] - (onealignment[-1][0] + onealignment[-1][3])
                if(looksize < 0):
                    print_log('error')
                    print_log(onealignment[-1], alignment_list[pre_onealignment_index][0])



                preitem = onealignment[-2]
                nowitem = onealignment[-1]
                if(debug): print_log('enlong end')
                if(debug): print_log(preitem)
                if(debug): print_log(nowitem)

                target, query, target_st, target_en, query_st, query_en = get_target_query_for_drop(preitem, nowitem, testseq, kmersize, pos2contig, contig2start, contig2seq, reverse = False, rc = False)

                testcontig = pos2contig(preitem[1], contig2start)

                if(preitem[2] == 1):
                    
                    #looksize = min(looksize, contig2start[testcontig] + len(contig2seq[testcontig]) - target_en)
                    looksize = min(looksize, get_contig2start(testcontig, contig2start) + get_reflen(testcontig, contig2seq) - target_en)
                    
                    if(looksize > max_extend_size):
                        looksize = max_extend_size
                        skip_loop = False
                    else:
                        skip_loop = True
                    if(looksize == 0):
                        if(debug): print_log('looksize == 0, skiped')
                        break
                    query = testseq[query_en: query_en + looksize]
                    
                    #target = contig2seq[testcontig][target_en-contig2start[testcontig]: target_en+len(query)-contig2start[testcontig]]
                    target = get_refseq(testcontig, target_en-get_contig2start(testcontig, contig2start), target_en+len(query)-get_contig2start(testcontig, contig2start), contig2seq)
                    
                    if(debug): print_log(len(target), len(query), target_st, target_en, query_st, query_en)

                    #zdroped, max_extend, cigarstring, Score = fast_globalms_align_extend(target, query, match = 2, mismatch = -4, gap_open = gap_open_extend, gap_extend = gap_extend_extend, zdrop = zdrop_value_extend, loc_only = True)
                    cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 4, gap_open_2 = 4, gap_extend_2 = 4, bw=bw, zdropvalue=zdrop_value_extend)
                    
                    max_extend = (t_e, q_e)
                    
                    onealignment[-1] = ((query_en + max_extend[1], target_en + max_extend[0], 1, 0))
                else:
                    
                    #looksize = min(looksize, target_st - contig2start[testcontig])
                    looksize = min(looksize, target_st - get_contig2start(testcontig, contig2start))
                    
                    if(looksize > max_extend_size):
                        looksize = max_extend_size
                        skip_loop = False
                    else:
                        skip_loop = True
                    if(looksize == 0):
                        if(debug): print_log('looksize == 0, skiped')
                        break
                    query = testseq[query_en: query_en + looksize]
                    
                    #target = get_reverse_complement(contig2seq[testcontig][target_st-contig2start[testcontig] - len(query): target_st-contig2start[testcontig]])
                    target = str(Seq(get_refseq(testcontig, target_st-get_contig2start(testcontig, contig2start) - len(query), target_st-get_contig2start(testcontig, contig2start), contig2seq)).reverse_complement())
                    
                    if(debug): print_log(len(target), len(query), target_st, target_en, query_st, query_en)

                    #zdroped, max_extend, cigarstring, Score = fast_globalms_align_extend(target, query, match = 2, mismatch = -4, gap_open = gap_open_extend, gap_extend = gap_extend_extend, zdrop = zdrop_value_extend, loc_only = True)
                    cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 4, gap_open_2 = 4, gap_extend_2 = 4, bw=bw, zdropvalue=zdrop_value_extend)
                    max_extend = (t_e, q_e)
                    onealignment[-1] = ((query_en + max_extend[1], target_st - max_extend[0], -1, 0))
                if(debug): print_log(max_extend)
                if(debug): print_log(onealignment[-1])
                if(debug): print_log(target)
                if(debug): print_log(query)

                break
        else:
            tmpitem = alignment_list[onealignment_index][-1]
            if(alignment_list[onealignment_index][-1][2] == 1):
                alignment_list[onealignment_index][-1] = ((tmpitem[0] + tmpitem[3], tmpitem[1] + tmpitem[3], 1, 0))
            else:
                alignment_list[onealignment_index][-1] = ((tmpitem[0] + tmpitem[3], tmpitem[1], -1, 0))


            
        if(debug): print_log()

    if(debug): print_log()


@njit
def testoverlap(alignment_list, iloc):
    refstart, refend = min(alignment_list[iloc][0][1], alignment_list[iloc][-1][1]), max(alignment_list[iloc][0][1], alignment_list[iloc][-1][1])
    overlapsize = 0
    for jloc in range(len(alignment_list)):
        if(jloc == iloc):
            continue
        else:
            test_refstart, test_refend = min(alignment_list[jloc][0][1], alignment_list[jloc][-1][1]), max(alignment_list[jloc][0][1], alignment_list[jloc][-1][1])

            if((test_refend > refstart) and (test_refstart < refend)):
                tmp = min(test_refend - refstart, refend - test_refstart)
                if(overlapsize < tmp):
                    overlapsize = tmp
                    
                
    return overlapsize
def getdupiloc(alignment_list):

    duplist = []
    if(len(alignment_list) >= 2):
            
        iloc = 0
        
        while((iloc + 1) < len(alignment_list)):
            readpos_1 = alignment_list[iloc][-1][0] + alignment_list[iloc][-1][3]
            if(alignment_list[iloc][-1][2] == 1):
                refpos_1 = alignment_list[iloc][-1][1] + alignment_list[iloc][-1][3]#highest position 
                strand_1 = 1
            else:
                refpos_1 = alignment_list[iloc][-1][1]#lowest position
                strand_1 = -1
            jloc = iloc
            hit = False
            dupsize = 0
            while((jloc + 1) < len(alignment_list)):
                jloc += 1
                
                if(alignment_list[jloc][-1][2] == 1):
                    refpos_2 = alignment_list[jloc][0][1]#lowest position
                    strand_2 = 1
                else:
                    refpos_2 = alignment_list[jloc][0][1] + alignment_list[jloc][0][2]#highest position
                    strand_2 = -1
                if(strand_1 != strand_2):
                    continue
                
                if(strand_1 == 1):
                    if((refpos_2 - refpos_1) < 50):
                        new_iloc = jloc
                        dupsize = refpos_2 - refpos_1
                        readpos_2 = alignment_list[jloc][0][0]
                        hit = True
          
                else:
                    if((refpos_1 - refpos_2) < 50):
                        new_iloc = jloc
                        dupsize = refpos_1 - refpos_2
                        readpos_2 = alignment_list[jloc][0][0]
                        hit = True

            if(hit == True):
                readgap = readpos_2 - readpos_1
                if(((iloc + 1) < new_iloc) or (((dupsize - readgap) < -30) and (readgap < 30))):
                    for skipiloc in range(iloc, new_iloc):
                        duplist.append(skipiloc)

                    
                iloc = new_iloc
            else:
                iloc += 1
    return duplist

@njit
def rebuild_chain_break(contig2start, raw_alignment_list, large_cost, small_alignment = 50, small_dup = -100):
    #rebuild chain break
    #step 1
    #remove anchor cause large cost and small alignment
    preitem = raw_alignment_list[0]
    alignment_list = List([List([preitem])])
    for nowitem in raw_alignment_list[1:]:
        if(preitem[2] == nowitem[2]):
            readgap = nowitem[0] - preitem[0] - preitem[3]


            if(preitem[2] == 1):
                refgap = nowitem[1] - preitem[1] - preitem[3]
            else:
                refgap = preitem[1]  - nowitem[1] - nowitem[3]

            if((refgap < 0) and (refgap > small_dup)):
                continue
            if((abs(readgap - refgap) <= large_cost) and (refgap >= 0) and (readgap < 100)):
                if(pos2contig(preitem[1], contig2start) == pos2contig(nowitem[1], contig2start)):
                    alignment_list[-1].append(nowitem)
                    preitem = nowitem
                    continue

        if(len(alignment_list[-1]) == 1):
            alignment_list.pop(-1)
        if((alignment_list[-1][-1][0] + alignment_list[-1][-1][3] - alignment_list[-1][0][0]) < small_alignment):
            alignment_list.pop(-1)
        if(len(alignment_list) > 0):
            if(len(alignment_list[-1]) >= 4):#
                alignment_list[-1] = alignment_list[-1][1:-1]
        alignment_list.append(List([nowitem]))
        preitem = nowitem
    if(len(alignment_list[-1]) == 1):
        alignment_list.pop(-1)
    if((alignment_list[-1][-1][0] + alignment_list[-1][-1][3] - alignment_list[-1][0][0]) < small_alignment):
        alignment_list.pop(-1)
    if(len(alignment_list[-1]) >= 4):
        alignment_list[-1] = alignment_list[-1][1:-1]

    return alignment_list

@njit
def rebuild_chain_break_H(contig2start, raw_alignment_list, large_cost, small_alignment = 20, small_dup = -20):
    #rebuild chain break
    #step 1
    #remove anchor cause large cost and small alignment
    small_alignment = 20
    small_dup = -20
    preitem = raw_alignment_list[0]
    alignment_list = List([List([preitem])])
    for nowitem in raw_alignment_list[1:]:
        if(preitem[2] == nowitem[2]):
            readgap = nowitem[0] - preitem[0] - preitem[3]


            if(preitem[2] == 1):
                refgap = nowitem[1] - preitem[1] - preitem[3]
            else:
                refgap = preitem[1]  - nowitem[1] - nowitem[3]

            if((refgap < 0) and (refgap > small_dup)):
                continue
            if((abs(readgap - refgap) <= large_cost) and (refgap >= 0)):
                if(pos2contig(preitem[1], contig2start) == pos2contig(nowitem[1], contig2start)):
                    alignment_list[-1].append(nowitem)
                    preitem = nowitem
                    continue

        if(len(alignment_list[-1]) == 1):
            alignment_list.pop(-1)
        if((alignment_list[-1][-1][0] + alignment_list[-1][-1][3] - alignment_list[-1][0][0]) < small_alignment):
            alignment_list.pop(-1)

        alignment_list.append(List([nowitem]))
        preitem = nowitem
    if(len(alignment_list[-1]) == 1):
        alignment_list.pop(-1)
    if((alignment_list[-1][-1][0] + alignment_list[-1][-1][3] - alignment_list[-1][0][0]) < small_alignment):
        alignment_list.pop(-1)


    return alignment_list


@njit
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = [0]
    for gapcost in range(1, maxdiff + 1):
        gapcost_list.append(0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))


    n = len(one_mapinfo)
    S = np.zeros(n, np.float64)
    P = np.zeros(n, np.int64)
    start = 0
    
    testdict = dict()
    testdict[-1] = 0.
    testdict.pop(-1)
    maxop = 2000*(one_mapinfo[-1][0]-one_mapinfo[0][0])
    #print_log(maxop)
    opcount = 0
    
    for i in range(n):
        P[i] = 9999999
        max_scores = one_mapinfo[i][3]
        pre_index = 9999999
        if(opcount > 2*maxop):
            path = []
            take_index = g_max_index
            path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
            while(True):
                if((P[take_index] == 9999999)):
                    break
                take_index = P[take_index]
                path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

            return g_max_scores, path
        if(opcount > maxop):
            if(max_scores < 15):
                S[i] = max_scores
                continue
        #while(((i - start) > 100) and ((one_mapinfo[i][0] - one_mapinfo[start][0]) > max_gap)):
            #start += 1
        j = i
        skipcount = 0
        #while(j > start):
            #j -= 1
        if(len(testdict) > 50):
            
            cutoff = g_max_scores - max(2.5*skipcost, 60.)
            for loc_in_one_mapinfo in testdict:
                if(S[loc_in_one_mapinfo] < cutoff and (one_mapinfo[i][0] - one_mapinfo[loc_in_one_mapinfo][0]>100)):
                    
                    testdict.pop(loc_in_one_mapinfo)
                    #else:
                        #break




               
        
        for j in testdict:
            
                    
            #(12240, 2791460820, 1, 21) (12251, 2791460711, 1, 108)
            nocost = False
            filtered = True
            if(one_mapinfo[i][0] == one_mapinfo[j][0]):
                continue

                      #      x_2 - x_1 - a_1
            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            if((readgap < 0)):
                continue
                #if((one_mapinfo[i][2] == one_mapinfo[j][2]) and ((readgap + one_mapinfo[i][3]) > 10)):
                    #pass
                #else:
                    #continue
            opcount += 1
            while(True):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    break

                if(one_mapinfo[i][2] == 1):
                              #    y_2 - y_1 - a_1
                    refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                else:
                             #     y_1 - y_2 - a_1
                    refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                
    
                if((one_mapinfo[i][2] == 1) and (refgap < 0)):
                    break
                if(one_mapinfo[i][2] == -1):
                    if(refgap > 0):
                        break
                    else:
                        refgap = abs(refgap)

                gapcost = abs(abs(readgap) - refgap)
                if(((readgap ) > maxgap) or (gapcost > maxdiff)):
                    break
                
                gapcost = gapcost_list[gapcost]
                #if(gapcost != 0):
                    #gapcost = 0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost)
                    #gapcost = 0.01 * kmersize * gapcost + 2.0 * np.log2(gapcost)
        
        
                test_scores = S[j] + (min(readgap, 0) + one_mapinfo[i][3]) - gapcost
                filtered = False
                break
            if(filtered == True):
                """if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(magic > abs(one_mapinfo[i][1] - one_mapinfo[j][1])):
                        continue"""
                test_scores = S[j] - skipcost + one_mapinfo[i][3]# - np.log2(readgap)# - np.log2(abs(abs(refgap) - readgap)+1)
                #continue
            if(test_scores >= max_scores):
                max_scores = test_scores
                pre_index = j

       
            
                
        S[i] = max_scores
        P[i] = pre_index 
        
        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
            
 
    
            

        testdict[i] = max_scores
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    while(True):
        if((P[take_index] == 9999999)):
            break
        take_index = P[take_index]
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

    #print_log(opcount)
    return g_max_scores, path#, testdict

@njit
def hit2wor1(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage = 0.5, hastra = False, H = False):
    def getsecond(x):
        return x[1]
    def getfirst(x):
        return x[0]
    def getlength(x):
        return len(x)

    def get_readloc_set_bin(one_mappos, bin_size):
        return set([i[0] // bin_size for i in one_mappos])
    def get_refloc_set_bin(one_mappos, bin_size):
        return set([i[1] // bin_size for i in one_mappos])
    def get_overlapsize(readloc_set_a, readloc_set_b):
        return len(readloc_set_a&readloc_set_b)/min(len(readloc_set_a), len(readloc_set_b))
    
    maxgap = 200
    
    one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:,1])] 

    cluster_list = List()
    iloc = 0
    preitem = one_mapinfo[iloc]
    st_iloc = 0
    for iloc in range(len(one_mapinfo)):
        nowitem = one_mapinfo[iloc]
        if(((nowitem[1] - preitem[1]) > c_bias)):
            if((iloc - st_iloc) < 3):
                continue
            cluster_list.append(one_mapinfo[st_iloc: iloc])
            st_iloc = iloc
        preitem = nowitem


    if((iloc - st_iloc) > 3):
        cluster_list.append(one_mapinfo[st_iloc: iloc + 1])
    cluster_list.sort(key = get_length)

    cluster_list = cluster_list[::-1][:check_num]

    
    hit = False
    minichain_scores = 40
    path_list = List()
    scores_list = []
    scores_list.append(0.)
    scores_list.pop()
    path_list.append([(0, 0, 0, 0)])
    path_list.pop()
    

    max_scores = 0
    from_repeat = False
    for one_mapinfo in cluster_list:
        repeat = False
        min_ref, max_ref = one_mapinfo[0][1], one_mapinfo[-1][1]
        testcontig = pos2contig(min_ref, contig2start)
        if(testcontig != pos2contig(max_ref, contig2start)):
            #print_log('hit2work: testcontig != pos2contig(max_ref, contig2start)')
            continue

        scores, path = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo[np.argsort(one_mapinfo[:,0])], kmersize = kmersize, skipcost = skipcost[0], maxdiff = maxdiff[0], maxgap = maxgap)
 
        
        if(scores < minichain_scores):
            continue
        if(H == False):
            readspan, refspan = max(abs(path[0][0] - path[-1][0]), 1), max(abs(path[0][1] - path[-1][1]), 1)

            if((max(readspan, refspan) / min(readspan, refspan)) > 4):
                repeat = True
                scores, path = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo[np.argsort(one_mapinfo[:,0])], kmersize = kmersize, skipcost = 1000, maxdiff = 200, maxgap = maxgap)
                #scores, path = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_linear(one_mapinfo[np.argsort(one_mapinfo[:,0])], kmersize = kmersize, skipcost = skipcost[0], maxdiff = 200, maxgap = 2000) 
            else:
                readspan_r, refspan_r = len(get_readloc_set_bin(path, 100))*100, len(get_refloc_set_bin(path, 100))*100

                if((readspan_r/readspan < 0.5) and (refspan_r/refspan < 0.5)):
                    repeat = True
                    scores, path = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo[np.argsort(one_mapinfo[:,0])], kmersize = kmersize, skipcost = 1000, maxdiff = 200, maxgap = maxgap)
                    #scores, path = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_linear(one_mapinfo[np.argsort(one_mapinfo[:,0])], kmersize = kmersize, skipcost = skipcost[0], maxdiff = 200, maxgap = 2000)
                else:
                    readspan, refspan = readspan_r, refspan_r

                    if((max(readspan, refspan) / min(readspan, refspan)) > 2):
                        repeat = True
                        scores, path = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo[np.argsort(one_mapinfo[:,0])], kmersize = kmersize, skipcost = 1000, maxdiff = 200, maxgap = maxgap)
                        #scores, path = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_linear(one_mapinfo[np.argsort(one_mapinfo[:,0])], kmersize = kmersize, skipcost = skipcost[0], maxdiff = 200, maxgap = 2000)



        hit = True

        path_list.append(path) 
        scores_list.append(scores)
        if(scores > max_scores):
            if(repeat == True):
                from_repeat = True
            else:
                from_repeat = False
            max_scores = scores
    
    if(hit == True):

        order = np.argsort(np.array(scores_list))[::-1]
        
        primary_rlocset_List = List()
        primary_scores_List = List()
        primary_index_List = List()
        
        iloc = order[0]
        primary_rlocset_List.append(get_readloc_set_bin(path_list[iloc], bin_size))
        primary_scores_List.append(List([scores_list[iloc]]))
        primary_index_List.append(iloc)

        for iloc in order[1:]:
            readloc_set_b = get_readloc_set_bin(path_list[iloc], bin_size)
            maxoverlapsize = 0.
            for p_loc in range(len(primary_rlocset_List)):
                tmp_overlapsize = get_overlapsize(primary_rlocset_List[p_loc], readloc_set_b)
                if(tmp_overlapsize > maxoverlapsize):
                    maxoverlapsize = tmp_overlapsize
                    prefer_p_loc = p_loc
            if(maxoverlapsize < overlapprecentage):
                primary_rlocset_List.append(readloc_set_b)
                primary_scores_List.append(List([scores_list[iloc]]))
                primary_index_List.append(iloc)
            else:
                primary_scores_List[prefer_p_loc].append(scores_list[iloc])
                
        m = len(path_list[order[0]])    
        if(len(primary_scores_List[0]) < 2):
            f1 = primary_scores_List[0][0]
            f2 = 0
        else:
            f1 = primary_scores_List[0][0]
            f2 = primary_scores_List[0][1]
        mapq = min(int(40*(1-f2/f1)*min(1, m/10)*np.log(f1)), 60)

        base_iloc = primary_index_List[0]

        for iloc in primary_index_List[1:4]:
            if(hastra == False):
                skip = True
                item = path_list[iloc][0]
                bias_1 = abs(item[1]-path_list[base_iloc][0][1])
                bias_2 = abs(item[1] - path_list[base_iloc][-1][1])
                if(min(bias_1, bias_2) < 200000):
                    if(pos2contig(item[1], contig2start) == pos2contig(path_list[base_iloc][0][1], contig2start)):
                        skip = False
                if(skip == True):
                    continue
            for item in path_list[iloc]:                            
                path_list[base_iloc].append(item)

        path_list[base_iloc].sort(key = getfirst)
        scores, path = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(np.array(path_list[base_iloc]), kmersize = kmersize, skipcost = skipcost[1], maxdiff = maxdiff[1], maxgap = maxgap)
        if(from_repeat == False):
            return mapq, scores, path
        else:
            return mapq, -scores, path

    else:
        return 0, 0., [(0, 0, 0, 0)]
import psutil
def check_func(raw_alignment_list, readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, debug = False):
    
    
    #step 1
    #rebuild chain break
    #remove anchor cause large cost and small alignment
    st = time.time()
    alignment_list = rebuild_chain_break(contig2start, raw_alignment_list, large_cost = setting_maxdiff, small_alignment = 50, small_dup = -100)
    #^print_log('step 1: rebuild chain break ', time.time() - st)
    max_alignment = 0
    size = 0
    for line in alignment_list:
        tempcontig = pos2contig(line[0][1], contig2start)
        temprefbias = contig2start[tempcontig]
        preitem, nowitem = line[0], line[-1]

        target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, setting_kmersize, contig2seq, contig2start)

        diffratio = edlib.align(query = query, target = target, task = 'distance')['editDistance']/min(len(target), len(query))
        size += len(query)
        if(diffratio < 0.2):
            if(max_alignment < min(len(target), len(query))):
                max_alignment = min(len(target), len(query))
    if(max_alignment > 400 and size > testseq_len/4):
        return True
    else:
        return False
    
def check_func_clean(raw_alignment_list, readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, debug = False):
    
    
    #step 1
    #rebuild chain break
    #remove anchor cause large cost and small alignment
    st = time.time()
    clean_alignment_list = []
    alignment_list = rebuild_chain_break(contig2start, raw_alignment_list, large_cost = setting_maxdiff, small_alignment = 50, small_dup = -100)
    #^print_log('step 1: rebuild chain break ', time.time() - st)
    max_alignment = 0
    size = 0
    for line in alignment_list:
        tempcontig = pos2contig(line[0][1], contig2start)
        temprefbias = contig2start[tempcontig]
        preitem, nowitem = line[0], line[-1]

        target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, setting_kmersize, contig2seq, contig2start)

        diffratio = edlib.align(query = query, target = target, task = 'distance')['editDistance']/min(len(target), len(query))
        size += len(query)
        if(diffratio < 0.2):
            for item in line:
                clean_alignment_list.append(item)
            if(max_alignment < min(len(target), len(query))):
                max_alignment = min(len(target), len(query))
    return clean_alignment_list[::-1]

def decode_hit(index_object, index2contig, testseq, testseq_len, contig2start, kmersize, skipcost = (50., 30.), maxdiff = (50, 30), maxgap = 2000, check_num = 20, c_bias = 5000, bin_size = 100, overlapprecentage = 0.5, hastra = True, H = False):            
    need_reverse, one_mapinfo = get_reversed_chain_numpy_rough(np.array(index_object.map(testseq)), testseq_len)
    mapq, scores, path =  hit2work(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage, hastra, H = H)
    if(need_reverse == True):
        return mapq, -scores, path
    else:
        return mapq, scores, path

def iftra(onemapinfolist, testseq_len, hitsize = 200):


    line = onemapinfolist[0]
    if(line[2] == '+'):
        readstart_1 = line[3]
        readend_1 = line[4]
    else:
        readstart_1 = testseq_len - line[4]
        readend_1 = testseq_len - line[3]
    refstart_1 = line[5]
    refend_1 = line[6]
    preline = line
    for line in onemapinfolist[1:]:
        if((preline[1] != line[1])):
            continue
        if(line[2] == '+'):
            readstart_2 = line[3]
            readend_2 = line[4]
        else:
            readstart_2 = testseq_len - line[4]
            readend_2 = testseq_len - line[3]
        refstart_2 = line[5]
        refend_2 = line[6]  
        
        readgap = readstart_2 - readend_1
        gap = min(abs(readgap - refstart_2 + refend_1), abs(readgap - refstart_1 + refend_2))
        
        if((readgap > hitsize) and (gap < 50)):
            return True
        readstart_1 = readstart_2
        readend_1 = readend_2
        refstart_1 = refstart_2
        refend_1 = refend_2
        preline = line
    return False
def getimtra(onemapinfolist):
    if(len(onemapinfolist) >= 3):
        for iloc in range(len(onemapinfolist)-2):

            strand_1 = onemapinfolist[iloc][2]
            if(strand_1 == '+'):
                readpos_1 = onemapinfolist[iloc][4]
                refpos_1 = onemapinfolist[iloc][6]
            else:
                readpos_1 = onemapinfolist[iloc][3]
                refpos_1 = onemapinfolist[iloc][5]
            for jloc in range(iloc + 2, len(onemapinfolist)):
                
                strand_2 = onemapinfolist[jloc][2]
                if((strand_1 == strand_2) and (onemapinfolist[iloc][1] == onemapinfolist[jloc][1]) and (abs(onemapinfolist[iloc][5] - onemapinfolist[jloc][5])>100000)):
 
                    if(strand_2 == '+'):
                        readpos_2 = onemapinfolist[jloc][4]
                        refpos_2 = onemapinfolist[jloc][6]
                        readgap = readpos_2 - readpos_1
                        if(readgap < 0):
                            continue
                        refgap = refpos_2 - refpos_1
                        if(refgap - readgap < -50):
                            return True
                    else:
                        readpos_2 = onemapinfolist[jloc][3]
                        refpos_2 = onemapinfolist[jloc][5]
                        readgap = readpos_1 - readpos_2
                        if(readgap < 0):
                            continue
                        refgap = refpos_1 - refpos_2
                        if(refgap - readgap < -50):
                            return True
        return False
    else:
        return False
def get_list_of_readmap1(raw_queue, savepath, index_object, contig2seq, hastra, H, header):
    st = time.time()
    
    contig2start = Dict()
    index2contig = List()
    contig2iloc = dict()
    iloc = -1
    for item in index_object.seq_offset:
        iloc += 1
        contig2start[item[0].decode()] = item[2]
        index2contig.append(item[0].decode())
        contig2iloc[item[0].decode()] = iloc

    
    iloc = 0
    unmapcountlist = []
    plotdata = []
    
    largeruntime = 0
    total_r_time = 0
    tmp_time = time.time()
    a_list = []
    large_rt = 0
    with pysam.AlignmentFile(savepath+'.bam', "wb", header=header) as outf:
        while(True):
            readidandseq = raw_queue.get()
            if(type(readidandseq) == int):
                break
            
            #print_log(iloc, readidandseq[0])
            r_st = time.time()
            try:

                if(H == True):
                    onemapinfolist, (alignment_list,raw_alignment_list), one_mapinfo = get_readmap_DP_test(readidandseq[0], readidandseq[1], contig2start, contig2seq, index_object, index2contig, hastra = hastra, H = True)
                else:
                    onemapinfolist, (alignment_list,raw_alignment_list), TRA_signal = get_readmap_DP_test(readidandseq[0], readidandseq[1], contig2start, contig2seq, index_object, index2contig, hastra = False)
                    if(TRA_signal == True):
                        onemapinfolist, (alignment_list,raw_alignment_list), TRA_signal = get_readmap_DP_test(readidandseq[0], readidandseq[1], contig2start, contig2seq, index_object, index2contig, hastra = True, check_num = 100)
            except:
                onemapinfolist = []
                unmapcountlist.append(readidandseq[0])
                total_r_time += time.time() - r_st
                if(time.time() - r_st > large_rt):
                    large_rt = time.time() - r_st
                    if(large_rt > 60):
                        print_log(readidandseq[0])
                continue
            if(len(onemapinfolist) == 0):
                onemapinfolist = []
                unmapcountlist.append(readidandseq[0])
                total_r_time += time.time() - r_st
                if(time.time() - r_st > large_rt):
                    large_rt = time.time() - r_st
                    if(large_rt > 60):
                        print_log(readidandseq[0])
                continue
            total_r_time += time.time() - r_st
            if(time.time() - r_st > large_rt):
                large_rt = time.time() - r_st
                if(large_rt > 60):
                    print_log(readidandseq[0])
            if(len(onemapinfolist) != 0):
                tmp_a_list = get_bam_dict(onemapinfolist, readidandseq[1], readidandseq[2], contig2iloc, contig2seq)
                if((tmp_a_list) == None):
                    print_log(readidandseq[0])
                    print_log()
                else:
                    for a in tmp_a_list:
                        outf.write(a)

            else:
                unmapcountlist.append(readidandseq[0])


    print_log(total_r_time, time.time() - st, large_rt, 'unmapsize', len(unmapcountlist))
    print_log('unmapsize', len(unmapcountlist))
    #print_log(unmapcountlist)
    #print_log(unmapcountlist)
    print_log()




    
    
    




@njit
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = [0]
    for gapcost in range(1, maxdiff + 1):
        gapcost_list.append(0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))


    n = len(one_mapinfo)
    S = np.zeros(n, np.float64)
    P = np.zeros(n, np.int64)
    start = 0
    
    testdict = dict()
    testdict[-1] = 0.
    testdict.pop(-1)
    maxop = 2000*(one_mapinfo[-1][0]-one_mapinfo[0][0])
    #print_log(maxop)
    opcount = 0
    
    coverage_dict = dict()
    for i in range(n):
        if(one_mapinfo[i][0] in coverage_dict):
            coverage_dict[one_mapinfo[i][0]] += 1
        else:
            coverage_dict[one_mapinfo[i][0]] = 1

    #coverage_array *= 2
    
    oskipcost = skipcost
    omaxdiff = maxdiff
    
    for i in range(n):
        P[i] = 9999999
        max_scores = one_mapinfo[i][3]
        pre_index = 9999999
        if(opcount > 2*maxop):
            path = []
            take_index = g_max_index
            path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
            while(True):
                if((P[take_index] == 9999999)):
                    break
                take_index = P[take_index]
                path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

            return g_max_scores, path
        if(opcount > maxop):
            if(max_scores < 15):
                S[i] = max_scores
                continue
        #while(((i - start) > 100) and ((one_mapinfo[i][0] - one_mapinfo[start][0]) > max_gap)):
            #start += 1
            
        skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
        maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)    
            
        j = i
        skipcount = 0
        #while(j > start):
            #j -= 1
        if(len(testdict) > 50):
            
            cutoff = g_max_scores - max(2.5*skipcost, 60.)

            tmp_poplist = List([0])
            for loc_in_one_mapinfo in testdict:
                if(S[loc_in_one_mapinfo] < cutoff and (one_mapinfo[i][0] - one_mapinfo[loc_in_one_mapinfo][0]>100)):
                    tmp_poplist.append(loc_in_one_mapinfo)

            if(len(tmp_poplist) != 1):
                tmp_poplist.pop(0)
                for loc_in_one_mapinfo in tmp_poplist:
                    testdict.pop(loc_in_one_mapinfo)




               
        
        for j in testdict:
            
                    
            #(12240, 2791460820, 1, 21) (12251, 2791460711, 1, 108)
            nocost = False
            filtered = True
            if(one_mapinfo[i][0] == one_mapinfo[j][0]):
                continue

                      #      x_2 - x_1 - a_1
            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            if((readgap < 0)):
                continue
                #if((one_mapinfo[i][2] == one_mapinfo[j][2]) and ((readgap + one_mapinfo[i][3]) > 10)):
                    #pass
                #else:
                    #continue
            opcount += 1
            while(True):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    break

                if(one_mapinfo[i][2] == 1):
                              #    y_2 - y_1 - a_1
                    refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                else:
                             #     y_1 - y_2 - a_1
                    refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                
    
                if((one_mapinfo[i][2] == 1) and (refgap < 0)):
                    break
                if(one_mapinfo[i][2] == -1):
                    if(refgap > 0):
                        break
                    else:
                        refgap = abs(refgap)

                gapcost = abs(abs(readgap) - refgap)
                if(((readgap ) > maxgap) or (gapcost > maxdiff)):
                    break
                
                gapcost = gapcost_list[gapcost]
                #if(gapcost != 0):
                    #gapcost = 0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost)
                    #gapcost = 0.01 * kmersize * gapcost + 2.0 * np.log2(gapcost)
        
        
                test_scores = S[j] + (min(readgap, 0) + one_mapinfo[i][3]) - gapcost
                filtered = False
                break
            if(filtered == True):
                """if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(magic > abs(one_mapinfo[i][1] - one_mapinfo[j][1])):
                        continue"""
                test_scores = S[j] - skipcost + one_mapinfo[i][3]# - np.log2(readgap)# - np.log2(abs(abs(refgap) - readgap)+1)
                #continue
            if(test_scores >= max_scores):
                max_scores = test_scores
                pre_index = j

       
            
                
        S[i] = max_scores
        P[i] = pre_index 
        
        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
            
 
    
            

        testdict[i] = max_scores
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    while(True):
        if((P[take_index] == 9999999)):
            break
        take_index = P[take_index]
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

    #print_log(opcount)
    return g_max_scores, path#, testdict

@njit
def hit2work(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage = 0.5, hastra = False, H = False):
    def getsecond(x):
        return x[1]
    def getfirst(x):
        return x[0]
    def getlength(x):
        return len(x)

    def get_readloc_set_bin(one_mappos, bin_size):
        return set([i[0] // bin_size for i in one_mappos])
    def get_refloc_set_bin(one_mappos, bin_size):
        return set([i[1] // bin_size for i in one_mappos])
    def get_overlapsize(readloc_set_a, readloc_set_b):
        return len(readloc_set_a&readloc_set_b)/min(len(readloc_set_a), len(readloc_set_b))
    
    maxgap = 200

    
    one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:,1])] 

    cluster_list = List()
    iloc = 0
    preitem = one_mapinfo[iloc]
    st_iloc = 0
    for iloc in range(len(one_mapinfo)):
        nowitem = one_mapinfo[iloc]
        if(((nowitem[1] - preitem[1]) > c_bias)):
            if((iloc - st_iloc) < 3):
                continue
            cluster_list.append(one_mapinfo[st_iloc: iloc])
            st_iloc = iloc
        preitem = nowitem


    if((iloc - st_iloc) > 3):
        cluster_list.append(one_mapinfo[st_iloc: iloc + 1])
    cluster_list.sort(key = get_length)

    cluster_list = cluster_list[::-1][:check_num]

    
    hit = False
    minichain_scores = 40
    path_list = List()
    scores_list = []
    scores_list.append(0.)
    scores_list.pop()
    path_list.append([(0, 0, 0, 0)])
    path_list.pop()
    

    max_scores = 0
    from_repeat = False
    for one_mapinfo in cluster_list:
        repeat = False
        min_ref, max_ref = one_mapinfo[0][1], one_mapinfo[-1][1]
        testcontig = pos2contig(min_ref, contig2start)
        if(testcontig != pos2contig(max_ref, contig2start)):
            #print_log('hit2work: testcontig != pos2contig(max_ref, contig2start)')
            continue

        scores, path = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo[np.argsort(one_mapinfo[:,0])], kmersize = kmersize, skipcost = skipcost[0], maxdiff = maxdiff[0], maxgap = maxgap)
 
        
        if(scores < minichain_scores):
            continue
        


        hit = True

        path_list.append(path) 
        scores_list.append(scores)
        if(scores > max_scores):
            if(repeat == True):
                from_repeat = True
            else:
                from_repeat = False
            max_scores = scores
    
    if(hit == True and max_scores > 100):

        order = np.argsort(np.array(scores_list))[::-1]
        
        primary_rlocset_List = List()
        primary_scores_List = List()
        primary_index_List = List()
        
        iloc = order[0]
        primary_rlocset_List.append(get_readloc_set_bin(path_list[iloc], bin_size))
        primary_scores_List.append(List([scores_list[iloc]]))
        primary_index_List.append(iloc)

        for iloc in order[1:]:
            readloc_set_b = get_readloc_set_bin(path_list[iloc], bin_size)
            maxoverlapsize = 0.
            for p_loc in range(len(primary_rlocset_List)):
                tmp_overlapsize = get_overlapsize(primary_rlocset_List[p_loc], readloc_set_b)
                if(tmp_overlapsize > maxoverlapsize):
                    maxoverlapsize = tmp_overlapsize
                    prefer_p_loc = p_loc
            if(maxoverlapsize < overlapprecentage):
                primary_rlocset_List.append(readloc_set_b)
                primary_scores_List.append(List([scores_list[iloc]]))
                primary_index_List.append(iloc)
            else:
                primary_scores_List[prefer_p_loc].append(scores_list[iloc])
                
        m = len(path_list[order[0]])    
        if(len(primary_scores_List[0]) < 2):
            f1 = primary_scores_List[0][0]
            f2 = 0
        else:
            f1 = primary_scores_List[0][0]
            f2 = primary_scores_List[0][1]
        mapq = min(int(40*(1-f2/f1)*min(1, m/10)*np.log(f1)), 60)

        base_iloc = primary_index_List[0]

        for iloc in primary_index_List[1:4]:
            if(hastra == False):
                skip = True
                item = path_list[iloc][0]
                bias_1 = abs(item[1]-path_list[base_iloc][0][1])
                bias_2 = abs(item[1] - path_list[base_iloc][-1][1])
                if(min(bias_1, bias_2) < 200000):
                    if(pos2contig(item[1], contig2start) == pos2contig(path_list[base_iloc][0][1], contig2start)):
                        skip = False
                if(skip == True):
                    continue
            for item in path_list[iloc]:                            
                path_list[base_iloc].append(item)

        path_list[base_iloc].sort(key = getfirst)
        scores, path = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(np.array(path_list[base_iloc]), kmersize = kmersize, skipcost = skipcost[1], maxdiff = maxdiff[1], maxgap = maxgap)

        return mapq, scores, path


    else:
        return 0, 0., [(0, 0, 0, 0)]
    
@njit
def get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start):#
    testcontig = pos2contig(preitem[1], contig2start)
    refbias = contig2start[testcontig]
    query_st, query_en = preitem[0], nowitem[0]
    if(preitem[2] == 1):
        target_st, target_en = preitem[1], nowitem[1]
        query = testseq[preitem[0]: nowitem[0]]
        target = contig2seq[testcontig][preitem[1] - refbias: nowitem[1] - refbias]

    else:
        target_st, target_en = nowitem[1] + nowitem[3], preitem[1] + preitem[3]
        query = rc_testseq[testseq_len - nowitem[0]: testseq_len - preitem[0]]
        target = contig2seq[testcontig][nowitem[1] + nowitem[3] - refbias: preitem[1] + preitem[3] - refbias]
    return target, query, target_st, target_en, query_st, query_en


def extend_edge_drop_test(testseq, testseq_len, alignment_list, kmersize, pos2contig, contig2start, contig2seq, san, debug):#



    zdrop_value = 400  

    mismatch_value = -4

    gap_open_extend = -4
    gap_extend_extend = -4
    zdrop_value_extend = 50
    bw = 100
    droppercentage = 0.5




    max_extend_size = 20000
    if(debug == True): print_log('Drop mismatch alignment by Extend')
    if(debug == True): print_log('len(alignment_list): ', len(alignment_list))
    onealignment_index = -1
    while(True):
        onealignment_index += 1
        if(debug == True): print_log('onealignment_index', onealignment_index)
        if((onealignment_index+2) >= len(alignment_list)):
            break
        top_alignment = alignment_list[onealignment_index]
        tail_alignment = alignment_list[onealignment_index+1]
        top_alignment_size = top_alignment[-1][0] + top_alignment[-1][3] - top_alignment[0][0]
        tail_alignment_size = tail_alignment[-1][0] + tail_alignment[-1][3] - tail_alignment[0][0]
        if(min(top_alignment_size, tail_alignment_size) > 1000):
            continue
        if(top_alignment_size > tail_alignment_size):
            if(debug): print_log('top_alignment_size large')

            looksize = tail_alignment[-1][0] + tail_alignment[-1][3] + 200 - (top_alignment[-1][0] + top_alignment[-1][3])
            
            looksize = min(max_extend_size, looksize)
            
            if(looksize < 0):
                print_log('extend_edge_drop_test: error looksize < 0')



            preitem = top_alignment[-2]
            nowitem = top_alignment[-1]
            if(debug): print_log('enlong top')
            if(debug): print_log(preitem)
            if(debug): print_log(nowitem)

            #target, query, target_st, target_en, query_st, query_en = get_target_query_for_drop(preitem, nowitem, testseq, kmersize, pos2contig, contig2start, contig2seq, reverse = False, rc = False)

            testcontig = pos2contig(preitem[1], contig2start)

            if(preitem[2] == 1):
                
                target_en = nowitem[1] + nowitem[3]
                query_en = nowitem[0] + nowitem[3]
                

                #looksize = min(looksize, contig2start[testcontig] + len(contig2seq[testcontig]) - target_en)
                looksize = min(looksize, get_contig2start(testcontig, contig2start) + get_reflen(testcontig, contig2seq) - (target_en - 1))


                if(looksize == 0):
                    if(debug): print_log('looksize == 0, skiped')
                    continue
                query = testseq[query_en: query_en + looksize]

                #target = contig2seq[testcontig][target_en-contig2start[testcontig]: target_en+len(query)-contig2start[testcontig]]
                target = get_refseq(testcontig, target_en-get_contig2start(testcontig, contig2start), target_en+len(query)-get_contig2start(testcontig, contig2start), contig2seq)

                if(debug): print_log(len(target), len(query), target_en, query_en)

                #zdroped, max_extend, cigarstring, Score = fast_globalms_align_extend(target, query, match = 2, mismatch = -4, gap_open = gap_open_extend, gap_extend = gap_extend_extend, zdrop = zdrop_value_extend, loc_only = True)
                cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 4, gap_open_2 = 4, gap_extend_2 = 4, bw=bw, zdropvalue=zdrop_value_extend)

                max_extend = (t_e, q_e)
                extend_size = query_en+q_e - tail_alignment[0][0]
                if(debug == True):
 
                    print_log('(extend_size > tail_alignment_size), (extend_size/tail_alignment_size>0.5)', (extend_size > tail_alignment_size), (extend_size/tail_alignment_size>droppercentage))
                    print_log('extend_size, tail_alignment_size', extend_size, tail_alignment_size)
                if((extend_size > tail_alignment_size) or ((extend_size/tail_alignment_size)>droppercentage)):
                    if(debug == True):
                        print_log('Alignment removed', alignment_list[onealignment_index+1][0][0], alignment_list[onealignment_index+1][-1][0]+alignment_list[onealignment_index+1][-1][3])
                    #alignment_list[onealignment_index][-1] = ((query_en + max_extend[1], target_en + max_extend[0], 1, 0))
                    alignment_list.pop(onealignment_index+1)
                    onealignment_index -= 1
            else:
                
                target_st = nowitem[1]
                query_en = nowitem[0] + nowitem[3]
                #looksize = min(looksize, target_st - contig2start[testcontig])
                looksize = min(looksize, target_st - get_contig2start(testcontig, contig2start))

                if(looksize == 0):
                    if(debug): print_log('looksize == 0, skiped')
                    continue
                query = testseq[query_en: query_en + looksize]

                #target = get_reverse_complement(contig2seq[testcontig][target_st-contig2start[testcontig] - len(query): target_st-contig2start[testcontig]])
                target = str(Seq(get_refseq(testcontig, target_st-get_contig2start(testcontig, contig2start) - len(query), target_st-get_contig2start(testcontig, contig2start), contig2seq)).reverse_complement())

                if(debug): print_log(len(target), len(query), target_st, query_en)

                #zdroped, max_extend, cigarstring, Score = fast_globalms_align_extend(target, query, match = 2, mismatch = -4, gap_open = gap_open_extend, gap_extend = gap_extend_extend, zdrop = zdrop_value_extend, loc_only = True)
                cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 4, gap_open_2 = 4, gap_extend_2 = 4, bw=bw, zdropvalue=zdrop_value_extend)
                max_extend = (t_e, q_e)
                extend_size = query_en+q_e - tail_alignment[0][0]
                if(debug == True):
 
                    print_log('(extend_size > tail_alignment_size), (extend_size/tail_alignment_size>0.5)', (extend_size > tail_alignment_size), (extend_size/tail_alignment_size>droppercentage))
                    print_log('extend_size, tail_alignment_size', extend_size, tail_alignment_size)
                if((extend_size > tail_alignment_size) or ((extend_size/tail_alignment_size)>droppercentage)):
                    if(debug == True):
                        print_log('Alignment removed', alignment_list[onealignment_index+1][0][0], alignment_list[onealignment_index+1][-1][0]+alignment_list[onealignment_index+1][-1][3])
                    #alignment_list[onealignment_index][-1] = ((query_en + max_extend[1], target_st - max_extend[0], 1, 0))
                    alignment_list.pop(onealignment_index+1)
                    onealignment_index -= 1

            if(debug): print_log(max_extend)
            if(debug): print_log(top_alignment[-1])
            #if(debug): print_log(target)
            #if(debug): print_log(query)
        else:
            if(debug): print_log('tail_alignment_size large')
            looksize = tail_alignment[0][0] - top_alignment[0][0] + 200

            looksize = min(max_extend_size, looksize)

            if(looksize < 0):

                print_log('extend_edge_drop_test: error looksize < 0')


            preitem = tail_alignment[0]
            nowitem = tail_alignment[1]
            if(debug): print_log('enlong tail')
            if(debug): print_log(preitem)
            if(debug): print_log(nowitem)

            #target, query, target_st, target_en, query_st, query_en = get_target_query_for_drop(preitem, nowitem, testseq, kmersize, pos2contig, contig2start, contig2seq, reverse = False, rc = False)

            testcontig = pos2contig(preitem[1], contig2start)

            if(preitem[2] == 1):
                
                target_st = preitem[1]
                query_st = preitem[0]
                
                #looksize = min(looksize, target_st - contig2start[testcontig])
                looksize = min(looksize, target_st - get_contig2start(testcontig, contig2start))


                if(looksize == 0):
                    if(debug): print_log('looksize == 0, skiped')
                    continue
                query = testseq[max(query_st - looksize, 0): query_st][::-1]

                #target = contig2seq[testcontig][target_st-contig2start[testcontig] - len(query): target_st-contig2start[testcontig]][::-1]
                target = get_refseq(testcontig, target_st-get_contig2start(testcontig, contig2start) - len(query), target_st-get_contig2start(testcontig, contig2start), contig2seq)[::-1]

                if(debug): print_log(len(target), len(query), target_st, query_st)


                #zdroped, max_extend, cigarstring, Score = fast_globalms_align_extend(target, query, match = 2, mismatch = -4, gap_open = gap_open_extend, gap_extend = gap_extend_extend, zdrop = zdrop_value_extend, loc_only = True)
                cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 4, gap_open_2 = 4, gap_extend_2 = 4, bw=bw, zdropvalue=zdrop_value_extend)


                max_extend = (t_e, q_e)
                extend_size = top_alignment[-1][0] + top_alignment[-1][3] - (query_st - q_e)
                if(debug == True):

                    print_log('(extend_size > top_alignment_size), (extend_size/top_alignment_size>0.5)', (extend_size > top_alignment_size), (extend_size/top_alignment_size>droppercentage))
                    print_log('extend_size, tail_alignment_size', extend_size, top_alignment_size)
                if((extend_size > top_alignment_size) or ((extend_size / top_alignment_size) > droppercentage)):
                    if(debug == True):
                        print_log('Alignment removed', alignment_list[onealignment_index][0][0], alignment_list[onealignment_index][-1][0]+alignment_list[onealignment_index][-1][3])
                    #alignment_list[onealignment_index+1][0] = ((query_st - max_extend[1], target_st - max_extend[0], 1, 0))
                    alignment_list.pop(onealignment_index)
                    onealignment_index -= 1
                    
            else:
                
                target_en = preitem[1] + preitem[3]
                query_st = preitem[0]
                
                #looksize = min(looksize, contig2start[testcontig] + len(contig2seq[testcontig]) - target_en)
                looksize = min(looksize, get_contig2start(testcontig, contig2start) + get_reflen(testcontig, contig2seq) - (target_en - 1))

                if(looksize == 0):
                    if(debug): print_log('looksize == 0, skiped')
                    continue
                query = testseq[max(query_st - looksize, 0): query_st][::-1]

                #target = get_reverse_complement(contig2seq[testcontig][target_en-contig2start[testcontig]: target_en + len(query)-contig2start[testcontig]])[::-1]
                target = str(Seq(get_refseq(testcontig, target_en-get_contig2start(testcontig, contig2start), target_en + len(query)-get_contig2start(testcontig, contig2start), contig2seq)).reverse_complement())[::-1]

                if(debug): print_log(len(target), len(query), target_en, query_st)


                #zdroped, max_extend, cigarstring, Score = fast_globalms_align_extend(target, query, match = 2, mismatch = -4, gap_open = gap_open_extend, gap_extend = gap_extend_extend, zdrop = zdrop_value_extend, loc_only = True)
                cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 4, gap_open_2 = 4, gap_extend_2 = 4, bw=bw, zdropvalue=zdrop_value_extend)


                max_extend = (t_e, q_e)
                extend_size = top_alignment[-1][0] + top_alignment[-1][3] - (query_st - q_e)
                if(debug == True):

                    print_log('(extend_size > top_alignment_size), (extend_size/top_alignment_size>0.5)', (extend_size > top_alignment_size), (extend_size/top_alignment_size>droppercentage))
                    print_log('extend_size, tail_alignment_size', extend_size, top_alignment_size)
                if((extend_size > top_alignment_size) or ((extend_size / top_alignment_size) > droppercentage)):
                    if(debug == True):
                        print_log('Alignment removed', alignment_list[onealignment_index][0][0], alignment_list[onealignment_index][-1][0]+alignment_list[onealignment_index][-1][3])
                    #alignment_list[onealignment_index+1][0] = ((query_st - max_extend[1], target_en + max_extend[0], -1, 0))
                    alignment_list.pop(onealignment_index)
                    onealignment_index -= 1
            if(debug): print_log(max_extend)
            if(debug): print_log(tail_alignment[0])
            #if(debug): print_log(target)
            #if(debug): print_log(query)

                    
      

            
        if(debug): print_log()

    if(debug): print_log()
        
def extend_edge_test(testseq, testseq_len, alignment_list, kmersize, pos2contig, contig2start, contig2seq, san, debug):#



    zdrop_value = 400  

    mismatch_value = -4

    gap_open_extend = -4
    gap_extend_extend = -4
    zdrop_value_extend = 50
    bw = 100




    max_extend_size = 20000
    if(debug == True): print_log('Extend')
    if(debug == True): print_log('len(alignment_list): ', len(alignment_list))
    onealignment_index = -1
    while(True):
        onealignment_index += 1
        if(debug == True): print_log('onealignment_index', onealignment_index)
        if(onealignment_index >= len(alignment_list)):
            break
        onealignment = alignment_list[onealignment_index]
        if(debug == True): print_log('Alignment: start end ', onealignment[0][0], onealignment[-1][0] + onealignment[-1][3])
        if(debug): print_log(onealignment[0], onealignment[-1])
        if(onealignment[0][0]>0):
            while(True):

                if(onealignment[0][0] == 0):
                    break
                pre_onealignment_index = max(onealignment_index - san, 0)
                if(onealignment_index == 0 or  onealignment_index - san < 0):
                    looksize = onealignment[0][0] - 0
                else:
                    looksize = onealignment[0][0] - (alignment_list[pre_onealignment_index][-1][0] + alignment_list[pre_onealignment_index][-1][3])

                if(looksize < 0):

                    print_log('error')
                    print_log(alignment_list[pre_onealignment_index][-1], onealignment[0])


                preitem = onealignment[0]
                nowitem = onealignment[1]
                if(debug): print_log('enlong start')
                if(debug): print_log(preitem)
                if(debug): print_log(nowitem)

                #target, query, target_st, target_en, query_st, query_en = get_target_query_for_drop(preitem, nowitem, testseq, kmersize, pos2contig, contig2start, contig2seq, reverse = False, rc = False)

                testcontig = pos2contig(preitem[1], contig2start)

                if(preitem[2] == 1):
                    
                    target_st = preitem[1]
                    query_st = preitem[0]
                    #looksize = min(looksize, target_st - contig2start[testcontig])
                    looksize = min(looksize, target_st - get_contig2start(testcontig, contig2start))
                    
                    if(looksize > max_extend_size):
                        looksize = max_extend_size
                        skip_loop = False
                    else:
                        skip_loop = True
                    if(looksize == 0):
                        if(debug): print_log('looksize == 0, skiped')
                        break
                    query = testseq[max(query_st - looksize, 0): query_st][::-1]
                    
                    #target = contig2seq[testcontig][target_st-contig2start[testcontig] - len(query): target_st-contig2start[testcontig]][::-1]
                    target = get_refseq(testcontig, target_st-get_contig2start(testcontig, contig2start) - len(query), target_st-get_contig2start(testcontig, contig2start), contig2seq)[::-1]
                    
                    if(debug): print_log(len(target), len(query), target_st,  query_st)


                    #zdroped, max_extend, cigarstring, Score = fast_globalms_align_extend(target, query, match = 2, mismatch = -4, gap_open = gap_open_extend, gap_extend = gap_extend_extend, zdrop = zdrop_value_extend, loc_only = True)
                    cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 4, gap_open_2 = 4, gap_extend_2 = 4, bw=bw, zdropvalue=zdrop_value_extend)
          
                    
                    max_extend = (t_e, q_e)
                    onealignment[0] = ((query_st - max_extend[1], target_st - max_extend[0], 1, 0))
                else:
                    
                    #looksize = min(looksize, contig2start[testcontig] + len(contig2seq[testcontig]) - target_en)
                    target_en = preitem[1] + preitem[3]
                    query_st = preitem[0]
                    looksize = min(looksize, get_contig2start(testcontig, contig2start) + get_reflen(testcontig, contig2seq) - (target_en - 1) )
                    
                    if(looksize > max_extend_size):
                        looksize = max_extend_size
                        skip_loop = False
                    else:
                        skip_loop = True
                    if(looksize == 0):
                        if(debug): print_log('looksize == 0, skiped')
                        break
                    query = testseq[max(query_st - looksize, 0): query_st][::-1]
                    
                    #target = get_reverse_complement(contig2seq[testcontig][target_en-contig2start[testcontig]: target_en + len(query)-contig2start[testcontig]])[::-1]
                    target = str(Seq(get_refseq(testcontig, target_en-get_contig2start(testcontig, contig2start), target_en + len(query)-get_contig2start(testcontig, contig2start), contig2seq)).reverse_complement())[::-1]
                    
                    if(debug): print_log(len(target), len(query),  target_en, query_st)
                        
  
                    #zdroped, max_extend, cigarstring, Score = fast_globalms_align_extend(target, query, match = 2, mismatch = -4, gap_open = gap_open_extend, gap_extend = gap_extend_extend, zdrop = zdrop_value_extend, loc_only = True)
                    cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 4, gap_open_2 = 4, gap_extend_2 = 4, bw=bw, zdropvalue=zdrop_value_extend)

                    
                    max_extend = (t_e, q_e)
                    onealignment[0] = ((query_st - max_extend[1], target_en + max_extend[0], -1, 0))
                if(debug): print_log(max_extend)
                if(debug): print_log(onealignment[0])
                if(debug): print_log(target)
                if(debug): print_log(query)


                break
        else:
            tmpitem = alignment_list[onealignment_index][0]
            if(alignment_list[onealignment_index][0][2] == 1):
                alignment_list[onealignment_index][0] = ((tmpitem[0], tmpitem[1], 1, 0))
            else:
                alignment_list[onealignment_index][0] = ((tmpitem[0], tmpitem[1] + tmpitem[3], -1, 0))
            
        if((onealignment[-1][0]+onealignment[-1][3])<len(testseq)):
            while(True):
                if((onealignment[-1][0]+onealignment[-1][3])>=len(testseq)):
                    break
                pre_onealignment_index = min(onealignment_index + san, len(alignment_list))
                if(pre_onealignment_index == len(alignment_list)):
                    looksize = testseq_len - (alignment_list[onealignment_index][-1][0] + alignment_list[onealignment_index][-1][3])
                else:
                    looksize = alignment_list[pre_onealignment_index][0][0] - (onealignment[-1][0] + onealignment[-1][3])
                if(looksize < 0):
                    print_log('error')
                    print_log(onealignment[-1], alignment_list[pre_onealignment_index][0])



                preitem = onealignment[-2]
                nowitem = onealignment[-1]
                if(debug): print_log('enlong end')
                if(debug): print_log(preitem)
                if(debug): print_log(nowitem)

                #target, query, target_st, target_en, query_st, query_en = get_target_query_for_drop(preitem, nowitem, testseq, kmersize, pos2contig, contig2start, contig2seq, reverse = False, rc = False)

                testcontig = pos2contig(preitem[1], contig2start)

                if(preitem[2] == 1):
                    
                    #looksize = min(looksize, contig2start[testcontig] + len(contig2seq[testcontig]) - target_en)
                    target_en = nowitem[1] + nowitem[3]
                    query_en = nowitem[0] + nowitem[3]
                    looksize = min(looksize, get_contig2start(testcontig, contig2start) + get_reflen(testcontig, contig2seq) - (target_en - 1))
                    
                    if(looksize > max_extend_size):
                        looksize = max_extend_size
                        skip_loop = False
                    else:
                        skip_loop = True
                    if(looksize == 0):
                        if(debug): print_log('looksize == 0, skiped')
                        break
                    query = testseq[query_en: query_en + looksize]
                    
                    #target = contig2seq[testcontig][target_en-contig2start[testcontig]: target_en+len(query)-contig2start[testcontig]]
                    target = get_refseq(testcontig, target_en-get_contig2start(testcontig, contig2start), target_en+len(query)-get_contig2start(testcontig, contig2start), contig2seq)
                    
                    if(debug): print_log(len(target), len(query),  target_en, query_en)

                    #zdroped, max_extend, cigarstring, Score = fast_globalms_align_extend(target, query, match = 2, mismatch = -4, gap_open = gap_open_extend, gap_extend = gap_extend_extend, zdrop = zdrop_value_extend, loc_only = True)
                    cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 4, gap_open_2 = 4, gap_extend_2 = 4, bw=bw, zdropvalue=zdrop_value_extend)
                    
                    max_extend = (t_e, q_e)
                    
                    onealignment[-1] = ((query_en + max_extend[1], target_en + max_extend[0], 1, 0))
                else:
                    
                    target_st = nowitem[1]
                    query_en = nowitem[0] + nowitem[3]
                    #looksize = min(looksize, target_st - contig2start[testcontig])
                    looksize = min(looksize, target_st - get_contig2start(testcontig, contig2start))
                    
                    if(looksize > max_extend_size):
                        looksize = max_extend_size
                        skip_loop = False
                    else:
                        skip_loop = True
                    if(looksize == 0):
                        if(debug): print_log('looksize == 0, skiped')
                        break
                    query = testseq[query_en: query_en + looksize]
                    
                    #target = get_reverse_complement(contig2seq[testcontig][target_st-contig2start[testcontig] - len(query): target_st-contig2start[testcontig]])
                    target = str(Seq(get_refseq(testcontig, target_st-get_contig2start(testcontig, contig2start) - len(query), target_st-get_contig2start(testcontig, contig2start), contig2seq)).reverse_complement())
                    
                    if(debug): print_log(len(target), len(query), target_st, query_en)

                    #zdroped, max_extend, cigarstring, Score = fast_globalms_align_extend(target, query, match = 2, mismatch = -4, gap_open = gap_open_extend, gap_extend = gap_extend_extend, zdrop = zdrop_value_extend, loc_only = True)
                    cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 4, gap_open_2 = 4, gap_extend_2 = 4, bw=bw, zdropvalue=zdrop_value_extend)
                    max_extend = (t_e, q_e)
                    onealignment[-1] = ((query_en + max_extend[1], target_st - max_extend[0], -1, 0))
                if(debug): print_log(max_extend)
                if(debug): print_log(onealignment[-1])
                if(debug): print_log(target)
                if(debug): print_log(query)

                break
        else:
            tmpitem = alignment_list[onealignment_index][-1]
            if(alignment_list[onealignment_index][-1][2] == 1):
                alignment_list[onealignment_index][-1] = ((tmpitem[0] + tmpitem[3], tmpitem[1] + tmpitem[3], 1, 0))
            else:
                alignment_list[onealignment_index][-1] = ((tmpitem[0] + tmpitem[3], tmpitem[1], -1, 0))


            
        if(debug): print_log()

    if(debug): print_log()
        
def check_drop_site(preitem, nowitem, testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start):#
    gap_open_extend = -4
    gap_extend_extend = -4
    zdrop_value_extend = 50
    if(preitem[2] == 1): #preitem readloc low
        target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start)
        
        cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 2, gap_open_2 = 24, gap_extend_2 = 1, bw=500, zdropvalue=zdrop_value_extend)
        max_extend = (t_e, q_e)
        preitem_extend = (query_st + max_extend[1], target_st + max_extend[0], 1, 0)
        
        cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target[::-1], query[::-1], match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 2, gap_open_2 = 24, gap_extend_2 = 1, bw=500, zdropvalue=zdrop_value_extend)
        max_extend = (t_e, q_e)
        nowitem_extend = (query_en - max_extend[1], target_en - max_extend[0], 1, 0)

        gap_size = nowitem_extend[0] - preitem_extend[0]



    else:#nowitem readloc low
        target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(nowitem, preitem, testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start)
        
        cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 2, gap_open_2 = 24, gap_extend_2 = 1, bw=500, zdropvalue=zdrop_value_extend)
        max_extend = (t_e, q_e)
        preitem_extend = (query_en - max_extend[1], target_st + max_extend[0], -1, 0)
        
        cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target[::-1], query[::-1], match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 2, gap_open_2 = 24, gap_extend_2 = 1, bw=500, zdropvalue=zdrop_value_extend)
        max_extend = (t_e, q_e)
        nowitem_extend = (query_st + max_extend[1], target_en - max_extend[0], -1, 0)


        gap_size = preitem_extend[0] - nowitem_extend[0]

    
    return [preitem_extend, nowitem_extend]

def split_alignment_test(alignment, testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start, debug, H = False, eqx = False):#no big cost >30
    #print_log(alignment_list)
    
    cigartime = 0
    droptime = 0
    testtime = 0
    splittime = 0
    
    cigarlist = []

    gap_open_extend = -4
    gap_extend_extend = -4
    zdrop_value_extend = 100


    gap_open_drop = -4
    gap_extend_drop = -2
    zdrop_value_drop = 400 #{clr,ccs}:200, {nano}:400
    merge_smallgap = 300
    
    min_gap_forcigar = 200

    new_alignment = List()
    cigarlist = []
    preDEL, preINS = False, False
    if(alignment[0][2] == 1):
        iloc = -1
        if(alignment[iloc][3] != 0):
            alignment[iloc] = (alignment[iloc][0]+alignment[iloc][3], alignment[iloc][1]+alignment[iloc][3], 1, 0)
        nowitem = alignment[0]
        new_alignment.append(List([nowitem]))
        preitem = nowitem
        cigarlist.append([])
        iloc = 1
        while(iloc < len(alignment)):
            nowitem = alignment[iloc]
            readgap = nowitem[0] - preitem[0] - preitem[3]
            refgap = nowitem[1] - preitem[1] - preitem[3]

            if((nowitem[3] < 19) or (min(readgap, refgap) < min_gap_forcigar)):
                if(iloc +1 != len(alignment)):
                    iloc += 1
                    continue
            
            #print_log(preitem, nowitem)#
            target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start)
            if(len(target) > 0 and len(query) >0):
                st = time.time()
                
                cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 2, gap_open_2 = 24, gap_extend_2 = 1, bw=-1, zdropvalue=zdrop_value_drop, eqx = eqx)
                if(tmpdelcount>0):
                    nowDEL = True
                else:
                    nowDEL = False
                if(tmpinscount>0):
                    nowINS = True
                else:
                    nowINS = False
                if((tmpdelcount * tmpinscount) >= 1):
                    indel_present = True
                else:
                    indel_present = False
                if(zdropedcode == 0):
                    zdroped = False
                else:
                    zdroped = True
                droptime += time.time() - st
                if(H == False):
                    if(min(readgap, refgap) < 500):
                        if(indel_present == True):
                            if(len(cigarlist[-1]) > 1):
                                new_alignment[-1].pop(-1)
                                cigarlist[-1].pop(-1)
                            preitem = new_alignment[-1][-1]
                            preDEL, preINS = False, False
                            iloc += 1
                            continue
                        if(preDEL == True or preINS == True):
                            if(nowDEL == True or nowINS == True):
                                if(len(cigarlist[-1]) > 1):
                                    new_alignment[-1].pop(-1)
                                    cigarlist[-1].pop(-1)
                                preitem = new_alignment[-1][-1]
                                preDEL, preINS = False, False
                                iloc += 1
                                continue
                    else:
                        indel_present = False
                        preDEL, preINS = False, False
                        
                preDEL, preINS = nowDEL, nowINS
                if(zdroped == True or indel_present == True):
                    if(debug): print_log('split anchor: preitem, nowitem, zdroped, indel_present', preitem, nowitem, zdroped, indel_present)
                    if(debug): print_log('min(readgap, refgap), merge_smallgap', min(readgap, refgap), merge_smallgap)
                    if(min(readgap, refgap) < 0):
                        print_log('error min(readgap, refgap) < 0', min(readgap, refgap), (readgap, refgap))
                    split = True
                    if(H == False):
                        if(indel_present == False):
                            if(min(readgap, refgap) < merge_smallgap):
                                preDEL, preINS = False, False
                                split = False

                    if(split == True):
                        preDEL, preINS = False, False
                        if(debug): print_log('split anchor: preitem, nowitem, zdroped, indel_present', preitem, nowitem, zdroped, indel_present)
                        st = time.time()
                        batch = check_drop_site(preitem, nowitem, testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start)
                        splittime += time.time() - st
                        if((((((batch[1][0] - batch[0][0]) > 10)) and (((batch[1][0] - batch[0][0]) > merge_smallgap) or indel_present == True)) and H == False) or ((H == True) and (((batch[1][0] - batch[0][0]) > 50)))):    
                            if(debug): print_log('split anchor: preitem, batch, nowitem', preitem, batch, nowitem)
                            target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, batch[0], testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start)

                            if(len(target) == 0 or len(query) == 0):
                                if(len(target) == 0 and len(query) == 0):
                                    new_alignment[-1][-1] = batch[0]
                                else:
                                    print_log('unexcept: split_alignment_test: preitem, nowitem ', preitem, nowitem)

                            else:
                                new_alignment[-1].append(batch[0])
                                cigarstring, _, _, _, _, _ = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 2, gap_open_2 = 24, gap_extend_2 = 1, bw=-1, zdropvalue=-1, eqx = eqx)
                                #cigarstring = wfa(target, query)[0]
                                cigarlist[-1].append(cigarstring)
                                #print_log(cigarstring)#


                            if(len(batch) > 2):
                                new_alignment.append(List([batch[3], batch[2]]))  
                                target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(batch[2], batch[3], testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start)
                                cigarstring, _, _, _, _, _ = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 2, gap_open_2 = 24, gap_extend_2 = 1, bw=-1, zdropvalue=-1, eqx = eqx)
                                cigarlist.append([cigarstring])
                                #print_log(cigarstring)#


                            target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(batch[1], nowitem, testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start)

                            if(len(target) == 0 or len(query) == 0):
                                if(len(target) == 0 and len(query) == 0):
                                    new_alignment.append(List([batch[1]]))
                                    cigarlist.append([])
                                    nowitem = batch[1]
                                else:
                                    print_log('unexcept: split_alignment_test: preitem, nowitem ', preitem, nowitem)
                            else:      
                                new_alignment.append(List([batch[1], nowitem]))
                                cigarstring, _, _, _, _, _ = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 2, gap_open_2 = 24, gap_extend_2 = 1, bw=-1, zdropvalue=-1, eqx = eqx)
                                cigarlist.append([cigarstring])
                                #print_log(cigarstring)#
                        else:
                            if(debug == True):
                                print_log('split anchor: extend overlaped')
                                print_log(preitem, nowitem)
                                print_log(batch)
                                print_log()
                            cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 2, gap_open_2 = 24, gap_extend_2 = 1, bw=-1, zdropvalue=-1, eqx = eqx)
                            new_alignment[-1].append(nowitem)
                            cigarlist[-1].append(cigarstring)
                            #print_log(cigarstring)#
                    else:
                        cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 2, gap_open_2 = 24, gap_extend_2 = 1, bw=-1, zdropvalue=-1, eqx = eqx)
                        new_alignment[-1].append(nowitem)
                        cigarlist[-1].append(cigarstring) 
                        #print_log(cigarstring)#
                else:
                    new_alignment[-1].append(nowitem)
                    cigarlist[-1].append(cigarstring)
                    #print_log(cigarstring)#
            else:
                if(len(query)==0 and len(target) == 0):
                    iloc += 1
                    continue
                if(len(target) == 0):
                    cigarstring = str(len(query))+'I'
                else:
                    cigarstring = str(len(target))+'D'
                new_alignment[-1].append(nowitem)
                cigarlist[-1].append(cigarstring)
                #print_log('len(target), len(query)', len(target), len(query), cigarstring)#

            preitem = nowitem
            iloc += 1
        if(cigarlist[-1] == []):
            new_alignment.pop(-1)
            cigarlist.pop(-1)

        return new_alignment, cigarlist
    else:
        
        if(alignment[0][3] != 0):
            alignment[0] = (alignment[0][0], alignment[0][1] + alignment[0][3], -1, 0)
        if(alignment[-1][3] != 0):
            alignment[-1] = (alignment[-1][0] + alignment[-1][3], alignment[-1][1], -1, 0)
        
        alignment = alignment[::-1]
        nowitem = alignment[0]
        new_alignment.append(List([nowitem]))
        preitem = nowitem
        cigarlist.append([])
        iloc = 1
        while(iloc < len(alignment)):
            nowitem = alignment[iloc]
            readgap = preitem[0] - nowitem[0] - nowitem[3]
            refgap = nowitem[1]  - preitem[1] - preitem[3]
            
            if((nowitem[3] < 19) or (min(readgap, refgap) < min_gap_forcigar)):
                if(iloc +1 != len(alignment)):
                    iloc += 1
                    continue
            #print_log(preitem, nowitem)
            target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(nowitem, preitem, testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start)
            if(len(target) > 0 and len(query) >0):
                cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 2, gap_open_2 = 24, gap_extend_2 = 1, bw=-1, zdropvalue=zdrop_value_drop, eqx = eqx)
                if(tmpdelcount>0):
                    nowDEL = True
                else:
                    nowDEL = False
                if(tmpinscount>0):
                    nowINS = True
                else:
                    nowINS = False
                if((tmpdelcount* tmpinscount) >= 1):
                    indel_present = True
                else:
                    indel_present = False
                if(zdropedcode == 0):
                    zdroped = False
                else:
                    zdroped = True
                if(H == False):
                    if(min(readgap, refgap) < 500):
                        if(indel_present == True):
                            if(len(cigarlist[-1]) > 1):
                                new_alignment[-1].pop(-1)
                                cigarlist[-1].pop(-1)
                            preitem = new_alignment[-1][-1]
                            preDEL, preINS = False, False
                            iloc += 1
                            continue
                        if(preDEL == True or preINS == True):
                            if(nowDEL == True or nowINS == True):
                                if(len(cigarlist[-1]) > 1):
                                    new_alignment[-1].pop(-1)
                                    cigarlist[-1].pop(-1)
                                preitem = new_alignment[-1][-1]
                                preDEL, preINS = False, False
                                iloc += 1
                                continue
                    else:
                        indel_present = False
                        preDEL, preINS = False, False
                preDEL, preINS = nowDEL, nowINS
                if(zdroped == True or indel_present == True):
                    if(debug): print_log('split anchor: preitem, nowitem, zdroped, indel_present', preitem, nowitem, zdroped, indel_present)
                    if(debug): print_log('min(readgap, refgap), merge_smallgap', min(readgap, refgap), merge_smallgap)
                    if(min(readgap, refgap) < 0):
                        print_log('error min(readgap, refgap) < 0', min(readgap, refgap), (readgap, refgap))
                    split = True
                    if(H == False):
                        if(indel_present == False):
                            if(min(readgap, refgap) < merge_smallgap):
                                preDEL, preINS = False, False
                                split = False

                    if(split == True):
                        preDEL, preINS = False, False
                        if(debug): print_log('split anchor: preitem, nowitem, zdroped, indel_present', preitem, nowitem, zdroped, indel_present)
                        batch = check_drop_site(preitem, nowitem, testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start)
                        
                        if((((((batch[1][0] - batch[0][0]) > 10)) and (((batch[1][0] - batch[0][0]) > merge_smallgap) or indel_present == True)) and H == False) or ((H == True) and (((batch[1][0] - batch[0][0]) > 50)))):
                        
                            if(debug): print_log('split anchor: preitem, batch, nowitem', preitem, batch, nowitem)
                            target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(batch[0], preitem, testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start)

                            if(len(target) == 0 or len(query) == 0):
                                if(len(target) == 0 and len(query) == 0):
                                    new_alignment[-1][-1] = batch[0]
                                else:
                                    print_log('unexcept: split_alignment_test: preitem, nowitem ', preitem, nowitem)

                            else:
                                new_alignment[-1].append(batch[0])
                                cigarstring, _, _, _, _, _ = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 2, gap_open_2 = 24, gap_extend_2 = 1, bw=-1, zdropvalue=-1, eqx = eqx)
                                #cigarstring = wfa(target, query)[0]
                                cigarlist[-1].append(cigarstring)
                                #print_log(cigarstring)#
                                

                            if(len(batch) > 2):
                                new_alignment.append(List([batch[2], batch[3]]))  
                                target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(batch[2], batch[3], testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start)
                                cigarstring, _, _, _, _, _ = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 2, gap_open_2 = 24, gap_extend_2 = 1, bw=-1, zdropvalue=-1, eqx = eqx)
                                #cigarstring = wfa(target, query)[0]
                                cigarlist.append([cigarstring])
                                #print_log(cigarstring)#

                            target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(nowitem, batch[1], testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start)
                            if(len(target) == 0 or len(query) == 0):
                                if(len(target) == 0 and len(query) == 0):
                                    new_alignment.append(List([batch[1]]))
                                    cigarlist.append([])
                                    nowitem = batch[1]
                                else:
                                    print_log('unexcept: split_alignment_test: preitem, nowitem ', preitem, nowitem)
                            else:      
                                new_alignment.append(List([batch[1], nowitem]))
                                cigarstring, _, _, _, _, _ = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 2, gap_open_2 = 24, gap_extend_2 = 1, bw=-1, zdropvalue=-1, eqx = eqx)
                                #cigarstring = wfa(target, query)[0]
                                cigarlist.append([cigarstring])
                                #print_log(cigarstring)#
                        else:
                            if(debug == True):
                                print_log('split anchor: extend overlaped')
                                print_log(preitem, nowitem)
                                print_log(batch)
                                print_log()
                            cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 2, gap_open_2 = 24, gap_extend_2 = 1, bw=-1, zdropvalue=-1, eqx = eqx)
                            new_alignment[-1].append(nowitem)
                            cigarlist[-1].append(cigarstring)
                            #print_log(cigarstring)#
                    else:
                        cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 2, gap_open_2 = 24, gap_extend_2 = 1, bw=-1, zdropvalue=-1, eqx = eqx)
                        new_alignment[-1].append(nowitem)
                        cigarlist[-1].append(cigarstring)
                        #print_log(cigarstring)#
                        
                else:
                    new_alignment[-1].append(nowitem)
                    cigarlist[-1].append(cigarstring)
                    #print_log(cigarstring)#
            else:

                if(len(query)==0 and len(target) == 0):
                    iloc += 1
                    continue
                if(len(target) == 0):
                    cigarstring = str(len(query))+'I'
                else:
                    cigarstring = str(len(target))+'D'
                new_alignment[-1].append(nowitem)
                cigarlist[-1].append(cigarstring)
                #print_log('len(target), len(query)', len(target), len(query), cigarstring)#


            preitem = nowitem
            iloc += 1

        return new_alignment, cigarlist


def get_readmap_DP_test(readid, testseq, contig2start, contig2seq, index_object, index2contig, kmersize = 15, local_kmersize = 9, local_maxdiff = 30, refine = True, local_skipcost = 50., golbal_skipcost = (30., 30.),  golbal_maxdiff = (50, 30), check_num = 20, bin_size = 100, hastra = False, debug = False, H = False):

    local_skipcost = 30.    
    mapq = 60
    setting_kmersize = kmersize
    setting_maxdiff = golbal_maxdiff[1]
    if(H == False):
        local_skipcost += local_kmersize
        golbal_skipcost = (golbal_skipcost[0] + kmersize, golbal_skipcost[1] + kmersize)
    else:
        local_skipcost = 30.
        golbal_skipcost = (30., 30.)
        hastra = True
        check_num = 100
    

    rc_testseq = str(Seq(testseq).reverse_complement())
    testseq_len = len(testseq)
    st = time.time()                            
    mapq, scores, raw_alignment_list = decode_hit(index_object, index2contig, testseq, testseq_len, contig2start, kmersize, skipcost = golbal_skipcost, maxdiff = golbal_maxdiff, maxgap = 2000, check_num = check_num, c_bias = 5000, bin_size = bin_size, overlapprecentage = 0.5, hastra = hastra, H = H)
    #raw_alignment_list = check_func_clean(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, 10, debug = debug)
    if(debug == True):
        plot_result = np.array(raw_alignment_list)
        plt.scatter(plot_result[:, 0], plot_result[:, 1])
        plt.show()
    if(scores == 0.):
        return
    
    need_reverse = False
    if((scores < 0.)):
        need_reverse = True

    if(debug == True): print_log('Normal time', time.time() - st, scores)
    if(refine == True):
        setting_maxdiff = local_maxdiff
        setting_kmersize = local_kmersize
        st = time.time()

        #need_reverse, raw_alignment_array = get_reversed_chain_numpy_rough(np.array(raw_alignment_list), testseq_len)
        raw_alignment_array = np.array(raw_alignment_list)
        if(need_reverse == False):
            scores, raw_alignment_list = get_localmap_multi_all_forDP_inv_guide(raw_alignment_array, testseq, rc_testseq, contig2start, contig2seq, kmersize = setting_kmersize, skipcost = local_skipcost, maxdiff = setting_maxdiff, maxgap = 100, shift = 1)

        else:
            scores, raw_alignment_list = get_localmap_multi_all_forDP_inv_guide(raw_alignment_array, rc_testseq, testseq, contig2start, contig2seq, kmersize = setting_kmersize, skipcost = local_skipcost, maxdiff = setting_maxdiff, maxgap = 100, shift = 1)
            testseq, rc_testseq = rc_testseq, testseq
            #raw_alignment_list = get_reversed_chain_numpy(np.array(raw_alignment_list), testseq_len)

        #raw_alignment_list = check_func_clean(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, 10, debug = debug)
        #for item in raw_alignment_list:
            #print_log(item)

        if(debug == True): print_log('Refine time', time.time() - st)


        
        

    st = time.time()
    alignment_list, onemapinfolist, TRA_signal = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, debug = debug, H = H)
    if(debug == True): print_log('extend_func', time.time() - st)
    if(len(onemapinfolist) > 1 and refine == 'auto'):
        setting_maxdiff = local_maxdiff
        setting_kmersize = local_kmersize
        st = time.time()                                                                         
        scores, raw_alignment_list = get_localmap_multi_all_forDP_inv_guide(np.array(raw_alignment_list), testseq, rc_testseq, contig2start, contig2seq, kmersize = setting_kmersize, skipcost = local_skipcost, maxdiff = setting_maxdiff, maxgap = 100, shift = 1)
        if(debug == True): print_log('Refine time', time.time() - st, len(raw_alignment_list))

        alignment_list, onemapinfolist, TRA_signal = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, need_reverse, setting_maxdiff, debug = debug)

    return onemapinfolist, (alignment_list, raw_alignment_list), TRA_signal


@njit
def get_reversed_chain_numpy(raw_alignment_array, testseq_len):
    raw_alignment_array[:,0] = testseq_len - raw_alignment_array[:,0] - raw_alignment_array[:,3]
    raw_alignment_array[:,2] *= -1 
    return [(item[0], item[1], item[2], item[3]) for item in raw_alignment_array[::-1]]
def get_reversed_chain_numpy_rough(raw_alignment_array, testseq_len):
    np_counts_dict = {-1: 0, 1: 0}
        
    np_counts = np.unique(raw_alignment_array[:,2], return_counts = True)
    for tmp_iloc_npcount in range(len(np_counts[0])):
        np_counts_dict[np_counts[0][tmp_iloc_npcount]] += np_counts[1][tmp_iloc_npcount]

    if(np_counts_dict[-1] > np_counts_dict[1]):

        raw_alignment_array[:,0] = testseq_len - raw_alignment_array[:,0] - raw_alignment_array[:,3]
        raw_alignment_array[:,2] *= -1 
        return True, raw_alignment_array[::-1]
    else:
        return False, raw_alignment_array
def get_onemapinfolist(new_alignment_list, cigarlist, readid, mapq, testseq_len, contig2start, need_reverse):
    onemapinfolist = []
    iloc = -1
    if(need_reverse == False):
        for alignment in new_alignment_list:
            contig = pos2contig(alignment[0][1], contig2start)
            iloc += 1
            refbias = contig2start[contig]
            if(alignment[0][2] == 1):
                query_st = alignment[0][0]
                query_en = alignment[-1][0] + alignment[-1][3]
                target_st = alignment[0][1]
                target_en = alignment[-1][1] + alignment[-1][3]
                if(query_st > 0):
                    topcigar = str(query_st) + 'S'
                else:

                    topcigar = ''
                if((testseq_len - query_en) > 0):
                    tailcigar = str(testseq_len - query_en) + 'S'
                else:
                    tailcigar = ''
                if(alignment[-1][3] > 0):
                    tailcigar = str(int(alignment[-1][3]))+'M'+tailcigar
                cigarstring = ''.join(cigarlist[iloc])
                onemapinfolist.append((readid, contig, '+', query_st, query_en, target_st-refbias, target_en-refbias, mapq, topcigar+cigarstring+tailcigar))
            else:
                query_st = testseq_len-alignment[0][0]-alignment[0][3]
                query_en = testseq_len-alignment[-1][0]
                target_st = alignment[0][1]
                target_en = alignment[-1][1] + alignment[-1][3]

                cigarstring = ''.join(cigarlist[iloc])
                if(query_st > 0):
                    topcigar = str(query_st) + 'S'
                else:
                    topcigar = ''
                if((testseq_len - query_en) > 0):
                    tailcigar = str(testseq_len - query_en) + 'S'
                else:
                    tailcigar = ''
                onemapinfolist.append((readid, contig, '-', query_st, query_en, target_st-refbias, target_en-refbias, mapq, topcigar+cigarstring+tailcigar))
        return new_alignment_list, onemapinfolist
    else:
        for alignment in new_alignment_list:
            contig = pos2contig(alignment[0][1], contig2start)
            iloc += 1
            refbias = contig2start[contig]
            if(alignment[0][2] == 1):
                query_st = alignment[0][0]
                query_en = alignment[-1][0] + alignment[-1][3]
                target_st = alignment[0][1]
                target_en = alignment[-1][1] + alignment[-1][3]
                if(query_st > 0):
                    topcigar = str(query_st) + 'S'
                else:

                    topcigar = ''
                if((testseq_len - query_en) > 0):
                    tailcigar = str(testseq_len - query_en) + 'S'
                else:
                    tailcigar = ''
                if(alignment[-1][3] > 0):
                    tailcigar = str(int(alignment[-1][3]))+'M'+tailcigar
                cigarstring = ''.join(cigarlist[iloc])
                onemapinfolist.append((readid, contig, '-', query_st, query_en, target_st-refbias, target_en-refbias, mapq, topcigar+cigarstring+tailcigar))
            else:
                query_st = testseq_len-alignment[0][0]-alignment[0][3]
                query_en = testseq_len-alignment[-1][0]
                target_st = alignment[0][1]
                target_en = alignment[-1][1] + alignment[-1][3]

                cigarstring = ''.join(cigarlist[iloc])
                if(query_st > 0):
                    topcigar = str(query_st) + 'S'
                else:
                    topcigar = ''
                if((testseq_len - query_en) > 0):
                    tailcigar = str(testseq_len - query_en) + 'S'
                else:
                    tailcigar = ''
                onemapinfolist.append((readid, contig, '+', query_st, query_en, target_st-refbias, target_en-refbias, mapq, topcigar+cigarstring+tailcigar))

        
        return new_alignment_list, onemapinfolist[::-1]

##############
#1117mapq
@njit
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list1(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = [0]
    for gapcost in range(1, maxdiff + 1):
        gapcost_list.append(0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))


    n = len(one_mapinfo)
    S = np.zeros(n, np.float64)
    P = np.zeros(n, np.int64)
    start = 0
    
    testdict = dict()
    testdict[-1] = 0.
    testdict.pop(-1)
    maxop = 2000*(one_mapinfo[-1][0]-one_mapinfo[0][0])
    #print_log(maxop)
    opcount = 0
    
    coverage_dict = dict()
    for i in range(n):
        if(one_mapinfo[i][0] in coverage_dict):
            coverage_dict[one_mapinfo[i][0]] += 1
        else:
            coverage_dict[one_mapinfo[i][0]] = 1

    #coverage_array *= 2
    
    oskipcost = skipcost
    omaxdiff = maxdiff
    
    for i in range(n):
        P[i] = 9999999
        max_scores = one_mapinfo[i][3]
        pre_index = 9999999
        if(opcount > 2*maxop):
            path = []
            take_index = g_max_index
            path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
            while(True):
                if((P[take_index] == 9999999)):
                    break
                take_index = P[take_index]
                path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

            return g_max_scores, path
        if(opcount > maxop):
            if(max_scores < 15):
                S[i] = max_scores
                continue
        #while(((i - start) > 100) and ((one_mapinfo[i][0] - one_mapinfo[start][0]) > max_gap)):
            #start += 1
            
        skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
        maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)    
            
        j = i
        skipcount = 0
        #while(j > start):
            #j -= 1
        if(len(testdict) > 50):
            
            cutoff = g_max_scores - max(2.5*skipcost, 60.)

            tmp_poplist = List([0])
            for loc_in_one_mapinfo in testdict:
                if(S[loc_in_one_mapinfo] < cutoff and (one_mapinfo[i][0] - one_mapinfo[loc_in_one_mapinfo][0]>100)):
                    tmp_poplist.append(loc_in_one_mapinfo)

            if(len(tmp_poplist) != 1):
                tmp_poplist.pop(0)
                for loc_in_one_mapinfo in tmp_poplist:
                    testdict.pop(loc_in_one_mapinfo)




               
        
        for j in testdict:
            
                    
            #(12240, 2791460820, 1, 21) (12251, 2791460711, 1, 108)
            nocost = False
            filtered = True
            if(one_mapinfo[i][0] == one_mapinfo[j][0]):
                continue

                      #      x_2 - x_1 - a_1
            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            if((readgap < 0)):
                continue
                #if((one_mapinfo[i][2] == one_mapinfo[j][2]) and ((readgap + one_mapinfo[i][3]) > 10)):
                    #pass
                #else:
                    #continue
            opcount += 1
            if(one_mapinfo[i][2] == 1):
                          #    y_2 - y_1 - a_1
                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
            else:
                         #     y_1 - y_2 - a_1
                refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
            while(True):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    break

              
                if((one_mapinfo[i][2] == 1) and (refgap < 0)):
                    break
                if(one_mapinfo[i][2] == -1):
                    if(refgap > 0):
                        break
                    else:
                        refgap = abs(refgap)

                gapcost = abs(abs(readgap) - refgap)
                if(((readgap ) > maxgap) or (gapcost > maxdiff)):
                    break
                
                gapcost = gapcost_list[gapcost]
                #if(gapcost != 0):
                    #gapcost = 0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost)
                    #gapcost = 0.01 * kmersize * gapcost + 2.0 * np.log2(gapcost)
        
        
                test_scores = S[j] + (min(readgap, 0) + one_mapinfo[i][3]) - gapcost
                filtered = False
                break
            if(filtered == True):
                """if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(magic > abs(one_mapinfo[i][1] - one_mapinfo[j][1])):
                        continue"""
                test_scores = S[j] - skipcost + one_mapinfo[i][3] - min(10., abs(refgap)/100)# - np.log2(readgap)# - np.log2(abs(abs(refgap) - readgap)+1)
                #continue
            if(test_scores >= max_scores):
                max_scores = test_scores
                pre_index = j

       
            
                
        S[i] = max_scores
        P[i] = pre_index 
        
        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
            
 
    
            

        testdict[i] = max_scores
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    while(True):
        if((P[take_index] == 9999999)):
            break
        take_index = P[take_index]
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

    #print_log(opcount)
    return g_max_scores, path#, testdict
########################
#1212
#print_log(1212, 'coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, 20)')
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    at = skipcost + 20.
    lastupdatereadloc = 0
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = [0]
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list.append(0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list.append(0.01 * kmersize * gapcost + 2 * np.log2(gapcost))


    n = len(one_mapinfo)
    S = np.zeros(n, np.float64)
    P = np.zeros(n, np.int64)

    
    testdict = dict()
    testdict[-1] = 0.
    testdict.pop(-1)
    
    maxop = 20000*(one_mapinfo[-1][0]-one_mapinfo[0][0])
    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    
    
    coverage_dict = dict()
    for i in range(n):
        if(one_mapinfo[i][0] in coverage_dict):
            coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, 20)
        else:
            coverage_dict[one_mapinfo[i][0]] = 1
    
    oskipcost = skipcost
    omaxdiff = maxdiff
    
    for i in range(n):
        
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        pre_index = -9999999
        
        if(prereadloc < one_mapinfo[i][0]):
            
            if(skipcost + 20. > at):
                at = skipcost + 20.
                lastupdatereadloc = i
            elif(one_mapinfo[i][0] - lastupdatereadloc > 500):
                at = skipcost + 20.
                lastupdatereadloc = one_mapinfo[i][0]
            
            prereadloc = one_mapinfo[i][0]
            cutoff = g_max_scores - at
            tmp_poplist = List([0])
            for loc_in_one_mapinfo in testdict:
                if(S[loc_in_one_mapinfo] < cutoff):
                    tmp_poplist.append(loc_in_one_mapinfo)

            if(len(tmp_poplist) != 1):
                tmp_poplist.pop(0)
                for loc_in_one_mapinfo in tmp_poplist:
                    testdict.pop(loc_in_one_mapinfo)     
        if(opcount > maxop):
            return 0., [(0, 0, 0, 0)]
        
        skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
        maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
        
        
        
        for j in testdict:
            
                    
            #(12240, 2791460820, 1, 21) (12251, 2791460711, 1, 108)
            if(S[j] < (max_scores - one_mapinfo[i][3])):
                continue
            nocost = False
            filtered = True
            if(one_mapinfo[i][0] == one_mapinfo[j][0]):
                continue

                      #      x_2 - x_1 - a_1
            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            if((readgap < 0)):
                continue
                
            

            if(one_mapinfo[i][2] == 1):
                          #    y_2 - y_1 - a_1
                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
            else:
                         #     y_1 - y_2 - a_1
                refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                
            opcount += 1
            while(True):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    break
                if((one_mapinfo[i][2] == 1) and (refgap < 0)):
                    break
                if(one_mapinfo[i][2] == -1):
                    if(refgap > 0):
                        break
                    else:
                        refgap = abs(refgap)

                gapcost = abs(readgap - refgap)
                if((readgap  > maxgap) or (gapcost > maxdiff)):
                    break
                
                gapcost = gapcost_list[gapcost]
                
                test_scores = S[j] + one_mapinfo[i][3] - gapcost
                filtered = False
                break
            if(filtered == True):
                test_scores = S[j] - skipcost + one_mapinfo[i][3]
                #gapcost = abs(readgap - refgap)
                #if(gapcost != 0):
                    #test_scores = S[j] - skipcost + one_mapinfo[i][3] - min(np.log(abs(gapcost)), abs(gapcost)/100, 12)
                #else:
                    #test_scores = S[j] - skipcost + one_mapinfo[i][3]
            if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
            elif((test_scores == max_scores) and (pre_index != -9999999)):
                pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                if(now_refdistance < pre_refdistance):
                    max_scores = test_scores
                    pre_index = j

        S[i] = max_scores
        P[i] = pre_index
        
        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
            
        testdict[i] = max_scores
        
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    while(True):
        if((P[take_index] == -9999999)):
            break
        take_index = P[take_index]
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

    return g_max_scores, path

@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    at = skipcost + 20.
    lastupdatereadloc = 0
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = [0]
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list.append(0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list.append(0.01 * kmersize * gapcost + 2 * np.log2(gapcost))


    n = len(one_mapinfo)
    S = np.zeros(n, np.float64)
    P = np.zeros(n, np.int64)

    
    testdict = dict()
    testdict[-1] = 0.
    testdict.pop(-1)
    
    maxop = 20000*(one_mapinfo[-1][0]-one_mapinfo[0][0])
    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    
    
    coverage_dict = dict()
    for i in range(n):
        if(one_mapinfo[i][0] in coverage_dict):
            coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, 20)
        else:
            coverage_dict[one_mapinfo[i][0]] = 1
    
    oskipcost = skipcost
    omaxdiff = maxdiff
    
    for i in range(n):
        
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        pre_index = -9999999
        
        if(prereadloc < one_mapinfo[i][0]):
            
            if(skipcost + 20. > at):
                at = skipcost + 20.
                lastupdatereadloc = i
            elif(one_mapinfo[i][0] - lastupdatereadloc > 500):
                at = skipcost + 20.
                lastupdatereadloc = one_mapinfo[i][0]
            
            prereadloc = one_mapinfo[i][0]
            cutoff = g_max_scores - at
            tmp_poplist = List([0])
            for loc_in_one_mapinfo in testdict:
                if(S[loc_in_one_mapinfo] < cutoff):
                    tmp_poplist.append(loc_in_one_mapinfo)

            if(len(tmp_poplist) != 1):
                tmp_poplist.pop(0)
                for loc_in_one_mapinfo in tmp_poplist:
                    testdict.pop(loc_in_one_mapinfo)     
        if(opcount > maxop):
            return 0., [(0, 0, 0, 0)]
        
        skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
        maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
        
        
        
        for j in testdict:
            
                    
            #(12240, 2791460820, 1, 21) (12251, 2791460711, 1, 108)
            if(S[j] < (max_scores - one_mapinfo[i][3])):
                continue
            nocost = False
            filtered = True
            if(one_mapinfo[i][0] == one_mapinfo[j][0]):
                continue

                      #      x_2 - x_1 - a_1
            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            if((readgap < 0)):
                continue
                
            

            if(one_mapinfo[i][2] == 1):
                          #    y_2 - y_1 - a_1
                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
            else:
                         #     y_1 - y_2 - a_1
                refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                
            opcount += 1
            while(True):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    break
                if((one_mapinfo[i][2] == 1) and (refgap < 0)):
                    break
                if(one_mapinfo[i][2] == -1):
                    if(refgap > 0):
                        break
                    else:
                        refgap = abs(refgap)

                gapcost = abs(readgap - refgap)
                if((readgap  > maxgap) or (gapcost > maxdiff)):
                    break
                
                gapcost = gapcost_list[gapcost]
                
                test_scores = S[j] + one_mapinfo[i][3] - gapcost
                filtered = False
                break
            if(filtered == True):
                #test_scores = S[j] - skipcost + one_mapinfo[i][3]
                gapcost = abs(readgap - refgap)
                if(gapcost != 0):
                    test_scores = S[j] - skipcost + one_mapinfo[i][3] - min(np.log(abs(gapcost)), abs(gapcost)/100, 12)
                else:
                    test_scores = S[j] - skipcost + one_mapinfo[i][3]
            if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
            elif((test_scores == max_scores) and (pre_index != -9999999)):
                pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                if(now_refdistance < pre_refdistance):
                    max_scores = test_scores
                    pre_index = j

        S[i] = max_scores
        P[i] = pre_index
        
        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
            
        testdict[i] = max_scores
        
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    while(True):
        if((P[take_index] == -9999999)):
            break
        take_index = P[take_index]
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

    return g_max_scores, path


@njit
def hit2work(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage = 0.5, hastra = False, H = False):
    def getsecond(x):
        return x[1]
    def getfirst(x):
        return x[0]
    def getlength(x):
        return len(x)

    def get_readloc_set_bin(one_mappos, bin_size):
        return set([i[0] // bin_size for i in one_mappos])
    def get_refloc_set_bin(one_mappos, bin_size):
        return set([i[1] // bin_size for i in one_mappos])
    def get_overlapsize(readloc_set_a, readloc_set_b):
        return len(readloc_set_a&readloc_set_b)/min(len(readloc_set_a), len(readloc_set_b))
    
    maxgap = 200

    
    one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:,1])] 

    cluster_list = List()
    iloc = 0
    preitem = one_mapinfo[iloc]
    st_iloc = 0
    for iloc in range(len(one_mapinfo)):
        nowitem = one_mapinfo[iloc]
        if(((nowitem[1] - preitem[1]) > c_bias)):
            if((iloc - st_iloc) < 3):
                continue
            cluster_list.append(one_mapinfo[st_iloc: iloc])
            st_iloc = iloc
        preitem = nowitem


    if((iloc - st_iloc) > 3):
        cluster_list.append(one_mapinfo[st_iloc: iloc + 1])
    cluster_list.sort(key = get_length)

    cluster_list = cluster_list[::-1][:check_num]


    
    hit = False
    minichain_scores = 40
    path_list = List()
    scores_list = []
    scores_list.append(0.)
    scores_list.pop()
    path_list.append([(0, 0, 0, 0)])
    path_list.pop()
    

    max_scores = 0
    from_repeat = False

    for one_mapinfo in cluster_list:
        repeat = False
        min_ref, max_ref = one_mapinfo[0][1], one_mapinfo[-1][1]
        testcontig = pos2contig(min_ref, contig2start)
        if(testcontig != pos2contig(max_ref, contig2start)):
            #print_log('hit2work: testcontig != pos2contig(max_ref, contig2start)')
            continue

        scores, path = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d(one_mapinfo[np.argsort(one_mapinfo[:,0])], kmersize = kmersize, skipcost = skipcost[0], maxdiff = maxdiff[0], maxgap = maxgap)
 
        
        if(scores < minichain_scores):
            if(hit == False):
                return 0, 0., [(0, 0, 0, 0)]
            continue
        
        

        hit = True

        path_list.append(path) 
        scores_list.append(scores)
        if(scores > max_scores):
            if(repeat == True):
                from_repeat = True
            else:
                from_repeat = False
            max_scores = scores
    
    if(hit == True and max_scores > 200):

        order = np.argsort(np.array(scores_list))[::-1]
        
        primary_rlocset_List = List()
        primary_scores_List = List()
        primary_index_List = List()
        
        iloc = order[0]
        primary_rlocset_List.append(get_readloc_set_bin(path_list[iloc], bin_size))
        primary_scores_List.append(List([scores_list[iloc]]))
        primary_index_List.append(iloc)

        for iloc in order[1:]:
            readloc_set_b = get_readloc_set_bin(path_list[iloc], bin_size)
            maxoverlapsize = 0.
            for p_loc in range(len(primary_rlocset_List)):
                tmp_overlapsize = get_overlapsize(primary_rlocset_List[p_loc], readloc_set_b)
                if(tmp_overlapsize > maxoverlapsize):
                    maxoverlapsize = tmp_overlapsize
                    prefer_p_loc = p_loc
            if(maxoverlapsize < overlapprecentage):
                primary_rlocset_List.append(readloc_set_b)
                primary_scores_List.append(List([scores_list[iloc]]))
                primary_index_List.append(iloc)
            else:
                primary_scores_List[prefer_p_loc].append(scores_list[iloc])
                
        m = len(path_list[order[0]])    
        if(len(primary_scores_List[0]) < 2):
            f1 = primary_scores_List[0][0]
            f2 = 0
        else:
            f1 = primary_scores_List[0][0]
            f2 = primary_scores_List[0][1]
        mapq = min(int(40*(1-f2/f1)*min(1, m/10)*np.log(f1)), 60)

        base_iloc = primary_index_List[0]

        for iloc in primary_index_List[1:4]:
            if(hastra == False):
                skip = True
                item = path_list[iloc][0]
                bias_1 = abs(item[1]-path_list[base_iloc][0][1])
                bias_2 = abs(item[1] - path_list[base_iloc][-1][1])
                if(min(bias_1, bias_2) < 200000):
                    if(pos2contig(item[1], contig2start) == pos2contig(path_list[base_iloc][0][1], contig2start)):
                        skip = False
                if(skip == True):
                    continue
            for item in path_list[iloc]:                            
                path_list[base_iloc].append(item)

        path_list[base_iloc].sort(key = getfirst)
        scores, path = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d(np.array(path_list[base_iloc]), kmersize = kmersize, skipcost = skipcost[1], maxdiff = maxdiff[1], maxgap = maxgap)

        return mapq, scores, path


    else:
        return 0, 0., [(0, 0, 0, 0)]

#print_log(1216)
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    at = skipcost + 20.
    lastupdatereadloc = 0
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = [0]
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list.append(0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list.append(0.01 * kmersize * gapcost + 2 * np.log2(gapcost))


    n = len(one_mapinfo)
    S = np.zeros(n, np.float64)
    P = np.zeros(n, np.int64)

    
    testdict = dict()
    testdict[-1] = 0.
    testdict.pop(-1)
    
    maxop = 20000*(one_mapinfo[-1][0]-one_mapinfo[0][0])
    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    
    
    coverage_dict = dict()
    for i in range(n):
        if(one_mapinfo[i][0] in coverage_dict):
            coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, 20)
        else:
            coverage_dict[one_mapinfo[i][0]] = 1
    
    oskipcost = skipcost
    omaxdiff = maxdiff
    
    gap_dict = dict()
    gap_dict[0] = 100
    
    for i in range(n):
        
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        pre_index = -9999999
        
        if(prereadloc < one_mapinfo[i][0]):
            
            gap_dict[prereadloc] = 200+ 2*(one_mapinfo[i][0] - prereadloc)#mark
            
            if(skipcost + 20. > at):
                at = skipcost + 20.
                lastupdatereadloc = i
            elif(one_mapinfo[i][0] - lastupdatereadloc > 500):
                at = skipcost + 20.
                lastupdatereadloc = one_mapinfo[i][0]
            
            prereadloc = one_mapinfo[i][0]
            cutoff = g_max_scores - at
            tmp_poplist = List([0])
            for loc_in_one_mapinfo in testdict:
                if(S[loc_in_one_mapinfo] < cutoff):
                    tmp_poplist.append(loc_in_one_mapinfo)

            if(len(tmp_poplist) != 1):
                tmp_poplist.pop(0)
                for loc_in_one_mapinfo in tmp_poplist:
                    testdict.pop(loc_in_one_mapinfo)     
        if(opcount > maxop):
            return 0., [(0, 0, 0, 0)]
        
        skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
        maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
        
        
        
        for j in testdict:
            
                    
            #(12240, 2791460820, 1, 21) (12251, 2791460711, 1, 108)
            if(S[j] < (max_scores - one_mapinfo[i][3])):
                continue
            nocost = False
            filtered = True
            if(one_mapinfo[i][0] == one_mapinfo[j][0]):
                continue

                      #      x_2 - x_1 - a_1
            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            if((readgap < 0)):
                continue
                
            

            if(one_mapinfo[i][2] == 1):
                          #    y_2 - y_1 - a_1
                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
            else:
                         #     y_1 - y_2 - a_1
                refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                
            opcount += 1
            while(True):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    break
                if((one_mapinfo[i][2] == 1) and (refgap < 0)):
                    break
                if(one_mapinfo[i][2] == -1):
                    if(refgap > 0):
                        break
                    else:
                        refgap = abs(refgap)

                gapcost = abs(readgap - refgap)
                if((readgap  > gap_dict[one_mapinfo[j][0]]) or (gapcost > maxdiff)):
                    break
                
                gapcost = gapcost_list[gapcost]
                
                test_scores = S[j] + one_mapinfo[i][3] - gapcost
                filtered = False
                break
            if(filtered == True):
                #test_scores = S[j] - skipcost + one_mapinfo[i][3]
                gapcost = abs(readgap - refgap)
                if(gapcost != 0):
                    test_scores = S[j] - skipcost + one_mapinfo[i][3] - min(np.log(abs(gapcost)), abs(gapcost)/100, 12)
                else:
                    test_scores = S[j] - skipcost + one_mapinfo[i][3]
            if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
            elif((test_scores == max_scores) and (pre_index != -9999999)):
                pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                if(now_refdistance < pre_refdistance):
                    max_scores = test_scores
                    pre_index = j

        S[i] = max_scores
        P[i] = pre_index
        
        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
            
        testdict[i] = max_scores
        
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    while(True):
        if((P[take_index] == -9999999)):
            break
        take_index = P[take_index]
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

    return g_max_scores, path

#print_log(1218)
@njit
def pairedindel(cigarlist, indelsize = 30):


    zero = ord('0')
    INSsyb = ord('I') - zero
    SOFTsyb = ord('S') - zero
    HARDsyb = ord('H') - zero
    PADsyb = ord('P') - zero
    DELsyb = ord('D') - zero
    delsyb = ord('^') - zero

    indel = []

    for cigar in cigarlist:
        
        
        typed_cigar = [ord(item) - zero for item in cigar] 
        
        number = 0.
        for item in typed_cigar:
            if(item < 10):
                number = number * 10. + item
  
            else:
                if(item != INSsyb and item != SOFTsyb and item != HARDsyb and item != PADsyb):
                    if(item == DELsyb and number > indelsize):
                        indel.append(number)
                    number = 0.

                else:

                    if(item == INSsyb and number > indelsize):
                        indel.append(number)
                        
                    number = 0.
    indel.sort()
    preitem = 0
    clustersize = 1
    for nowitem in indel:
        if((min(preitem, nowitem)/max(preitem, nowitem))>0.7):
            clustersize += 1
            if(clustersize > 1):
                return True
        else:
            clustersize = 1
        preitem = nowitem
    return False
def extend_func1(raw_alignment_list, readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, debug = False, H = False, nofilter = False, eqx = False):
    
    TRA_signal = False
    maxdiffratio = 0.1
    #step 1
    #rebuild chain break
    #remove anchor cause large cost and small alignment
    st = time.time()
    if(H == False):
        alignment_list = rebuild_chain_break(contig2start, raw_alignment_list, large_cost = setting_maxdiff, small_alignment = 50, small_dup = -100)
    else:
        alignment_list = rebuild_chain_break_H(contig2start, raw_alignment_list, large_cost = setting_maxdiff, small_alignment = 30, small_dup = -30)
        nofilter = True
    
        
    tmpiloc = -1
    while((tmpiloc + 1) < len(alignment_list)):
        tmpiloc += 1
        preitem, nowitem = alignment_list[tmpiloc][0], alignment_list[tmpiloc][-1]
        target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, setting_kmersize, contig2seq, contig2start)
        diffratio = edlib.align(query = query, target = target, task = 'distance')['editDistance']/min(len(target), len(query))
        if((diffratio>maxdiffratio)):
            alignment_list.pop(tmpiloc)
            tmpiloc -= 1
            
    #^print_log('step 1: rebuild chain break ', time.time() - st)
    if(debug == True):
        print_log('step 1: rebuild chain break')
        for line in alignment_list:
            tempcontig = pos2contig(line[0][1], contig2start)
            temprefbias = contig2start[tempcontig]
            preitem, nowitem = line[0], line[-1]

            target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, setting_kmersize, contig2seq, contig2start)

            print_log(preitem, nowitem)
            diffratio = edlib.align(query = query, target = target, task = 'distance')['editDistance']/min(len(target), len(query))
            print_log(line[0][0], line[-1][0], line[0][1] - temprefbias, line[-1][1] - temprefbias, diffratio)
            plot_result = np.array(line)
            plt.scatter(plot_result[:,0], plot_result[:,1])
        plt.show()
        
    #step 2
    #extend edge to recover small dup and misplaced alignment
    # and also to merge alignment gaped with deletion or insertion
    st = time.time()
    extend_edge_test(testseq, testseq_len, alignment_list, setting_kmersize, pos2contig, contig2start, contig2seq, san = 1, debug = debug)
    if(debug == True):print_log('step 2: extend edge ', time.time() - st)
    if(debug == True):
        print_log('After extend edge')
        for line in alignment_list:
            tempcontig = pos2contig(line[0][1], contig2start)
            temprefbias = contig2start[tempcontig]
            print_log(line[0][0], line[-1][0], line[0][1] - temprefbias, line[-1][1] - temprefbias)
        
    #step 3
    #remove miss placed alignment which cause del/ins or ins/del in start and end
    st = time.time()
    
    o_alignment_list_len = len(alignment_list)
    filtered = False
    if((len(alignment_list) > 2) and (nofilter == False)):    
        iloc = 0
        while(iloc < (len(alignment_list) - 2)):
            removed = drop_misplaced_alignment_test(alignment_list, iloc, debug = debug)
            if(removed == True):
                continue
            else:
                iloc += 1
    
        extend_edge_drop_test(testseq, testseq_len, alignment_list, setting_kmersize, pos2contig, contig2start, contig2seq, san = 1, debug = debug)
        
    if(debug == True):print_log('step 3: remove miss placed alignment ', time.time() - st)
    if(debug == True):
        print_log('After remove miss placed alignment')
        for line in alignment_list:
            tempcontig = pos2contig(line[0][1], contig2start)
            temprefbias = contig2start[tempcontig]
            print_log(line[0][0], line[-1][0], line[0][1] - temprefbias, line[-1][1] - temprefbias)
            plot_result = np.array(line)
            plt.scatter(plot_result[:,0], plot_result[:,1])
        plt.show()
            
    tmpiloc = -1
    while(H == False and (tmpiloc + 1) < len(alignment_list)):
        tmpiloc += 1
        preitem, nowitem = alignment_list[tmpiloc][0], alignment_list[tmpiloc][-1]
        target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, setting_kmersize, contig2seq, contig2start)
        diffratio = edlib.align(query = query, target = target, task = 'distance')['editDistance']/min(len(target), len(query))
        if(diffratio>maxdiffratio):
            alignment_list.pop(tmpiloc)
            tmpiloc -= 1
            
    
    if(len(alignment_list) <  o_alignment_list_len):#fill the gap
        filtered = True
        st = time.time()
        extend_edge_test(testseq, testseq_len, alignment_list, setting_kmersize, pos2contig, contig2start, contig2seq, san = 1, debug = debug)
        if(debug == True):print_log('step 4: fill the gap by extend edge ', time.time() - st)
        if(debug == True):
            print_log('After extend edge')
            for line in alignment_list:
                tempcontig = pos2contig(line[0][1], contig2start)
                temprefbias = contig2start[tempcontig]
                print_log(line[0][0], line[-1][0], line[0][1] - temprefbias, line[-1][1] - temprefbias)

        
    #step 4
    #merge conjacent alignment with small readgap or refgap
    st = time.time()
    merge_smallgap = 2000
    too_large_gap = 5000
    if(len(alignment_list) >= 2):
        iloc = 0
        duplist = getdupiloc(alignment_list)
        while((iloc + 1) < len(alignment_list)):
            skiped = False
            if(iloc in duplist):
                iloc += 1
                continue
            while(True):
                preitem = alignment_list[iloc][-1]
                nowitem = alignment_list[iloc + 1][0]
                presize = alignment_list[iloc][-1][0] + alignment_list[iloc][-1][3] - alignment_list[iloc][0][0]
                nowsize = alignment_list[iloc + 1][-1][0] + alignment_list[iloc + 1][-1][3] - alignment_list[iloc + 1][0][0]
                if(alignment_list[iloc][-1][2] == 1):
                    ref_presize = alignment_list[iloc][-1][1] + alignment_list[iloc][-1][3] - alignment_list[iloc][0][1]
                else:
                    ref_presize = alignment_list[iloc][0][1] + alignment_list[iloc][0][3] - alignment_list[iloc][-1][1]
                if(alignment_list[iloc+1][-1][2] == 1):
                    ref_nowsize = alignment_list[iloc+1][-1][1] + alignment_list[iloc+1][-1][3] - alignment_list[iloc+1][0][1]
                else:
                    ref_nowsize = alignment_list[iloc+1][0][1] + alignment_list[iloc+1][0][3] - alignment_list[iloc+1][-1][1]

                if(preitem[2] != nowitem[2] or (pos2contig(preitem[1], contig2start) != pos2contig(nowitem[1], contig2start))):
                    iloc += 1
                    skiped = True
                    break
                readgap = nowitem[0] - preitem[0] - preitem[3]
                if(preitem[2] == 1):
                    refgap = nowitem[1] - preitem[1] - preitem[3]

                else:
                    refgap = preitem[1]  - nowitem[1] - nowitem[3]
                if(refgap < 0):
                    if(refgap > -100):
                        if(len(alignment_list[iloc]) == 2):
                            alignment_list.pop(iloc)
                            skiped = True
                            break
                        if(len(alignment_list[iloc+1]) == 2):
                            alignment_list.pop(iloc + 1)
                            skiped = True
                            break
                        alignment_list[iloc].pop(-1)
                        alignment_list[iloc+1].pop(0)
                        continue
                if(presize < nowsize):
                    bias_ratio = nowsize/ref_nowsize  
                else:
                    bias_ratio = presize/ref_presize 
                if(min(refgap, readgap) > 100 and abs(readgap/bias_ratio - refgap) < 50):
                    target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, setting_kmersize, contig2seq, contig2start)
                    diffratio = edlib.align(query = query, target = target, task = 'distance')['editDistance']/min(len(target), len(query))
                    if(diffratio>maxdiffratio):
                        if(refgap > 200):
                            TRA_signal = True
                        iloc += 1
                        if(debug):
                            print_log('Large difference: readgap, refgap, diffratio', readgap, refgap, diffratio)
                        skiped = True
                        break
                break
            if(skiped == True):
                continue
                        
            if(refgap < 0):
                iloc += 1
                continue
            #if((min(readgap, refgap) < merge_smallgap) and (max(readgap, refgap) < too_large_gap)):
            if((H == False and ((min(presize, nowsize) > 500) or ((abs(readgap - refgap) / min(presize, nowsize)) < 0.5)) and (max(readgap, refgap) < 20000)) or (H == True and (((testoverlap(alignment_list, iloc) < 50) and (min(readgap, refgap)<100) and (abs(readgap - refgap) < 2000))))):
                alignment_list[iloc] = List_merge((alignment_list[iloc], alignment_list[iloc + 1]))
                alignment_list.pop(iloc+1)
            else:
                iloc += 1
    if(debug == True):print_log('step 4: merge conjacent alignment ', time.time() - st)
    if(debug == True):
        print_log('After merge conjacent alignment')
        for line in alignment_list:
            tempcontig = pos2contig(line[0][1], contig2start)
            temprefbias = contig2start[tempcontig]
            print_log(line[0][0], line[-1][0], line[0][1] - temprefbias, line[-1][1] - temprefbias)
    
    #step 5
    #split unrelate read
    #note inside alignment there is no way two conjacent anchors
    #has cost large than 30
    #so there is no big gap
    st = time.time()
    new_alignment_list = List()
    cigarlist = []
    for alignment in alignment_list: 

        tmp_alignment_list, tmp_cigarlist = split_alignment_test(alignment, testseq, rc_testseq, testseq_len, kmersize=setting_kmersize , contig2seq = contig2seq, contig2start = contig2start, debug = debug, H = H, eqx = eqx)
        if(debug): print_log(len(tmp_alignment_list), len(tmp_cigarlist))
        iloc = -1
        for alignment in tmp_alignment_list:
            iloc += 1
            new_alignment_list.append(alignment)
            cigarlist.append(tmp_cigarlist[iloc])


    if(debug == True):print_log('step 5: split unrelate read ', time.time() - st)
    if(debug == True):
        print_log('After split unrelate read')
        for line in new_alignment_list:
            tempcontig = pos2contig(line[0][1], contig2start)
            temprefbias = contig2start[tempcontig]
            print_log(line[0][0], line[-1][0], line[0][1] - temprefbias, line[-1][1] - temprefbias)
    
    alignment_list, onemapinfolist = get_onemapinfolist(new_alignment_list, cigarlist, readid, mapq, testseq_len, contig2start, need_reverse)
    return alignment_list, onemapinfolist, TRA_signal, filtered

def get_readmap_DP_test(readid, testseq, contig2start, contig2seq, index_object, index2contig, kmersize = 15, local_kmersize = 9, local_maxdiff = 30, refine = True, local_skipcost = 50., golbal_skipcost = (30., 30.),  golbal_maxdiff = (50, 30), check_num = 20, bin_size = 100, hastra = False, debug = False, H = False):
    local_skipcost = 30.    
    mapq = 60
    setting_kmersize = kmersize
    setting_maxdiff = golbal_maxdiff[1]
    if(H == False):
        local_skipcost += local_kmersize
        golbal_skipcost = (golbal_skipcost[0] + kmersize, golbal_skipcost[1] + kmersize)
    else:
        local_skipcost = 30.
        golbal_skipcost = (30., 30.)
        hastra = True
        check_num = 100
    

    rc_testseq = str(Seq(testseq).reverse_complement())
    testseq_len = len(testseq)
    st = time.time()                            
    mapq, scores, raw_alignment_list = decode_hit(index_object, index2contig, testseq, testseq_len, contig2start, kmersize, skipcost = golbal_skipcost, maxdiff = golbal_maxdiff, maxgap = 2000, check_num = check_num, c_bias = 5000, bin_size = bin_size, overlapprecentage = 0.5, hastra = hastra, H = H)
    #raw_alignment_list = check_func_clean(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, 10, debug = debug)
    if(debug == True):
        plot_result = np.array(raw_alignment_list)
        plt.scatter(plot_result[:, 0], plot_result[:, 1])
        plt.show()
    if(scores == 0.):
        return
    
    need_reverse = False
    if((scores < 0.)):
        need_reverse = True

    if(debug == True): print_log('Normal time', time.time() - st, scores)
    if(refine == True):
        setting_maxdiff = local_maxdiff
        setting_kmersize = local_kmersize
        st = time.time()

        #need_reverse, raw_alignment_array = get_reversed_chain_numpy_rough(np.array(raw_alignment_list), testseq_len)
        raw_alignment_array = np.array(raw_alignment_list)
        if(need_reverse == False):
            scores, raw_alignment_list = get_localmap_multi_all_forDP_inv_guide(raw_alignment_array, testseq, rc_testseq, contig2start, contig2seq, kmersize = setting_kmersize, skipcost = local_skipcost, maxdiff = setting_maxdiff, maxgap = 100, shift = 1)

        else:
            scores, raw_alignment_list = get_localmap_multi_all_forDP_inv_guide(raw_alignment_array, rc_testseq, testseq, contig2start, contig2seq, kmersize = setting_kmersize, skipcost = local_skipcost, maxdiff = setting_maxdiff, maxgap = 100, shift = 1)
            testseq, rc_testseq = rc_testseq, testseq
            #raw_alignment_list = get_reversed_chain_numpy(np.array(raw_alignment_list), testseq_len)

        #raw_alignment_list = check_func_clean(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, 10, debug = debug)
        #for item in raw_alignment_list:
            #print_log(item)

        if(debug == True): print_log('Refine time', time.time() - st)


        
        

    st = time.time()
    alignment_list, onemapinfolist, TRA_signal, filtered = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, debug = debug, H = H, nofilter = False)
    if(filtered == True and pairedindel(List([line[-1] for line in onemapinfolist]), indelsize = 30) == True):
        alignment_list, onemapinfolist, TRA_signal, filtered = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, debug = debug, H = H, nofilter = True)
        
    if(debug == True): print_log('extend_func', time.time() - st)
    if(len(onemapinfolist) > 1 and refine == 'auto'):
        setting_maxdiff = local_maxdiff
        setting_kmersize = local_kmersize
        st = time.time()                                                                         
        scores, raw_alignment_list = get_localmap_multi_all_forDP_inv_guide(np.array(raw_alignment_list), testseq, rc_testseq, contig2start, contig2seq, kmersize = setting_kmersize, skipcost = local_skipcost, maxdiff = setting_maxdiff, maxgap = 100, shift = 1)
        if(debug == True): print_log('Refine time', time.time() - st, len(raw_alignment_list))

        alignment_list, onemapinfolist, TRA_signal = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, need_reverse, setting_maxdiff, debug = debug)

    return onemapinfolist, (alignment_list, raw_alignment_list), TRA_signal

#print_log(1224)

@njit
def smallorequal2target(arr, target):
    n = len(arr)
    if(target < arr[0][0]):
        return -1

    if(target >= arr[n - 1][0]):
        return n-1

    i = 0
    j = n
    mid = 0
    while(i < j):
        mid = (i + j) // 2

        if(target == arr[mid][0]):
            if(arr[mid+1][0] > target):
                return mid
            else:
                i = mid + 1

        elif(target < arr[mid][0]) :

            if(mid > 0 and target >= arr[mid - 1][0]):
                return mid-1

            j = mid

        else: # target>arr[mid][0]
            if(mid < n - 1 and target < arr[mid + 1][0]):
                return mid

            i = mid + 1

    return mid

@njit
def smallorequal2target_1d(arr, target, n):
    if(target < arr[0]):
        return -1

    if(target >= arr[n - 1]):
        return n-1

    i = 0
    j = n
    mid = 0
    while(i < j):
        mid = (i + j) // 2

        if(target == arr[mid]):
            if(arr[mid+1] > target):
                return mid
            else:
                i = mid + 1

        elif(target < arr[mid]) :

            if(mid > 0 and target >= arr[mid - 1]):
                return mid-1

            j = mid

        else: # target>arr[mid][0]
            if(mid < n - 1 and target < arr[mid + 1]):
                return mid

            i = mid + 1

    return mid




@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20
    max_kmersize = 0
    
    #gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    #gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))


    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)

    

    

    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    readend_arr = one_mapinfo[:, 0] + one_mapinfo[:, 3]

    coverage_dict = np.zeros(one_mapinfo[-1][0], np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_st = 0
    testspace_en = 0
    
    ms_arr = np.empty(one_mapinfo[-1][0] + 1)
    
    test_S = np.empty(0, np.float64)
    
    for i in range(n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        
        max_kmersize = max(one_mapinfo[i][3], max_kmersize)
        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            if((opcount/one_mapinfo[i][0] > 2000) and one_mapinfo[i][0] > 1000):
                return 0., [(0, 0, 0, 0)]
            
            testspace_en = i
            
            #gap_arr[prereadloc] = 200+ 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            ms_arr[prereadloc: one_mapinfo[i][0]] = g_max_scores
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
            low_bound = 0
            loc_in_one_mapinfo = smallorequal2target(one_mapinfo, prereadloc - max_kmersize)#loose
            
            if(loc_in_one_mapinfo != -1):

                low_bound = ms_arr[one_mapinfo[loc_in_one_mapinfo][0]]-skipcost
                
                ms_arr_st = smallorequal2target_1d(ms_arr, low_bound-1e-4, prereadloc)

                if(ms_arr_st == -1):
                    
                    testspace_st = 0
                    
                else:
                    
                    testspace_st = smallorequal2target(one_mapinfo, ms_arr_st)+1
                       
            else:
                
                testspace_st = 0
                
            uesable = (prereadloc - readend_arr[testspace_st: testspace_en]) >= 0
            #print_log(one_mapinfo[testspace_st: testspace_en, 0].shape, uesable.shape)    
            test_S = S[testspace_st: testspace_en][uesable]
            #print_log(test_S.shape)    
            testspace = np.arange(testspace_st, testspace_en, dtype = np.int64)[uesable][np.argsort(test_S)][::-1]    
                      


        
        early_stop_score = 0.
        for j in testspace:
            
            
            opcount += 1     
            
            if(S[j] < max(max_scores - one_mapinfo[i][3], early_stop_score)):
                
                break

            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            #if((readgap < 0)):
                #print_log('error')
                #continue
                
            nocost = False
            filtered = True
            
            if(one_mapinfo[i][2] == 1):
                
                        #    y_2 - y_1 - a_1
                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                
            else:
                
                        #     y_1 - y_2 - a_2
                refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                
            
            while(True):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    
                    break
                    
                if((one_mapinfo[i][2] == 1) and (refgap < 0)):
                    
                    break
                    
                if(one_mapinfo[i][2] == -1):
                    
                    if(refgap > 0):
                        
                        break
                        
                    else:
                        
                        refgap = abs(refgap)

                gapcost = abs(readgap - refgap)
                
                if((readgap  > maxgap) or (gapcost > maxdiff)):
                #if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                    break
                    
                if(gapcost == 0):
                    
                    early_stop_score = S[j]
                    
                gapcost = gapcost_list[gapcost]
                
                test_scores = S[j] + one_mapinfo[i][3] - gapcost
                
                filtered = False
                
                break
            if(filtered == True):

                test_scores = S[j] - skipcost + one_mapinfo[i][3]
                '''gapcost = abs(readgap - abs(refgap))
                if(gapcost != 0):
                    test_scores = S[j] - skipcost + one_mapinfo[i][3] - min(np.log(abs(gapcost)), abs(gapcost)/100, 12)
                else:
                    test_scores = S[j] - skipcost + one_mapinfo[i][3]'''

            if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
            elif((test_scores == max_scores) and (pre_index != -9999999)):
                pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                if(now_refdistance < pre_refdistance):
                    max_scores = test_scores
                    pre_index = j

        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
            
        
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    while(True):
        if((P[take_index] == -9999999)):
            break
        take_index = P[take_index]
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

    return g_max_scores, path

@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20
    max_kmersize = 0
    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))

    var_gap = np.arange(1, 10001)
    var_gap = np.minimum(var_gap, np.log(var_gap))
    
    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)

    

    

    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    readend_arr = one_mapinfo[:, 0] + one_mapinfo[:, 3]

    coverage_dict = np.zeros(one_mapinfo[-1][0], np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_st = 0
    testspace_en = 0
    
    ms_arr = np.empty(one_mapinfo[-1][0] + 1)
    
    test_S = np.empty(0, np.float64)
    
    for i in range(n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        
        max_kmersize = max(one_mapinfo[i][3], max_kmersize)
        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            if((opcount/one_mapinfo[i][0] > 2000) and one_mapinfo[i][0] > 1000):
                return 0., [(0, 0, 0, 0)]
            
            testspace_en = i
            
            gap_arr[prereadloc] = 200+ 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            ms_arr[prereadloc: one_mapinfo[i][0]] = g_max_scores
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
            low_bound = 0
            loc_in_one_mapinfo = smallorequal2target(one_mapinfo, prereadloc - max_kmersize)#loose
            
            if(loc_in_one_mapinfo != -1):

                low_bound = ms_arr[one_mapinfo[loc_in_one_mapinfo][0]]-skipcost-10.
                
                ms_arr_st = smallorequal2target_1d(ms_arr, low_bound-1e-4, prereadloc)

                if(ms_arr_st == -1):
                    
                    testspace_st = 0
                    
                else:
                    
                    testspace_st = smallorequal2target(one_mapinfo, ms_arr_st)+1
                       
            else:
                
                testspace_st = 0
                
            uesable = (prereadloc - readend_arr[testspace_st: testspace_en]) >= 0
            #print_log(one_mapinfo[testspace_st: testspace_en, 0].shape, uesable.shape)    
            test_S = S[testspace_st: testspace_en][uesable]
            #print_log(test_S.shape)    
            testspace = np.arange(testspace_st, testspace_en, dtype = np.int64)[uesable][np.argsort(test_S)][::-1]    
                      


        
        early_stop_score = 0.
        for j in testspace:
            
            
            opcount += 1     
            
            if(S[j] < max(max_scores - one_mapinfo[i][3], early_stop_score)):
                
                break

            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            #if((readgap < 0)):
                #print_log('error')
                #continue
                
            nocost = False
            filtered = True
            
            if(one_mapinfo[i][2] == 1):
                
                        #    y_2 - y_1 - a_1
                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                
            else:
                
                        #     y_1 - y_2 - a_2
                refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                
            
            while(True):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    
                    break
                    
                if((one_mapinfo[i][2] == 1) and (refgap < 0)):
                    
                    break
                    
                if(one_mapinfo[i][2] == -1):
                    
                    if(refgap > 0):
                        
                        break
                        
                    else:
                        
                        refgap = abs(refgap)

                gapcost = abs(readgap - refgap)
                
                #if((readgap  > maxgap) or (gapcost > maxdiff)):
                if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                    break
                    
                if(gapcost == 0):
                    
                    early_stop_score = S[j]
                    
                gapcost = gapcost_list[gapcost]
                
                test_scores = S[j] + one_mapinfo[i][3] - gapcost
                
                filtered = False
                
                break
            if(filtered == True):

                #test_scores = S[j] - skipcost + one_mapinfo[i][3]
                gapcost = abs(readgap - abs(refgap))
                if(gapcost < 10000):

                    test_scores = S[j] - skipcost + one_mapinfo[i][3] - var_gap[gapcost]
                else:
                    test_scores = S[j] - skipcost + one_mapinfo[i][3] - 10.

            if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
            elif((test_scores == max_scores) and (pre_index != -9999999)):
                pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                if(now_refdistance < pre_refdistance):
                    max_scores = test_scores
                    pre_index = j

        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
            
        
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    while(True):
        if((P[take_index] == -9999999)):
            break
        take_index = P[take_index]
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

    return g_max_scores, path

#1226
@njit
def mergecigar_(cigarstring):

    numberrange = (ord('0'), ord('9'))
    
    oplist = []
    num = 0
    preop = '0'
    prenum = 0
    for ch in cigarstring:
        c = ord(ch)
        if(numberrange[0] <= c and numberrange[1] >= c):
            num = num * 10 + c - numberrange[0]
        else:
            if(preop == ch):
                prenum = prenum + num
                oplist[-2] = str(prenum)
                num = 0
            else:
                prenum = num
                oplist.append(str(num))
                oplist.append(ch)
                preop = ch
                num = 0
    return oplist
                
def mergecigar(cigarstring):  
    return ''.join(mergecigar_(cigarstring))
def fastq_q2b(qstring):
    b = bytearray()
    b.extend(map(ord, qstring))

    a = bytearray()
    a.extend(np.array(b)-33)
    a2 = array.array("B")
    a2.frombytes(a)
    return a2
def get_bam_dict(mapinfo, query, qual, contig2iloc, contig2seq):
    #'hhk',         ,  '1', '+', 11, 9192, 2767041, 2776138, 60
    #      0            1    2   3    4      5         6      7
    #'18_19897150_+', '18', '+', 0, 4776, 19832244, 19837393, 1]
    for iloc in range(len(mapinfo)):
        mapinfo[iloc] = list(mapinfo[iloc])
    '''for line in mapinfo:
        line[-1] = mergecigar(line[-1])
        if(len(Cigar(line[-1])) != len(query)):
            print_log()
            print_log(line[0], len(Cigar(line[-1])), len(query))
            print_log()
            return 
            break'''
    mq = mapinfo[-1][7]
    if(mq != 0):
        mq = 60
    else:
        mq = 1
    rc_query = str(Seq(query).reverse_complement())
    mapinfo.sort(key = sort_by_length)
    mapinfo = mapinfo[::-1]
    iloc2nm = dict()
    tmpiloc = -1
    for item in mapinfo:
        tmpiloc += 1
        if(item[2] == '+'):
            nm = compute_NM_tag(query[item[3]: item[4]], get_refseq(item[1], item[5], item[6], contig2seq))
        else:
            nm = compute_NM_tag(rc_query[item[3]: item[4]], get_refseq(item[1], item[5], item[6], contig2seq))
        iloc2nm[tmpiloc] = nm
    if((qual != None) and (len(qual) == len(query))):
        query_qualities = fastq_q2b(qual)
        rc_query_qualities = query_qualities[::-1]
    a_list = []
    for iloc in range(len(mapinfo)):
        bam_dict = dict()
        primary = mapinfo[iloc]
        bam_dict['readid'] = primary[0]
        bam_dict['contig'] = primary[1]
        if(iloc == 0):
            base_value = 0
        else:
            base_value = 2048
        if(primary[2] == '+'):
            bam_dict['flag'] = 0 + base_value

        else:
            bam_dict['flag'] = 16 + base_value


        bam_dict['refstart'] = primary[5]


        bam_dict['cigar'] = primary[8]

        if(len(mapinfo) > 1):
            salist = []
            tmpiloc = -1
            for item in mapinfo:
                tmpiloc += 1
                if(tmpiloc == iloc):
                    continue
                nm = iloc2nm[tmpiloc]
                salist.append(''.join((item[1], ',', str(item[5]+1), ',', item[2], ',', item[8], ',', str(mq), ',', str(nm)+';')))

            bam_dict['sa'] = ''.join(salist)


        a = pysam.AlignedSegment()
        a.query_name = bam_dict['readid']

        a.flag = bam_dict['flag']
        a.reference_id = contig2iloc[bam_dict['contig']]
        a.reference_start = bam_dict['refstart']
        item = primary

        a.mapping_quality = mq
        a.cigarstring = bam_dict['cigar']
        if(item[2] == '+'):
            a.query_sequence = query
            a.template_length = len(query)
            if((qual != None) and (len(qual) == len(query))):
                a.query_qualities = query_qualities
        else:
            a.query_sequence = rc_query
            a.template_length = len(rc_query)
            if((qual != None) and (len(qual) == len(query))):
                a.query_qualities = rc_query_qualities
        if('sa' in bam_dict):
            a.tags = [('NM', iloc2nm[iloc]), ("SA", bam_dict['sa'])]
        else:
            a.tags = [('NM', iloc2nm[iloc])]
        a_list.append(a)
    return a_list

def get_list_of_readmap_asm_1(raw_queue, savepath, minimap, contig2seq, hastra, H, header):
    st = time.time()
    
    contig2start = Dict()
    index2contig = List()
    contig2iloc = dict()
    iloc = -1
    for item in minimap.seq_offset:
        iloc += 1
        contig2start[item[0].decode()] = item[2]
        index2contig.append(item[0].decode())
        contig2iloc[item[0].decode()] = iloc

    
    iloc = 0
    unmapcountlist = []
    plotdata = []
    
    rt_list = []
    f_redo_ratio_list = []
    


    with open(savepath+'.temp', "w") as outf:
        while(True):
            readidandseq = raw_queue.get()
            if(type(readidandseq) == int):
                break
            

            tmp_st = time.time()
            try:
                if(H == True):
                    onemapinfolist, (alignment_list,raw_alignment_list), TRA_signal, f_redo_ratio = get_readmap_DP_test(readidandseq[0], readidandseq[1], contig2start, contig2seq, minimap, index2contig, hastra = hastra, H = True)
                else:
                    onemapinfolist, (alignment_list,raw_alignment_list), TRA_signal, f_redo_ratio = get_readmap_DP_test(readidandseq[0], readidandseq[1], contig2start, contig2seq, minimap, index2contig, hastra = False)
                    if(TRA_signal == True):
                        onemapinfolist, (alignment_list,raw_alignment_list), TRA_signal, f_redo_ratio = get_readmap_DP_test(readidandseq[0], readidandseq[1], contig2start, contig2seq, minimap, index2contig, hastra = True, check_num = 100)

            except:
                onemapinfolist = []
            if(len(onemapinfolist) == 0):
                onemapinfolist = []
                unmapcountlist.append(readidandseq[0])
                continue

            if(len(onemapinfolist) != 0):
                rt_list.append(time.time() - tmp_st)
                f_redo_ratio_list.append(f_redo_ratio)

                for a in onemapinfolist:
                    outf.write(str(a)+'\n')

            else:
                unmapcountlist.append(readidandseq[0])

    rt_list.sort()
    rt_list_len = len(rt_list)
    print_log(rt_list[0], rt_list[rt_list_len//2], rt_list[rt_list_len*99//100], rt_list[-1])
    rt_list = f_redo_ratio_list
    rt_list.sort()
    rt_list_len = len(rt_list)
    for iloc in list(range(rt_list_len))[::-1]:
        if(rt_list[iloc] < 5):
            break
    print_log(rt_list[0], rt_list[rt_list_len//2], rt_list[rt_list_len*99//100], rt_list[-1], iloc / rt_list_len)
    print_log(time.time() - st, 'unmapsize', len(unmapcountlist))

    print_log()

def get_list_of_readmap(raw_queue, savepath, minimap, contig2seq, hastra, H, header):
    st = time.time()
    redo_ratio = 5
    
    contig2start = Dict()
    index2contig = List()
    contig2iloc = dict()
    iloc = -1
    for item in minimap.seq_offset:
        iloc += 1
        contig2start[item[0].decode()] = item[2]
        index2contig.append(item[0].decode())
        contig2iloc[item[0].decode()] = iloc

    
    iloc = 0
    unmapcountlist = []
    plotdata = []
    
    rt_list = []
    f_redo_ratio_list = []
    


    rev_read_count = 0
    with pysam.AlignmentFile(savepath+'.bam', "wb", header=header) as outf:
        while(True):
            readidandseq = raw_queue.get()
            if(type(readidandseq) == int):
                break
            
            rev_read_count += 1
            tmp_st = time.time()
            try:
                if(H == True):
                    onemapinfolist, (alignment_list,raw_alignment_list), TRA_signal, f_redo_ratio = get_readmap_DP_test(readidandseq[0], readidandseq[1], contig2start, contig2seq, minimap, index2contig, hastra = hastra, H = True)
                else:
                    onemapinfolist, (alignment_list,raw_alignment_list), TRA_signal, f_redo_ratio = get_readmap_DP_test(readidandseq[0], readidandseq[1], contig2start, contig2seq, minimap, index2contig, hastra = False, redo_ratio = redo_ratio)
                    if(TRA_signal == True):
                        onemapinfolist, (alignment_list,raw_alignment_list), TRA_signal, f_redo_ratio = get_readmap_DP_test(readidandseq[0], readidandseq[1], contig2start, contig2seq, minimap, index2contig, hastra = True, check_num = 100, redo_ratio = redo_ratio)
            except:
                onemapinfolist = []
                


            if(len(onemapinfolist) == 0):
                onemapinfolist = []
                unmapcountlist.append(readidandseq[0])
                unmapcountlist.append(len(readidandseq[1]))
                continue

            if(len(onemapinfolist) != 0):
                rt_list.append(time.time() - tmp_st)
                if(rt_list[-1] > 10):
                    print_log(readidandseq[0], rt_list[-1])
                f_redo_ratio_list.append(f_redo_ratio)
                tmp_a_list = get_bam_dict(onemapinfolist, readidandseq[1], readidandseq[2], contig2iloc, contig2seq)
                if((tmp_a_list) == None):
                    print_log(readidandseq[0])
                    print_log()
                else:
                    for a in tmp_a_list:
                        outf.write(a)

            else:
                unmapcountlist.append(readidandseq[0])
    if(len(rt_list) > 0):
        rt_list.sort()
        rt_list_len = len(rt_list)
        print_log(rt_list[0], rt_list[rt_list_len//2], rt_list[rt_list_len*99//100], rt_list[-1])
        rt_list = f_redo_ratio_list
        rt_list.sort()
        rt_list_len = len(rt_list)
        for iloc in list(range(rt_list_len))[::-1]:
            if(rt_list[iloc] < redo_ratio):
                break
        print_log(rt_list[0], rt_list[rt_list_len//2], rt_list[rt_list_len*99//100], rt_list[-1], iloc / rt_list_len)

    print_log('mapped ', str(len(rt_list)) + '/' + str(rev_read_count))
    print_log(time.time() - st, 'unmapsize', len(unmapcountlist))
    print_log(unmapcountlist)
    print_log()
def get_bam_dict_stdout(mapinfo, query, qual, contig2iloc, contig2seq):
    #'hhk',         ,  '1', '+', 11, 9192, 2767041, 2776138, 60
    #      0            1    2   3    4      5         6      7
    #'18_19897150_+', '18', '+', 0, 4776, 19832244, 19837393, 1]
    for iloc in range(len(mapinfo)):
        mapinfo[iloc] = list(mapinfo[iloc])
    for line in mapinfo:
        line[-1] = mergecigar(line[-1])

    mq = mapinfo[-1][7]
    if(mq != 0):
        mq = 60
    else:
        mq = 1
    rc_query = str(Seq(query).reverse_complement())
    mapinfo.sort(key = sort_by_length)
    mapinfo = mapinfo[::-1]
    iloc2nm = dict()
    tmpiloc = -1
    for item in mapinfo:
        tmpiloc += 1
        if(item[2] == '+'):
            nm = compute_NM_tag(query[item[3]: item[4]], get_refseq(item[1], item[5], item[6], contig2seq))
        else:
            nm = compute_NM_tag(rc_query[item[3]: item[4]], get_refseq(item[1], item[5], item[6], contig2seq))
        iloc2nm[tmpiloc] = nm
    if((qual != None) and (len(qual) == len(query))):
        query_qualities = fastq_q2b(qual)
        rc_query_qualities = query_qualities[::-1]
    a_list = []
    for iloc in range(len(mapinfo)):
                    #  0     1     2     3     4     5     6      7      8     9     10   11  12
                    #QNAME FLAG  RNAME  POS  MAPQ  CIGAR RNEXT  PNEXT  TLEN   SEQ   QUAL  NM  SA
        infolist = [None,    None,   None,   None,  None,  None,   None,    None,     None,   None,   None, None, None]
        primary = mapinfo[iloc]
        infolist[0] = primary[0]
        infolist[2] = contig2iloc[primary[1]]
        infolist[4] = mq
        if(iloc == 0):
            base_value = 0
        else:
            base_value = 2048
        if(primary[2] == '+'):
            infolist[1] = base_value
            infolist[9] = query
            if((qual != None) and (len(qual) == len(query))):
                infolist[10] = query_qualities

        else:
            infolist[1] = 16 + base_value
            infolist[9] = rc_query
            if((qual != None) and (len(qual) == len(query))):
                infolist[10] = rc_query_qualities


        infolist[3] = primary[5]


        infolist[5] = primary[8]
        infolist[11] = iloc2nm[iloc]

        if(len(mapinfo) > 1):
            salist = []
            tmpiloc = -1
            for item in mapinfo:
                tmpiloc += 1
                if(tmpiloc == iloc):
                    continue
                nm = iloc2nm[tmpiloc]
                salist.append(''.join((item[1], ',', str(item[5]+1), ',', item[2], ',', item[8], ',', str(mq), ',', str(nm)+';')))

            infolist[12] = ''.join(salist)


        
        a_list.append(infolist)
    return a_list
def bam_dict_2_AlignedSegment(infolist):
    #  0     1     2     3     4     5     6      7      8     9     10   11  12
    #QNAME FLAG  RNAME  POS  MAPQ  CIGAR RNEXT  PNEXT  TLEN   SEQ   QUAL  NM  SA
    a = pysam.AlignedSegment()
    a.query_name = infolist[0]

    a.flag = infolist[1]
    a.reference_id = infolist[2]
    a.reference_start = infolist[3]

    a.mapping_quality = infolist[4]

    a.query_sequence = infolist[9]
    if(infolist[10] != None):
        a.query_qualities = infolist[10]

    if(infolist[12] != None):
        a.tags = [('NM', infolist[11]), ("SA", infolist[12])]
    else:
        a.tags = [('NM', infolist[11])]
    if(len(infolist[5]) > 65530):
        a.tags.append(('CG', infolist[5]))
    else:
        a.cigarstring = infolist[5]
        
    return a
def get_list_of_readmap_stdout(raw_queue, cooked_queue, minimap, contig2seq, hastra, H, header):
    st = time.time()
    redo_ratio = 5
    
    contig2start = Dict()
    index2contig = List()
    contig2iloc = dict()
    iloc = -1
    for item in minimap.seq_offset:
        iloc += 1
        contig2start[item[0].decode()] = item[2]
        index2contig.append(item[0].decode())
        contig2iloc[item[0].decode()] = iloc

    
    iloc = 0
    unmapcountlist = []
    plotdata = []
    
    rt_list = []
    f_redo_ratio_list = []
    


    rev_read_count = 0

    while(True):
        readidandseq = raw_queue.get()
        if(type(readidandseq) == int):
            break


        try:
            if(H == True):
                onemapinfolist, (alignment_list,raw_alignment_list), TRA_signal, f_redo_ratio = get_readmap_DP_test(readidandseq[0], readidandseq[1], contig2start, contig2seq, minimap, index2contig, hastra = hastra, H = True)
            else:
                onemapinfolist, (alignment_list,raw_alignment_list), TRA_signal, f_redo_ratio = get_readmap_DP_test(readidandseq[0], readidandseq[1], contig2start, contig2seq, minimap, index2contig, hastra = False, redo_ratio = redo_ratio)
                if(TRA_signal == True):
                    onemapinfolist, (alignment_list,raw_alignment_list), TRA_signal, f_redo_ratio = get_readmap_DP_test(readidandseq[0], readidandseq[1], contig2start, contig2seq, minimap, index2contig, hastra = True, check_num = 100, redo_ratio = redo_ratio)
        except:
            onemapinfolist = []



        if(len(onemapinfolist) == 0):

            continue

        if(len(onemapinfolist) != 0):

            tmp_a_list = get_bam_dict_stdout(onemapinfolist, readidandseq[1], readidandseq[2], contig2iloc, contig2seq)
            if((tmp_a_list) == None):
                continue
            else:
                for a in tmp_a_list:
                    #print_log(str(a))
                    cooked_queue.put(a)





    
def stdout_writer(cooked_queue, header):


    with pysam.AlignmentFile('-', "w", header=header) as outf:
        while(True):
            bam_dict = cooked_queue.get()
            if(type(bam_dict) == int):
                break
            try:
                outf.write(bam_dict_2_AlignedSegment(bam_dict))
            except:
                print_log(str(bam_dict))


#########

@njit
def string2int(x):
    num = 0
    zero = ord('0')
    for c in x:
        num = num * 10 + ord(c) - zero    
    return num
@njit
def get_MD_CSshort(info, target, query):
    mdstring = ''
    csstring = ''
    refloc = 0
    readloc = 0
    preop = ''
    equal_value = 0
    for iloc in range(1, len(info), 2):
        value = string2int(info[iloc - 1])
        if(info[iloc] == 'X'):
            if(equal_value > 0):
                mdstring += str(equal_value)
            elif(preop == 'D'):
                mdstring += '0'
            mdstring += target[refloc]
            csstring += '*' + (target[refloc] + query[readloc]).lower()
            for jloc in range(1, value):
                mdstring += '0' + target[refloc + jloc]
                csstring += '*' + (target[refloc + jloc] + query[readloc + jloc]).lower()
            refloc += value
            readloc += value
            equal_value = 0
        elif(info[iloc] == '='):
            #mdstring += info[iloc - 1]
            refloc += value
            readloc += value
            equal_value += value
            csstring += ':' + info[iloc - 1]
        elif(info[iloc] == 'D'):
            if(equal_value > 0):
                mdstring += str(equal_value)
            elif(preop == 'X'):
                mdstring += '0'
            mdstring += '^' + target[refloc: refloc + value]
            csstring += '-' + target[refloc: refloc + value].lower()
            refloc += value
            equal_value = 0
        elif(info[iloc] == 'I'):
            csstring += '+' + query[readloc: readloc + value].lower()
            readloc += value
            continue
        elif(info[iloc] == 'S'):
            continue
        else:
            return '', ''
        preop = info[iloc]
    if(equal_value > 0):
        mdstring += str(equal_value)
    return mdstring, csstring
@njit
def get_MD_CSlong(info, target, query):
    mdstring = ''
    csstring = ''
    refloc = 0
    readloc = 0
    preop = ''
    equal_value = 0
    for iloc in range(1, len(info), 2):
        value = string2int(info[iloc - 1])
        if(info[iloc] == 'X'):
            if(equal_value > 0):
                mdstring += str(equal_value)
            elif(preop == 'D'):
                mdstring += '0'
            mdstring += target[refloc]
            csstring += '*' + (target[refloc] + query[readloc]).lower()
            for jloc in range(1, value):
                mdstring += '0' + target[refloc + jloc]
                csstring += '*' + (target[refloc + jloc] + query[readloc + jloc]).lower()
            refloc += value
            readloc += value
            equal_value = 0
        elif(info[iloc] == '='):

            csstring += '=' + target[refloc: refloc + value].upper()
            refloc += value
            readloc += value
            equal_value += value
            
        elif(info[iloc] == 'D'):
            if(equal_value > 0):
                mdstring += str(equal_value)
            elif(preop == 'X'):
                mdstring += '0'
            mdstring += '^' + target[refloc: refloc + value]
            csstring += '-' + target[refloc: refloc + value].lower()
            refloc += value
            equal_value = 0
        elif(info[iloc] == 'I'):
            csstring += '+' + query[readloc: readloc + value].lower()
            readloc += value
            continue
        elif(info[iloc] == 'S'):
            continue
        else:
            return '', ''
        preop = info[iloc]
    if(equal_value > 0):
        mdstring += str(equal_value)
    return mdstring, csstring
@njit
def mergecigar_md_(cigarstring, target, query, shortcs = True):

    numberrange = (ord('0'), ord('9'))
    
    oplist = []
    num = 0
    preop = '0'
    prenum = 0
    for ch in cigarstring:
        c = ord(ch)
        if(numberrange[0] <= c and numberrange[1] >= c):
            num = num * 10 + c - numberrange[0]
        else:
            if(preop == ch):
                prenum = prenum + num
                oplist[-2] = str(prenum)
                num = 0
            else:
                prenum = num
                oplist.append(str(num))
                oplist.append(ch)
                preop = ch
                num = 0
    top = ''
    if(oplist[1] == 'S'):
        top = query[:string2int(oplist[0])].lower()
    tail = ''
    if(oplist[-1] == 'S'):
        tail = query[-string2int(oplist[-2]):].lower()
    if(shortcs == True):
        mdstring, csstring = get_MD_CSshort(oplist, target, query)
    else:
        mdstring, csstring = get_MD_CSlong(oplist, target, query)
    return oplist, mdstring, top+csstring+tail
def mergecigar_md_cs(cigarstring, target, query, shortcs = True):  
    oplist, mdstring, csstring = mergecigar_md_(cigarstring, target, query, shortcs)
    return ''.join(oplist), mdstring, csstring
def mergecigar_md_cs(cigarstring, target, query, shortcs = True):  

    oplist, mdstring, csstring = mergecigar_md_(cigarstring, target, query, shortcs)

    n_cigar = len(oplist)
    return ''.join(oplist), mdstring, csstring, n_cigar
def P_alignmentstring(infodict):
                #  0     1     2     3     4     5     6      7      8     9     10
                #QNAME FLAG  RNAME  POS  MAPQ  CIGAR RNEXT  PNEXT  TLEN   SEQ   QUAL
    infolist = ['*',   '4',   '*',  '0', '255', '*',  '*',   '0',  '0',   '*',   '*']
    name2iloc = {
        'QNAME': 0,
        'FLAG': 1,
        'RNAME': 2,
        'POS': 3,
        'MAPQ': 4,
        'CIGAR': 5,
        'RNEXT': 6,
        'PNEXT': 7,
        'TLEN': 8,
        'SEQ': 9,
        'QUAL': 10
    }
    def p_other_tag(tag, value):
        if(type(value) == int):
            code = 'i'
        elif(type(value) == float):
            code = 'f'
        elif(type(value) == str):
            code = 'Z'
        else:
            code = 'Z'
        return tag+':'+code+':'+str(value)

    for key in infodict:
        if(key in name2iloc):
            infolist[name2iloc[key]] = infodict[key]
        else:
            infolist.append(p_other_tag(key, infodict[key]))
    return '\t'.join(infolist)
def mergecigar_n(cigarstring):  
    oplist = mergecigar_(cigarstring)
    n_cigar = len(oplist)
    return ''.join(oplist), n_cigar
def get_bam_dict_str(mapinfo, query, qual, contig2iloc, contig2seq, md, shortcs, cigar2cg):
    #'hhk',         ,  '1', '+', 11, 9192, 2767041, 2776138, 60
    #      0            1    2   3    4      5         6      7
    #'18_19897150_+', '18', '+', 0, 4776, 19832244, 19837393, 1]
    for iloc in range(len(mapinfo)):
        mapinfo[iloc] = list(mapinfo[iloc])
    

    mq = mapinfo[-1][7]
    if(mq != 0):
        mq = 60
    else:
        mq = 1
    rc_query = str(Seq(query).reverse_complement())
    mapinfo.sort(key = sort_by_length)
    mapinfo = mapinfo[::-1]
    iloc2nm = dict()
    iloc2md = dict()
    iloc2cs = dict()
    iloc2n_cigar = dict()
    tmpiloc = -1
    if(md == False):
        for item in mapinfo:
            item[-1], n_cigar = mergecigar_n(item[-1])
            tmpiloc += 1
            if(item[2] == '+'):
                nm = compute_NM_tag(query[item[3]: item[4]], get_refseq(item[1], item[5], item[6], contig2seq))
            else:
                nm = compute_NM_tag(rc_query[item[3]: item[4]], get_refseq(item[1], item[5], item[6], contig2seq))
            iloc2nm[tmpiloc] = nm
            iloc2n_cigar[tmpiloc] = n_cigar
    else:
        for item in mapinfo:
            tmpiloc += 1
            if(item[2] == '+'):
                tmp_query = query[item[3]: item[4]]
                tmp_target = get_refseq(item[1], item[5], item[6], contig2seq)
            else:
                tmp_query = rc_query[item[3]: item[4]]
                tmp_target = get_refseq(item[1], item[5], item[6], contig2seq)
            cigarstring, mdstring, csstring, n_cigar = mergecigar_md_cs(item[-1], tmp_target, tmp_query, shortcs)
            nm = compute_NM_tag(tmp_query, tmp_target)
            item[-1] = cigarstring
            iloc2nm[tmpiloc] = nm
            iloc2md[tmpiloc] = mdstring
            iloc2cs[tmpiloc] = csstring
            iloc2n_cigar[tmpiloc] = n_cigar

    if((qual != None) and (len(qual) == len(query))):
        query_qualities = qual
        rc_query_qualities = query_qualities[::-1]
    a_list = []

    #QNAME FLAG  RNAME  POS  MAPQ  CIGAR RNEXT  PNEXT  TLEN   SEQ   QUAL
    for iloc in range(len(mapinfo)):
        bam_dict = dict()
        primary = mapinfo[iloc]
        bam_dict['QNAME'] = primary[0]
        bam_dict['RNAME'] = primary[1]
        if(iloc == 0):
            base_value = 0
        else:
            base_value = 2048
        if(primary[2] == '+'):
            bam_dict['FLAG'] = str(base_value)

        else:
            bam_dict['FLAG'] = str(16 + base_value)


        bam_dict['POS'] = str(primary[5] + 1)# SAM Format

        if(iloc2n_cigar[iloc] > 65535):
            if(cigar2cg == True):
                bam_dict['CG'] = primary[8]
                logging.info('Write long CIGAR to CG tag.')
            else:
                bam_dict['CIGAR'] = primary[8]
        else:
            bam_dict['CIGAR'] = primary[8]

        if(len(mapinfo) > 1):
            salist = []
            tmpiloc = -1
            for item in mapinfo:
                tmpiloc += 1
                if(tmpiloc == iloc):
                    continue
                nm = iloc2nm[tmpiloc]
                salist.append(''.join((item[1], ',', str(item[5]+1), ',', item[2], ',', item[8], ',', str(mq), ',', str(nm)+';')))

            bam_dict['SA'] = ''.join(salist)



        item = primary

        bam_dict['MAPQ'] = str(mq)

        if(item[2] == '+'):
            bam_dict['SEQ'] = query

            if((qual != None) and (len(qual) == len(query))):
                bam_dict['QUAL'] = query_qualities
        else:
            bam_dict['SEQ'] = rc_query

            if((qual != None) and (len(qual) == len(query))):
                bam_dict['QUAL'] = rc_query_qualities
        bam_dict['NM'] = iloc2nm[iloc]
        if(md == True):
            bam_dict['MD'] = iloc2md[iloc]
            bam_dict['cs'] = iloc2cs[iloc]

        a_list.append(P_alignmentstring(bam_dict))
    return a_list

def get_list_of_readmap_stdout(raw_queue, cooked_queue, minimap, contig2seq, hastra, H, header, option):
    cache_size = 100
    a_list = []
    st = time.time()
    redo_ratio = 5
    
    contig2start = Dict()
    index2contig = List()
    contig2iloc = dict()
    iloc = -1
    for item in minimap.seq_offset:
        iloc += 1
        contig2start[item[0].decode()] = item[2]
        index2contig.append(item[0].decode())
        contig2iloc[item[0].decode()] = iloc

    
    iloc = 0
    unmapcountlist = []
    plotdata = []
    
    rt_list = []
    f_redo_ratio_list = []
    


    rev_read_count = 0

    while(True):
        readidandseq = raw_queue.get()
        if(type(readidandseq) == int):
            break


        try:
            if(H == True):
                onemapinfolist, (alignment_list,raw_alignment_list), TRA_signal, f_redo_ratio = get_readmap_DP_test(readidandseq[0], readidandseq[1], contig2start, contig2seq, minimap, index2contig, hastra = hastra, H = True, eqx = option['eqx'])
            else:
                onemapinfolist, (alignment_list,raw_alignment_list), TRA_signal, f_redo_ratio = get_readmap_DP_test(readidandseq[0], readidandseq[1], contig2start, contig2seq, minimap, index2contig, hastra = False, redo_ratio = redo_ratio, eqx = option['eqx'])
                if(TRA_signal == True):
                    onemapinfolist, (alignment_list,raw_alignment_list), TRA_signal, f_redo_ratio = get_readmap_DP_test(readidandseq[0], readidandseq[1], contig2start, contig2seq, minimap, index2contig, hastra = True, check_num = 100, redo_ratio = redo_ratio, eqx = option['eqx'])
        except:
            onemapinfolist = []
            continue





        if(len(onemapinfolist) != 0):
            try:
                tmp_a_list = get_bam_dict_str(onemapinfolist, readidandseq[1], readidandseq[2], contig2iloc, contig2seq, option['md'], option['shortcs'], option['cigar2cg'])
            except:
                logging.info('Failed to convert to sam record')
                logging.info(readidandseq[0])
                logging.info(readidandseq[1])
                continue
            if((tmp_a_list) == None):
                continue
            else:
                a_list += tmp_a_list
                if(len(a_list) > cache_size):
                    cooked_queue.put(a_list)
                    a_list = []
                else:
                    continue
    if(len(a_list) > 0):
        cooked_queue.put(a_list)

def stdout_writer(cooked_queue, header):


    with pysam.AlignmentFile('-', "w", header=header, threads = 16) as outf:
        while(True):
            a_list = cooked_queue.get()
            if(type(a_list) == int):
                break
            for a in a_list:
                print(a)
########


def decode_hit(index_object, index2contig, testseq, testseq_len, contig2start, kmersize, skipcost = (50., 30.), maxdiff = (50, 30), maxgap = 2000, check_num = 20, c_bias = 5000, bin_size = 100, overlapprecentage = 0.5, hastra = True, H = False):            
    need_reverse, one_mapinfo = get_reversed_chain_numpy_rough(np.array(index_object.map(testseq)), testseq_len)
    mapq, scores, path =  hit2work(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage, hastra, H = H)
    if(need_reverse == True):
        return mapq, -scores, path
    else:
        return mapq, scores, path
    
def get_readmap_DP_test1(readid, testseq, contig2start, contig2seq, index_object, index2contig, kmersize = 15, local_kmersize = 9, local_maxdiff = 30, refine = True, local_skipcost = 50., golbal_skipcost = (30., 30.),  golbal_maxdiff = (50, 30), check_num = 20, bin_size = 100, hastra = False, debug = False, H = False):
    local_skipcost = 30.    
    mapq = 60
    setting_kmersize = kmersize
    setting_maxdiff = golbal_maxdiff[1]
    if(H == False):
        local_skipcost += local_kmersize
        golbal_skipcost = (golbal_skipcost[0] + kmersize, golbal_skipcost[1] + kmersize)
    else:
        local_skipcost = 30.
        golbal_skipcost = (30., 30.)
        hastra = True
        check_num = 100
    

    rc_testseq = str(Seq(testseq).reverse_complement())
    testseq_len = len(testseq)
    st = time.time()                            
    mapq, scores, raw_alignment_list = decode_hit(index_object, index2contig, testseq, testseq_len, contig2start, kmersize, skipcost = golbal_skipcost, maxdiff = golbal_maxdiff, maxgap = 2000, check_num = check_num, c_bias = 5000, bin_size = bin_size, overlapprecentage = 0.5, hastra = hastra, H = H)
    #raw_alignment_list = check_func_clean(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, 10, debug = debug)
    if(debug == True):
        plot_result = np.array(raw_alignment_list)
        plt.scatter(plot_result[:, 0], plot_result[:, 1])
        plt.show()
    if(scores == 0.):
        return [], ([], []), 0
    
    need_reverse = False
    if((scores < 0.)):
        need_reverse = True

    if(debug == True): print_log('Normal time', time.time() - st, scores)
    if(refine == True):
        setting_maxdiff = local_maxdiff
        setting_kmersize = local_kmersize
        st = time.time()

        #need_reverse, raw_alignment_array = get_reversed_chain_numpy_rough(np.array(raw_alignment_list), testseq_len)
        raw_alignment_array = np.array(raw_alignment_list)
        if(need_reverse == False):
            scores, raw_alignment_list = get_localmap_multi_all_forDP_inv_guide(raw_alignment_array, testseq, rc_testseq, contig2start, contig2seq, kmersize = setting_kmersize, skipcost = local_skipcost, maxdiff = setting_maxdiff, maxgap = 100, shift = 1)

        else:
            scores, raw_alignment_list = get_localmap_multi_all_forDP_inv_guide(raw_alignment_array, rc_testseq, testseq, contig2start, contig2seq, kmersize = setting_kmersize, skipcost = local_skipcost, maxdiff = setting_maxdiff, maxgap = 100, shift = 1)
            testseq, rc_testseq = rc_testseq, testseq
            #raw_alignment_list = get_reversed_chain_numpy(np.array(raw_alignment_list), testseq_len)

        #raw_alignment_list = check_func_clean(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, 10, debug = debug)
        #for item in raw_alignment_list:
            #print_log(item)

        if(debug == True): print_log('Refine time', time.time() - st)


        
        

    st = time.time()
    try:
        filtered = False
        #alignment_list, onemapinfolist, TRA_signal = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, debug = debug, H = H)
        alignment_list, onemapinfolist, TRA_signal, filtered = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, debug = debug, H = H, nofilter = False)
        if(len(onemapinfolist) == 0):
            return onemapinfolist, (alignment_list, raw_alignment_list), TRA_signal
    except:
        print_log(testseq)
    
    

    if(filtered == True and pairedindel(List([line[-1] for line in onemapinfolist]), indelsize = 30) == True):
        alignment_list, onemapinfolist, TRA_signal, filtered = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, debug = debug, H = H, nofilter = True)
        
    if(debug == True): print_log('extend_func', time.time() - st)
    if(len(onemapinfolist) > 1 and refine == 'auto'):
        setting_maxdiff = local_maxdiff
        setting_kmersize = local_kmersize
        st = time.time()                                                                         
        scores, raw_alignment_list = get_localmap_multi_all_forDP_inv_guide(np.array(raw_alignment_list), testseq, rc_testseq, contig2start, contig2seq, kmersize = setting_kmersize, skipcost = local_skipcost, maxdiff = setting_maxdiff, maxgap = 100, shift = 1)
        if(debug == True): print_log('Refine time', time.time() - st, len(raw_alignment_list))

        alignment_list, onemapinfolist, TRA_signal = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, need_reverse, setting_maxdiff, debug = debug)

    return onemapinfolist, (alignment_list, raw_alignment_list), TRA_signal
@njit
def pairedindel(cigarlist, indelsize = 30):


    zero = ord('0')
    INSsyb = ord('I') - zero
    SOFTsyb = ord('S') - zero
    HARDsyb = ord('H') - zero
    PADsyb = ord('P') - zero
    DELsyb = ord('D') - zero
    delsyb = ord('^') - zero

    indel = []

    for cigar in cigarlist:
        
        
        typed_cigar = [ord(item) - zero for item in cigar] 
        
        number = 0.
        for item in typed_cigar:
            if(item < 10):
                number = number * 10. + item
  
            else:
                if(item != INSsyb and item != SOFTsyb and item != HARDsyb and item != PADsyb):
                    if(item == DELsyb and number > indelsize):
                        indel.append(number)
                    number = 0.

                else:

                    if(item == INSsyb and number > indelsize):
                        indel.append(number)
                        
                    number = 0.
    indel.sort()
    preitem = 0
    clustersize = 1
    for nowitem in indel:
        if((min(preitem, nowitem)/max(preitem, nowitem))>0.7):
            clustersize += 1
            if(clustersize > 1):
                return True
        else:
            clustersize = 1
        preitem = nowitem
    return False


@njit
def hit2work(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage = 0.5, hastra = False, H = False):
    def getsecond(x):
        return x[1]
    def getfirst(x):
        return x[0]
    def getlength(x):
        return len(x)

    def get_readloc_set_bin(one_mappos, bin_size):
        return set([i[0] // bin_size for i in one_mappos])
    def get_refloc_set_bin(one_mappos, bin_size):
        return set([i[1] // bin_size for i in one_mappos])
    def get_overlapsize(readloc_set_a, readloc_set_b):
        return len(readloc_set_a&readloc_set_b)/min(len(readloc_set_a), len(readloc_set_b))
    
    maxgap = 200

    
    one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:,1])] 

    cluster_list = List()
    iloc = 0
    preitem = one_mapinfo[iloc]
    st_iloc = 0
    for iloc in range(len(one_mapinfo)):
        nowitem = one_mapinfo[iloc]
        if(((nowitem[1] - preitem[1]) > c_bias)):
            if((iloc - st_iloc) < 3):
                continue
            cluster_list.append(one_mapinfo[st_iloc: iloc])
            st_iloc = iloc
        preitem = nowitem


    if((iloc - st_iloc) > 3):
        cluster_list.append(one_mapinfo[st_iloc: iloc + 1])
    cluster_list.sort(key = get_length)

    cluster_list = cluster_list[::-1][:check_num]


    
    hit = False
    minichain_scores = 40
    path_list = List()
    scores_list = []
    scores_list.append(0.)
    scores_list.pop()
    path_list.append([(0, 0, 0, 0)])
    path_list.pop()
    

    max_scores = 0
    from_repeat = False

    for one_mapinfo in cluster_list:
        repeat = False
        min_ref, max_ref = one_mapinfo[0][1], one_mapinfo[-1][1]
        testcontig = pos2contig(min_ref, contig2start)
        if(testcontig != pos2contig(max_ref, contig2start)):
            #print_log('hit2work: testcontig != pos2contig(max_ref, contig2start)')
            continue

        scores, path = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d(one_mapinfo[np.argsort(one_mapinfo[:,0])], kmersize = kmersize, skipcost = skipcost[0], maxdiff = maxdiff[0], maxgap = maxgap)
         
        if(scores == 0.):
            scores, path = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_fast(one_mapinfo[np.argsort(one_mapinfo[:,0])], kmersize = kmersize, skipcost = skipcost[0], maxdiff = maxdiff[0], maxgap = maxgap)
        if(scores < minichain_scores):
            if(hit == False):
                return 0, 0., [(0, 0, 0, 0)]
            continue
        
        

        hit = True

        path_list.append(path) 
        scores_list.append(scores)
        if(scores > max_scores):
            if(repeat == True):
                from_repeat = True
            else:
                from_repeat = False
            max_scores = scores
    
    if(hit == True and max_scores > 200):

        order = np.argsort(np.array(scores_list))[::-1]
        
        primary_rlocset_List = List()
        primary_scores_List = List()
        primary_index_List = List()
        
        iloc = order[0]
        primary_rlocset_List.append(get_readloc_set_bin(path_list[iloc], bin_size))
        primary_scores_List.append(List([scores_list[iloc]]))
        primary_index_List.append(iloc)

        for iloc in order[1:]:
            readloc_set_b = get_readloc_set_bin(path_list[iloc], bin_size)
            maxoverlapsize = 0.
            for p_loc in range(len(primary_rlocset_List)):
                tmp_overlapsize = get_overlapsize(primary_rlocset_List[p_loc], readloc_set_b)
                if(tmp_overlapsize > maxoverlapsize):
                    maxoverlapsize = tmp_overlapsize
                    prefer_p_loc = p_loc
            if(maxoverlapsize < overlapprecentage):
                primary_rlocset_List.append(readloc_set_b)
                primary_scores_List.append(List([scores_list[iloc]]))
                primary_index_List.append(iloc)
            else:
                primary_scores_List[prefer_p_loc].append(scores_list[iloc])
                
        m = len(path_list[order[0]])    
        if(len(primary_scores_List[0]) < 2):
            f1 = primary_scores_List[0][0]
            f2 = 0
        else:
            f1 = primary_scores_List[0][0]
            f2 = primary_scores_List[0][1]
        mapq = min(int(40*(1-f2/f1)*min(1, m/10)*np.log(f1)), 60)

        base_iloc = primary_index_List[0]

        for iloc in primary_index_List[1:4]:
            if(hastra == False):
                skip = True
                item = path_list[iloc][0]
                bias_1 = abs(item[1]-path_list[base_iloc][0][1])
                bias_2 = abs(item[1] - path_list[base_iloc][-1][1])
                if(min(bias_1, bias_2) < 200000):
                    if(pos2contig(item[1], contig2start) == pos2contig(path_list[base_iloc][0][1], contig2start)):
                        skip = False
                if(skip == True):
                    continue
            for item in path_list[iloc]:                            
                path_list[base_iloc].append(item)

        path_list[base_iloc].sort(key = getfirst)
        scores, path = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d(np.array(path_list[base_iloc]), kmersize = kmersize, skipcost = skipcost[1], maxdiff = maxdiff[1], maxgap = maxgap)

        return mapq, scores, path


    else:
        return 0, 0., [(0, 0, 0, 0)]
    
@njit
def get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start):#
    
    query_st, query_en = preitem[0], nowitem[0]
    if(preitem[2] == 1):
        testcontig = pos2contig(preitem[1], contig2start)
        refbias = contig2start[testcontig]
        target_st, target_en = preitem[1], nowitem[1]
        query = testseq[preitem[0]: nowitem[0]]
        target = contig2seq[testcontig][preitem[1] - refbias: nowitem[1] - refbias]

    else:
        testcontig = pos2contig(nowitem[1], contig2start)
        refbias = contig2start[testcontig]
        target_st, target_en = nowitem[1] + nowitem[3], preitem[1] + preitem[3]
        query = rc_testseq[testseq_len - nowitem[0]: testseq_len - preitem[0]]
        target = contig2seq[testcontig][nowitem[1] + nowitem[3] - refbias: preitem[1] + preitem[3] - refbias]
    return target, query, target_st, target_en, query_st, query_en
@njit
def rebuild_chain_break(contig2start, raw_alignment_list, large_cost, small_alignment = 50, small_dup = -100):
    #rebuild chain break
    #step 1
    #remove anchor cause large cost and small alignment
    preitem = raw_alignment_list[0]
    alignment_list = List([List([preitem])])
    for nowitem in raw_alignment_list[1:]:
        if(preitem[2] == nowitem[2]):
            readgap = nowitem[0] - preitem[0] - preitem[3]


            if(preitem[2] == 1):
                refgap = nowitem[1] - preitem[1] - preitem[3]
            else:
                refgap = preitem[1]  - nowitem[1] - nowitem[3]

            if((refgap < 0) and (refgap > small_dup)):
                continue
            if((abs(readgap - refgap) <= large_cost) and (refgap >= 0) and (readgap < 100)):
                if(pos2contig(preitem[1], contig2start) == pos2contig(nowitem[1], contig2start)):
                    alignment_list[-1].append(nowitem)
                    preitem = nowitem
                    continue

        if(len(alignment_list[-1]) == 1):
            alignment_list.pop(-1)
        if(len(alignment_list) > 0):
            if((alignment_list[-1][-1][0] + alignment_list[-1][-1][3] - alignment_list[-1][0][0]) < small_alignment):
                alignment_list.pop(-1)
        if(len(alignment_list) > 0):
            if(len(alignment_list[-1]) >= 4):#
                alignment_list[-1] = alignment_list[-1][1:-1]
        alignment_list.append(List([nowitem]))
        preitem = nowitem
    if(len(alignment_list[-1]) == 1):
        alignment_list.pop(-1)
    if((alignment_list[-1][-1][0] + alignment_list[-1][-1][3] - alignment_list[-1][0][0]) < small_alignment):
        alignment_list.pop(-1)
    if(len(alignment_list[-1]) >= 4):
        alignment_list[-1] = alignment_list[-1][1:-1]

    return alignment_list

@njit
def smallorequal2target(arr, target):
    n = len(arr)
    if(target < arr[0][0]):
        return -1

    if(target >= arr[n - 1][0]):
        return n-1

    i = 0
    j = n
    mid = 0
    while(i < j):
        mid = (i + j) // 2

        if(target == arr[mid][0]):
            if(arr[mid+1][0] > target):
                return mid
            else:
                i = mid + 1

        elif(target < arr[mid][0]) :

            if(mid > 0 and target >= arr[mid - 1][0]):
                return mid-1

            j = mid

        else: # target>arr[mid][0]
            if(mid < n - 1 and target < arr[mid + 1][0]):
                return mid

            i = mid + 1

    return mid

@njit
def smallorequal2target_1d(arr, target, n):
    if(target < arr[0]):
        return -1

    if(target >= arr[n - 1]):
        return n-1

    i = 0
    j = n
    mid = 0
    while(i < j):
        mid = (i + j) // 2

        if(target == arr[mid]):
            if(arr[mid+1] > target):
                return mid
            else:
                i = mid + 1

        elif(target < arr[mid]) :

            if(mid > 0 and target >= arr[mid - 1]):
                return mid-1

            j = mid

        else: # target>arr[mid][0]
            if(mid < n - 1 and target < arr[mid + 1]):
                return mid

            i = mid + 1

    return mid




@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20
    max_kmersize = 0
    
    #gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    #gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))


    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)

    

    

    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    readend_arr = one_mapinfo[:, 0] + one_mapinfo[:, 3]

    coverage_dict = np.zeros(one_mapinfo[-1][0], np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_st = 0
    testspace_en = 0
    
    ms_arr = np.empty(one_mapinfo[-1][0] + 1)
    
    test_S = np.empty(0, np.float64)
    
    for i in range(n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        
        max_kmersize = max(one_mapinfo[i][3], max_kmersize)
        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            #if((opcount/one_mapinfo[i][0] > 2000) and one_mapinfo[i][0] > 1000):
                #return 0., [(0, 0, 0, 0)]
            
            testspace_en = i
            
            #gap_arr[prereadloc] = 200+ 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            ms_arr[prereadloc: one_mapinfo[i][0]] = g_max_scores
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
            low_bound = 0
            loc_in_one_mapinfo = smallorequal2target(one_mapinfo, prereadloc - max_kmersize)#loose
            
            if(loc_in_one_mapinfo != -1):

                low_bound = ms_arr[one_mapinfo[loc_in_one_mapinfo][0]]-skipcost
                
                ms_arr_st = smallorequal2target_1d(ms_arr, low_bound-1e-4, prereadloc)

                if(ms_arr_st == -1):
                    
                    testspace_st = 0
                    
                else:
                    
                    testspace_st = smallorequal2target(one_mapinfo, ms_arr_st)+1
                       
            else:
                
                testspace_st = 0
                
            uesable = (prereadloc - readend_arr[testspace_st: testspace_en]) >= 0
            #print_log(one_mapinfo[testspace_st: testspace_en, 0].shape, uesable.shape)    
            test_S = S[testspace_st: testspace_en][uesable]
            #print_log(test_S.shape)    
            testspace = np.arange(testspace_st, testspace_en, dtype = np.int64)[uesable][np.argsort(test_S)][::-1]    
                      


        
        early_stop_score = 0.
        for j in testspace:
            
            
            #opcount += 1     
            
            if(S[j] < max(max_scores - one_mapinfo[i][3], early_stop_score)):
                
                break

            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            #if((readgap < 0)):
                #print_log('error')
                #continue
                
            nocost = False
            filtered = True
            
            if(one_mapinfo[i][2] == 1):
                
                        #    y_2 - y_1 - a_1
                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                
            else:
                
                        #     y_1 - y_2 - a_2
                refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                
            
            while(True):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    
                    break
                    
                if((one_mapinfo[i][2] == 1) and (refgap < 0)):
                    
                    break
                    
                if(one_mapinfo[i][2] == -1):
                    
                    if(refgap > 0):
                        
                        break
                        
                    else:
                        
                        refgap = abs(refgap)

                gapcost = abs(readgap - refgap)
                
                if((readgap  > maxgap) or (gapcost > maxdiff)):
                #if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                    break
                    
                if(gapcost == 0):
                    
                    early_stop_score = S[j]
                    
                gapcost = gapcost_list[gapcost]
                
                test_scores = S[j] + one_mapinfo[i][3] - gapcost
                
                filtered = False
                
                break
            if(filtered == True):

                test_scores = S[j] - skipcost + one_mapinfo[i][3]
                '''gapcost = abs(readgap - abs(refgap))
                if(gapcost != 0):
                    test_scores = S[j] - skipcost + one_mapinfo[i][3] - min(np.log(abs(gapcost)), abs(gapcost)/100, 12)
                else:
                    test_scores = S[j] - skipcost + one_mapinfo[i][3]'''

            if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
            elif((test_scores == max_scores) and (pre_index != -9999999)):
                pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                if(now_refdistance < pre_refdistance):
                    max_scores = test_scores
                    pre_index = j

        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
            
        
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    while(True):
        if((P[take_index] == -9999999)):
            break
        take_index = P[take_index]
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

    return g_max_scores, path

@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20
    max_kmersize = 0
    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))

    var_gap = np.arange(1, 10001)
    var_gap = np.minimum(var_gap/100, np.log(var_gap))
    
    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)

    

    

    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    readend_arr = one_mapinfo[:, 0] + one_mapinfo[:, 3]

    coverage_dict = np.zeros(one_mapinfo[-1][0], np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_st = 0
    testspace_en = 0
    
    ms_arr = np.empty(one_mapinfo[-1][0] + 1)
    
    test_S = np.empty(0, np.float64)
    
    for i in range(n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        
        max_kmersize = max(one_mapinfo[i][3], max_kmersize)
        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            if((opcount/one_mapinfo[i][0] > 2000) and one_mapinfo[i][0] > 1000):
                return 0., [(0, 0, 0, 0)]
            
            testspace_en = i
            
            gap_arr[prereadloc] = 200+ 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            ms_arr[prereadloc: one_mapinfo[i][0]] = g_max_scores
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
            low_bound = 0
            loc_in_one_mapinfo = smallorequal2target(one_mapinfo, prereadloc - max_kmersize)#loose
            
            if(loc_in_one_mapinfo != -1):

                low_bound = ms_arr[one_mapinfo[loc_in_one_mapinfo][0]]-skipcost-10.
                
                ms_arr_st = smallorequal2target_1d(ms_arr, low_bound-1e-4, prereadloc)

                if(ms_arr_st == -1):
                    
                    testspace_st = 0
                    
                else:
                    
                    testspace_st = smallorequal2target(one_mapinfo, ms_arr_st)+1
                       
            else:
                
                testspace_st = 0
                
            uesable = (prereadloc - readend_arr[testspace_st: testspace_en]) >= 0
            #print_log(one_mapinfo[testspace_st: testspace_en, 0].shape, uesable.shape)    
            test_S = S[testspace_st: testspace_en][uesable]
            #print_log(test_S.shape)    
            testspace = np.arange(testspace_st, testspace_en, dtype = np.int64)[uesable][np.argsort(test_S)][::-1]    
                      


        
        early_stop_score = 0.
        for j in testspace:
            
            
            opcount += 1     
            
            if(S[j] < max(max_scores - one_mapinfo[i][3], early_stop_score)):
                
                break

            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            #if((readgap < 0)):
                #print_log('error')
                #continue
                
            nocost = False
            filtered = True
            
            if(one_mapinfo[i][2] == 1):
                
                        #    y_2 - y_1 - a_1
                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                
            else:
                
                        #     y_1 - y_2 - a_2
                refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                
            
            while(True):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    
                    break
                    
                if((one_mapinfo[i][2] == 1) and (refgap < 0)):
                    
                    break
                    
                if(one_mapinfo[i][2] == -1):
                    
                    if(refgap > 0):
                        
                        break
                        
                    else:
                        
                        refgap = abs(refgap)

                gapcost = abs(readgap - refgap)
                
                #if((readgap  > maxgap) or (gapcost > maxdiff)):
                if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                    break
                    
                if(gapcost == 0):
                    
                    early_stop_score = S[j]
                    
                gapcost = gapcost_list[gapcost]
                
                test_scores = S[j] + one_mapinfo[i][3] - gapcost
                
                filtered = False
                
                break
            if(filtered == True):

                #test_scores = S[j] - skipcost + one_mapinfo[i][3]
                gapcost = abs(readgap - abs(refgap))
                if(gapcost < 10000):

                    test_scores = S[j] - skipcost + one_mapinfo[i][3] - var_gap[gapcost]
                else:
                    test_scores = S[j] - skipcost + one_mapinfo[i][3] - 10.

            if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
            elif((test_scores == max_scores) and (pre_index != -9999999)):
                pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                if(now_refdistance < pre_refdistance):
                    max_scores = test_scores
                    pre_index = j

        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
            
        
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    while(True):
        if((P[take_index] == -9999999)):
            break
        take_index = P[take_index]
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

    return g_max_scores, path


@njit
def smallorequal2target_1d_point(arr, target, n, point):
    if(target < arr[point[0]]):
        return -1

    if(target >= arr[point[n - 1]]):
        return n-1

    i = 0
    j = n
    mid = 0
    while(i < j):
        mid = (i + j) // 2

        if(target == arr[point[mid]]):
            if(mid < n - 1):
                if(arr[point[mid+1]] > target):
                    return mid
                else:
                    i = mid + 1
            else:
                return mid

        elif(target < arr[point[mid]]) :

            if(mid > 0 and target >= arr[point[mid - 1]]):
                return mid-1

            j = mid

        else: # target>arr[mid][0]
            if(mid < n - 1 and target < arr[point[mid + 1]]):
                return mid

            i = mid + 1

    return mid
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20
    max_kmersize = 0
    
    #gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    #gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
    
    #var_gap = np.arange(1, 10001)
    #var_gap = np.minimum(var_gap/100, np.log(var_gap))


    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)

    

    

    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    readend_arr = one_mapinfo[:, 0] + one_mapinfo[:, 3]

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_st = 0
    testspace_en = 1
    
    ms_arr = np.empty(one_mapinfo[-1][0] + 1)
    
    test_S = np.empty(0, np.float64)
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0
    max_kmersize = one_mapinfo[i][3]
    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        
        max_kmersize = max(one_mapinfo[i][3], max_kmersize)
        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            #if((opcount/one_mapinfo[i][0] > 2000) and one_mapinfo[i][0] > 1000):
                #return 0., [(0, 0, 0, 0)]
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i
            
            #gap_arr[prereadloc] = 200+ 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
               
                      


        

        for j in S_arg[:testspace_en][::-1]:
            
            
            #opcount += 1     
            
            if(S[j] < (max_scores - one_mapinfo[i][3])):
                
                break


            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            if((readgap < 0)):
                #print_log('error')
                continue
                
            nocost = False
            filtered = True
            
            if(one_mapinfo[i][2] == 1):
                
                        #    y_2 - y_1 - a_1
                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                
            else:
                
                        #     y_1 - y_2 - a_2
                refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                
            
            while(True):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    
                    break
                    
                if((one_mapinfo[i][2] == 1) and (refgap < 0)):
                    
                    break
                    
                if(one_mapinfo[i][2] == -1):
                    
                    if(refgap > 0):
                        
                        break
                        
                    else:
                        
                        refgap = abs(refgap)

                gapcost = abs(readgap - refgap)
                
                if((readgap  > maxgap) or (gapcost > maxdiff)):
                #if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                    break
                    

                    
                gapcost = gapcost_list[gapcost]
                
                test_scores = S[j] + one_mapinfo[i][3] - gapcost
                
                filtered = False
                
                break
            if(filtered == True):

                test_scores = S[j] - skipcost + one_mapinfo[i][3]
                '''gapcost = abs(readgap - abs(refgap))
                if(gapcost < 10000):

                    test_scores = S[j] - skipcost + one_mapinfo[i][3] - var_gap[gapcost]
                else:
                    test_scores = S[j] - skipcost + one_mapinfo[i][3] - 10.'''

            if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
            elif((test_scores == max_scores) and (pre_index != -9999999)):
                pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                if(now_refdistance < pre_refdistance):
                    max_scores = test_scores
                    pre_index = j

        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i

        
            
        
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    while(True):
        if((P[take_index] == -9999999)):
            break
        take_index = P[take_index]
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

    return g_max_scores, path

@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20
    max_kmersize = 0
    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
    
    var_gap = np.arange(1, 10001)#mark 1
    var_gap = np.minimum(var_gap/100, np.log(var_gap))#mark 1


    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)

    

    

    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    readend_arr = one_mapinfo[:, 0] + one_mapinfo[:, 3]

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_st = 0
    testspace_en = 0
    
    ms_arr = np.empty(one_mapinfo[-1][0] + 1)
    
    test_S = np.empty(0, np.float64)
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0
    max_kmersize = one_mapinfo[i][3]
    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        
        max_kmersize = max(one_mapinfo[i][3], max_kmersize)
        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            if((opcount/one_mapinfo[i][0] > 2000) and one_mapinfo[i][0] > 1000):#mark 1
                return 0., [(0, 0, 0, 0)]#mark 1
            
            k = testspace_en
            while(one_mapinfo[k][0] + kmersize <= one_mapinfo[i][0]):
                if(k == 0):
                    k += 1
                    continue
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1

            testspace_en = k
            
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
               
                      


        
        for j in S_arg[:testspace_en][::-1]:
            
            
            opcount += 1     #mark 1
            
            if(S[j] < (max_scores - one_mapinfo[i][3])):
                
                break


            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            if((readgap < 0)):
                #print_log('error')
                continue
                
            nocost = False
            filtered = True
            
            if(one_mapinfo[i][2] == 1):
                
                        #    y_2 - y_1 - a_1
                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                
            else:
                
                        #     y_1 - y_2 - a_2
                refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                
            
            while(True):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    
                    break
                    
                if((one_mapinfo[i][2] == 1) and (refgap < 0)):
                    
                    break
                    
                if(one_mapinfo[i][2] == -1):
                    
                    if(refgap > 0):
                        
                        break
                        
                    else:
                        
                        refgap = abs(refgap)

                gapcost = abs(readgap - refgap)
                
                #if((readgap  > maxgap) or (gapcost > maxdiff)):
                if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                    break
                    
                    
                gapcost = gapcost_list[gapcost]
                
                test_scores = S[j] + one_mapinfo[i][3] - gapcost
                
                filtered = False
                
                break
            if(filtered == True):

                #test_scores = S[j] - skipcost + one_mapinfo[i][3]
                gapcost = abs(readgap - abs(refgap))#mark 1
                if(gapcost < 10000):#mark 1

                    test_scores = S[j] - skipcost + one_mapinfo[i][3] - var_gap[gapcost]#mark 1
                else:
                    test_scores = S[j] - skipcost + one_mapinfo[i][3] - 10.#mark 1

            if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
            elif((test_scores == max_scores) and (pre_index != -9999999)):
                pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                if(now_refdistance < pre_refdistance):
                    max_scores = test_scores
                    pre_index = j

        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i

        
            
        
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    while(True):
        if((P[take_index] == -9999999)):
            break
        take_index = P[take_index]
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

    return g_max_scores, path

@njit
def smallorequal2target_1d_point(arr, target, n, point):
    if(target < arr[point[0]]):
        return -1

    if(target >= arr[point[n - 1]]):
        return n-1

    i = 0
    j = n
    mid = 0
    while(i < j):
        mid = (i + j) // 2

        if(target == arr[point[mid]]):
            if(mid < n - 1):
                if(arr[point[mid+1]] > target):
                    return mid
                else:
                    i = mid + 1
            else:
                return mid

        elif(target < arr[point[mid]]) :
            


            if(target >= arr[point[mid - 1]]):
                return mid-1

            j = mid

        else: # target>arr[point[mid]]

            if(target < arr[point[mid + 1]]):
                return mid

            i = mid + 1

    return mid


@njit
def smallorequal2target_1d_point_target(arr, target, n, point, fixpos, target_arr):
    
    if(target < arr[point[0]]):
        
        return -1
    if(target == arr[point[0]]):
        if(fixpos < target_arr[point[0]]):
            return -1

    if(target > arr[point[n - 1]]):
        
        return n-1
    
    if(target == arr[point[n - 1]]):
        if(fixpos >= target_arr[point[n - 1]]):
            return n - 1
    

    i = 0
    j = n
    mid = 0
    while(i < j):
        
        mid = (i + j) // 2

        if(target == arr[point[mid]]):
            if(fixpos > target_arr[point[mid]]):
                if(mid < n - 1):
                    if(target < arr[point[mid + 1]]):
                        return mid
                    else:#target == arr[point[mid + 1]
                        if(fixpos < target_arr[point[mid + 1]]):
                            #print_log('fixpos <= target_arr[point[mid + 1]]')
                            return mid
                        else:
                            i = mid + 1
                else:
                    #print_log('mid == n - 1')
                    return mid
            elif(fixpos < target_arr[point[mid]]):
                if(mid > 0):
                    if(arr[point[mid - 1]] < target):
                        #print_log('arr[point[mid - 1]] < target')
                        return mid-1
                    else:#arr[point[mid - 1]] == target
                        if(fixpos >= target_arr[point[mid - 1]]):
                            #print_log('fixpos >= target_arr[point[mid - 1]]')
                            return mid - 1
                        else:#fixpos < target_arr[point[mid - 1]]
                            j = mid

                else:
                    #print_log('mid == 0')
                    return -1
            else:
                if(mid < n - 1):
                    if(target < arr[point[mid + 1]]):
                        return mid
                    else:#target == arr[point[mid + 1]
                        if(fixpos < target_arr[point[mid + 1]]):
                            #print_log('fixpos <= target_arr[point[mid + 1]]')
                            return mid
                        else:
                            i = mid + 1
                #print_log('fixpos < target_arr[point[mid]]')
                else:
                    return mid
                    


        elif(target < arr[point[mid]]) :

            if(mid > 0):
                if(target > arr[point[mid - 1]]):
                    return mid-1



            j = mid

        else: # target > arr[point[mid]]
            if(mid < n - 1):
                if(target < arr[point[mid + 1]]):
                    return mid


            i = mid + 1

    return mid

@njit
def closest2target_1d_point_pos(arr, target, st, en, point, pos, pos_arr):
    if(target < arr[point[st]]):
        return st

    if(target > arr[point[en - 1]]):
        
        return en-1

    i = st
    j = en
    mid = st
    while(i < j):
        mid = (i + j) // 2

        if(target == arr[point[mid]]):
            if(pos_arr[point[mid]] < pos):
                if(mid + 1 < en):
                    if(target == arr[point[mid + 1]]):
                        if(pos_arr[point[mid + 1]] > pos):
                            if(pos - pos_arr[point[mid]] < pos_arr[point[mid + 1]] - pos):
                                return mid
                            else:
                                return mid + 1
                        else:#pos > pos_arr[point[mid + 1]]
                            i = mid + 1
                    else:
                        return mid
                else:
                    return mid
            else:#pos_arr[point[mid]] > pos
                if(mid > st):
                    if(target == arr[point[mid - 1]]):
                        if(pos > pos_arr[point[mid - 1]]):
                            if(pos - pos_arr[point[mid - 1]] < pos_arr[point[mid]] - pos):
                                return mid - 1
                            else:
                                return mid
                        else:#pos < pos_arr[point[mid - 1]]
                            j = mid
                    else:
                        return mid
                else:
                    return mid
                
                
            
            
            

        elif(target < arr[point[mid]]) :
            


            if(target >= arr[point[mid - 1]]):
                if(target - arr[point[mid - 1]] < arr[point[mid]] - target):
                    return mid-1
                else:
                    return mid

            j = mid

        else: # target>arr[point[mid]]

            if(target <= arr[point[mid + 1]]):
                if(arr[point[mid + 1]] - target > target - arr[point[mid]]):
                    return mid
                else:
                    return mid + 1

            i = mid + 1

    return mid

@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_fast(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20
    max_kmersize = 0
    
    #gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    #gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))//1
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))//1
    
    #var_gap = np.arange(1, 10001)#mark 1
    #var_gap = np.minimum(var_gap//100, np.log(var_gap)//1)#mark 1

    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)

    

    

    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    readend_arr = one_mapinfo[:, 0] + one_mapinfo[:, 3]

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    readlength = one_mapinfo[-1][0] + 1000
    target_arr = np.zeros(n, dtype = np.float64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
        if(one_mapinfo[i][2] == 1):
            target_arr[i] = one_mapinfo[i][1] - one_mapinfo[i][0] + readlength 
        else:
            target_arr[i] = -(one_mapinfo[i][1] + one_mapinfo[i][0] + readlength) 
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_st = 0
    testspace_en = 0
    
    ms_arr = np.empty(one_mapinfo[-1][0] + 1)
    
    test_S = np.empty(0, np.float64)
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0
    max_kmersize = one_mapinfo[i][3]
    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    
    cache_dict = dict()
    cache_dict[0] = 0
    cache_dict.pop(0)
    
    fast = False
    
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        
        max_kmersize = max(one_mapinfo[i][3], max_kmersize)
        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            #if((opcount/one_mapinfo[i][0] > 2000) and one_mapinfo[i][0] > 1000):
                #return 0., [(0, 0, 0, 0)]
                
            k = testspace_en
            while(one_mapinfo[k][0] + kmersize <= one_mapinfo[i][0]):
                if(k == 0):
                    k += 1
                    continue
                loc_in_sorted_S = smallorequal2target_1d_point_target(S, S[k], k, S_arg, target_arr[k], target_arr) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1

            testspace_en = k
            
            
            #gap_arr[prereadloc] = 200+ 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
            
            
            fast = False
            if(testspace_en > 0):
                temp_endpos = (testspace_en - max(testspace_en - 1000, 0))

                scorerange = S[S_arg[testspace_en - 1]] - S[S_arg[max(testspace_en - 1000, 0)]] + 1
                
                if(temp_endpos / scorerange > 20):
                    fast = True
                    cache_dict = dict()
                    cache_dict[0] = 0
                    cache_dict.pop(0)
            
            
               
                      
        if(fast == True):

            #print_log(i, prereadloc)

            st_loc = testspace_en
            en_loc = testspace_en

            while(en_loc > 0):
                if(st_loc not in cache_dict):
                    pre_st_loc = st_loc
                    st_loc = smallorequal2target_1d_point(S, S[S_arg[st_loc-1]]-1e-7, en_loc, S_arg) + 1
                    cache_dict[pre_st_loc] = st_loc
                else:
                    st_loc = cache_dict[st_loc]

                if(S[S_arg[st_loc]] < (max_scores - one_mapinfo[i][3])):

                    break
                #print_log(S[S_arg[st_loc:en_loc]])
                #print_log(target_arr[S_arg[st_loc:en_loc]])
                #print_log(st_loc, en_loc)
                #print_log()

                #S_arg_j = closest2target_1d_point(target_arr, target_arr[i], st_loc, en_loc, S_arg)
                S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg, one_mapinfo[i][0], one_mapinfo[:,0])




                j = S_arg[S_arg_j]

                opcount += 1     



                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])



                nocost = False
                filtered = True

                if(one_mapinfo[i][2] == 1):

                            #    y_2 - y_1 - a_1
                    refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])

                else:

                            #     y_1 - y_2 - a_2
                    refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])


                while(True):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):

                        break

                    if((one_mapinfo[i][2] == 1) and (refgap < 0)):

                        break

                    if(one_mapinfo[i][2] == -1):

                        if(refgap > 0):

                            break

                        else:

                            refgap = abs(refgap)

                    gapcost = abs(readgap - refgap)

                    if(gapcost > maxdiff):

                        break

                    if(readgap  > maxgap):
                    #if(readgap  > gap_arr[one_mapinfo[j][0]]):#mark 1                   
                        break



                    gapcost = gapcost_list[gapcost]

                    test_scores = S[j] + one_mapinfo[i][3] - gapcost

                    filtered = False


                    break
                if(filtered == True):

                    test_scores = S[j] - skipcost + one_mapinfo[i][3]
                    '''gapcost = abs(readgap - abs(refgap))
                    if(gapcost < 10000):

                        test_scores = S[j] - skipcost + one_mapinfo[i][3] - var_gap[gapcost]
                    else:
                        test_scores = S[j] - skipcost + one_mapinfo[i][3] - 10.'''

                if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j
                elif((test_scores == max_scores) and (pre_index != -9999999)):
                    pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                    now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                    if(now_refdistance < pre_refdistance):
                        max_scores = test_scores
                        pre_index = j






                en_loc = st_loc
        
        else:
            
            for j in S_arg[:testspace_en][::-1]:
            
            
                opcount += 1     

                if(S[j] < (max_scores - one_mapinfo[i][3])):

                    break


                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                if((readgap < 0)):
                    #print_log('error')
                    continue

                nocost = False
                filtered = True

                if(one_mapinfo[i][2] == 1):

                            #    y_2 - y_1 - a_1
                    refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])

                else:

                            #     y_1 - y_2 - a_2
                    refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])


                while(True):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):

                        break

                    if((one_mapinfo[i][2] == 1) and (refgap < 0)):

                        break

                    if(one_mapinfo[i][2] == -1):

                        if(refgap > 0):

                            break

                        else:

                            refgap = abs(refgap)

                    gapcost = abs(readgap - refgap)

                    if((readgap  > maxgap) or (gapcost > maxdiff)):
                    #if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                        break


                    gapcost = gapcost_list[gapcost]

                    test_scores = S[j] + one_mapinfo[i][3] - gapcost

                    filtered = False

                    break
                if(filtered == True):

                    test_scores = S[j] - skipcost + one_mapinfo[i][3]
                    '''gapcost = abs(readgap - abs(refgap))
                    if(gapcost < 10000):

                        test_scores = S[j] - skipcost + one_mapinfo[i][3] - var_gap[gapcost]
                    else:
                        test_scores = S[j] - skipcost + one_mapinfo[i][3] - 10.'''

                if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j
                elif((test_scores == max_scores) and (pre_index != -9999999)):
                    pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                    now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                    if(now_refdistance < pre_refdistance):
                        max_scores = test_scores
                        pre_index = j

            S[i] = max_scores
            P[i] = pre_index

            if(max_scores > g_max_scores):

                g_max_scores = max_scores
                g_max_index = i
            
            
            
            
            
                


        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
        

    k = testspace_en
    while(k < (i+1)):
        if(k == 0):
            k += 1
            continue
        loc_in_sorted_S = smallorequal2target_1d_point_target(S, S[k], k, S_arg, target_arr[k], target_arr) + 1

        S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
        S_arg[loc_in_sorted_S] = k

        k += 1     
        
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    while(True):
        if((P[take_index] == -9999999)):
            break
        take_index = P[take_index]
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    #print_log(opcount)
    return g_max_scores, path

@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_fast(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20
    max_kmersize = 0
    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))//1
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))//1
    
    var_gap = np.arange(1, 10001)#mark 1
    var_gap = np.minimum(var_gap//100, np.log(var_gap)//1)#mark 1

    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)

    

    

    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    readend_arr = one_mapinfo[:, 0] + one_mapinfo[:, 3]

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    readlength = one_mapinfo[-1][0] + 1000
    target_arr = np.zeros(n, dtype = np.float64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
        if(one_mapinfo[i][2] == 1):
            target_arr[i] = one_mapinfo[i][1] - one_mapinfo[i][0] + readlength 
        else:
            target_arr[i] = -(one_mapinfo[i][1] + one_mapinfo[i][0] + readlength) 
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_st = 0
    testspace_en = 0
    
    ms_arr = np.empty(one_mapinfo[-1][0] + 1)
    
    test_S = np.empty(0, np.float64)
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0
    max_kmersize = one_mapinfo[i][3]
    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    
    cache_dict = dict()
    cache_dict[0] = 0
    cache_dict.pop(0)
    
    fast = False
    
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        
        max_kmersize = max(one_mapinfo[i][3], max_kmersize)
        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            #if((opcount/one_mapinfo[i][0] > 2000) and one_mapinfo[i][0] > 1000):
                #return 0., [(0, 0, 0, 0)]
                
            k = testspace_en
            while(one_mapinfo[k][0] + kmersize <= one_mapinfo[i][0]):
                if(k == 0):
                    k += 1
                    continue
                loc_in_sorted_S = smallorequal2target_1d_point_target(S, S[k], k, S_arg, target_arr[k], target_arr) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1

            testspace_en = k
            
            
            gap_arr[prereadloc] = 200+ 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
            
            
            fast = False
            if(testspace_en > 0):
                temp_endpos = (testspace_en - max(testspace_en - 1000, 0))

                scorerange = S[S_arg[testspace_en - 1]] - S[S_arg[max(testspace_en - 1000, 0)]] + 1
                
                if(temp_endpos / scorerange > 20):
                    fast = True
                    cache_dict = dict()
                    cache_dict[0] = 0
                    cache_dict.pop(0)
            
            
               
                      
        if(fast == True):

            #print_log(i, prereadloc)

            st_loc = testspace_en
            en_loc = testspace_en

            while(en_loc > 0):
                if(st_loc not in cache_dict):
                    pre_st_loc = st_loc
                    st_loc = smallorequal2target_1d_point(S, S[S_arg[st_loc-1]]-1e-7, en_loc, S_arg) + 1
                    cache_dict[pre_st_loc] = st_loc
                else:
                    st_loc = cache_dict[st_loc]

                if(S[S_arg[st_loc]] < (max_scores - one_mapinfo[i][3])):

                    break

                S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg, one_mapinfo[i][0], one_mapinfo[:,0])




                j = S_arg[S_arg_j]

                opcount += 1     



                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])



                nocost = False
                filtered = True

                if(one_mapinfo[i][2] == 1):

                            #    y_2 - y_1 - a_1
                    refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])

                else:

                            #     y_1 - y_2 - a_2
                    refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])


                while(True):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):

                        break

                    if((one_mapinfo[i][2] == 1) and (refgap < 0)):

                        break

                    if(one_mapinfo[i][2] == -1):

                        if(refgap > 0):

                            break

                        else:

                            refgap = abs(refgap)

                    gapcost = abs(readgap - refgap)

                    if(gapcost > maxdiff):

                        break

                    #if(readgap  > maxgap):
                    if(readgap  > gap_arr[one_mapinfo[j][0]]):#mark 1                   
                        break



                    gapcost = gapcost_list[gapcost]

                    test_scores = S[j] + one_mapinfo[i][3] - gapcost

                    filtered = False


                    break
                if(filtered == True):

                    #test_scores = S[j] - skipcost + one_mapinfo[i][3]
                    gapcost = abs(readgap - abs(refgap))
                    if(gapcost < 10000):

                        test_scores = S[j] - skipcost + one_mapinfo[i][3] - var_gap[gapcost]
                    else:
                        test_scores = S[j] - skipcost + one_mapinfo[i][3] - 10.

                if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j
                elif((test_scores == max_scores) and (pre_index != -9999999)):
                    pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                    now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                    if(now_refdistance < pre_refdistance):
                        max_scores = test_scores
                        pre_index = j






                en_loc = st_loc
        
        else:
            
            for j in S_arg[:testspace_en][::-1]:
            
            
                opcount += 1     

                if(S[j] < (max_scores - one_mapinfo[i][3])):

                    break


                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                if((readgap < 0)):
                    #print_log('error')
                    continue

                nocost = False
                filtered = True

                if(one_mapinfo[i][2] == 1):

                            #    y_2 - y_1 - a_1
                    refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])

                else:

                            #     y_1 - y_2 - a_2
                    refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])


                while(True):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):

                        break

                    if((one_mapinfo[i][2] == 1) and (refgap < 0)):

                        break

                    if(one_mapinfo[i][2] == -1):

                        if(refgap > 0):

                            break

                        else:

                            refgap = abs(refgap)

                    gapcost = abs(readgap - refgap)

                    #if((readgap  > maxgap) or (gapcost > maxdiff)):
                    if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                        break


                    gapcost = gapcost_list[gapcost]

                    test_scores = S[j] + one_mapinfo[i][3] - gapcost

                    filtered = False

                    break
                if(filtered == True):

                    #test_scores = S[j] - skipcost + one_mapinfo[i][3]
                    gapcost = abs(readgap - abs(refgap))
                    if(gapcost < 10000):

                        test_scores = S[j] - skipcost + one_mapinfo[i][3] - var_gap[gapcost]
                    else:
                        test_scores = S[j] - skipcost + one_mapinfo[i][3] - 10.

                if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j
                elif((test_scores == max_scores) and (pre_index != -9999999)):
                    pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                    now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                    if(now_refdistance < pre_refdistance):
                        max_scores = test_scores
                        pre_index = j

            S[i] = max_scores
            P[i] = pre_index

            if(max_scores > g_max_scores):

                g_max_scores = max_scores
                g_max_index = i
            
            
            
            
            
                


        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
        

    k = testspace_en
    while(k < (i+1)):
        if(k == 0):
            k += 1
            continue
        loc_in_sorted_S = smallorequal2target_1d_point_target(S, S[k], k, S_arg, target_arr[k], target_arr) + 1

        S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
        S_arg[loc_in_sorted_S] = k

        k += 1     
        
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    while(True):
        if((P[take_index] == -9999999)):
            break
        take_index = P[take_index]
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    #print_log(opcount)
    return g_max_scores, path
#20240127
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_all(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20
    max_kmersize = 0
    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
    
    var_gap = np.arange(1, 10001)#mark 1
    var_gap = np.minimum(var_gap/100, np.log(var_gap))#mark 1


    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)

    

    

    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    readend_arr = one_mapinfo[:, 0] + one_mapinfo[:, 3]

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_st = 0
    testspace_en = 1
    
    ms_arr = np.empty(one_mapinfo[-1][0] + 1)
    
    test_S = np.empty(0, np.float64)
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0
    max_kmersize = one_mapinfo[i][3]
    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        
        max_kmersize = max(one_mapinfo[i][3], max_kmersize)
        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            if((opcount/one_mapinfo[i][0] > 2000) and one_mapinfo[i][0] > 1000):#mark 1
                return -1, S, P, S_arg#mark 1
            
            k = testspace_en
            while(one_mapinfo[k][0] + kmersize <= one_mapinfo[i][0]):
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1

            testspace_en = k
            
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
               
                      


        
        for j in S_arg[:testspace_en][::-1]:
            
            
            opcount += 1     #mark 1
            
            if(S[j] < (max_scores - one_mapinfo[i][3])):
                
                break


            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            if((readgap < 0)):
                #print_log('error')
                continue
                
            nocost = False
            filtered = True
            
            if(one_mapinfo[i][2] == 1):
                
                        #    y_2 - y_1 - a_1
                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                
            else:
                
                        #     y_1 - y_2 - a_2
                refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                
            
            while(True):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    
                    break
                    
                if((one_mapinfo[i][2] == 1) and (refgap < 0)):
                    
                    break
                    
                if(one_mapinfo[i][2] == -1):
                    
                    if(refgap > 0):
                        
                        break
                        
                    else:
                        
                        refgap = abs(refgap)

                gapcost = abs(readgap - refgap)
                
                #if((readgap  > maxgap) or (gapcost > maxdiff)):
                if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                    break
                    
                    
                gapcost = gapcost_list[gapcost]
                
                test_scores = S[j] + one_mapinfo[i][3] - gapcost
                
                filtered = False
                
                break
            if(filtered == True):

                #test_scores = S[j] - skipcost + one_mapinfo[i][3]
                gapcost = abs(readgap - abs(refgap))#mark 1
                if(gapcost < 10000):#mark 1

                    test_scores = S[j] - skipcost + one_mapinfo[i][3] - var_gap[gapcost]#mark 1
                else:
                    test_scores = S[j] - skipcost + one_mapinfo[i][3] - 10.#mark 1

            if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
            elif((test_scores == max_scores) and (pre_index != -9999999)):
                pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                if(now_refdistance < pre_refdistance):
                    max_scores = test_scores
                    pre_index = j

        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i

        
    k = testspace_en
    while(k < n):

        loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

        S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
        S_arg[loc_in_sorted_S] = k

        k += 1

    testspace_en = k  
    
    

        
    
    return g_max_index, S, P, S_arg
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_fast_all(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20
    max_kmersize = 0
    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))//1
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))//1
            
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = gapcost
        else:
            gapcost_list[gapcost] = gapcost
    
    var_gap = np.arange(1, 10001)#mark 1
    var_gap = np.minimum(var_gap//100, np.log(var_gap)//1)#mark 1

    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)

    

    

    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    readend_arr = one_mapinfo[:, 0] + one_mapinfo[:, 3]

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    readlength = one_mapinfo[-1][0] + 1000
    target_arr = np.zeros(n, dtype = np.float64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
        if(one_mapinfo[i][2] == 1):
            target_arr[i] = one_mapinfo[i][1] - one_mapinfo[i][0] + readlength 
        else:
            target_arr[i] = -(one_mapinfo[i][1] + one_mapinfo[i][0] + readlength) 
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_st = 0
    testspace_en = 1
    
    ms_arr = np.empty(one_mapinfo[-1][0] + 1)
    
    test_S = np.empty(0, np.float64)
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0
    max_kmersize = one_mapinfo[i][3]
    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    
    cache_dict = dict()
    cache_dict[0] = 0
    cache_dict.pop(0)
    
    fast = False
    
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        
        max_kmersize = max(one_mapinfo[i][3], max_kmersize)
        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            #if((opcount/one_mapinfo[i][0] > 2000) and one_mapinfo[i][0] > 1000):
                #return 0., [(0, 0, 0, 0)]
                
            k = testspace_en
            while(one_mapinfo[k][0] + kmersize <= one_mapinfo[i][0]):

                loc_in_sorted_S = smallorequal2target_1d_point_target(S, S[k], k, S_arg, target_arr[k], target_arr) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1

            testspace_en = k
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
            
            
            fast = False
            if(testspace_en > 0):
                temp_endpos = (testspace_en - max(testspace_en - 1000, 0))

                scorerange = S[S_arg[testspace_en - 1]] - S[S_arg[max(testspace_en - 1000, 0)]] + 1
                
                if(temp_endpos / scorerange > 20):
                    fast = True
                    cache_dict = dict()
                    cache_dict[0] = 0
                    cache_dict.pop(0)
            
            
               
                      
        if(fast == True):

            #print_log(i, prereadloc)

            st_loc = testspace_en
            en_loc = testspace_en

            while(en_loc > 0):
                if(st_loc not in cache_dict):
                    pre_st_loc = st_loc
                    st_loc = smallorequal2target_1d_point(S, S[S_arg[st_loc-1]]-1e-7, en_loc, S_arg) + 1
                    cache_dict[pre_st_loc] = st_loc
                else:
                    st_loc = cache_dict[st_loc]

                if(S[S_arg[st_loc]] < (max_scores - one_mapinfo[i][3])):

                    break

                S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg, one_mapinfo[i][0], one_mapinfo[:,0])




                j = S_arg[S_arg_j]

                opcount += 1     



                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])



                nocost = False
                filtered = True

                if(one_mapinfo[i][2] == 1):

                            #    y_2 - y_1 - a_1
                    refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])

                else:

                            #     y_1 - y_2 - a_2
                    refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])


                while(True):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):

                        break

                    if((one_mapinfo[i][2] == 1) and (refgap < 0)):

                        break

                    if(one_mapinfo[i][2] == -1):

                        if(refgap > 0):

                            break

                        else:

                            refgap = abs(refgap)

                    gapcost = abs(readgap - refgap)

                    if(gapcost > maxdiff):

                        break

                    #if(readgap  > maxgap):
                    if(readgap  > gap_arr[one_mapinfo[j][0]]):#mark 1                   
                        break



                    gapcost = gapcost_list[gapcost]

                    test_scores = S[j] + one_mapinfo[i][3] - gapcost

                    filtered = False


                    break
                if(filtered == True):

                    #test_scores = S[j] - skipcost + one_mapinfo[i][3]
                    gapcost = abs(readgap - abs(refgap))
                    if(gapcost < 10000):

                        test_scores = S[j] - skipcost + one_mapinfo[i][3] - var_gap[gapcost]
                    else:
                        test_scores = S[j] - skipcost + one_mapinfo[i][3] - 10.

                if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j
                elif((test_scores == max_scores) and (pre_index != -9999999)):
                    pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                    now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                    if(now_refdistance < pre_refdistance):
                        max_scores = test_scores
                        pre_index = j






                en_loc = st_loc
        
        else:
            
            for j in S_arg[:testspace_en][::-1]:
            
            
                opcount += 1     

                if(S[j] < (max_scores - one_mapinfo[i][3])):

                    break


                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                if((readgap < 0)):
                    #print_log('error')
                    continue

                nocost = False
                filtered = True

                if(one_mapinfo[i][2] == 1):

                            #    y_2 - y_1 - a_1
                    refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])

                else:

                            #     y_1 - y_2 - a_2
                    refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])


                while(True):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):

                        break

                    if((one_mapinfo[i][2] == 1) and (refgap < 0)):

                        break

                    if(one_mapinfo[i][2] == -1):

                        if(refgap > 0):

                            break

                        else:

                            refgap = abs(refgap)

                    gapcost = abs(readgap - refgap)

                    #if((readgap  > maxgap) or (gapcost > maxdiff)):
                    if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                        break


                    gapcost = gapcost_list[gapcost]

                    test_scores = S[j] + one_mapinfo[i][3] - gapcost

                    filtered = False

                    break
                if(filtered == True):

                    #test_scores = S[j] - skipcost + one_mapinfo[i][3]
                    gapcost = abs(readgap - abs(refgap))
                    if(gapcost < 10000):

                        test_scores = S[j] - skipcost + one_mapinfo[i][3] - var_gap[gapcost]
                    else:
                        test_scores = S[j] - skipcost + one_mapinfo[i][3] - 10.

                if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j
                elif((test_scores == max_scores) and (pre_index != -9999999)):
                    pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                    now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                    if(now_refdistance < pre_refdistance):
                        max_scores = test_scores
                        pre_index = j

            S[i] = max_scores
            P[i] = pre_index

            if(max_scores > g_max_scores):

                g_max_scores = max_scores
                g_max_index = i
            
            
            
            
            
                


        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
        

    k = testspace_en
    while(k < n):

        loc_in_sorted_S = smallorequal2target_1d_point_target(S, S[k], k, S_arg, target_arr[k], target_arr) + 1

        S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
        S_arg[loc_in_sorted_S] = k

        k += 1     
        
    return g_max_index, S, P, S_arg
@njit
def hit2work(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage = 0.5, hastra = False, H = False):
    def getsecond(x):
        return x[1]
    def getfirst(x):
        return x[0]
    def getlength(x):
        return len(x)

    def get_readloc_set_bin(one_mappos, bin_size):
        return set([i[0] // bin_size for i in one_mappos])
    def get_refloc_set_bin(one_mappos, bin_size):
        return set([i[1] // bin_size for i in one_mappos])
    def get_overlapsize(readloc_set_a, readloc_set_b):
        return len(readloc_set_a&readloc_set_b)/min(len(readloc_set_a), len(readloc_set_b))
    #print_log(len(one_mapinfo))


    
    one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:,1])] 

    cluster_list = List()
    iloc = 0
    preitem = one_mapinfo[iloc]
    st_iloc = 0
    for iloc in range(len(one_mapinfo)):
        nowitem = one_mapinfo[iloc]
        if(((nowitem[1] - preitem[1]) > c_bias)):
            if((iloc - st_iloc) < 3):
                continue
            cluster_list.append(one_mapinfo[st_iloc: iloc])
            #print_log(one_mapinfo[st_iloc][1], one_mapinfo[iloc - 1][1])
            st_iloc = iloc
        preitem = nowitem


    if((iloc - st_iloc) > 3):
        cluster_list.append(one_mapinfo[st_iloc: iloc + 1])
    cluster_list.sort(key = get_length)

    cluster_list = cluster_list[::-1][:check_num]


    
    hit = False
    minichain_scores = 40
    path_list = List()
    scores_list = []
    scores_list.append(0.)
    scores_list.pop()
    path_list.append([(0, 0, 0, 0)])
    path_list.pop()
    

    max_scores = 0
    from_repeat = False
    
    fast_enable = False
    g_max_index = 0

    for one_mapinfo in cluster_list:
        repeat = False
        min_ref, max_ref = one_mapinfo[0][1], one_mapinfo[-1][1]
        testcontig = pos2contig(min_ref, contig2start)
        if(testcontig != pos2contig(max_ref, contig2start)):
            #print_log('hit2work: testcontig != pos2contig(max_ref, contig2start)')
            continue
        
        one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:,0])]
        if(fast_enable == False):
            g_max_index, S, P, S_arg = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_all(one_mapinfo, kmersize = kmersize, skipcost = skipcost[0], maxdiff = maxdiff[0], maxgap = maxgap)
        

        if(fast_enable == True or g_max_index == -1):
            fast_enable = True
            g_max_index, S, P, S_arg = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_fast_all(one_mapinfo, kmersize = kmersize, skipcost = skipcost[0], maxdiff = maxdiff[0], maxgap = maxgap)

        scores = S[g_max_index]
        #if(scores < minichain_scores):
            #if(hit == False):
                #return 0, 0., [(0, 0, 0, 0)]
            #continue
        
        ##################
        usedindex = set()

        path = [(0, 0, 0, 0)]
        path.pop()
        take_index = g_max_index
        usedindex.add(take_index)
        score = S[take_index]
        while(True):
            if((P[take_index] == -9999999)):
                path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
                break
            path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
            take_index = P[take_index]
            usedindex.add(take_index)
        if(score > 40):
            hit = True
            scores_list.append(score)
            path_list.append(path)
        if(scores > max_scores):
            if(repeat == True):
                from_repeat = True
            else:
                from_repeat = False
            max_scores = scores
        #print_log(score, len(path))

        for take_index in S_arg[::-1]:
            if(take_index in usedindex):
                continue
            path = [(0, 0, 0, 0)]
            path.pop()
            usedindex.add(take_index)
            score = S[take_index]
            while(True):
                if((P[take_index] == -9999999)):
                    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
                    break
                path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
                take_index = P[take_index]
                if(take_index in usedindex):
                    score = score - S[take_index]
                    break
                usedindex.add(take_index)
            if(score > 40):
                scores_list.append(score)
                path_list.append(path)
            if(scores > max_scores):
                if(repeat == True):
                    from_repeat = True
                else:
                    from_repeat = False
                max_scores = scores
        ##################
        
        

        

    #print_log(max_scores)
    #print_log()
    
    if(hit == True and max_scores > 150):

        order = np.argsort(np.array(scores_list))[::-1]
        #print_log('scores_list[order[0]], len(path_list[order[0]])', scores_list[order[0]], len(path_list[order[0]]))
        primary_rlocset_List = List()
        primary_scores_List = List()
        primary_index_List = List()
        
        iloc = order[0]
        primary_rlocset_List.append(get_readloc_set_bin(path_list[iloc], bin_size))
        primary_scores_List.append(List([scores_list[iloc]]))
        primary_index_List.append(iloc)
        #print_log(scores_list[order[0]], path_list[order[0]][-1][1], path_list[order[0]][0][1])
        for iloc in order[1:]:
            
            #print_log(scores_list[iloc], path_list[iloc][-1][1], path_list[iloc][0][1])
            readloc_set_b = get_readloc_set_bin(path_list[iloc], bin_size)
            maxoverlapsize = 0.
            for p_loc in range(len(primary_rlocset_List)):
                tmp_overlapsize = get_overlapsize(primary_rlocset_List[p_loc], readloc_set_b)
                if(tmp_overlapsize > maxoverlapsize):
                    maxoverlapsize = tmp_overlapsize
                    prefer_p_loc = p_loc
            if(maxoverlapsize < overlapprecentage):
                primary_rlocset_List.append(readloc_set_b)
                primary_scores_List.append(List([scores_list[iloc]]))
                primary_index_List.append(iloc)
            else:
                primary_scores_List[prefer_p_loc].append(scores_list[iloc])
                
        m = len(path_list[order[0]])    
        if(len(primary_scores_List[0]) < 2):
            f1 = primary_scores_List[0][0]
            f2 = 0
        else:
            f1 = primary_scores_List[0][0]
            f2 = primary_scores_List[0][1]
        mapq = min(int(40*(1-f2/f1)*min(1, m/10)*np.log(f1)), 60)

        base_iloc = primary_index_List[0]

        for iloc in primary_index_List[1:]:
            if(hastra == False):
                skip = True
                item = path_list[iloc][0]
                bias_1 = abs(item[1]-path_list[base_iloc][0][1])
                bias_2 = abs(item[1] - path_list[base_iloc][-1][1])
                if(min(bias_1, bias_2) < 200000):
                    if(pos2contig(item[1], contig2start) == pos2contig(path_list[base_iloc][0][1], contig2start)):
                        skip = False
                if(skip == True):
                    continue
            for item in path_list[iloc]:                            
                path_list[base_iloc].append(item)

        path_list[base_iloc].sort(key = getfirst)
        scores, path = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d(np.array(path_list[base_iloc]), kmersize = kmersize, skipcost = skipcost[1], maxdiff = maxdiff[1], maxgap = maxgap)
        #print_log(scores, path[-1][1], path[0][1]) 
        #print_log(path)
        #print_log()
        #print_log(scores, len(path))
        return mapq, scores, path


    else:
        return 0, 0., [(0, 0, 0, 0)]
def decode_hit(index_object, index2contig, testseq, testseq_len, contig2start, kmersize, skipcost = (50., 30.), maxdiff = (50, 30), maxgap = 200, check_num = 20, c_bias = 5000, bin_size = 100, overlapprecentage = 0.5, hastra = True, H = False, mid_occ = -1):            
    need_reverse, one_mapinfo = get_reversed_chain_numpy_rough(np.array(index_object.map(testseq, check_num = check_num,mid_occ = mid_occ)), testseq_len)
    mapq, scores, path =  hit2work(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage, hastra, H = H)
    if(need_reverse == True):
        return mapq, -scores, path
    else:
        return mapq, scores, path
def get_readmap_DP_test(readid, testseq, contig2start, contig2seq, index_object, index2contig, kmersize = 15, local_kmersize = 9, local_maxdiff = 30, refine = True, local_skipcost = 50., golbal_skipcost = (30., 30.),  golbal_maxdiff = (50, 30), check_num = 20, bin_size = 100, hastra = False, debug = False, H = False, mid_occ = -1, redo_ratio = 5):
    local_skipcost = 30.    
    mapq = 60
    setting_kmersize = kmersize
    setting_maxdiff = golbal_maxdiff[1]
    if(H == False):
        local_skipcost += local_kmersize
        golbal_skipcost = (golbal_skipcost[0] + kmersize, golbal_skipcost[1] + kmersize)
    else:
        local_skipcost = 30.
        golbal_skipcost = (30., 30.)
        hastra = True
        check_num = 100
    

    rc_testseq = str(Seq(testseq).reverse_complement())
    testseq_len = len(testseq)
    st = time.time()                            
    mapq, scores, raw_alignment_list = decode_hit(index_object, index2contig, testseq, testseq_len, contig2start, kmersize, contig2seq, skipcost = golbal_skipcost, maxdiff = golbal_maxdiff, maxgap = 200, check_num = check_num, c_bias = 5000, bin_size = bin_size, overlapprecentage = 0.5, hastra = hastra, H = H, mid_occ = mid_occ)
    f_redo_ratio = min(abs(testseq_len/(scores+1e-7)), 20)
    if(f_redo_ratio > redo_ratio or mapq == 0):
        mapq, scores, raw_alignment_list = decode_hit(index_object, index2contig, testseq, testseq_len, contig2start, kmersize, contig2seq, skipcost = golbal_skipcost, maxdiff = golbal_maxdiff, maxgap = 200, check_num = 100, c_bias = 5000, bin_size = bin_size, overlapprecentage = 0.5, hastra = hastra, H = H, mid_occ = mid_occ)


    if(debug == True):
        plot_result = np.array(raw_alignment_list)
        plt.scatter(plot_result[:, 0], plot_result[:, 1])
        plt.show()
    if(scores == 0.):
        return [], ([], []), 0, 0
    
    need_reverse = False
    if((scores < 0.)):
        need_reverse = True

    if(debug == True): print_log('Normal time', time.time() - st, scores)
    if(refine == True):
        setting_maxdiff = local_maxdiff
        setting_kmersize = local_kmersize
        st = time.time()

        #need_reverse, raw_alignment_array = get_reversed_chain_numpy_rough(np.array(raw_alignment_list), testseq_len)
        raw_alignment_array = np.array(raw_alignment_list)
        if(need_reverse == False):
            scores, raw_alignment_list = get_localmap_multi_all_forDP_inv_guide(raw_alignment_array, testseq, rc_testseq, contig2start, contig2seq, kmersize = setting_kmersize, skipcost = local_skipcost, maxdiff = setting_maxdiff, maxgap = 100, shift = 1)

        else:
            scores, raw_alignment_list = get_localmap_multi_all_forDP_inv_guide(raw_alignment_array, rc_testseq, testseq, contig2start, contig2seq, kmersize = setting_kmersize, skipcost = local_skipcost, maxdiff = setting_maxdiff, maxgap = 100, shift = 1)
            testseq, rc_testseq = rc_testseq, testseq
            #raw_alignment_list = get_reversed_chain_numpy(np.array(raw_alignment_list), testseq_len)

        #raw_alignment_list = check_func_clean(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, 10, debug = debug)
        #for item in raw_alignment_list:
            #print_log(item)
        if(len(raw_alignment_list) <= 1):
            return [], ([], []), 0, 0
        if(debug == True): print_log('Refine time', time.time() - st)


        
        

    st = time.time()
    try:
        filtered = False
        #alignment_list, onemapinfolist, TRA_signal = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, debug = debug, H = H)
        alignment_list, onemapinfolist, TRA_signal, filtered = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, debug = debug, H = H, nofilter = False)
        if(len(onemapinfolist) == 0):
            return [], ([], []), 0, 0
    except:
        print_log(testseq)
    
    

    if(filtered == True and pairedindel(List([line[-1] for line in onemapinfolist]), indelsize = 30) == True):
        alignment_list, onemapinfolist, TRA_signal, filtered = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, debug = debug, H = H, nofilter = True)
        
    if(debug == True): print_log('extend_func', time.time() - st)
    if(len(onemapinfolist) > 1 and refine == 'auto'):
        setting_maxdiff = local_maxdiff
        setting_kmersize = local_kmersize
        st = time.time()                                                                         
        scores, raw_alignment_list = get_localmap_multi_all_forDP_inv_guide(np.array(raw_alignment_list), testseq, rc_testseq, contig2start, contig2seq, kmersize = setting_kmersize, skipcost = local_skipcost, maxdiff = setting_maxdiff, maxgap = 100, shift = 1)
        if(debug == True): print_log('Refine time', time.time() - st, len(raw_alignment_list))

        alignment_list, onemapinfolist, TRA_signal = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, need_reverse, setting_maxdiff, debug = debug)

    return onemapinfolist, (alignment_list, raw_alignment_list), TRA_signal, f_redo_ratio

@njit
def hit2work_1(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage = 0.5, hastra = False, H = False):
    def getsecond(x):
        return x[1]
    def getfirst(x):
        return x[0]
    def getlength(x):
        return len(x)

    def get_readloc_set_bin(one_mappos, bin_size):
        return set([i[0] // bin_size for i in one_mappos])
    def get_refloc_set_bin(one_mappos, bin_size):
        return set([i[1] // bin_size for i in one_mappos])
    def get_overlapsize(readloc_set_a, readloc_set_b):
        return len(readloc_set_a&readloc_set_b)/min(len(readloc_set_a), len(readloc_set_b))
    #print_log(len(one_mapinfo))


    
    one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:,1])] 

    cluster_list = List()
    iloc = 0
    preitem = one_mapinfo[iloc]
    st_iloc = 0
    for iloc in range(len(one_mapinfo)):
        nowitem = one_mapinfo[iloc]
        if(((nowitem[1] - preitem[1]) > c_bias)):
            if((iloc - st_iloc) < 3):
                continue
            cluster_list.append(one_mapinfo[st_iloc: iloc])
            #print_log(one_mapinfo[st_iloc][1], one_mapinfo[iloc - 1][1])
            st_iloc = iloc
        preitem = nowitem


    if((iloc - st_iloc) > 3):
        cluster_list.append(one_mapinfo[st_iloc: iloc + 1])
    cluster_list.sort(key = get_length)

    cluster_list = cluster_list[::-1][:check_num]


    
    hit = False
    minichain_scores = 40
    path_list = List()
    scores_list = []
    scores_list.append(0.)
    scores_list.pop()
    path_list.append([(0, 0, 0, 0)])
    path_list.pop()
    

    max_scores = 0
    from_repeat = False
    
    fast_enable = False
    g_max_index = 0

    for one_mapinfo in cluster_list:
        repeat = False
        min_ref, max_ref = one_mapinfo[0][1], one_mapinfo[-1][1]
        testcontig = pos2contig(min_ref, contig2start)
        if(testcontig != pos2contig(max_ref, contig2start)):
            #print_log('hit2work: testcontig != pos2contig(max_ref, contig2start)')
            continue
        
        one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:,0])]
        if(fast_enable == False):
            g_max_index, S, P, S_arg = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_all(one_mapinfo, kmersize = kmersize, skipcost = skipcost[0], maxdiff = maxdiff[0], maxgap = maxgap)
        

        if(fast_enable == True or g_max_index == -1):
            fast_enable = True
            g_max_index, S, P, S_arg = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_fast_all(one_mapinfo, kmersize = kmersize, skipcost = skipcost[0], maxdiff = maxdiff[0], maxgap = maxgap)

        scores = S[g_max_index]
        #if(scores < minichain_scores):
            #if(hit == False):
                #return 0, 0., [(0, 0, 0, 0)]
            #continue
        
        ##################
        usedindex = set()

        path = [(0, 0, 0, 0)]
        path.pop()
        take_index = g_max_index
        usedindex.add(take_index)
        score = S[take_index]
        while(True):
            if((P[take_index] == -9999999)):
                path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
                break
            path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
            take_index = P[take_index]
            usedindex.add(take_index)
        if(score > 40):
            hit = True
            scores_list.append(score)
            path_list.append(path)
        if(scores > max_scores):
            if(repeat == True):
                from_repeat = True
            else:
                from_repeat = False
            max_scores = scores
        #print_log(score, len(path))

        for take_index in S_arg[::-1]:
            if(take_index in usedindex):
                continue
            path = [(0, 0, 0, 0)]
            path.pop()
            usedindex.add(take_index)
            score = S[take_index]
            while(True):
                if((P[take_index] == -9999999)):
                    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
                    break
                path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
                take_index = P[take_index]
                if(take_index in usedindex):
                    score = score - S[take_index]
                    break
                usedindex.add(take_index)
            if(score > 40):
                scores_list.append(score)
                path_list.append(path)
            if(scores > max_scores):
                if(repeat == True):
                    from_repeat = True
                else:
                    from_repeat = False
                max_scores = scores
        ##################
        
        

        

    #print_log(max_scores)
    #print_log()
    
    if(hit == True and max_scores > 50):

        order = np.argsort(np.array(scores_list))[::-1]
        #print_log('scores_list[order[0]], len(path_list[order[0]])', scores_list[order[0]], len(path_list[order[0]]))
        primary_rlocset_List = List()
        
        primary_scores_List = List()

        
        primary_index_List = List()

        
        all_index_List = List()

        
        iloc = order[0]
        primary_rlocset_List.append(get_readloc_set_bin(path_list[iloc], bin_size))
        primary_scores_List.append(List([scores_list[iloc]]))
        primary_index_List.append(iloc)
        all_index_List.append(List([iloc]))
        #print_log(scores_list[order[0]], path_list[order[0]][-1][1], path_list[order[0]][0][1])
        for iloc in order[1:]:
            
            #print_log(scores_list[iloc], path_list[iloc][-1][1], path_list[iloc][0][1])
            readloc_set_b = get_readloc_set_bin(path_list[iloc], bin_size)
            maxoverlapsize = 0.
            for p_loc in range(len(primary_rlocset_List)):
                tmp_overlapsize = get_overlapsize(primary_rlocset_List[p_loc], readloc_set_b)
                if(tmp_overlapsize > maxoverlapsize):
                    maxoverlapsize = tmp_overlapsize
                    prefer_p_loc = p_loc
            if(maxoverlapsize < overlapprecentage):
                primary_rlocset_List.append(readloc_set_b)
                primary_scores_List.append(List([scores_list[iloc]]))
                all_index_List.append(List([iloc]))
                primary_index_List.append(iloc)
            else:
                primary_scores_List[prefer_p_loc].append(scores_list[iloc])
                all_index_List[prefer_p_loc].append(iloc)
        
        
        m = len(path_list[order[0]])    
        if(len(primary_scores_List[0]) < 2):
            f1 = primary_scores_List[0][0]
            f2 = 0
        else:
            f1 = primary_scores_List[0][0]
            f2 = primary_scores_List[0][1]
        mapq = min(int(40*(1-f2/f1)*min(1, m/10)*np.log(f1)), 60)
        
        return path_list, primary_index_List, primary_scores_List, all_index_List, mapq



    else:
        primary_scores_List = List()
        primary_scores_List.append(List([0.]))
        primary_scores_List.pop(0)
        
        primary_index_List = List()
        primary_index_List.append(0)
        primary_index_List.pop(0)
        
        all_index_List = List()
        all_index_List.append(List([0]))
        all_index_List.pop(0)
        return path_list, primary_index_List, primary_scores_List, all_index_List, 0

@njit
def hit2work_2(path_list, primary_index_List, primary_scores_List, all_index_List, kmersize, skipcost, maxdiff, maxgap, hastra, contig2start, base_iloc):


    for iloc in primary_index_List:
        if(iloc == base_iloc):
            continue
        if(hastra == False):
            skip = True
            item = path_list[iloc][0]
            bias_1 = abs(item[1]-path_list[base_iloc][0][1])
            bias_2 = abs(item[1] - path_list[base_iloc][-1][1])
            if(min(bias_1, bias_2) < 200000):
                if(pos2contig(item[1], contig2start) == pos2contig(path_list[base_iloc][0][1], contig2start)):
                    skip = False
            if(skip == True):
                continue
        for item in path_list[iloc]:                            
            path_list[base_iloc].append(item)

    path_list[base_iloc].sort(key = getfirst)
    return get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d(np.array(path_list[base_iloc]), kmersize = kmersize, skipcost = skipcost[1], maxdiff = maxdiff[1], maxgap = maxgap)
    
@njit
def return_main_alignment_size(contig2start, raw_alignment_list):

    preitem = raw_alignment_list[0]
    
    size = 0
     
    st_item = preitem
    for nowitem in raw_alignment_list[1:]:
        if(preitem[2] == nowitem[2]):
            readgap = nowitem[0] - preitem[0] - preitem[3]


            if(preitem[2] == 1):
                refgap = nowitem[1] - preitem[1] - preitem[3]
            else:
                refgap = preitem[1]  - nowitem[1] - nowitem[3]

            if((abs(readgap - refgap) <= 30) and (refgap >= 0)):
                if(pos2contig(preitem[1], contig2start) == pos2contig(nowitem[1], contig2start)):
                    preitem = nowitem
                    continue

        if(preitem[0] - st_item[0] > size):
            size = preitem[0] - st_item[0]
            return_pack = (st_item, preitem)
        preitem = nowitem
        
        st_item = preitem
    if(preitem[0] - st_item[0] > size):
        size = preitem[0] - st_item[0]
        return_pack = (st_item, preitem)

    return return_pack    
def decode_hit(index_object, index2contig, testseq, testseq_len, contig2start, kmersize, contig2seq, skipcost = (50., 30.), maxdiff = (50, 30), maxgap = 200, check_num = 20, c_bias = 5000, bin_size = 100, overlapprecentage = 0.5, hastra = True, H = False, mid_occ = -1):            
    need_reverse, one_mapinfo = get_reversed_chain_numpy_rough(np.array(index_object.map(testseq, check_num = check_num,mid_occ = mid_occ)), testseq_len)
    maxgap = 1000
    path_list, primary_index_List, primary_scores_List, all_index_List, mapq =  hit2work_1(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage, hastra, H = H)
    if(len(primary_index_List) != 0):
        base_iloc = primary_index_List[0]
        
        if(mapq == 0):
            #print_log('len(primary_scores_List)', len(primary_scores_List))
            base_score = primary_scores_List[0][0]
            rc_testseq = str(Seq(testseq).reverse_complement())
            min_diff = 10
            for tmpiloc in range(len(primary_scores_List[0])):
                if(primary_scores_List[0][tmpiloc] / base_score < 0.999):
                    break
                #print_log(all_index_List[0][tmpiloc])
                #print_log(primary_scores_List[0][tmpiloc])
                #print_log(path_list[all_index_List[0][tmpiloc]][-1], ',', path_list[all_index_List[0][tmpiloc]][0])
                #preitem, nowitem = path_list[all_index_List[0][tmpiloc]][-1], path_list[all_index_List[0][tmpiloc]][0]
                preitem, nowitem = return_main_alignment_size(contig2start, np.array(path_list[all_index_List[0][tmpiloc]][::-1]))
                #print_log(preitem, nowitem)
                if(preitem[2] != nowitem[2]):
                    continue
                if(need_reverse == False):
                    target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start)
                else:
                    target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, rc_testseq, testseq, testseq_len, kmersize, contig2seq, contig2start)
                if(min(len(target), len(query)) == 0):
                    print_log(preitem, nowitem)
                    print_log(testseq)
                diff = edlib.align(query = query, target = target, task = 'distance')['editDistance']/min(len(target), len(query))
                #print_log(diff)
                if(diff <= min_diff):
                    min_diff = diff
                    #print_log('change', base_iloc, '->', all_index_List[0][tmpiloc])
                    base_iloc = all_index_List[0][tmpiloc]
                    primary_index_List[0] = base_iloc
 
                    
                #print_log()
        
        scores, path = hit2work_2(path_list, primary_index_List, primary_scores_List, all_index_List, kmersize, skipcost, maxdiff, maxgap, hastra, contig2start, base_iloc)
    else:
        scores, path = 0., [(0, 0, 0, 0)]
    if(need_reverse == True):
        return mapq, -scores, path
    else:
        return mapq, scores, path
    
#2024312
@njit
def get_localmap_multi_all_forDP_inv_guide(raw_alignment_array, testseq, rc_testseq, contig2start, contig2seq, kmersize, skipcost, maxdiff, maxgap, shift = 1):#
    




    def seq2hashtable_multi_test(onelookuptable, seq, start, kmersize):
        skiphash = hash('N'*kmersize)
        for iloc in range(0, len(seq) - kmersize + 1, 1):
            hashedkmer = hash(seq[iloc:iloc+kmersize])
            if(skiphash == hashedkmer):
                continue
            if(hashedkmer in onelookuptable):

                onelookuptable[hashedkmer].append(start + iloc)
            else:
                onelookuptable[hashedkmer] = List([start + iloc])

    readgap = 0
    pre = raw_alignment_array[0]
    for now in raw_alignment_array[1:]:
        if(abs(now[0] - pre[0]) > readgap):
            readgap = abs(now[0] - pre[0])
        pre = now
    readgap += 1000  
    readgap = max(readgap, 5000)
    raw_alignment_array = raw_alignment_array[np.argsort(raw_alignment_array[:, 1])]
    startandend = List([(raw_alignment_array[0][1], raw_alignment_array[0][1])])
    for item in raw_alignment_array[1:]:
        if(((item[1] - startandend[-1][1]) < readgap)):
            startandend[-1] = (startandend[-1][0], item[1])
        else:
            if(startandend[-1][0] == startandend[-1][1]):
                startandend.pop(-1)
            startandend.append((item[1], item[1]))
    if(startandend[-1][0] == startandend[-1][1]):
        startandend.pop(-1)
    local_lookuptable = Dict()
    local_lookuptable[0] = List([0])

    for item in startandend:
        min_ref, max_ref = item[0], item[1]
        testcontig = pos2contig(min_ref, contig2start)
        if(testcontig != pos2contig(max_ref, contig2start)):
            continue
        lookfurther = min(2000, min_ref-contig2start[testcontig])
        min_ref -= lookfurther
        max_ref += 2000
        refseq = contig2seq[testcontig][min_ref-contig2start[testcontig]: max_ref-contig2start[testcontig]]
        seq2hashtable_multi_test(local_lookuptable, refseq, min_ref, kmersize)
    local_lookuptable.pop(0)
    raw_alignment_array = raw_alignment_array[np.argsort(raw_alignment_array[:, 0])]
    readstart = max(0, raw_alignment_array[0][0]-500)
    readend = min(len(testseq)-kmersize+1, raw_alignment_array[-1][0]+500)
    
    
    iloc = readstart
    iloc -= shift
    
    maxop = 0
    pointdict = Dict()
    
    while(True):
   
        iloc += shift
        if(iloc >= readend):
            break

        if(hash(testseq[iloc: iloc + kmersize]) == hash(rc_testseq[-(iloc + kmersize): -iloc])):
            continue
        biasvalue, biasvalue_1, closest_index, closest_index_1 = findClosest(raw_alignment_array, target = iloc)


        interval = min(biasvalue + biasvalue_1 + 500, 2000)
        upperrefloc =  (raw_alignment_array[closest_index][1] + interval, raw_alignment_array[closest_index_1][1] + interval)
        lowerrefloc =  (raw_alignment_array[closest_index][1] - interval, raw_alignment_array[closest_index_1][1] - interval)
        readgap = abs(iloc - raw_alignment_array[closest_index][0])
   
        hashedkmer = hash(testseq[iloc: iloc + kmersize])  
        if(hashedkmer in local_lookuptable):
            for refloc in local_lookuptable[hashedkmer]:
                refgap = abs(refloc - raw_alignment_array[closest_index][1])
                diff = abs(readgap - refgap)
                if((diff < 500) or (upperrefloc[0] >= refloc and lowerrefloc[0] <= refloc) or (upperrefloc[1] >= refloc and lowerrefloc[1] <= refloc)):
                    item = (iloc, refloc, 1)
                    point = item[1] - item[0]

                    if(point in pointdict):

                        if((pointdict[point][-1][0] + pointdict[point][-1][3]) >= item[0]):
                            bouns = item[0] - (pointdict[point][-1][0] + pointdict[point][-1][3]) + kmersize
                            if(bouns > 0):
                                if((pointdict[point][-1][3] + bouns < 20)):
                                    pointdict[point][-1] = (pointdict[point][-1][0], pointdict[point][-1][1], 1, pointdict[point][-1][3] + bouns)
                                else:
                                    pointdict[point].append((pointdict[point][-1][0] + pointdict[point][-1][3], pointdict[point][-1][1] + pointdict[point][-1][3], 1, bouns))
                          


                        else:    
                            pointdict[point].append((item[0], item[1], item[2], kmersize))
                    else:
                        pointdict[point] = List([(item[0], item[1], item[2], kmersize)])


        hashedkmer = hash(rc_testseq[-(iloc + kmersize): -iloc]) 
        if(hashedkmer in local_lookuptable):

            for refloc in local_lookuptable[hashedkmer]:
                refgap = abs(refloc - raw_alignment_array[closest_index][1])
                diff = abs(readgap - refgap)
                if((diff < 500) or (upperrefloc[0] >= refloc and lowerrefloc[0] <= refloc) or (upperrefloc[1] >= refloc and lowerrefloc[1] <= refloc)):
                    item = (iloc, refloc, -1)
                    point = -(item[1] + item[0])
                    if(point in pointdict):

                        if((pointdict[point][-1][0] + pointdict[point][-1][3]) >= item[0]):
                            bouns = item[0] - (pointdict[point][-1][0] + pointdict[point][-1][3]) + kmersize
                            if(bouns > 0):
                                if((pointdict[point][-1][3] + bouns < 20)):
                                    pointdict[point][-1] = (pointdict[point][-1][0], item[1], -1, pointdict[point][-1][3] + bouns)
                                else:

                                    pointdict[point].append((pointdict[point][-1][0] + pointdict[point][-1][3], item[1], -1, bouns))
                                        #print_log(pointdict[point][-1])


                        else:    
                            pointdict[point].append((item[0], item[1], item[2], kmersize))
                    else:
                        pointdict[point] = List([(item[0], item[1], item[2], kmersize)])
    
    one_mapinfo = [(-1, -1, -1, -1)]
    one_mapinfo.pop(0)
    for key in pointdict:
        for item in pointdict[key]:
            one_mapinfo.append(item)
    one_mapinfo = np.array(one_mapinfo)
    one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:, 0])]
  


    return get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo, kmersize = kmersize, skipcost = skipcost, maxdiff = maxdiff, maxgap = maxgap)

#20240313
@njit
def findClosest_point_inguide(arr, target):
    def getClosest(val1, val2, target):
 
        if (target - val1 >= val2 - target):
            return val2
        else:
            return val1
    n = len(arr)
    # Corner cases
    if (target <= arr[0]):
        return arr[0]
    if (target >= arr[n - 1]):
        return arr[n - 1]
 
    # Doing binary search
    i = 0; j = n; mid = 0
    while (i < j): 
        mid = (i + j) // 2
 
        if (arr[mid] == target):
            return arr[mid]
 
        # If target is less than array 
        # element, then search in left
        if (target < arr[mid]) :
 
            # If target is greater than previous
            # to mid, return closest of two
            if (mid > 0 and target > arr[mid - 1]):
                return getClosest(arr[mid - 1], arr[mid], target)
 
            # Repeat for left half 
            j = mid
         
        # If target is greater than mid
        else :
            if (mid < n - 1 and target < arr[mid + 1]):
                return getClosest(arr[mid], arr[mid + 1], target)
                 
            # update i
            i = mid + 1
         
    # Only single element left after search
    return arr[mid]

#20240401

@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20
    max_kmersize = 0
    
    #gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    #gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
    
    #var_gap = np.arange(1, 10001)
    #var_gap = np.minimum(var_gap/100, np.log(var_gap))


    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)
    T = np.empty(n, np.int64)


    

    

    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    readend_arr = one_mapinfo[:, 0] + one_mapinfo[:, 3]

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_st = 0
    testspace_en = 1
    
    ms_arr = np.empty(one_mapinfo[-1][0] + 1)
    
    test_S = np.empty(0, np.float64)
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0
    max_kmersize = one_mapinfo[i][3]
    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    T[0] = 0
    

    
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        
        T[i] = 0
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        
        max_kmersize = max(one_mapinfo[i][3], max_kmersize)
        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            #if((opcount/one_mapinfo[i][0] > 2000) and one_mapinfo[i][0] > 1000):
                #return 0., [(0, 0, 0, 0)]
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i
            
            #gap_arr[prereadloc] = 200+ 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
               
                      


        

        for j in S_arg[:testspace_en][::-1]:
            
            
            #opcount += 1     
            
            if(S[j] < (max_scores - one_mapinfo[i][3])):
                
                break


            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            if((readgap < 0)):
                #print_log('error')
                continue
                
            nocost = False
            filtered = True
            
            if(one_mapinfo[i][2] == 1):
                
                        #    y_2 - y_1 - a_1
                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                
            else:
                
                        #     y_1 - y_2 - a_2
                refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                
            
            while(True):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    
                    break
                    
                if((one_mapinfo[i][2] == 1) and (refgap < 0)):
                    
                    break
                    
                if(one_mapinfo[i][2] == -1):
                    
                    if(refgap > 0):
                        
                        break
                        
                    else:
                        
                        refgap = abs(refgap)

                gapcost = abs(readgap - refgap)
                
                if((readgap  > maxgap) or (gapcost > maxdiff)):
                #if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                    break
                    

                    
                gapcost = gapcost_list[gapcost]
                
                test_scores = S[j] + one_mapinfo[i][3] - gapcost
                
                filtered = False
                
                break
            if(filtered == True):

                #test_scores = S[j] - skipcost + one_mapinfo[i][3]
                gapcost = abs(readgap - abs(refgap))#mark 1

                #test_scores = S[j] - skipcost + one_mapinfo[i][3] - min(gapcost/100, 2 * np.log(gapcost))
                test_scores = S[j] - skipcost + one_mapinfo[i][3] -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))

                '''gapcost = abs(readgap - abs(refgap))
                if(gapcost < 10000):

                    test_scores = S[j] - skipcost + one_mapinfo[i][3] - var_gap[gapcost]
                else:
                    test_scores = S[j] - skipcost + one_mapinfo[i][3] - 10.'''

            '''if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
            elif((test_scores == max_scores) and (pre_index != -9999999)):
                pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                if(now_refdistance < pre_refdistance):
                    max_scores = test_scores
                    pre_index = j'''
            if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    T[i] = abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1]))
                else:
                    T[i] = T[j]
                    
            elif((test_scores == max_scores) and (pre_index != -9999999)):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    if(abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1])) < T[i]):
                        T[i] = abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1]))
                        pre_index = j
                else:
                    if(T[j] < T[i]):
                        T[i] = T[j]
                        pre_index = j

        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
        elif(max_scores >= g_max_scores):
            if(T[g_max_index] > T[i]):
                g_max_index = i
      

        
            
        
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    while(True):
        if((P[take_index] == -9999999)):
            break
        take_index = P[take_index]
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

    return g_max_scores, path
@njit
def rebuild_chain_break(contig2start, raw_alignment_list, large_cost, small_alignment = 50, small_dup = -100):
    #rebuild chain break
    #step 1
    #remove anchor cause large cost and small alignment
    preitem = raw_alignment_list[0]
    alignment_list = List([List([preitem])])
    for nowitem in raw_alignment_list[1:]:
        if(preitem[2] == nowitem[2]):
            readgap = nowitem[0] - preitem[0] - preitem[3]


            if(preitem[2] == 1):
                refgap = nowitem[1] - preitem[1] - preitem[3]
            else:
                refgap = preitem[1]  - nowitem[1] - nowitem[3]

            if((refgap < 0) and (refgap > small_dup)):
                continue
            if((abs(readgap - refgap) <= large_cost) and (refgap >= 0) and (readgap < 100)):
                if(pos2contig(preitem[1], contig2start) == pos2contig(nowitem[1], contig2start)):
                    alignment_list[-1].append(nowitem)
                    preitem = nowitem
                    continue

        if(len(alignment_list[-1]) == 1):
            alignment_list.pop(-1)
        if(len(alignment_list) > 0):
            if((alignment_list[-1][-1][0] + alignment_list[-1][-1][3] - alignment_list[-1][0][0]) < small_alignment):
                alignment_list.pop(-1)
  
        alignment_list.append(List([nowitem]))
        preitem = nowitem
    if(len(alignment_list[-1]) == 1):
        alignment_list.pop(-1)
    if((alignment_list[-1][-1][0] + alignment_list[-1][-1][3] - alignment_list[-1][0][0]) < small_alignment):
        alignment_list.pop(-1)


    return alignment_list
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_all(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20
    max_kmersize = 0
    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
    
    var_gap = np.arange(1, 10001)#mark 1
    var_gap = np.minimum(var_gap/100, np.log(var_gap))#mark 1


    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)
    T = np.empty(n, np.int64)


    

    

    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    readend_arr = one_mapinfo[:, 0] + one_mapinfo[:, 3]

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_st = 0
    testspace_en = 1
    
    ms_arr = np.empty(one_mapinfo[-1][0] + 1)
    
    test_S = np.empty(0, np.float64)
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0
    max_kmersize = one_mapinfo[i][3]
    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    
    T[0] = 0
    

    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        
        T[i] = 0
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        
        max_kmersize = max(one_mapinfo[i][3], max_kmersize)
        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            if((opcount/one_mapinfo[i][0] > 2000) and one_mapinfo[i][0] > 1000):#mark 1
                return -1, S, P, S_arg#mark 1
            
            
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i
            
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
               
                      


        
        for j in S_arg[:testspace_en][::-1]:
            
            
            opcount += 1     #mark 1
            
            if(S[j] < (max_scores - one_mapinfo[i][3])):
                
                break

            #print_log(one_mapinfo[j])
            #print_log(one_mapinfo[i])
            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])
            #print_log(readgap)
            goodsignal = False
            if((readgap < 0)):
                
                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):

                                #    y_2 - y_1 - a_1
                        refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])


                    else:
                         #1645 17384    -1    15]
                        #[ 1647 17382    -1    15

                                #     y_1 - y_2 - a_2
                        refgap = one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3]
                    if(refgap == readgap):
                        goodsignal = True
                if(goodsignal == False):    
                #print_log('error')
                    continue
            #print_log(goodsignal)  
            #print_log()
            nocost = False
            filtered = True
            
            if(one_mapinfo[i][2] == 1):
                
                        #    y_2 - y_1 - a_1
                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                
            else:
                
                        #     y_1 - y_2 - a_2
                refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                
            if(goodsignal == False):
                while(True):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):

                        break

                    if((one_mapinfo[i][2] == 1) and (refgap < 0)):

                        break

                    if(one_mapinfo[i][2] == -1):

                        if(refgap > 0):

                            break

                        else:

                            refgap = abs(refgap)

                    gapcost = abs(readgap - refgap)

                    #if((readgap  > maxgap) or (gapcost > maxdiff)):
                    if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                        break


                    gapcost = gapcost_list[gapcost]

                    test_scores = S[j] + one_mapinfo[i][3] - gapcost
                    

                    filtered = False

                    break
            else:
                test_scores = S[j] + one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                filtered = False

            if(filtered == True):
                #test_scores = S[j] - skipcost + one_mapinfo[i][3]
                gapcost = abs(readgap - abs(refgap))#mark 1

                #test_scores = S[j] - skipcost + one_mapinfo[i][3] - min(gapcost/100, 2 * np.log(gapcost))
                test_scores = S[j] - skipcost + one_mapinfo[i][3] -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))

                '''if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    #test_scores = S[j] - skipcost + one_mapinfo[i][3]
                    gapcost = abs(readgap - abs(refgap))#mark 1
                    if(gapcost < 10000):#mark 1

                        test_scores = S[j] - skipcost + one_mapinfo[i][3] - var_gap[gapcost]#mark 1
                    else:
                        test_scores = S[j] - skipcost + one_mapinfo[i][3] - 10.#mark 1
                else:
                    test_scores = S[j] - skipcost + one_mapinfo[i][3]'''

            '''if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
            elif((test_scores == max_scores) and (pre_index != -9999999)):
                pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                if(now_refdistance < pre_refdistance):
                    max_scores = test_scores
                    pre_index = j'''
            if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    T[i] = abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1]))
                else:
                    T[i] = T[j]
                    
            elif((test_scores == max_scores) and (pre_index != -9999999)):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    if(abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1])) < T[i]):
                        T[i] = abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1]))
                        pre_index = j
                else:
                    if(T[j] < T[i]):
                        T[i] = T[j]
                        pre_index = j

        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
        elif(max_scores >= g_max_scores):
            if(T[g_max_index] > T[i]):
                g_max_index = i


        
    k = testspace_en
    while(k < n):

        loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

        S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
        S_arg[loc_in_sorted_S] = k

        k += 1

    testspace_en = k  
    
    

        
    
    return g_max_index, S, P, S_arg
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_fast_all(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20
    max_kmersize = 0
    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))//1
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))//1
            
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = gapcost
        else:
            gapcost_list[gapcost] = gapcost
    
    var_gap = np.arange(1, 10001)#mark 1
    var_gap = np.minimum(var_gap//100, np.log(var_gap)//1)#mark 1

    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)
    T = np.empty(n, np.int64)

    

    

    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    readend_arr = one_mapinfo[:, 0] + one_mapinfo[:, 3]

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    readlength = one_mapinfo[-1][0] + 1000
    target_arr = np.zeros(n, dtype = np.float64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
        if(one_mapinfo[i][2] == 1):
            target_arr[i] = one_mapinfo[i][1] - one_mapinfo[i][0] + readlength 
        else:
            target_arr[i] = -(one_mapinfo[i][1] + one_mapinfo[i][0] + readlength) 
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_st = 0
    testspace_en = 1
    
    ms_arr = np.empty(one_mapinfo[-1][0] + 1)
    
    test_S = np.empty(0, np.float64)
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0
    max_kmersize = one_mapinfo[i][3]
    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    
    cache_dict = dict()
    cache_dict[0] = 0
    cache_dict.pop(0)
    
    fast = False
    
    T[0] = 0
    
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        T[i] = 0
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        
        max_kmersize = max(one_mapinfo[i][3], max_kmersize)
        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            #if((opcount/one_mapinfo[i][0] > 2000) and one_mapinfo[i][0] > 1000):
                #return 0., [(0, 0, 0, 0)]
                
            k = testspace_en
            while(one_mapinfo[k][0] + kmersize <= one_mapinfo[i][0]):

                loc_in_sorted_S = smallorequal2target_1d_point_target(S, S[k], k, S_arg, target_arr[k], target_arr) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1

            testspace_en = k
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
            
            
            fast = False
            if(testspace_en > 0):
                temp_endpos = (testspace_en - max(testspace_en - 1000, 0))

                scorerange = S[S_arg[testspace_en - 1]] - S[S_arg[max(testspace_en - 1000, 0)]] + 1
                
                if(temp_endpos / scorerange > 20):
                    fast = True
                    cache_dict = dict()
                    cache_dict[0] = 0
                    cache_dict.pop(0)
            
            
               
                      
        if(fast == True):

            #print_log(i, prereadloc)

            st_loc = testspace_en
            en_loc = testspace_en

            while(en_loc > 0):
                if(st_loc not in cache_dict):
                    pre_st_loc = st_loc
                    st_loc = smallorequal2target_1d_point(S, S[S_arg[st_loc-1]]-1e-7, en_loc, S_arg) + 1
                    cache_dict[pre_st_loc] = st_loc
                else:
                    st_loc = cache_dict[st_loc]

                if(S[S_arg[st_loc]] < (max_scores - one_mapinfo[i][3])):

                    break

                S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg, one_mapinfo[i][0], one_mapinfo[:,0])




                j = S_arg[S_arg_j]

                opcount += 1     



                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])



                nocost = False
                filtered = True

                if(one_mapinfo[i][2] == 1):

                            #    y_2 - y_1 - a_1
                    refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])

                else:

                            #     y_1 - y_2 - a_2
                    refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])


                while(True):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):

                        break

                    if((one_mapinfo[i][2] == 1) and (refgap < 0)):

                        break

                    if(one_mapinfo[i][2] == -1):

                        if(refgap > 0):

                            break

                        else:

                            refgap = abs(refgap)

                    gapcost = abs(readgap - refgap)

                    if(gapcost > maxdiff):

                        break

                    #if(readgap  > maxgap):
                    if(readgap  > gap_arr[one_mapinfo[j][0]]):#mark 1                   
                        break



                    gapcost = gapcost_list[gapcost]

                    test_scores = S[j] + one_mapinfo[i][3] - gapcost

                    filtered = False


                    break
                if(filtered == True):

                    gapcost = abs(readgap - abs(refgap))#mark 1

                    #test_scores = S[j] - skipcost + one_mapinfo[i][3] - min(gapcost/100, 2 * np.log(gapcost))
                    test_scores = S[j] - skipcost + one_mapinfo[i][3] -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))

                '''if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j
                elif((test_scores == max_scores) and (pre_index != -9999999)):
                    pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                    now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                    if(now_refdistance < pre_refdistance):
                        max_scores = test_scores
                        pre_index = j'''
                if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                        T[i] = abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1]))
                    else:
                        T[i] = T[j]

                elif((test_scores == max_scores) and (pre_index != -9999999)):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                        if(abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1])) < T[i]):
                            T[i] = abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1]))
                            pre_index = j
                    else:
                        if(T[j] < T[i]):
                            T[i] = T[j]
                            pre_index = j






                en_loc = st_loc
        
        else:
            
            for j in S_arg[:testspace_en][::-1]:
            
            
                opcount += 1     

                if(S[j] < (max_scores - one_mapinfo[i][3])):

                    break


                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                if((readgap < 0)):
                    #print_log('error')
                    continue

                nocost = False
                filtered = True

                if(one_mapinfo[i][2] == 1):

                            #    y_2 - y_1 - a_1
                    refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])

                else:

                            #     y_1 - y_2 - a_2
                    refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])


                while(True):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):

                        break

                    if((one_mapinfo[i][2] == 1) and (refgap < 0)):

                        break

                    if(one_mapinfo[i][2] == -1):

                        if(refgap > 0):

                            break

                        else:

                            refgap = abs(refgap)

                    gapcost = abs(readgap - refgap)

                    #if((readgap  > maxgap) or (gapcost > maxdiff)):
                    if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                        break


                    gapcost = gapcost_list[gapcost]

                    test_scores = S[j] + one_mapinfo[i][3] - gapcost

                    filtered = False

                    break
                if(filtered == True):

                    #test_scores = S[j] - skipcost + one_mapinfo[i][3]
                    gapcost = abs(readgap - abs(refgap))#mark 1

                    #test_scores = S[j] - skipcost + one_mapinfo[i][3] - min(gapcost/100, 2 * np.log(gapcost))
                    test_scores = S[j] - skipcost + one_mapinfo[i][3] -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))

                '''if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j
                elif((test_scores == max_scores) and (pre_index != -9999999)):
                    pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                    now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                    if(now_refdistance < pre_refdistance):
                        max_scores = test_scores
                        pre_index = j'''
                if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                        T[i] = abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1]))
                    else:
                        T[i] = T[j]

                elif((test_scores == max_scores) and (pre_index != -9999999)):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                        if(abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1])) < T[i]):
                            T[i] = abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1]))
                            pre_index = j
                    else:
                        if(T[j] < T[i]):
                            T[i] = T[j]
                            pre_index = j


            
            
            
            
            
                


        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
        elif(max_scores >= g_max_scores):
            if(T[g_max_index] > T[i]):
                g_max_index = i
        

    k = testspace_en
    while(k < n):

        loc_in_sorted_S = smallorequal2target_1d_point_target(S, S[k], k, S_arg, target_arr[k], target_arr) + 1

        S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
        S_arg[loc_in_sorted_S] = k

        k += 1     
        
    return g_max_index, S, P, S_arg
@njit
def return_main_alignment_size(contig2start, raw_alignment_list):

    preitem = raw_alignment_list[0]
    
    size = 0
     
    st_item = preitem
    for nowitem in raw_alignment_list[1:]:
        if(preitem[2] == nowitem[2]):
            readgap = nowitem[0] - preitem[0] - preitem[3]
            if(readgap < 0):
                continue


            if(preitem[2] == 1):
                refgap = nowitem[1] - preitem[1] - preitem[3]
            else:
                refgap = preitem[1]  - nowitem[1] - nowitem[3]

            if((abs(readgap - refgap) <= 30) and (refgap >= 0)):
                if(pos2contig(preitem[1], contig2start) == pos2contig(nowitem[1], contig2start)):
                    preitem = nowitem
                    continue

        if(preitem[0] - st_item[0] > size):
            size = preitem[0] - st_item[0]
            return_pack = (st_item, preitem)
        preitem = nowitem
        
        st_item = preitem
    if(preitem[0] - st_item[0] > size):
        size = preitem[0] - st_item[0]
        return_pack = (st_item, preitem)

    return return_pack 
#20240415
def decode_hit(index_object, index2contig, testseq, testseq_len, contig2start, kmersize, contig2seq, skipcost = (50., 30.), maxdiff = (50, 30), maxgap = 200, check_num = 20, c_bias = 5000, bin_size = 100, overlapprecentage = 0.5, hastra = True, H = False, mid_occ = -1):            
    redo_ratio = 5
    st = time.time()
    need_reverse, one_mapinfo = get_reversed_chain_numpy_rough(np.array(index_object.map(testseq, check_num = check_num,mid_occ = mid_occ)), testseq_len)
    #print_log(need_reverse)
    #print_log('index_object.map', time.time() - st)
    maxgap = 1000
    #print_log('len(one_mapinfo)', len(one_mapinfo))
    st = time.time()
    path_list, primary_index_List, primary_scores_List, all_index_List, mapq, scores_list =  hit2work_1(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage, hastra, H = H)
    f_redo_ratio = min(abs(testseq_len/(primary_scores_List[0][0]+1e-7)), 20)
    redo_flag = 0
    if(f_redo_ratio > redo_ratio or mapq == 0):
        
        need_reverse, one_mapinfo = get_reversed_chain_numpy_force(need_reverse, np.array(index_object.map(testseq, check_num = 100, mid_occ = mid_occ))[len(one_mapinfo): ], testseq_len)
        #print_log('len(one_mapinfo)', len(one_mapinfo))
        path_list, primary_index_List, primary_scores_List, all_index_List, mapq, scores_list =  hit2work_1_add(path_list, np.array(scores_list), one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage, hastra, H = H)
        redo_flag = 1
    #print_log(scores_list)
    #print_log('hit2work_1', time.time() - st)
    #print_log('mapq', mapq)
    st = time.time()
    if(len(primary_index_List) != 0):
        base_iloc = primary_index_List[0]
        
        if(mapq == 0):
            #print_log('len(primary_scores_List)', len(primary_scores_List))
            base_score = primary_scores_List[0][0]
            rc_testseq = str(Seq(testseq).reverse_complement())
            min_diff = 10
            for tmpiloc in range(len(primary_scores_List[0])):
                if(primary_scores_List[0][tmpiloc] / base_score < 0.999):
                    break
                #print_log(all_index_List[0][tmpiloc])
                #print_log(primary_scores_List[0][tmpiloc])
                #print_log(path_list[all_index_List[0][tmpiloc]][-1], ',', path_list[all_index_List[0][tmpiloc]][0])
                #preitem, nowitem = path_list[all_index_List[0][tmpiloc]][-1], path_list[all_index_List[0][tmpiloc]][0]
                preitem, nowitem = return_main_alignment_size(contig2start, np.array(path_list[all_index_List[0][tmpiloc]][::-1]))
                #print_log(preitem, nowitem)
                if(preitem[2] != nowitem[2]):
                    continue
                if(need_reverse == False):
                    target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start)
                else:
                    target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, rc_testseq, testseq, testseq_len, kmersize, contig2seq, contig2start)
                if(min(len(target), len(query)) == 0):
                    print_log(preitem, nowitem)
                    print_log(testseq)
                diff = edlib.align(query = query, target = target, task = 'distance')['editDistance']/min(len(target), len(query))
                #print_log(diff)
                if(diff <= min_diff):
                    min_diff = diff
                    #print_log('change', base_iloc, '->', all_index_List[0][tmpiloc])
                    base_iloc = all_index_List[0][tmpiloc]
                    primary_index_List[0] = base_iloc
 
                    
                #print_log()
        
        scores, path = hit2work_2(path_list, primary_index_List, primary_scores_List, all_index_List, kmersize, skipcost, maxdiff, maxgap, hastra, contig2start, base_iloc)
    else:
        scores, path = 0., [(0, 0, 0, 0)]
    #print_log('edlib', time.time() - st)
    if(need_reverse == True):
        return mapq, -scores, path, redo_flag
    else:
        return mapq, scores, path, redo_flag
def get_reversed_chain_numpy_force(need_reverse, raw_alignment_array, testseq_len):

    if(need_reverse):

        raw_alignment_array[:,0] = testseq_len - raw_alignment_array[:,0] - raw_alignment_array[:,3]
        raw_alignment_array[:,2] *= -1 
        return True, raw_alignment_array[::-1]
    else:
        return False, raw_alignment_array
@njit
def hit2work_1(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage = 0.5, hastra = False, H = False):
    def getsecond(x):
        return x[1]
    def getfirst(x):
        return x[0]
    def getlength(x):
        return len(x)

    def get_readloc_set_bin(one_mappos, bin_size):
        return set([i[0] // bin_size for i in one_mappos])
    def get_refloc_set_bin(one_mappos, bin_size):
        return set([i[1] // bin_size for i in one_mappos])
    def get_overlapsize(readloc_set_a, readloc_set_b):
        return len(readloc_set_a&readloc_set_b)/min(len(readloc_set_a), len(readloc_set_b))
    #print_log(len(one_mapinfo))


    
    #one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:,1])] 

    cluster_list = List()
    iloc = 0
    preitem = one_mapinfo[iloc]
    st_iloc = 0
    for iloc in range(len(one_mapinfo)):
        nowitem = one_mapinfo[iloc]
        if((abs(nowitem[1] - preitem[1]) > c_bias)):
            if((iloc - st_iloc) < 3):
                continue
            cluster_list.append(one_mapinfo[st_iloc: iloc])
            #print_log(one_mapinfo[st_iloc][1], one_mapinfo[iloc - 1][1])
            st_iloc = iloc
        preitem = nowitem


    if((iloc - st_iloc) > 3):
        cluster_list.append(one_mapinfo[st_iloc: iloc + 1])
    cluster_list.sort(key = get_length)

    cluster_list = cluster_list[::-1][:check_num]
    #print_log(len(cluster_list[0]), cluster_list[0][:,2].sum())


    
    hit = False
    minichain_scores = 40
    path_list = List()
    scores_list = []
    scores_list.append(0.)
    scores_list.pop()
    path_list.append([(0, 0, 0, 0)])
    path_list.pop()
    

    max_scores = 0
    from_repeat = False
    
    fast_enable = False
    g_max_index = 0

    for one_mapinfo in cluster_list:
        repeat = False
        min_ref, max_ref = one_mapinfo[0][1], one_mapinfo[-1][1]
        testcontig = pos2contig(min_ref, contig2start)
        if(testcontig != pos2contig(max_ref, contig2start)):
            #print_log('hit2work: testcontig != pos2contig(max_ref, contig2start)')
            continue
        if(len(one_mapinfo) / testseq_len > 5):
            fast_enable = True
        one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:,0])]
        if(fast_enable == False):
            g_max_index, S, P, S_arg = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_all(one_mapinfo, kmersize = kmersize, skipcost = skipcost[0], maxdiff = maxdiff[0], maxgap = maxgap)
        

        if(fast_enable == True or g_max_index == -1):
            fast_enable = True
            g_max_index, S, P, S_arg = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_fast_all(one_mapinfo, kmersize = kmersize, skipcost = skipcost[0], maxdiff = maxdiff[0], maxgap = maxgap)

        scores = S[g_max_index]
        #if(scores < minichain_scores):
            #if(hit == False):
                #return 0, 0., [(0, 0, 0, 0)]
            #continue
        
        ##################
        usedindex = set()

        path = [(0, 0, 0, 0)]
        path.pop()
        take_index = g_max_index
        usedindex.add(take_index)
        score = S[take_index]
        while(True):
            if((P[take_index] == -9999999)):
                path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
                break
            path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
            take_index = P[take_index]
            usedindex.add(take_index)

        if(score > 40):
            hit = True
            scores_list.append(score)
            path_list.append(path)
        if(scores > max_scores):
            if(repeat == True):
                from_repeat = True
            else:
                from_repeat = False
            max_scores = scores
        #print_log(score, len(path))

        for take_index in S_arg[::-1]:
            if(take_index in usedindex):
                continue
            path = [(0, 0, 0, 0)]
            path.pop()
            usedindex.add(take_index)
            score = S[take_index]
            while(True):
                if((P[take_index] == -9999999)):
                    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
                    break
                path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
                take_index = P[take_index]
                if(take_index in usedindex):
                    score = score - S[take_index]
                    break
                usedindex.add(take_index)
            if(score > 40):
                scores_list.append(score)
                path_list.append(path)
            if(scores > max_scores):
                if(repeat == True):
                    from_repeat = True
                else:
                    from_repeat = False
                max_scores = scores
        ##################
        
        

        

    #print_log(max_scores)
    #print_log()

    if(hit == True and max_scores > 50):

        order = np.argsort(np.array(scores_list))[::-1]
        #print_log('scores_list[order[0]], len(path_list[order[0]])', scores_list[order[0]], len(path_list[order[0]]))
        primary_rlocset_List = List()
        
        primary_scores_List = List()

        
        primary_index_List = List()

        
        all_index_List = List()

        
        iloc = order[0]
        primary_rlocset_List.append(get_readloc_set_bin(path_list[iloc], bin_size))
        primary_scores_List.append(List([scores_list[iloc]]))
        primary_index_List.append(iloc)
        all_index_List.append(List([iloc]))
        #print_log(scores_list[order[0]], path_list[order[0]][-1][1], path_list[order[0]][0][1])
        for iloc in order[1:]:
            
            #print_log(scores_list[iloc], path_list[iloc][-1][1], path_list[iloc][0][1])
            readloc_set_b = get_readloc_set_bin(path_list[iloc], bin_size)
            maxoverlapsize = 0.
            for p_loc in range(len(primary_rlocset_List)):
                tmp_overlapsize = get_overlapsize(primary_rlocset_List[p_loc], readloc_set_b)
                if(tmp_overlapsize > maxoverlapsize):
                    maxoverlapsize = tmp_overlapsize
                    prefer_p_loc = p_loc
            if(maxoverlapsize < overlapprecentage):
                primary_rlocset_List.append(readloc_set_b)
                primary_scores_List.append(List([scores_list[iloc]]))
                all_index_List.append(List([iloc]))
                primary_index_List.append(iloc)
            else:
                primary_scores_List[prefer_p_loc].append(scores_list[iloc])
                all_index_List[prefer_p_loc].append(iloc)
        
        
        m = len(path_list[order[0]])    
        if(len(primary_scores_List[0]) < 2):
            f1 = primary_scores_List[0][0]
            f2 = 0
        else:
            f1 = primary_scores_List[0][0]
            f2 = primary_scores_List[0][1]
        mapq = min(int(40*(1-f2/f1)*min(1, m/10)*np.log(f1)), 60)
        
        return path_list, primary_index_List, primary_scores_List, all_index_List, mapq, scores_list



    else:
        primary_scores_List = List()
        primary_scores_List.append(List([0.]))
        primary_scores_List.pop(0)
        
        primary_index_List = List()
        primary_index_List.append(0)
        primary_index_List.pop(0)
        
        all_index_List = List()
        all_index_List.append(List([0]))
        all_index_List.pop(0)
        return path_list, primary_index_List, primary_scores_List, all_index_List, 0, scores_list

def get_readmap_DP_test(readid, testseq, contig2start, contig2seq, index_object, index2contig, kmersize = 15, local_kmersize = 9, local_maxdiff = 30, refine = True, local_skipcost = 50., golbal_skipcost = (30., 30.),  golbal_maxdiff = (50, 30), check_num = 20, bin_size = 100, hastra = False, debug = False, H = False, mid_occ = -1, redo_ratio = 5, eqx = False):
    kmersize = index_object.k
    local_skipcost = 30.    
    mapq = 60
    setting_kmersize = kmersize
    setting_maxdiff = golbal_maxdiff[1]
    if(H == False):
        local_skipcost += local_kmersize
        golbal_skipcost = (golbal_skipcost[0] + kmersize, golbal_skipcost[1] + kmersize)
    else:
        local_skipcost = 30.
        golbal_skipcost = (30., 30.)
        hastra = True
        check_num = 100
    

    rc_testseq = str(Seq(testseq).reverse_complement())
    testseq_len = len(testseq)
    st = time.time()                            
    mapq, scores, raw_alignment_list, redo_flag = decode_hit(index_object, index2contig, testseq, testseq_len, contig2start, kmersize, contig2seq, skipcost = golbal_skipcost, maxdiff = golbal_maxdiff, maxgap = 200, check_num = check_num, c_bias = 5000, bin_size = bin_size, overlapprecentage = 0.5, hastra = hastra, H = H, mid_occ = mid_occ)



    if(debug == True):
        plot_result = np.array(raw_alignment_list)
        plt.scatter(plot_result[:, 0], plot_result[:, 1])
        plt.show()
    if(scores == 0.):
        return [], ([], []), 0, redo_flag
    
    need_reverse = False
    if((scores < 0.)):
        need_reverse = True

    if(debug == True): print_log('Normal time', time.time() - st, scores)
    if(refine == True):
        setting_maxdiff = local_maxdiff
        setting_kmersize = local_kmersize
        st = time.time()

        #need_reverse, raw_alignment_array = get_reversed_chain_numpy_rough(np.array(raw_alignment_list), testseq_len)
        raw_alignment_array = np.array(raw_alignment_list)
        if(need_reverse == False):
            scores, raw_alignment_list = get_localmap_multi_all_forDP_inv_guide(raw_alignment_array, testseq, rc_testseq, contig2start, contig2seq, kmersize = setting_kmersize, skipcost = local_skipcost, maxdiff = setting_maxdiff, maxgap = 100, shift = 1)

        else:
            scores, raw_alignment_list = get_localmap_multi_all_forDP_inv_guide(raw_alignment_array, rc_testseq, testseq, contig2start, contig2seq, kmersize = setting_kmersize, skipcost = local_skipcost, maxdiff = setting_maxdiff, maxgap = 100, shift = 1)
            testseq, rc_testseq = rc_testseq, testseq
            #raw_alignment_list = get_reversed_chain_numpy(np.array(raw_alignment_list), testseq_len)

        #raw_alignment_list = check_func_clean(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, 10, debug = debug)
        #for item in raw_alignment_list:
            #print_log(item)
        if(len(raw_alignment_list) <= 1):
            return [], ([], []), 0, redo_flag
        if(debug == True): print_log('Refine time', time.time() - st)


        
        

    st = time.time()
    try:
        filtered = False
        #alignment_list, onemapinfolist, TRA_signal = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, debug = debug, H = H)
        alignment_list, onemapinfolist, TRA_signal, filtered = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, debug = debug, H = H, nofilter = False, eqx = eqx)
        if(len(onemapinfolist) == 0):
            return [], ([], []), 0, redo_flag
    except:
        print_log(testseq)
    
    

    if(filtered == True and pairedindel(List([line[-1] for line in onemapinfolist]), indelsize = 30) == True):
        alignment_list, onemapinfolist, TRA_signal, filtered = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, debug = debug, H = H, nofilter = True, eqx = eqx)
        
    if(debug == True): print_log('extend_func', time.time() - st)
    if(len(onemapinfolist) > 1 and refine == 'auto'):
        setting_maxdiff = local_maxdiff
        setting_kmersize = local_kmersize
        st = time.time()                                                                         
        scores, raw_alignment_list = get_localmap_multi_all_forDP_inv_guide(np.array(raw_alignment_list), testseq, rc_testseq, contig2start, contig2seq, kmersize = setting_kmersize, skipcost = local_skipcost, maxdiff = setting_maxdiff, maxgap = 100, shift = 1)
        if(debug == True): print_log('Refine time', time.time() - st, len(raw_alignment_list))

        alignment_list, onemapinfolist, TRA_signal = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, need_reverse, setting_maxdiff, debug = debug, eqx = eqx)

    return onemapinfolist, (alignment_list, raw_alignment_list), TRA_signal, redo_flag
@njit
def hit2work_1_add(path_list, scores_array, one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage = 0.5, hastra = False, H = False):
    def getsecond(x):
        return x[1]
    def getfirst(x):
        return x[0]
    def getlength(x):
        return len(x)

    def get_readloc_set_bin(one_mappos, bin_size):
        return set([i[0] // bin_size for i in one_mappos])
    def get_refloc_set_bin(one_mappos, bin_size):
        return set([i[1] // bin_size for i in one_mappos])
    def get_overlapsize(readloc_set_a, readloc_set_b):
        return len(readloc_set_a&readloc_set_b)/min(len(readloc_set_a), len(readloc_set_b))
    #print_log(len(one_mapinfo))


    
    #one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:,1])] 

    cluster_list = List()
    iloc = 0
    preitem = one_mapinfo[iloc]
    st_iloc = 0
    for iloc in range(len(one_mapinfo)):
        nowitem = one_mapinfo[iloc]
        if((abs(nowitem[1] - preitem[1]) > c_bias)):
            if((iloc - st_iloc) < 3):
                continue
            cluster_list.append(one_mapinfo[st_iloc: iloc])
            #print_log(one_mapinfo[st_iloc][1], one_mapinfo[iloc - 1][1])
            st_iloc = iloc
        preitem = nowitem


    if((iloc - st_iloc) > 3):
        cluster_list.append(one_mapinfo[st_iloc: iloc + 1])
    cluster_list.sort(key = get_length)

    cluster_list = cluster_list[::-1][:check_num]


    
    hit = False
    minichain_scores = 40

    scores_list = [tmp_score for tmp_score in scores_array]


    

    max_scores = 0
    from_repeat = False
    
    fast_enable = False
    g_max_index = 0

    for one_mapinfo in cluster_list:
        repeat = False
        min_ref, max_ref = one_mapinfo[0][1], one_mapinfo[-1][1]
        testcontig = pos2contig(min_ref, contig2start)
        if(testcontig != pos2contig(max_ref, contig2start)):
            #print_log('hit2work: testcontig != pos2contig(max_ref, contig2start)')
            continue
        
        one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:,0])]
        if(fast_enable == False):
            g_max_index, S, P, S_arg = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_all(one_mapinfo, kmersize = kmersize, skipcost = skipcost[0], maxdiff = maxdiff[0], maxgap = maxgap)
        

        if(fast_enable == True or g_max_index == -1):
            fast_enable = True
            g_max_index, S, P, S_arg = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_fast_all(one_mapinfo, kmersize = kmersize, skipcost = skipcost[0], maxdiff = maxdiff[0], maxgap = maxgap)

        scores = S[g_max_index]
        #if(scores < minichain_scores):
            #if(hit == False):
                #return 0, 0., [(0, 0, 0, 0)]
            #continue
        
        ##################
        usedindex = set()

        path = [(0, 0, 0, 0)]
        path.pop()
        take_index = g_max_index
        usedindex.add(take_index)
        score = S[take_index]
        while(True):
            if((P[take_index] == -9999999)):
                path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
                break
            path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
            take_index = P[take_index]
            usedindex.add(take_index)
        if(score > 40):
            hit = True
            scores_list.append(score)
            path_list.append(path)
        if(scores > max_scores):
            if(repeat == True):
                from_repeat = True
            else:
                from_repeat = False
            max_scores = scores
        #print_log(score, len(path))

        for take_index in S_arg[::-1]:
            if(take_index in usedindex):
                continue
            path = [(0, 0, 0, 0)]
            path.pop()
            usedindex.add(take_index)
            score = S[take_index]
            while(True):
                if((P[take_index] == -9999999)):
                    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
                    break
                path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
                take_index = P[take_index]
                if(take_index in usedindex):
                    score = score - S[take_index]
                    break
                usedindex.add(take_index)
            if(score > 40):
                scores_list.append(score)
                path_list.append(path)
            if(scores > max_scores):
                if(repeat == True):
                    from_repeat = True
                else:
                    from_repeat = False
                max_scores = scores
        ##################
        
        

        



    order = np.argsort(np.array(scores_list))[::-1]
    #print_log('scores_list[order[0]], len(path_list[order[0]])', scores_list[order[0]], len(path_list[order[0]]))
    primary_rlocset_List = List()

    primary_scores_List = List()


    primary_index_List = List()


    all_index_List = List()


    iloc = order[0]
    primary_rlocset_List.append(get_readloc_set_bin(path_list[iloc], bin_size))
    primary_scores_List.append(List([scores_list[iloc]]))
    primary_index_List.append(iloc)
    all_index_List.append(List([iloc]))
    #print_log(scores_list[order[0]], path_list[order[0]][-1][1], path_list[order[0]][0][1])
    for iloc in order[1:]:

        #print_log(scores_list[iloc], path_list[iloc][-1][1], path_list[iloc][0][1])
        readloc_set_b = get_readloc_set_bin(path_list[iloc], bin_size)
        maxoverlapsize = 0.
        for p_loc in range(len(primary_rlocset_List)):
            tmp_overlapsize = get_overlapsize(primary_rlocset_List[p_loc], readloc_set_b)
            if(tmp_overlapsize > maxoverlapsize):
                maxoverlapsize = tmp_overlapsize
                prefer_p_loc = p_loc
        if(maxoverlapsize < overlapprecentage):
            primary_rlocset_List.append(readloc_set_b)
            primary_scores_List.append(List([scores_list[iloc]]))
            all_index_List.append(List([iloc]))
            primary_index_List.append(iloc)
        else:
            primary_scores_List[prefer_p_loc].append(scores_list[iloc])
            all_index_List[prefer_p_loc].append(iloc)


    m = len(path_list[order[0]])    
    if(len(primary_scores_List[0]) < 2):
        f1 = primary_scores_List[0][0]
        f2 = 0
    else:
        f1 = primary_scores_List[0][0]
        f2 = primary_scores_List[0][1]
    mapq = min(int(40*(1-f2/f1)*min(1, m/10)*np.log(f1)), 60)

    return path_list, primary_index_List, primary_scores_List, all_index_List, mapq, scores_list

@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_all(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20
    max_kmersize = 0
    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
    
    var_gap = np.arange(1, 10001)#mark 1
    var_gap = np.minimum(var_gap/100, np.log(var_gap))#mark 1


    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)
    T = np.empty(n, np.int64)


    

    

    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    readend_arr = one_mapinfo[:, 0] + one_mapinfo[:, 3]

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_st = 0
    testspace_en = 1
    
    ms_arr = np.empty(one_mapinfo[-1][0] + 1)
    
    test_S = np.empty(0, np.float64)
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0
    max_kmersize = one_mapinfo[i][3]
    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    
    T[0] = 0
    

    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        
        T[i] = 0
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        
        max_kmersize = max(one_mapinfo[i][3], max_kmersize)
        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            if((opcount/one_mapinfo[i][0] > 1000) and one_mapinfo[i][0] > 500):#mark 1
                return -1, S, P, S_arg#mark 1
            
            
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i
            
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
               
                      


        
        for j in S_arg[:testspace_en][::-1]:
            
            
            opcount += 1     #mark 1
            
            if(S[j] < (max_scores - one_mapinfo[i][3])):
                
                break

            #print_log(one_mapinfo[j])
            #print_log(one_mapinfo[i])
            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])
            #print_log(readgap)
            goodsignal = False
            if((readgap < 0)):

                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):

                                #    y_2 - y_1 - a_1
                        refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])


                    else:
                         #1645 17384    -1    15]
                        #[ 1647 17382    -1    15

                                #     y_1 - y_2 - a_2
                        refgap = one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3]
                    if(refgap == readgap):
                        goodsignal = True
                if(goodsignal == False):    
                #print_log('error')
                    continue
            #print_log(goodsignal)  
            #print_log()
            nocost = False
            filtered = True
            
            if(one_mapinfo[i][2] == 1):
                
                        #    y_2 - y_1 - a_1
                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                
            else:
                
                        #     y_1 - y_2 - a_2
                refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                
            if(goodsignal == False):
                while(True):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):

                        break

                    if((one_mapinfo[i][2] == 1) and (refgap < 0)):

                        break

                    if(one_mapinfo[i][2] == -1):

                        if(refgap > 0):

                            break

                        else:

                            refgap = abs(refgap)

                    gapcost = abs(readgap - refgap)

                    #if((readgap  > maxgap) or (gapcost > maxdiff)):
                    if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                        break


                    gapcost = gapcost_list[gapcost]

                    test_scores = S[j] + one_mapinfo[i][3] - gapcost
                    

                    filtered = False

                    break
            else:
                test_scores = S[j] + one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                filtered = False

            if(filtered == True):
                #test_scores = S[j] - skipcost + one_mapinfo[i][3]
                gapcost = abs(readgap - abs(refgap))#mark 1

                #test_scores = S[j] - skipcost + one_mapinfo[i][3] - min(gapcost/100, 2 * np.log(gapcost))
                test_scores = S[j] - skipcost + one_mapinfo[i][3] -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))

                '''if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    #test_scores = S[j] - skipcost + one_mapinfo[i][3]
                    gapcost = abs(readgap - abs(refgap))#mark 1
                    if(gapcost < 10000):#mark 1

                        test_scores = S[j] - skipcost + one_mapinfo[i][3] - var_gap[gapcost]#mark 1
                    else:
                        test_scores = S[j] - skipcost + one_mapinfo[i][3] - 10.#mark 1
                else:
                    test_scores = S[j] - skipcost + one_mapinfo[i][3]'''

            '''if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
            elif((test_scores == max_scores) and (pre_index != -9999999)):
                pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                if(now_refdistance < pre_refdistance):
                    max_scores = test_scores
                    pre_index = j'''
            if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    T[i] = abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1]))
                else:
                    T[i] = T[j]
                    
            elif((test_scores == max_scores) and (pre_index != -9999999)):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    if(abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1])) < T[i]):
                        T[i] = abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1]))
                        pre_index = j
                else:
                    if(T[j] < T[i]):
                        T[i] = T[j]
                        pre_index = j

        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
        elif(max_scores >= g_max_scores):
            if(T[g_max_index] > T[i]):
                g_max_index = i


        
    k = testspace_en
    while(k < n):

        loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

        S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
        S_arg[loc_in_sorted_S] = k

        k += 1

    testspace_en = k  
    

    
    return g_max_index, S, P, S_arg
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_fast_all(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500, fast_t = 5):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20
    max_kmersize = 0
    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
            

    
    var_gap = np.arange(1, 10001)#mark 1
    var_gap = np.minimum(var_gap//100, np.log(var_gap)//1)#mark 1

    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    S_i = np.empty(n, np.int64)
    P = np.empty(n, np.int64)
    T = np.empty(n, np.int64)

    

    

    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    readend_arr = one_mapinfo[:, 0] + one_mapinfo[:, 3]

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    readlength = one_mapinfo[-1][0] + 1000
    target_arr = np.zeros(n, dtype = np.float64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
        if(one_mapinfo[i][2] == 1):
            target_arr[i] = one_mapinfo[i][1] - one_mapinfo[i][0] + readlength 
        else:
            target_arr[i] = -(one_mapinfo[i][1] + one_mapinfo[i][0] + readlength) 
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_st = 0
    testspace_en = 1
    testspace_en_i = 1
    
    ms_arr = np.empty(one_mapinfo[-1][0] + 1)
    
    test_S = np.empty(0, np.float64)
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    S_arg_i = np.empty(n, np.int64)
    S_arg_i[0] = 0
    
    i = 0
    max_kmersize = one_mapinfo[i][3]
    S[i] = one_mapinfo[i][3]
    S_i[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    
    cache_dict = dict()
    cache_dict[0] = 0
    cache_dict.pop(0)
    
    fast = False
    
    T[0] = 0
    
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        T[i] = 0
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        
        max_kmersize = max(one_mapinfo[i][3], max_kmersize)
        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            #if((opcount/one_mapinfo[i][0] > 2000) and one_mapinfo[i][0] > 1000):
                #return 0., [(0, 0, 0, 0)]
                
                
                
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i
                
            k = testspace_en_i
            while(one_mapinfo[k][0] + kmersize <= one_mapinfo[i][0]):

                loc_in_sorted_S = smallorequal2target_1d_point_target(S_i, S_i[k], k, S_arg_i, target_arr[k], target_arr) + 1

                S_arg_i[loc_in_sorted_S + 1: k + 1] = S_arg_i[loc_in_sorted_S: k]
                S_arg_i[loc_in_sorted_S] = k
                
                k += 1

            testspace_en_i = k
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
            
            
            fast = False
            if(testspace_en_i > 0):
                temp_endpos = (testspace_en_i - max(testspace_en_i - 1000, 0))

                scorerange = S_i[S_arg_i[testspace_en_i - 1]] - S_i[S_arg_i[max(testspace_en_i - 1000, 0)]] + 1
                
                if(temp_endpos / scorerange > fast_t):
                    fast = True
                    cache_dict = dict()
                    cache_dict[0] = 0
                    cache_dict.pop(0)
            
            
               
                      
        if(fast == True):

            #print_log(i, prereadloc)

            st_loc = testspace_en_i
            en_loc = testspace_en_i

            while(en_loc > 0):
                if(st_loc not in cache_dict):
                    pre_st_loc = st_loc
                    st_loc = smallorequal2target_1d_point(S_i, S_i[S_arg_i[st_loc-1]]-1e-7, en_loc, S_arg_i) + 1
                    cache_dict[pre_st_loc] = st_loc
                else:
                    st_loc = cache_dict[st_loc]

                if(S_i[S_arg_i[st_loc]] < round(max_scores - one_mapinfo[i][3])):

                    break

                S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg_i, one_mapinfo[i][0], one_mapinfo[:,0])




                j = S_arg_i[S_arg_j]

                opcount += 1     



                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])



                nocost = False
                filtered = True

                if(one_mapinfo[i][2] == 1):

                            #    y_2 - y_1 - a_1
                    refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])

                else:

                            #     y_1 - y_2 - a_2
                    refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])


                while(True):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):

                        break

                    if((one_mapinfo[i][2] == 1) and (refgap < 0)):

                        break

                    if(one_mapinfo[i][2] == -1):

                        if(refgap > 0):

                            break

                        else:

                            refgap = abs(refgap)

                    gapcost = abs(readgap - refgap)

                    if(gapcost > maxdiff):

                        break

                    #if(readgap  > maxgap):
                    if(readgap  > gap_arr[one_mapinfo[j][0]]):#mark 1                   
                        break



                    gapcost = gapcost_list[gapcost]

                    test_scores = S[j] + one_mapinfo[i][3] - gapcost

                    filtered = False


                    break
                if(filtered == True):

                    gapcost = abs(readgap - abs(refgap))#mark 1

                    #test_scores = S[j] - skipcost + one_mapinfo[i][3] - min(gapcost/100, 2 * np.log(gapcost))
                    test_scores = S[j] - skipcost + one_mapinfo[i][3] -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))

                '''if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j
                elif((test_scores == max_scores) and (pre_index != -9999999)):
                    pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                    now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                    if(now_refdistance < pre_refdistance):
                        max_scores = test_scores
                        pre_index = j'''
                if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                        T[i] = abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1]))
                    else:
                        T[i] = T[j]

                elif((test_scores == max_scores) and (pre_index != -9999999)):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                        if(abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1])) < T[i]):
                            T[i] = abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1]))
                            pre_index = j
                    else:
                        if(T[j] < T[i]):
                            T[i] = T[j]
                            pre_index = j






                en_loc = st_loc
        
        else:
            
            for j in S_arg[:testspace_en][::-1]:
            
            
                opcount += 1     #mark 1

                if(S[j] < (max_scores - one_mapinfo[i][3])):

                    break

                #print_log(one_mapinfo[j])
                #print_log(one_mapinfo[i])
                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])
                #print_log(readgap)
                goodsignal = False
                if((readgap < 0)):

                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):

                                    #    y_2 - y_1 - a_1
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])


                        else:
                             #1645 17384    -1    15]
                            #[ 1647 17382    -1    15

                                    #     y_1 - y_2 - a_2
                            refgap = one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3]
                        if(refgap == readgap):
                            goodsignal = True
                    if(goodsignal == False):    
                    #print_log('error')
                        continue
                #print_log(goodsignal)  
                #print_log()
                nocost = False
                filtered = True

                if(one_mapinfo[i][2] == 1):

                            #    y_2 - y_1 - a_1
                    refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])

                else:

                            #     y_1 - y_2 - a_2
                    refgap = -(one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])

                if(goodsignal == False):
                    while(True):
                        if(one_mapinfo[i][2] != one_mapinfo[j][2]):

                            break

                        if((one_mapinfo[i][2] == 1) and (refgap < 0)):

                            break

                        if(one_mapinfo[i][2] == -1):

                            if(refgap > 0):

                                break

                            else:

                                refgap = abs(refgap)

                        gapcost = abs(readgap - refgap)

                        #if((readgap  > maxgap) or (gapcost > maxdiff)):
                        if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                            break


                        gapcost = gapcost_list[gapcost]

                        test_scores = S[j] + one_mapinfo[i][3] - gapcost


                        filtered = False

                        break
                else:
                    test_scores = S[j] + one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                    filtered = False

                if(filtered == True):
                    #test_scores = S[j] - skipcost + one_mapinfo[i][3]
                    gapcost = abs(readgap - abs(refgap))#mark 1

                    #test_scores = S[j] - skipcost + one_mapinfo[i][3] - min(gapcost/100, 2 * np.log(gapcost))
                    test_scores = S[j] - skipcost + one_mapinfo[i][3] -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))

                    '''if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        #test_scores = S[j] - skipcost + one_mapinfo[i][3]
                        gapcost = abs(readgap - abs(refgap))#mark 1
                        if(gapcost < 10000):#mark 1

                            test_scores = S[j] - skipcost + one_mapinfo[i][3] - var_gap[gapcost]#mark 1
                        else:
                            test_scores = S[j] - skipcost + one_mapinfo[i][3] - 10.#mark 1
                    else:
                        test_scores = S[j] - skipcost + one_mapinfo[i][3]'''

                '''if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j
                elif((test_scores == max_scores) and (pre_index != -9999999)):
                    pre_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[pre_index][1])
                    now_refdistance = abs(one_mapinfo[i][1] - one_mapinfo[j][1])
                    if(now_refdistance < pre_refdistance):
                        max_scores = test_scores
                        pre_index = j'''
                if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                        T[i] = abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1]))
                    else:
                        T[i] = T[j]

                elif((test_scores == max_scores) and (pre_index != -9999999)):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                        if(abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1])) < T[i]):
                            T[i] = abs(T[j] - abs(one_mapinfo[i][1] - one_mapinfo[j][1]))
                            pre_index = j
                    else:
                        if(T[j] < T[i]):
                            T[i] = T[j]
                            pre_index = j


            
            
            
            
            
                


        S[i] = max_scores
        S_i[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
        elif(max_scores >= g_max_scores):
            if(T[g_max_index] > T[i]):
                g_max_index = i
        

    k = testspace_en
    while(k < n):

        loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

        S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
        S_arg[loc_in_sorted_S] = k

        k += 1

    
    return g_max_index, S, P, S_arg

#20240602
def get_list_of_readmap_stdout(raw_queue, cooked_queue, minimap, contig2seq, hastra, H, header, option):
    cache_size = 100
    a_list = []
    st = time.time()
    redo_ratio = 5
    
    contig2start = Dict()
    index2contig = List()
    contig2iloc = dict()
    iloc = -1
    for item in minimap.seq_offset:
        iloc += 1
        contig2start[item[0].decode()] = item[2]
        index2contig.append(item[0].decode())
        contig2iloc[item[0].decode()] = iloc

    
    iloc = 0
    unmapcountlist = []
    plotdata = []
    
    rt_list = []
    f_redo_ratio_list = []
    


    rev_read_count = 0

    while(True):
        readidandseq = raw_queue.get()
        if(type(readidandseq) == int):
            break


        try:
            onemapinfolist, (alignment_list,raw_alignment_list), TRA_signal, f_redo_ratio = get_readmap_DP_test(readidandseq[0], readidandseq[1], contig2start, contig2seq, minimap, index2contig, option, hastra = False, redo_ratio = redo_ratio, eqx = option['eqx'], check_num = option['c'])
        except Exception as error:
            #logging.info('Failed to convert to sam record')
            if(option['debug'] == True):
                logging.error(error)
                logging.info(readidandseq[0])
                logging.info(readidandseq[1])
                #break
            continue






        if(len(onemapinfolist) != 0):
            try:
                tmp_a_list = get_bam_dict_str(onemapinfolist, readidandseq[1], readidandseq[2], contig2iloc, contig2seq, option['md'], option['shortcs'], option['cigar2cg'], option['markunbalancetra'])
            except:
                if(option['debug'] == True):
                    logging.info('Failed to convert to sam record')
                    logging.info(readidandseq[0])
                    logging.info(readidandseq[1])
                continue
            if((tmp_a_list) == None):
                continue
            else:
                #if(f_redo_ratio != -1):
                    #f_redo_ratio_list.append(f_redo_ratio)
                a_list += tmp_a_list
                if(len(a_list) > cache_size):
                    cooked_queue.put(a_list)
                    a_list = []
                else:
                    continue
    if(len(a_list) > 0):
        cooked_queue.put(a_list)
    #f_redo_ratio_list = np.array(f_redo_ratio_list)
    #info = 'median: '+ str(np.median(f_redo_ratio_list)) + ' avg: ' + str(f_redo_ratio_list.mean())
    #logging.info(info)
def sortbycontig(x):
    return x[0]
def get_bam_dict_str(mapinfo, query, qual, contig2iloc, contig2seq, md, shortcs, cigar2cg, markunbalancetra):
    #'hhk',         ,  '1', '+', 11, 9192, 2767041, 2776138, 60
    #      0            1    2   3    4      5         6      7
    #'18_19897150_+', '18', '+', 0, 4776, 19832244, 19837393, 1]
    if(markunbalancetra == True):
        mapinfo = reassign_mapq(mapinfo)
    else:
        for iloc in range(len(mapinfo)):
            mapinfo[iloc] = list(mapinfo[iloc])

    

    
    rc_query = str(Seq(query).reverse_complement())
    #mapinfo.sort(key = sort_by_length)
    #mapinfo = mapinfo[::-1]
    mapinfo.sort(key = sortbycontig)
    iloc2nm = dict()
    iloc2md = dict()
    iloc2cs = dict()
    iloc2n_cigar = dict()
    tmpiloc = -1
    if(md == False):
        for item in mapinfo:
            item[-1], n_cigar = mergecigar_n(item[-1])
            tmpiloc += 1
            if(item[2] == '+'):
                nm = compute_NM_tag(query[item[3]: item[4]], get_refseq(item[1], item[5], item[6], contig2seq))
            else:
                nm = compute_NM_tag(rc_query[item[3]: item[4]], get_refseq(item[1], item[5], item[6], contig2seq))
            iloc2nm[tmpiloc] = nm
            iloc2n_cigar[tmpiloc] = n_cigar
    else:
        for item in mapinfo:
            tmpiloc += 1
            if(item[2] == '+'):
                tmp_query = query[item[3]: item[4]]
                tmp_target = get_refseq(item[1], item[5], item[6], contig2seq)
            else:
                tmp_query = rc_query[item[3]: item[4]]
                tmp_target = get_refseq(item[1], item[5], item[6], contig2seq)
            cigarstring, mdstring, csstring, n_cigar = mergecigar_md_cs(item[-1], tmp_target, tmp_query, shortcs)
            nm = compute_NM_tag(tmp_query, tmp_target)
            item[-1] = cigarstring
            iloc2nm[tmpiloc] = nm
            iloc2md[tmpiloc] = mdstring
            iloc2cs[tmpiloc] = csstring
            iloc2n_cigar[tmpiloc] = n_cigar

    if((qual != None) and (len(qual) == len(query))):
        query_qualities = qual
        rc_query_qualities = query_qualities[::-1]
    a_list = []
    primary_iloc = 0
    if(len(mapinfo) > 1):
        if(mapinfo[0][7] == 1 and mapinfo[1][7] != 1):
            primary_iloc = 1

    #QNAME FLAG  RNAME  POS  MAPQ  CIGAR RNEXT  PNEXT  TLEN   SEQ   QUAL
    for iloc in range(len(mapinfo)):
        
        bam_dict = dict()
        primary = mapinfo[iloc]
        bam_dict['QNAME'] = primary[0]
        bam_dict['RNAME'] = primary[1]
        if(iloc == primary_iloc):
            base_value = 0
        else:
            base_value = 2048
        if(primary[2] == '+'):
            bam_dict['FLAG'] = str(base_value)

        else:
            bam_dict['FLAG'] = str(16 + base_value)


        bam_dict['POS'] = str(primary[5] + 1)# SAM Format

        if(iloc2n_cigar[iloc] > 65535):
            if(cigar2cg == True):
                bam_dict['CG'] = primary[8]
                logging.info('Write long CIGAR to CG tag.')
            else:
                bam_dict['CIGAR'] = primary[8]
        else:
            bam_dict['CIGAR'] = primary[8]

        if(len(mapinfo) > 1):
            salist = []
            tmpiloc = -1
            for item in mapinfo:
                tmpiloc += 1
                if(tmpiloc == iloc):
                    continue
                mq = mapinfo[tmpiloc][7]
                if(mq != 0):
                    mq = 60
                else:
                    mq = 1
                nm = iloc2nm[tmpiloc]
                salist.append(''.join((item[1], ',', str(item[5]+1), ',', item[2], ',', item[8], ',', str(mq), ',', str(nm)+';')))

            bam_dict['SA'] = ''.join(salist)
        mq = mapinfo[iloc][7]
        if(mq != 0):
            mq = 60
        else:
            mq = 1
        item = primary

        bam_dict['MAPQ'] = str(mq)

        if(item[2] == '+'):
            bam_dict['SEQ'] = query

            if((qual != None) and (len(qual) == len(query))):
                bam_dict['QUAL'] = query_qualities
        else:
            bam_dict['SEQ'] = rc_query

            if((qual != None) and (len(qual) == len(query))):
                bam_dict['QUAL'] = rc_query_qualities
        bam_dict['NM'] = iloc2nm[iloc]
        if(md == True):
            bam_dict['MD'] = iloc2md[iloc]
            bam_dict['cs'] = iloc2cs[iloc]

        a_list.append(P_alignmentstring(bam_dict))
    return a_list
def reassign_mapq(onemapinfolist):
    iloc = 0
    g_list = [iloc]
    while(g_list[-1] < len(onemapinfolist)-1):
        iloc = g_list[-1]
        test_iloc = iloc
        b_contig = onemapinfolist[iloc][1]
        b_strand = onemapinfolist[iloc][2]
        b_q_st, b_q_en, b_r_st, b_r_en = onemapinfolist[iloc][3], onemapinfolist[iloc][4], onemapinfolist[iloc][5], onemapinfolist[iloc][6]
        hit = False
        #print(b_r_st, b_r_en, b_q_st, b_q_en)
        while(test_iloc+1 < len(onemapinfolist)):
            test_iloc += 1
            #print(test_iloc)
            contig = onemapinfolist[test_iloc][1]
            strand = onemapinfolist[test_iloc][2]
            if(contig != b_contig):
                continue
            q_st, q_en, r_st, r_en = onemapinfolist[test_iloc][3], onemapinfolist[test_iloc][4], onemapinfolist[test_iloc][5], onemapinfolist[test_iloc][6]
            #print(b_r_st, b_r_en, b_q_st, b_q_en)
            #print(r_st, r_en, q_st, q_en)
            if(strand == '+'):
                refgap = r_st - b_r_en
                readgap = q_st - b_q_en

            else:
                refgap = b_r_st - r_en
                readgap = b_q_st - q_en
            if(abs(refgap) > 100000):
                continue
            if(refgap < 10):
                g_list.append(test_iloc)
                hit = True
                break

        #print()
        if(hit == True):
            continue
        else:
            g_list.append(iloc + 1)
    for iloc in range(len(onemapinfolist)):
        onemapinfolist[iloc] = list(onemapinfolist[iloc])
        if(iloc not in g_list):
            
            onemapinfolist[iloc][7] = 0

    return onemapinfolist
            

def decode_hit(index_object, index2contig, testseq, testseq_len, contig2start, kmersize, contig2seq, skipcost = (50., 30.), maxdiff = (50, 30), maxgap = 200, check_num = 20, c_bias = 5000, bin_size = 100, overlapprecentage = 0.5, hastra = True, H = False, mid_occ = -1):            
    redo_ratio = 5
    st = time.time()
    need_reverse, one_mapinfo = get_reversed_chain_numpy_rough(np.array(index_object.map(testseq, check_num = check_num,mid_occ = mid_occ)), testseq_len)
    #print_log(need_reverse)
    #print_log('index_object.map', time.time() - st)
    maxgap = 1000
    #print_log('len(one_mapinfo)', len(one_mapinfo))
    st = time.time()
    path_list, primary_index_List, primary_scores_List, all_index_List, mapq, scores_list, factor =  hit2work_1(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage, hastra, H = H)
    if(len(path_list) == 0):
        scores, path = 0., [(0, 0, 0, 0)]
        return 0, scores, path, factor
    f_redo_ratio = min(abs(testseq_len/(primary_scores_List[0][0]+1e-7)), 20)
    redo_flag = 0
    if(check_num < 100 or check_num != -1):
        if(f_redo_ratio > redo_ratio or mapq == 0):

            need_reverse, one_mapinfo = get_reversed_chain_numpy_rough(np.array(index_object.map(testseq, check_num = 100, mid_occ = mid_occ)), testseq_len)
            #print_log('len(one_mapinfo)', len(one_mapinfo))
            path_list, primary_index_List, primary_scores_List, all_index_List, mapq, scores_list, factor = hit2work_1(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage, hastra, H = H)
            redo_flag = 1
    #print_log(scores_list)
    #print_log('hit2work_1', time.time() - st)
    #print_log('mapq', mapq)
    st = time.time()
    if(len(primary_index_List) != 0):
        base_iloc = primary_index_List[0]
        
        if(mapq == 0):
            #print_log('len(primary_scores_List)', len(primary_scores_List))
            base_score = primary_scores_List[0][0]
            rc_testseq = str(Seq(testseq).reverse_complement())
            min_diff = 10
            for tmpiloc in range(len(primary_scores_List[0])):
                if(primary_scores_List[0][tmpiloc] / base_score < 0.999):
                    break
                #print_log(all_index_List[0][tmpiloc])
                #print_log(primary_scores_List[0][tmpiloc])
                #print_log(path_list[all_index_List[0][tmpiloc]][-1], ',', path_list[all_index_List[0][tmpiloc]][0])
                #preitem, nowitem = path_list[all_index_List[0][tmpiloc]][-1], path_list[all_index_List[0][tmpiloc]][0]
                preitem, nowitem = return_main_alignment_size(contig2start, np.array(path_list[all_index_List[0][tmpiloc]][::-1]))
                #print_log(preitem, nowitem)
                if(preitem[2] != nowitem[2]):
                    continue
                if(need_reverse == False):
                    target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start)
                else:
                    target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, rc_testseq, testseq, testseq_len, kmersize, contig2seq, contig2start)
                if(min(len(target), len(query)) == 0):
                    print_log(preitem, nowitem)
                    print_log(testseq)
                diff = edlib.align(query = query, target = target, task = 'distance')['editDistance']/min(len(target), len(query))
                #print_log(diff)
                if(diff <= min_diff):
                    min_diff = diff
                    #print_log('change', base_iloc, '->', all_index_List[0][tmpiloc])
                    base_iloc = all_index_List[0][tmpiloc]
                    primary_index_List[0] = base_iloc
 
                    
                #print_log()
        
        scores, path = scores_list[base_iloc], path_list[base_iloc]
    else:
        scores, path = 0., [(0, 0, 0, 0)]
    #print_log('edlib', time.time() - st)
    if(need_reverse == True):
        return mapq, -scores, path, factor
    else:
        return mapq, scores, path, factor
def get_reversed_chain_numpy_force(need_reverse, raw_alignment_array, testseq_len):

    if(need_reverse):

        raw_alignment_array[:,0] = testseq_len - raw_alignment_array[:,0] - raw_alignment_array[:,3]
        raw_alignment_array[:,2] *= -1 
        return True, raw_alignment_array[::-1]
    else:
        return False, raw_alignment_array
@njit
def hit2work_1(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage = 0.5, hastra = False, H = False):
    def getsecond(x):
        return x[1]
    def getfirst(x):
        return x[0]
    def getlength(x):
        return len(x)

    def get_readloc_set_bin(one_mappos, bin_size):
        return set([i[0] // bin_size for i in one_mappos])
    def get_refloc_set_bin(one_mappos, bin_size):
        return set([i[1] // bin_size for i in one_mappos])
    def get_overlapsize(readloc_set_a, readloc_set_b):
        return len(readloc_set_a&readloc_set_b)/min(len(readloc_set_a), len(readloc_set_b))
    #print_log(len(one_mapinfo))


    
    #one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:,1])] 

    cluster_list = List([one_mapinfo])

    #print_log(len(cluster_list[0]), cluster_list[0][:,2].sum())


    
    hit = False
    minichain_scores = 40
    path_list = List()
    scores_list = []
    scores_list.append(0.)
    scores_list.pop()
    path_list.append([(0, 0, 0, 0)])
    path_list.pop()
    

    max_scores = 0
    from_repeat = False
    
    fast_enable = False
    g_max_index = 0
    factor = -1.
    for one_mapinfo in cluster_list:
        repeat = False

        if(len(one_mapinfo) / testseq_len > 5):
            fast_enable = True
        one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:,0])]
        if(fast_enable == False):
            g_max_index, S, P, S_arg, factor = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_all(one_mapinfo, kmersize = kmersize, skipcost = skipcost[0], maxdiff = maxdiff[0], maxgap = maxgap)
        

        if(fast_enable == True or g_max_index == -1):
            fast_enable = True
            g_max_index, S, P, S_arg = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_fast_all(one_mapinfo, kmersize = kmersize, skipcost = skipcost[0], maxdiff = maxdiff[0], maxgap = maxgap)

        scores = S[g_max_index]
        #if(scores < minichain_scores):
            #if(hit == False):
                #return 0, 0., [(0, 0, 0, 0)]
            #continue
        
        ##################
        usedindex = set()

        path = [(0, 0, 0, 0)]
        path.pop()
        take_index = g_max_index
        usedindex.add(take_index)
        score = S[take_index]
        while(True):
            if((P[take_index] == -9999999)):
                path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
                break
            path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
            take_index = P[take_index]
            usedindex.add(take_index)

        if(score > 40):
            hit = True
            scores_list.append(score)
            path_list.append(path)
        if(scores > max_scores):
            if(repeat == True):
                from_repeat = True
            else:
                from_repeat = False
            max_scores = scores
        #print_log(score, len(path))

        for take_index in S_arg[::-1]:
            if(take_index in usedindex):
                continue
            path = [(0, 0, 0, 0)]
            path.pop()
            usedindex.add(take_index)
            score = S[take_index]
            while(True):
                if((P[take_index] == -9999999)):
                    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
                    break
                path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
                take_index = P[take_index]
                if(take_index in usedindex):
                    score = score - S[take_index]
                    break
                usedindex.add(take_index)
            if(score > 40):
                scores_list.append(score)
                path_list.append(path)
            if(scores > max_scores):
                if(repeat == True):
                    from_repeat = True
                else:
                    from_repeat = False
                max_scores = scores
        ##################
        
        

        

    #print_log(max_scores)
    #print_log()

    if(hit == True and max_scores > 50):

        order = np.argsort(np.array(scores_list))[::-1]
        #print_log('scores_list[order[0]], len(path_list[order[0]])', scores_list[order[0]], len(path_list[order[0]]))
        primary_rlocset_List = List()
        
        primary_scores_List = List()

        
        primary_index_List = List()

        
        all_index_List = List()

        
        iloc = order[0]
        primary_rlocset_List.append(get_readloc_set_bin(path_list[iloc], bin_size))
        primary_scores_List.append(List([scores_list[iloc]]))
        primary_index_List.append(iloc)
        all_index_List.append(List([iloc]))
        #print_log(scores_list[order[0]], path_list[order[0]][-1][1], path_list[order[0]][0][1])
        for iloc in order[1:]:
            
            #print_log(scores_list[iloc], path_list[iloc][-1][1], path_list[iloc][0][1])
            readloc_set_b = get_readloc_set_bin(path_list[iloc], bin_size)
            maxoverlapsize = 0.
            for p_loc in range(len(primary_rlocset_List)):
                tmp_overlapsize = get_overlapsize(primary_rlocset_List[p_loc], readloc_set_b)
                if(tmp_overlapsize > maxoverlapsize):
                    maxoverlapsize = tmp_overlapsize
                    prefer_p_loc = p_loc
            if(maxoverlapsize < overlapprecentage):
                primary_rlocset_List.append(readloc_set_b)
                primary_scores_List.append(List([scores_list[iloc]]))
                all_index_List.append(List([iloc]))
                primary_index_List.append(iloc)
            else:
                primary_scores_List[prefer_p_loc].append(scores_list[iloc])
                all_index_List[prefer_p_loc].append(iloc)
        
        
        m = len(path_list[order[0]])    
        if(len(primary_scores_List[0]) < 2):
            f1 = primary_scores_List[0][0]
            f2 = 0
        else:
            f1 = primary_scores_List[0][0]
            f2 = primary_scores_List[0][1]
        mapq = min(int(40*(1-f2/f1)*min(1, m/10)*np.log(f1)), 60)
        
        return path_list, primary_index_List, primary_scores_List, all_index_List, mapq, scores_list, factor



    else:
        path_list = List()
        scores_list = []
        scores_list.append(0.)
        scores_list.pop()
        path_list.append([(0, 0, 0, 0)])
        path_list.pop()
        
        primary_scores_List = List()
        primary_scores_List.append(List([0.]))
        primary_scores_List.pop(0)
        
        primary_index_List = List()
        primary_index_List.append(0)
        primary_index_List.pop(0)
        
        all_index_List = List()
        all_index_List.append(List([0]))
        all_index_List.pop(0)
        return path_list, primary_index_List, primary_scores_List, all_index_List, 0, scores_list, factor

def get_readmap_DP_test(readid, testseq, contig2start, contig2seq, index_object, index2contig, option, kmersize = 15, local_kmersize = 9, local_maxdiff = 30, refine = True, local_skipcost = 59., golbal_skipcost = (45., 45.),  golbal_maxdiff = (30, 30), check_num = 20, bin_size = 100, hastra = False, debug = False, H = False, mid_occ = -1, redo_ratio = 5, eqx = False):
    kmersize = index_object.k
    golbal_skipcost = (option['golbal_skipcost'], option['golbal_skipcost'])
    golbal_maxdiff = (option['golbal_maxdiff'], option['golbal_maxdiff'])
    local_skipcost = option['local_skipcost']
    local_maxdiff = option['local_maxdiff']
    mapq = 60
    setting_kmersize = kmersize
    setting_maxdiff = golbal_maxdiff[1]

    

    rc_testseq = str(Seq(testseq).reverse_complement())
    testseq_len = len(testseq)
    st = time.time()                            
    mapq, scores, raw_alignment_list, redo_flag = decode_hit(index_object, index2contig, testseq, testseq_len, contig2start, kmersize, contig2seq, skipcost = golbal_skipcost, maxdiff = golbal_maxdiff, maxgap = 200, check_num = check_num, c_bias = 5000, bin_size = bin_size, overlapprecentage = 0.5, hastra = hastra, H = H, mid_occ = mid_occ)



    if(debug == True):
        plot_result = np.array(raw_alignment_list)
        plt.scatter(plot_result[:, 0], plot_result[:, 1])
        plt.show()
    if(scores == 0.):
        return [], ([], []), 0, redo_flag
    
    need_reverse = False
    if((scores < 0.)):
        need_reverse = True

    if(debug == True): print_log('Normal time', time.time() - st, scores)
    if(refine == True):
        setting_maxdiff = local_maxdiff
        setting_kmersize = local_kmersize
        st = time.time()

        #need_reverse, raw_alignment_array = get_reversed_chain_numpy_rough(np.array(raw_alignment_list), testseq_len)
        raw_alignment_array = np.array(raw_alignment_list)
        if(need_reverse == False):
            scores, raw_alignment_list = get_localmap_multi_all_forDP_inv_guide(raw_alignment_array, testseq, rc_testseq, contig2start, contig2seq, kmersize = setting_kmersize, skipcost = local_skipcost, maxdiff = setting_maxdiff, maxgap = 50, shift = 1)

        else:
            scores, raw_alignment_list = get_localmap_multi_all_forDP_inv_guide(raw_alignment_array, rc_testseq, testseq, contig2start, contig2seq, kmersize = setting_kmersize, skipcost = local_skipcost, maxdiff = setting_maxdiff, maxgap = 50, shift = 1)
            testseq, rc_testseq = rc_testseq, testseq
            #raw_alignment_list = get_reversed_chain_numpy(np.array(raw_alignment_list), testseq_len)

        #raw_alignment_list = check_func_clean(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, 10, debug = debug)
        #for item in raw_alignment_list:
            #print_log(item)
        if(len(raw_alignment_list) <= 1):
            return [], ([], []), 0, redo_flag
        if(debug == True): print_log('Refine time', time.time() - st)


        
        

    st = time.time()

    filtered = False
    #alignment_list, onemapinfolist, TRA_signal = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, debug = debug, H = H)
    alignment_list, onemapinfolist, TRA_signal, filtered = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, maxdiffratio = option['maxdivergence'], debug = debug, H = H, nofilter = option['nodiscard'], eqx = eqx)
    if(len(onemapinfolist) == 0):
        return [], ([], []), 0, redo_flag

    
    

    if(option['nodiscard'] == False and filtered == True and pairedindel(List([line[-1] for line in onemapinfolist]), indelsize = 30) == True):
        alignment_list, onemapinfolist, TRA_signal, filtered = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, maxdiffratio = option['maxdivergence'], debug = debug, H = H, nofilter = True, eqx = eqx)
        
    if(debug == True): print_log('extend_func', time.time() - st)
    

    return onemapinfolist, (alignment_list, raw_alignment_list), TRA_signal, redo_flag

@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20
    max_kmersize = 0
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
    

    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)
    T = np.empty(n, np.int64)
    

    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    readend_arr = one_mapinfo[:, 0] + one_mapinfo[:, 3]

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_st = 0
    testspace_en = 1
    
    ms_arr = np.empty(one_mapinfo[-1][0] + 1)
    
    test_S = np.empty(0, np.float64)
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0
    max_kmersize = one_mapinfo[i][3]
    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    T[0] = 0
    c_repeat_weight = 0

    
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        
        T[i] = 0
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        
        max_kmersize = max(one_mapinfo[i][3], max_kmersize)
        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i            
            c_repeat_weight = coverage_dict[one_mapinfo[i][0]]
            skipcost = oskipcost + c_repeat_weight
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
               
                      


        

        for j in S_arg[:testspace_en][::-1]:
            
              
            
            if(S[j] < (max_scores - one_mapinfo[i][3])):
                
                break


            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            if((readgap < 0)):
                continue
                
            nocost = False
            filtered = True
            if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                if(one_mapinfo[i][2] == 1):
                    refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                else:
                    refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
            else:
                if(one_mapinfo[j][2] == -1):
                    refgap = one_mapinfo[i][1] - one_mapinfo[j][1]
                else:
                    refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - one_mapinfo[j][3]

   
            
            while(refgap >= 0):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    
                    break
                    
                gapcost = abs(readgap - refgap)
                
                if((readgap  > maxgap) or (gapcost > maxdiff)):                 
                    break
                                        
                gapcost = gapcost_list[gapcost]
                
                test_scores = S[j] + one_mapinfo[i][3] - gapcost
                
                filtered = False
                
                break
            if(filtered == True):
                gapcost = abs(readgap - refgap)#mark 1
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    invweight = 0
                    if(T[j] != 0):
                        invweight = min(abs(T[j] - gapcost)*0.3, 20)
                    test_scores = S[j] - invweight - skipcost + one_mapinfo[i][3] -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))
                else:
                    test_scores = S[j] - skipcost + one_mapinfo[i][3] -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))

            if(test_scores > max_scores):

                max_scores = test_scores
                pre_index = j
                
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    if(T[j] != 0):
                        T[i] = 0
                    else:
                        T[i] = gapcost
                else:
                    T[i] = T[j]
                    
            elif((test_scores == max_scores) and (pre_index != -9999999)):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    if(T[j] != 0):
                        if(abs(T[j] - gapcost) < abs(T[pre_index] - gapcost)):
                            pre_index = j
                            


        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i

      

        
            
        
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    while(True):
        if((P[take_index] == -9999999)):
            break
        take_index = P[take_index]
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

    return g_max_scores, path

@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_all(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20
    max_kmersize = 0
    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
    
    var_gap = np.arange(1, 10001)#mark 1
    var_gap = np.minimum(var_gap/100, np.log(var_gap))#mark 1


    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)
    T = np.empty(n, np.int64)


    

    

    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    readend_arr = one_mapinfo[:, 0] + one_mapinfo[:, 3]

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_st = 0
    testspace_en = 1
    
    ms_arr = np.empty(one_mapinfo[-1][0] + 1)
    
    test_S = np.empty(0, np.float64)
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0
    max_kmersize = one_mapinfo[i][3]
    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    
    T[0] = 0
    

    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        
        T[i] = 0
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        
        max_kmersize = max(one_mapinfo[i][3], max_kmersize)
        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            if((opcount/one_mapinfo[i][0] > 1000) and one_mapinfo[i][0] > 500):#mark 1
                return -1, S, P, S_arg, -1.#mark 1
            
            
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i
            
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
               
                      


        
        for j in S_arg[:testspace_en][::-1]:
            
            
            opcount += 1     #mark 1
            
            if(S[j] < (max_scores - one_mapinfo[i][3])):
                
                break


            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            goodsignal = False
            if((readgap < 0)):

                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):
                        refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                    else:
                        refgap = one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3]
                    if(refgap == readgap):
                        goodsignal = True
                if(goodsignal == False):    
                    continue
            nocost = False
            filtered = True
            
            if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                if(one_mapinfo[i][2] == 1):

                    refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                else:
                    refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
            else:
                if(one_mapinfo[j][2] == -1):
                    refgap = one_mapinfo[i][1] - one_mapinfo[j][1]
                else:
                    refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - one_mapinfo[j][3]

                
            if(goodsignal == False):
                while(refgap >= 0):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                        break                   
                    gapcost = abs(readgap - refgap)

                    if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                        break


                    gapcost = gapcost_list[gapcost]

                    test_scores = S[j] + one_mapinfo[i][3] - gapcost
                    

                    filtered = False

                    break
            else:
                test_scores = S[j] + one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                filtered = False

            if(filtered == True):

                gapcost = abs(readgap - refgap)#mark 1

                test_scores = S[j] - skipcost + one_mapinfo[i][3] -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))


            if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
                

        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
       


        
    k = testspace_en
    while(k < n):

        loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

        S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
        S_arg[loc_in_sorted_S] = k

        k += 1

    testspace_en = k  
    

    
    return g_max_index, S, P, S_arg, opcount/len(one_mapinfo)
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_fast_all(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500, fast_t = 5):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20
    max_kmersize = 0
    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
            

    
    var_gap = np.arange(1, 10001)#mark 1
    var_gap = np.minimum(var_gap//100, np.log(var_gap)//1)#mark 1

    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    S_i = np.empty(n, np.int64)
    P = np.empty(n, np.int64)
    T = np.empty(n, np.int64)

    

    

    opcount = 0
    
    prereadloc = one_mapinfo[0][0]
    
    readend_arr = one_mapinfo[:, 0] + one_mapinfo[:, 3]

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    readlength = one_mapinfo[-1][0] + 1000
    target_arr = np.zeros(n, dtype = np.float64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
        if(one_mapinfo[i][2] == 1):
            target_arr[i] = one_mapinfo[i][1] - one_mapinfo[i][0] + readlength 
        else:
            target_arr[i] = -(one_mapinfo[i][1] + one_mapinfo[i][0] + readlength) 
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_st = 0
    testspace_en = 1
    testspace_en_i = 1
    
    ms_arr = np.empty(one_mapinfo[-1][0] + 1)
    
    test_S = np.empty(0, np.float64)
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    S_arg_i = np.empty(n, np.int64)
    S_arg_i[0] = 0
    
    i = 0
    max_kmersize = one_mapinfo[i][3]
    S[i] = one_mapinfo[i][3]
    S_i[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    
    cache_dict = dict()
    cache_dict[0] = 0
    cache_dict.pop(0)
    
    fast = False
    
    T[0] = 0
    
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        T[i] = 0
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        
        max_kmersize = max(one_mapinfo[i][3], max_kmersize)
        
        
        if(prereadloc < one_mapinfo[i][0]):                
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i
                
            k = testspace_en_i
            while(one_mapinfo[k][0] + kmersize <= one_mapinfo[i][0]):

                loc_in_sorted_S = smallorequal2target_1d_point_target(S_i, S_i[k], k, S_arg_i, target_arr[k], target_arr) + 1

                S_arg_i[loc_in_sorted_S + 1: k + 1] = S_arg_i[loc_in_sorted_S: k]
                S_arg_i[loc_in_sorted_S] = k
                
                k += 1

            testspace_en_i = k
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
            
            
            fast = False
            if(testspace_en_i > 0):
                temp_endpos = (testspace_en_i - max(testspace_en_i - 1000, 0))

                scorerange = S_i[S_arg_i[testspace_en_i - 1]] - S_i[S_arg_i[max(testspace_en_i - 1000, 0)]] + 1
                
                if(temp_endpos / scorerange > fast_t):
                    fast = True
                    cache_dict = dict()
                    cache_dict[0] = 0
                    cache_dict.pop(0)
            
            
               
                      
        if(fast == True):

            #print_log(i, prereadloc)

            st_loc = testspace_en_i
            en_loc = testspace_en_i

            while(en_loc > 0):
                if(st_loc not in cache_dict):
                    pre_st_loc = st_loc
                    st_loc = smallorequal2target_1d_point(S_i, S_i[S_arg_i[st_loc-1]]-1e-7, en_loc, S_arg_i) + 1
                    cache_dict[pre_st_loc] = st_loc
                else:
                    st_loc = cache_dict[st_loc]

                if(S_i[S_arg_i[st_loc]] < round(max_scores - one_mapinfo[i][3])):

                    break

                S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg_i, one_mapinfo[i][0], one_mapinfo[:,0])




                j = S_arg_i[S_arg_j]

                opcount += 1     



                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])



                nocost = False
                filtered = True

                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):

                        refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                    else:
                        refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                else:
                    if(one_mapinfo[j][2] == -1):
                        refgap = one_mapinfo[i][1] - one_mapinfo[j][1]
                    else:
                        refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - one_mapinfo[j][3]



                while(refgap >= 0):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                        break                  
                    gapcost = abs(readgap - refgap)

                    if(gapcost > maxdiff):
                        break

                    if(readgap  > gap_arr[one_mapinfo[j][0]]):#mark 1                   
                        break

                    gapcost = gapcost_list[gapcost]

                    test_scores = S[j] + one_mapinfo[i][3] - gapcost

                    filtered = False


                    break
                if(filtered == True):

                    gapcost = abs(readgap - refgap)#mark 1

                    test_scores = S[j] - skipcost + one_mapinfo[i][3] -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))

                if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j
                    


                en_loc = st_loc
        
        else:
            
            for j in S_arg[:testspace_en][::-1]:
            
            
                opcount += 1     #mark 1

                if(S[j] < (max_scores - one_mapinfo[i][3])):

                    break


                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])
                goodsignal = False
                if((readgap < 0)):

                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                        else:
                            refgap = one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3]
                        if(refgap == readgap):
                            goodsignal = True
                    if(goodsignal == False):    
                        continue
                nocost = False
                filtered = True

                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):

                        refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                    else:
                        refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                else:
                    if(one_mapinfo[j][2] == -1):
                        refgap = one_mapinfo[i][1] - one_mapinfo[j][1]
                    else:
                        refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - one_mapinfo[j][3]

                if(goodsignal == False):
                    while(refgap >= 0):
                        if(one_mapinfo[i][2] != one_mapinfo[j][2]):

                            break

                        
                        gapcost = abs(readgap - refgap)

                        if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                            break


                        gapcost = gapcost_list[gapcost]

                        test_scores = S[j] + one_mapinfo[i][3] - gapcost


                        filtered = False

                        break
                else:
                    test_scores = S[j] + one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                    filtered = False

                if(filtered == True):
                    gapcost = abs(readgap - refgap)#mark 1
                    test_scores = S[j] - skipcost + one_mapinfo[i][3] -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))
                if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j

        S[i] = max_scores
        S_i[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i

        

    k = testspace_en
    while(k < n):

        loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

        S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
        S_arg[loc_in_sorted_S] = k

        k += 1

    
    return g_max_index, S, P, S_arg

def extend_func(raw_alignment_list, readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, maxdiffratio, debug = False, H = False, nofilter = False, eqx = False):

    TRA_signal = False
    #step 1
    #rebuild chain break
    #remove anchor cause large cost and small alignment
    st = time.time()
    if(H == False):
        alignment_list = rebuild_chain_break(contig2start, raw_alignment_list, large_cost = setting_maxdiff, small_alignment = 50, small_dup = -100)
    else:
        alignment_list = rebuild_chain_break_H(contig2start, raw_alignment_list, large_cost = setting_maxdiff, small_alignment = 30, small_dup = -30)
        nofilter = True
    
        
    tmpiloc = -1
    while((tmpiloc + 1) < len(alignment_list)):
        tmpiloc += 1
        preitem, nowitem = alignment_list[tmpiloc][0], alignment_list[tmpiloc][-1]
        target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, setting_kmersize, contig2seq, contig2start)
        diffratio = edlib.align(query = query, target = target, task = 'distance')['editDistance']/min(len(target), len(query))
        if((diffratio>maxdiffratio)):
            alignment_list.pop(tmpiloc)
            tmpiloc -= 1
            
    #^print_log('step 1: rebuild chain break ', time.time() - st)
    if(debug == True):
        print_log('step 1: rebuild chain break')
        for line in alignment_list:
            tempcontig = pos2contig(line[0][1], contig2start)
            temprefbias = contig2start[tempcontig]
            preitem, nowitem = line[0], line[-1]

            target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, setting_kmersize, contig2seq, contig2start)

            print_log(preitem, nowitem)
            diffratio = edlib.align(query = query, target = target, task = 'distance')['editDistance']/min(len(target), len(query))
            print_log(line[0][0], line[-1][0], line[0][1] - temprefbias, line[-1][1] - temprefbias, diffratio)
            plot_result = np.array(line)
            plt.scatter(plot_result[:,0], plot_result[:,1])
        plt.show()
        
    #step 2
    #extend edge to recover small dup and misplaced alignment
    # and also to merge alignment gaped with deletion or insertion
    st = time.time()
    extend_edge_test(testseq, testseq_len, alignment_list, setting_kmersize, pos2contig, contig2start, contig2seq, san = 1, debug = debug)
    if(debug == True):print_log('step 2: extend edge ', time.time() - st)
    if(debug == True):
        print_log('After extend edge')
        for line in alignment_list:
            tempcontig = pos2contig(line[0][1], contig2start)
            temprefbias = contig2start[tempcontig]
            print_log(line[0][0], line[-1][0], line[0][1] - temprefbias, line[-1][1] - temprefbias)
        
    #step 3
    #remove miss placed alignment which cause del/ins or ins/del in start and end
    st = time.time()
    
    o_alignment_list_len = len(alignment_list)
    filtered = False
    if((len(alignment_list) > 2) and (nofilter == False)):    
        iloc = 0
        while(iloc < (len(alignment_list) - 2)):
            removed = drop_misplaced_alignment_test(alignment_list, iloc, debug = debug)
            if(removed == True):
                continue
            else:
                iloc += 1
                
    
        #extend_edge_drop_test(testseq, testseq_len, alignment_list, setting_kmersize, pos2contig, contig2start, contig2seq, san = 1, debug = debug)
        
    if(debug == True):print_log('step 3: remove miss placed alignment ', time.time() - st)
    if(debug == True):
        print_log('After remove miss placed alignment')
        for line in alignment_list:
            tempcontig = pos2contig(line[0][1], contig2start)
            temprefbias = contig2start[tempcontig]
            print_log(line[0][0], line[-1][0], line[0][1] - temprefbias, line[-1][1] - temprefbias)
            plot_result = np.array(line)
            plt.scatter(plot_result[:,0], plot_result[:,1])
        plt.show()
            

            
    
    if(len(alignment_list) <  o_alignment_list_len):#fill the gap
        filtered = True
        st = time.time()
        extend_edge_test(testseq, testseq_len, alignment_list, setting_kmersize, pos2contig, contig2start, contig2seq, san = 1, debug = debug)
        if(debug == True):print_log('step 4: fill the gap by extend edge ', time.time() - st)
        if(debug == True):
            print_log('After extend edge')
            for line in alignment_list:
                tempcontig = pos2contig(line[0][1], contig2start)
                temprefbias = contig2start[tempcontig]
                print_log(line[0][0], line[-1][0], line[0][1] - temprefbias, line[-1][1] - temprefbias)

        
    #step 4
    #merge conjacent alignment with small readgap or refgap
    st = time.time()
    merge_smallgap = 2000
    too_large_gap = 5000
    if(len(alignment_list) >= 2):
        iloc = 0
        duplist = getdupiloc(alignment_list)
        while((iloc + 1) < len(alignment_list)):
            if(iloc in duplist):
                iloc += 1
                continue
            
            preitem = alignment_list[iloc][-1]
            nowitem = alignment_list[iloc + 1][0]
            if(preitem[2] != nowitem[2] or (pos2contig(preitem[1], contig2start) != pos2contig(nowitem[1], contig2start))):
                iloc += 1
                continue
            presize = alignment_list[iloc][-1][0] + alignment_list[iloc][-1][3] - alignment_list[iloc][0][0]
            nowsize = alignment_list[iloc + 1][-1][0] + alignment_list[iloc + 1][-1][3] - alignment_list[iloc + 1][0][0]
            
            if(alignment_list[iloc][-1][2] == 1):
                ref_presize = alignment_list[iloc][-1][1] + alignment_list[iloc][-1][3] - alignment_list[iloc][0][1]
            else:
                ref_presize = alignment_list[iloc][0][1] + alignment_list[iloc][0][3] - alignment_list[iloc][-1][1]
            if(alignment_list[iloc+1][-1][2] == 1):
                ref_nowsize = alignment_list[iloc+1][-1][1] + alignment_list[iloc+1][-1][3] - alignment_list[iloc+1][0][1]
            else:
                ref_nowsize = alignment_list[iloc+1][0][1] + alignment_list[iloc+1][0][3] - alignment_list[iloc+1][-1][1]
                
            readgap = nowitem[0] - preitem[0] - preitem[3]
            if(preitem[2] == 1):
                refgap = nowitem[1] - preitem[1] - preitem[3]

            else:
                refgap = preitem[1]  - nowitem[1] - nowitem[3]                        
            if(refgap < 0):
                iloc += 1
                continue
            #if((min(readgap, refgap) < merge_smallgap) and (max(readgap, refgap) < too_large_gap)):
            #if((((min(presize, nowsize) > 500) or ((abs(readgap - refgap) / min(presize, nowsize)) < 0.5)) and (max(readgap, refgap) < 20000))):
            if((min(readgap, refgap) < 50) and (abs(readgap - refgap) < 2000)):
                alignment_list[iloc] = List_merge((alignment_list[iloc], alignment_list[iloc + 1]))
                alignment_list.pop(iloc+1)
            else:
                iloc += 1
    if(debug == True):print_log('step 4: merge conjacent alignment ', time.time() - st)
    if(debug == True):
        print_log('After merge conjacent alignment')
        for line in alignment_list:
            tempcontig = pos2contig(line[0][1], contig2start)
            temprefbias = contig2start[tempcontig]
            print_log(line[0][0], line[-1][0], line[0][1] - temprefbias, line[-1][1] - temprefbias)
    
    #step 5
    #split unrelate read
    #note inside alignment there is no way two conjacent anchors
    #has cost large than 30
    #so there is no big gap
    st = time.time()
    new_alignment_list = List()
    cigarlist = []
    for alignment in alignment_list: 

        tmp_alignment_list, tmp_cigarlist = split_alignment_test(alignment, testseq, rc_testseq, testseq_len, kmersize=setting_kmersize , contig2seq = contig2seq, contig2start = contig2start, debug = debug, H = H, eqx = eqx)
        if(debug): print_log(len(tmp_alignment_list), len(tmp_cigarlist))
        iloc = -1
        for alignment in tmp_alignment_list:
            iloc += 1
            new_alignment_list.append(alignment)
            cigarlist.append(tmp_cigarlist[iloc])


    if(debug == True):print_log('step 5: split unrelate read ', time.time() - st)
    if(debug == True):
        print_log('After split unrelate read')
        for line in new_alignment_list:
            tempcontig = pos2contig(line[0][1], contig2start)
            temprefbias = contig2start[tempcontig]
            print_log(line[0][0], line[-1][0], line[0][1] - temprefbias, line[-1][1] - temprefbias)
    
    alignment_list, onemapinfolist = get_onemapinfolist(new_alignment_list, cigarlist, readid, mapq, testseq_len, contig2start, need_reverse)
    return alignment_list, onemapinfolist, TRA_signal, filtered
@njit
def rebuild_chain_break(contig2start, raw_alignment_list, large_cost, small_alignment = 50, small_dup = -100):
    #rebuild chain break
    #step 1
    #remove anchor cause large cost and small alignment
    preitem = raw_alignment_list[0]
    alignment_list = List([List([preitem])])
    for nowitem in raw_alignment_list[1:]:
        if(preitem[2] == nowitem[2]):
            readgap = nowitem[0] - preitem[0] - preitem[3]


            if(preitem[2] == 1):
                refgap = nowitem[1] - preitem[1] - preitem[3]
            else:
                refgap = preitem[1]  - nowitem[1] - nowitem[3]


            if((abs(readgap - refgap) <= large_cost) and (refgap >= 0) and (readgap < 100)):
                if(pos2contig(preitem[1], contig2start) == pos2contig(nowitem[1], contig2start)):
                    alignment_list[-1].append(nowitem)
                    preitem = nowitem
                    continue

        if(len(alignment_list[-1]) == 1):
            alignment_list.pop(-1)
        if(len(alignment_list) > 0):
            if((alignment_list[-1][-1][0] + alignment_list[-1][-1][3] - alignment_list[-1][0][0]) < small_alignment):
                alignment_list.pop(-1)
  
        alignment_list.append(List([nowitem]))
        preitem = nowitem
    if(len(alignment_list[-1]) == 1):
        alignment_list.pop(-1)
    if((alignment_list[-1][-1][0] + alignment_list[-1][-1][3] - alignment_list[-1][0][0]) < small_alignment):
        alignment_list.pop(-1)


    return alignment_list

@njit
def get_localmap_multi_all_forDP_inv_guide(raw_alignment_array, testseq, rc_testseq, contig2start, contig2seq, kmersize, skipcost, maxdiff, maxgap, shift = 1):#
    




    def seq2hashtable_multi_test(onelookuptable, seq, start, kmersize):
        skiphash = hash('N'*kmersize)
        for iloc in range(0, len(seq) - kmersize + 1, 1):
            hashedkmer = hash(seq[iloc:iloc+kmersize])
            if(skiphash == hashedkmer):
                continue
            if(hashedkmer in onelookuptable):

                onelookuptable[hashedkmer].append(start + iloc)
            else:
                onelookuptable[hashedkmer] = List([start + iloc])

    readgap = 0
    pre = raw_alignment_array[0]
    for now in raw_alignment_array[1:]:
        if(abs(now[0] - pre[0]) > readgap):
            readgap = abs(now[0] - pre[0])
        pre = now
    readgap += 1000  
    readgap = max(readgap, 5000)
    raw_alignment_array = raw_alignment_array[np.argsort(raw_alignment_array[:, 1])]
    startandend = List([(raw_alignment_array[0][1], raw_alignment_array[0][1])])
    for item in raw_alignment_array[1:]:
        if(((item[1] - startandend[-1][1]) < readgap)):
            startandend[-1] = (startandend[-1][0], item[1])
        else:
            if(startandend[-1][0] == startandend[-1][1]):
                startandend.pop(-1)
            startandend.append((item[1], item[1]))
    if(startandend[-1][0] == startandend[-1][1]):
        startandend.pop(-1)
    local_lookuptable = Dict()
    local_lookuptable[0] = List([0])
    
    retry_diffcontig = False
    for item in startandend:
        min_ref, max_ref = item[0], item[1]
        testcontig = pos2contig(min_ref, contig2start)
        if(testcontig != pos2contig(max_ref, contig2start)):
            retry_diffcontig = True
            break
        lookfurther = min(2000, min_ref-contig2start[testcontig])
        min_ref -= lookfurther
        max_ref += 2000
        refseq = contig2seq[testcontig][min_ref-contig2start[testcontig]: max_ref-contig2start[testcontig]]
        seq2hashtable_multi_test(local_lookuptable, refseq, min_ref, kmersize)
    local_lookuptable.pop(0)
    if(retry_diffcontig == True):
        startandend = List([(raw_alignment_array[0][1], raw_alignment_array[0][1])])
        current_contig = pos2contig(raw_alignment_array[0][1], contig2start)
        for item in raw_alignment_array[1:]:
            if(((item[1] - startandend[-1][1]) < readgap) and (current_contig == pos2contig(item[1], contig2start))):
                startandend[-1] = (startandend[-1][0], item[1])
            else:
                if(startandend[-1][0] == startandend[-1][1]):
                    startandend.pop(-1)
                startandend.append((item[1], item[1]))
                current_contig = pos2contig(item[1], contig2start)
        if(startandend[-1][0] == startandend[-1][1]):
            startandend.pop(-1)
        local_lookuptable = Dict()
        local_lookuptable[0] = List([0])

        retry_diffcontig = False
        for item in startandend:
            min_ref, max_ref = item[0], item[1]
            testcontig = pos2contig(min_ref, contig2start)
            if(testcontig != pos2contig(max_ref, contig2start)):
                retry_diffcontig = True
                break
            lookfurther = min(2000, min_ref-contig2start[testcontig])
            min_ref -= lookfurther
            max_ref += 2000
            refseq = contig2seq[testcontig][min_ref-contig2start[testcontig]: max_ref-contig2start[testcontig]]
            seq2hashtable_multi_test(local_lookuptable, refseq, min_ref, kmersize)
        local_lookuptable.pop(0)
        
    raw_alignment_array = raw_alignment_array[np.argsort(raw_alignment_array[:, 0])]
    readstart = max(0, raw_alignment_array[0][0]-500)
    readend = min(len(testseq)-kmersize+1, raw_alignment_array[-1][0]+500)

    
    
    iloc = readstart
    iloc -= shift
    
    maxop = 0
    pointdict = Dict()
    
    while(True):
   
        iloc += shift
        if(iloc >= readend):
            break

        if(hash(testseq[iloc: iloc + kmersize]) == hash(rc_testseq[-(iloc + kmersize): -iloc])):
            continue
        biasvalue, biasvalue_1, closest_index, closest_index_1 = findClosest(raw_alignment_array, target = iloc)


        interval = min(biasvalue + biasvalue_1 + 500, 2000)
        upperrefloc =  (raw_alignment_array[closest_index][1] + interval, raw_alignment_array[closest_index_1][1] + interval)
        lowerrefloc =  (raw_alignment_array[closest_index][1] - interval, raw_alignment_array[closest_index_1][1] - interval)
        readgap = abs(iloc - raw_alignment_array[closest_index][0])
   
        hashedkmer = hash(testseq[iloc: iloc + kmersize])  
        if(hashedkmer in local_lookuptable):
            for refloc in local_lookuptable[hashedkmer]:
                refgap = abs(refloc - raw_alignment_array[closest_index][1])
                diff = abs(readgap - refgap)
                if((diff < 500) or (upperrefloc[0] >= refloc and lowerrefloc[0] <= refloc) or (upperrefloc[1] >= refloc and lowerrefloc[1] <= refloc)):
                    item = (iloc, refloc, 1)
                    point = item[1] - item[0]

                    if(point in pointdict):

                        if((pointdict[point][-1][0] + pointdict[point][-1][3]) >= item[0]):
                            bouns = item[0] - (pointdict[point][-1][0] + pointdict[point][-1][3]) + kmersize
                            if(bouns > 0):
                                if((pointdict[point][-1][3] + bouns < 20)):
                                    pointdict[point][-1] = (pointdict[point][-1][0], pointdict[point][-1][1], 1, pointdict[point][-1][3] + bouns)
                                else:
                                    pointdict[point].append((pointdict[point][-1][0] + pointdict[point][-1][3], pointdict[point][-1][1] + pointdict[point][-1][3], 1, bouns))
                          


                        else:    
                            pointdict[point].append((item[0], item[1], item[2], kmersize))
                    else:
                        pointdict[point] = List([(item[0], item[1], item[2], kmersize)])


        hashedkmer = hash(rc_testseq[-(iloc + kmersize): -iloc]) 
        if(hashedkmer in local_lookuptable):

            for refloc in local_lookuptable[hashedkmer]:
                refgap = abs(refloc - raw_alignment_array[closest_index][1])
                diff = abs(readgap - refgap)
                if((diff < 500) or (upperrefloc[0] >= refloc and lowerrefloc[0] <= refloc) or (upperrefloc[1] >= refloc and lowerrefloc[1] <= refloc)):
                    item = (iloc, refloc, -1)
                    point = -(item[1] + item[0])
                    if(point in pointdict):

                        if((pointdict[point][-1][0] + pointdict[point][-1][3]) >= item[0]):
                            bouns = item[0] - (pointdict[point][-1][0] + pointdict[point][-1][3]) + kmersize
                            if(bouns > 0):
                                if((pointdict[point][-1][3] + bouns < 20)):
                                    pointdict[point][-1] = (pointdict[point][-1][0], item[1], -1, pointdict[point][-1][3] + bouns)
                                else:

                                    pointdict[point].append((pointdict[point][-1][0] + pointdict[point][-1][3], item[1], -1, bouns))
                                        #print_log(pointdict[point][-1])


                        else:    
                            pointdict[point].append((item[0], item[1], item[2], kmersize))
                    else:
                        pointdict[point] = List([(item[0], item[1], item[2], kmersize)])
    
    one_mapinfo = [(-1, -1, -1, -1)]
    one_mapinfo.pop(0)
    for key in pointdict:
        for item in pointdict[key]:
            one_mapinfo.append(item)
    one_mapinfo = np.array(one_mapinfo)
    one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:, 0])]
  


    return get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo, kmersize = kmersize, skipcost = skipcost, maxdiff = maxdiff, maxgap = maxgap)


@njit
def smallorequal2target_1d_point(arr, target, n, point):
    if(target < arr[point[0]]):
        return -1

    if(target >= arr[point[n - 1]]):
        return n-1

    i = 0
    j = n
    mid = 0
    while(i < j):
        mid = (i + j) // 2

        if(target == arr[point[mid]]):
            if(mid < n - 1):
                if(arr[point[mid+1]] > target):
                    return mid
                else:
                    i = mid + 1
            else:
                return mid

        elif(target < arr[point[mid]]) :

            if(mid > 0 and target >= arr[point[mid - 1]]):
                return mid-1

            j = mid

        else: # target>arr[mid][0]
            if(mid < n - 1 and target < arr[point[mid + 1]]):
                return mid

            i = mid + 1

    return mid
#20240618
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20

    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
    

    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)
    T = np.empty(n, np.int64)
    

    opcount = 0
    
    

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)

    testspace_en = 1
    

    
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0

    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    T[0] = 0
    c_repeat_weight = 0

    
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        
        T[i] = 0
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        

        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i            
            c_repeat_weight = coverage_dict[one_mapinfo[i][0]]
            skipcost = oskipcost + c_repeat_weight
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
               
                      


        

        for j in S_arg[:testspace_en][::-1]:
            
              
            
            if(S[j] < (max_scores - one_mapinfo[i][3])):
                
                break


            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            if((readgap < 0)):
                continue
                
            nocost = False
            filtered = True
            if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                if(one_mapinfo[i][2] == 1):
                    refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                else:
                    refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
            else:
                if(one_mapinfo[j][2] == -1):
                    refgap = one_mapinfo[i][1] - one_mapinfo[j][1]
                else:
                    refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - one_mapinfo[j][3]

   
            
            while(refgap >= 0):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    
                    break
                    
                gapcost = abs(readgap - refgap)
                
                if((readgap  > maxgap) or (gapcost > maxdiff)):                 
                    break
                                        
                gapcost = gapcost_list[gapcost]
                
                test_scores = S[j] + one_mapinfo[i][3] - gapcost
                
                filtered = False
                
                break
            if(filtered == True):
                gapcost = abs(readgap - refgap)#mark 1
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    invweight = 0
                    if(T[j] != 0):
                        invweight = min(abs(T[j] - gapcost)*0.3, 20)
                    test_scores = S[j] - invweight - skipcost + one_mapinfo[i][3] -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))
                else:
                    test_scores = S[j] - skipcost + one_mapinfo[i][3] -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))

            if(test_scores > max_scores):

                max_scores = test_scores
                pre_index = j
                
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    if(T[j] != 0):
                        T[i] = 0
                    else:
                        T[i] = gapcost
                else:
                    T[i] = T[j]
                    
            elif((test_scores == max_scores) and (pre_index != -9999999)):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    if(T[j] != 0):
                        if(abs(T[j] - gapcost) < abs(T[pre_index] - gapcost)):
                            pre_index = j
                            


        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i

      

        
            
        
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    while(True):
        if((P[take_index] == -9999999)):
            break
        take_index = P[take_index]
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))

    return g_max_scores, path

@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_all(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20

    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
    



    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)



    

    

    opcount = 0
    
    

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)

    testspace_en = 1
    
    
    test_S = np.empty(0, np.float64)
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0

    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    

    

    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        


        max_scores = one_mapinfo[i][3]

        
        pre_index = -9999999
        

        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            if((opcount/one_mapinfo[i][0] > 1000) and one_mapinfo[i][0] > 500):#mark 1
                return -1, S, P, S_arg, -1.#mark 1
            
            
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i
            
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
               
                      


        
        for j in S_arg[:testspace_en][::-1]:
            
            
            opcount += 1     #mark 1
            
            if(S[j] < (max_scores - one_mapinfo[i][3])):
                
                break


            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            goodsignal = False
            if((readgap < 0)):

                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):
                        refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                    else:
                        refgap = one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3]
                    if(refgap == readgap):
                        goodsignal = True
                if(goodsignal == False):    
                    continue
            nocost = False
            filtered = True
            
            if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                if(one_mapinfo[i][2] == 1):

                    refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                else:
                    refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
            else:
                if(one_mapinfo[j][2] == -1):
                    refgap = one_mapinfo[i][1] - one_mapinfo[j][1]
                else:
                    refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - one_mapinfo[j][3]

                
            if(goodsignal == False):
                while(refgap >= 0):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                        break                   
                    gapcost = abs(readgap - refgap)

                    if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                        break


                    gapcost = gapcost_list[gapcost]

                    test_scores = S[j] + one_mapinfo[i][3] - gapcost
                    

                    filtered = False

                    break
            else:
                test_scores = S[j] + one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                filtered = False

            if(filtered == True):

                gapcost = abs(readgap - refgap)#mark 1

                test_scores = S[j] - skipcost + one_mapinfo[i][3] -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))


            if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
                

        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
       


        
    k = testspace_en
    while(k < n):

        loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

        S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
        S_arg[loc_in_sorted_S] = k

        k += 1

    testspace_en = k  
    

    
    return g_max_index, S, P, S_arg, opcount/len(one_mapinfo)
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_fast_all(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500, fast_t = 5):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20

    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
            

    

    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    S_i = np.empty(n, np.int64)
    P = np.empty(n, np.int64)

    

    

    
    

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    readlength = one_mapinfo[-1][0] + 1000
    target_arr = np.zeros(n, dtype = np.float64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
        if(one_mapinfo[i][2] == 1):
            target_arr[i] = one_mapinfo[i][1] - one_mapinfo[i][0] + readlength 
        else:
            target_arr[i] = -(one_mapinfo[i][1] + one_mapinfo[i][0] + readlength) 
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_en = 1
    testspace_en_i = 1
    
    
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    S_arg_i = np.empty(n, np.int64)
    S_arg_i[0] = 0
    
    i = 0

    S[i] = one_mapinfo[i][3]
    S_i[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    
    cache_dict = dict()
    cache_dict[0] = 0
    cache_dict.pop(0)
    
    fast = False
    
    
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])

        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        

        
        
        if(prereadloc < one_mapinfo[i][0]):                
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i
                
            k = testspace_en_i
            while(k < i):

                loc_in_sorted_S = smallorequal2target_1d_point_target(S_i, S_i[k], k, S_arg_i, target_arr[k], target_arr) + 1

                S_arg_i[loc_in_sorted_S + 1: k + 1] = S_arg_i[loc_in_sorted_S: k]
                S_arg_i[loc_in_sorted_S] = k
                
                k += 1

            testspace_en_i = k
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
            
            
            fast = False
            if(testspace_en_i > 0):
                temp_endpos = (testspace_en_i - max(testspace_en_i - 1000, 0))

                scorerange = S_i[S_arg_i[testspace_en_i - 1]] - S_i[S_arg_i[max(testspace_en_i - 1000, 0)]] + 1
                
                if(temp_endpos / scorerange > fast_t):
                    fast = True
                    cache_dict = dict()
                    cache_dict[0] = 0
                    cache_dict.pop(0)
            
            
               
                      
        if(fast == True):

            #print_log(i, prereadloc)

            st_loc = testspace_en_i
            en_loc = testspace_en_i

            while(en_loc > 0):
                if(st_loc not in cache_dict):
                    pre_st_loc = st_loc
                    st_loc = smallorequal2target_1d_point(S_i, S_i[S_arg_i[st_loc-1]]-1e-7, en_loc, S_arg_i) + 1
                    cache_dict[pre_st_loc] = st_loc
                else:
                    st_loc = cache_dict[st_loc]

                if(S_i[S_arg_i[st_loc]] < round(max_scores - one_mapinfo[i][3])):

                    break

                S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg_i, one_mapinfo[i][0], one_mapinfo[:,0])




                j = S_arg_i[S_arg_j]
    


                #########
                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])
                goodsignal = False
                if((readgap < 0)):

                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                        else:
                            refgap = one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3]
                        if(refgap == readgap):
                            goodsignal = True
                    if(goodsignal == False):  
                        en_loc = st_loc
                        continue
                nocost = False
                filtered = True

                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):

                        refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                    else:
                        refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                else:
                    if(one_mapinfo[j][2] == -1):
                        refgap = one_mapinfo[i][1] - one_mapinfo[j][1]
                    else:
                        refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - one_mapinfo[j][3]

                if(goodsignal == False):
                    while(refgap >= 0):
                        if(one_mapinfo[i][2] != one_mapinfo[j][2]):

                            break

                        
                        gapcost = abs(readgap - refgap)

                        if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                            break


                        gapcost = gapcost_list[gapcost]

                        test_scores = S[j] + one_mapinfo[i][3] - gapcost


                        filtered = False

                        break
                else:
                    test_scores = S[j] + one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                    filtered = False

                if(filtered == True):
                    gapcost = abs(readgap - refgap)#mark 1
                    test_scores = S[j] - skipcost + one_mapinfo[i][3] -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))
                if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j
                    


                en_loc = st_loc
        
        else:
            
            for j in S_arg[:testspace_en][::-1]:
            
            

                if(S[j] < (max_scores - one_mapinfo[i][3])):

                    break


                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])
                goodsignal = False
                if((readgap < 0)):

                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                        else:
                            refgap = one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3]
                        if(refgap == readgap):
                            goodsignal = True
                    if(goodsignal == False):    
                        continue
                nocost = False
                filtered = True

                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):

                        refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                    else:
                        refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                else:
                    if(one_mapinfo[j][2] == -1):
                        refgap = one_mapinfo[i][1] - one_mapinfo[j][1]
                    else:
                        refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - one_mapinfo[j][3]

                if(goodsignal == False):
                    while(refgap >= 0):
                        if(one_mapinfo[i][2] != one_mapinfo[j][2]):

                            break

                        
                        gapcost = abs(readgap - refgap)

                        if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                            break


                        gapcost = gapcost_list[gapcost]

                        test_scores = S[j] + one_mapinfo[i][3] - gapcost


                        filtered = False

                        break
                else:
                    test_scores = S[j] + one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                    filtered = False

                if(filtered == True):
                    gapcost = abs(readgap - refgap)#mark 1
                    test_scores = S[j] - skipcost + one_mapinfo[i][3] -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))
                if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j

        S[i] = max_scores
        S_i[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i

        

    k = testspace_en
    while(k < n):

        loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

        S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
        S_arg[loc_in_sorted_S] = k

        k += 1

    
    return g_max_index, S, P, S_arg

@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_all(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20

    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
    



    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)



    

    

    opcount = 0
    
    

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)

    testspace_en = 1
    
    
    test_S = np.empty(0, np.float64)
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0

    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    

    

    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        


        max_scores = one_mapinfo[i][3]

        
        pre_index = -9999999
        

        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            if((opcount/one_mapinfo[i][0] > 1000) and one_mapinfo[i][0] > 500):#mark 1
                return -1, S, P, S_arg, -1.#mark 1
            
            
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i
            
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
               
                      


        
        for j in S_arg[:testspace_en][::-1]:
            
            
            opcount += 1     #mark 1
            
            if(S[j] < (max_scores - one_mapinfo[i][3])):
                
                break


            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            goodsignal = False
            if((readgap < 0)):

                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):
                        refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                    else:
                        refgap = one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3]
                    if(refgap == readgap):
                        goodsignal = True
                if(goodsignal == False):    
                    continue
            nocost = False
            filtered = True
            
            if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                if(one_mapinfo[i][2] == 1):

                    refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                else:
                    refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
            else:
                if(one_mapinfo[j][2] == -1):
                    refgap = one_mapinfo[i][1] - one_mapinfo[j][1]
                else:
                    refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - one_mapinfo[j][3]

                
            if(goodsignal == False):
                while(refgap >= 0):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                        break                   
                    gapcost = abs(readgap - refgap)

                    if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                        break


                    gapcost = gapcost_list[gapcost]

                    test_scores = S[j] + one_mapinfo[i][3] - gapcost
                    

                    filtered = False

                    break
            else:
                test_scores = S[j] + one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                filtered = False

            if(filtered == True):

                gapcost = abs(readgap - refgap)#mark 1

                test_scores = S[j] - min(max(30, gapcost), skipcost) + one_mapinfo[i][3] -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))


            if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
                

        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
       


        
    k = testspace_en
    while(k < n):

        loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

        S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
        S_arg[loc_in_sorted_S] = k

        k += 1

    testspace_en = k  
    

    
    return g_max_index, S, P, S_arg, opcount/len(one_mapinfo)
#20240629
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_fast_all(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500, fast_t = 5):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20

    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
            

    

    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    S_i = np.empty(n, np.int64)
    P = np.empty(n, np.int64)

    

    

    
    

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    readlength = one_mapinfo[-1][0] + 1000
    target_arr = np.zeros(n, dtype = np.float64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
        if(one_mapinfo[i][2] == 1):
            target_arr[i] = one_mapinfo[i][1] - one_mapinfo[i][0] + readlength 
        else:
            target_arr[i] = -(one_mapinfo[i][1] + one_mapinfo[i][0] + readlength) 
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_en = 1
    testspace_en_i = 1
    
    
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    S_arg_i = np.empty(n, np.int64)
    S_arg_i[0] = 0
    
    i = 0

    S[i] = one_mapinfo[i][3]
    S_i[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    
    cache_dict = dict()
    cache_dict[0] = 0
    cache_dict.pop(0)
    
    fast = False
    
    
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])

        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        

        
        
        if(prereadloc < one_mapinfo[i][0]):                
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i
                
            k = testspace_en_i
            while(k < i):

                loc_in_sorted_S = smallorequal2target_1d_point_target(S_i, S_i[k], k, S_arg_i, target_arr[k], target_arr) + 1

                S_arg_i[loc_in_sorted_S + 1: k + 1] = S_arg_i[loc_in_sorted_S: k]
                S_arg_i[loc_in_sorted_S] = k
                
                k += 1

            testspace_en_i = k
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
            
            
            fast = False
            if(testspace_en_i > 0):
                temp_endpos = (testspace_en_i - max(testspace_en_i - 1000, 0))

                scorerange = S_i[S_arg_i[testspace_en_i - 1]] - S_i[S_arg_i[max(testspace_en_i - 1000, 0)]] + 1
                
                if(temp_endpos / scorerange > fast_t):
                    fast = True
                    cache_dict = dict()
                    cache_dict[0] = 0
                    cache_dict.pop(0)
            
            
               
                      
        if(fast == True):

            #print_log(i, prereadloc)

            st_loc = testspace_en_i
            en_loc = testspace_en_i

            while(en_loc > 0):
                if(st_loc not in cache_dict):
                    pre_st_loc = st_loc
                    st_loc = smallorequal2target_1d_point(S_i, S_i[S_arg_i[st_loc-1]]-1e-7, en_loc, S_arg_i) + 1
                    cache_dict[pre_st_loc] = st_loc
                else:
                    st_loc = cache_dict[st_loc]

                if(S_i[S_arg_i[st_loc]] < round(max_scores - one_mapinfo[i][3])):

                    break

                S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg_i, one_mapinfo[i][0], one_mapinfo[:,0])




                j = S_arg_i[S_arg_j]
    


                #########
                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                if((readgap < 0)):
                    bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                    readgap = 0
                    non_overlap_size = one_mapinfo[i][0] - one_mapinfo[j][0]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - non_overlap_size)
                        else:
                            refgap = (one_mapinfo[j][1] + one_mapinfo[j][3] - non_overlap_size - one_mapinfo[i][1] - one_mapinfo[i][3])
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] + one_mapinfo[j][3] - non_overlap_size - one_mapinfo[j][1]
                        else:
                            refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - non_overlap_size

                else:
                    bonus = one_mapinfo[i][3]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                        else:
                            refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] - one_mapinfo[j][1]
                        else:
                            refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - one_mapinfo[j][3]

                filtered = True

  
                while(refgap >= 0):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):

                        break


                    gapcost = abs(readgap - refgap)

                    if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                        break


                    gapcost = gapcost_list[gapcost]

                    test_scores = S[j] + bonus - gapcost


                    filtered = False

                    break


                if(filtered == True):
                    gapcost = abs(readgap - refgap)#mark 1
                    test_scores = S[j] - skipcost + bonus -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))
                if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j
                    


                en_loc = st_loc
        
        else:
            
            for j in S_arg[:testspace_en][::-1]:
            
            

                if(S[j] < (max_scores - one_mapinfo[i][3])):

                    break


                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                if((readgap < 0)):
                    bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                    readgap = 0
                    non_overlap_size = one_mapinfo[i][0] - one_mapinfo[j][0]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - non_overlap_size)
                        else:
                            refgap = (one_mapinfo[j][1] + one_mapinfo[j][3] - non_overlap_size - one_mapinfo[i][1] - one_mapinfo[i][3])
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] + one_mapinfo[j][3] - non_overlap_size - one_mapinfo[j][1]
                        else:
                            refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - non_overlap_size

                else:
                    bonus = one_mapinfo[i][3]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                        else:
                            refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] - one_mapinfo[j][1]
                        else:
                            refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - one_mapinfo[j][3]

                filtered = True

  
                while(refgap >= 0):
                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):

                        break


                    gapcost = abs(readgap - refgap)

                    if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                        break


                    gapcost = gapcost_list[gapcost]

                    test_scores = S[j] + bonus - gapcost


                    filtered = False

                    break


                if(filtered == True):
                    gapcost = abs(readgap - refgap)#mark 1
                    test_scores = S[j] - skipcost + bonus -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))
                if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j

        S[i] = max_scores
        S_i[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i

        

    k = testspace_en
    while(k < n):

        loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

        S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
        S_arg[loc_in_sorted_S] = k

        k += 1

    
    return g_max_index, S, P, S_arg

@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_all(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20

    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
    



    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)



    

    

    opcount = 0
    
    

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)

    testspace_en = 1
    
    
    test_S = np.empty(0, np.float64)
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0

    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    

    

    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        


        max_scores = one_mapinfo[i][3]

        
        pre_index = -9999999
        

        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            if((opcount/one_mapinfo[i][0] > 1000) and one_mapinfo[i][0] > 500):#mark 1
                return -1, S, P, S_arg, -1.#mark 1
            
            
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i
            
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
               
                      


        
        for j in S_arg[:testspace_en][::-1]:
            
            
            opcount += 1     #mark 1
            
            if(S[j] < (max_scores - one_mapinfo[i][3])):
                
                break


            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            if((readgap < 0)):
                bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                readgap = 0
                non_overlap_size = one_mapinfo[i][0] - one_mapinfo[j][0]
                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):
                        refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - non_overlap_size)
                    else:
                        refgap = (one_mapinfo[j][1] + one_mapinfo[j][3] - non_overlap_size - one_mapinfo[i][1] - one_mapinfo[i][3])
                else:
                    if(one_mapinfo[j][2] == -1):
                        refgap = one_mapinfo[i][1] + one_mapinfo[j][3] - non_overlap_size - one_mapinfo[j][1]
                    else:
                        refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - non_overlap_size

            else:
                bonus = one_mapinfo[i][3]
                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):
                        refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                    else:
                        refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                else:
                    if(one_mapinfo[j][2] == -1):
                        refgap = one_mapinfo[i][1] - one_mapinfo[j][1]
                    else:
                        refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - one_mapinfo[j][3]
                
            filtered = True

                

            while(refgap >= 0):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    break                   
                gapcost = abs(readgap - refgap)

                if((readgap  > gap_arr[one_mapinfo[j][0]]) or (gapcost > maxdiff)):#mark 1                   
                    break


                gapcost = gapcost_list[gapcost]

                test_scores = S[j] + bonus - gapcost


                filtered = False

                break


            if(filtered == True):

                gapcost = abs(readgap - refgap)#mark 1

                test_scores = S[j] - skipcost + bonus -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))


            if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
                

        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
       


        
    k = testspace_en
    while(k < n):

        loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

        S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
        S_arg[loc_in_sorted_S] = k

        k += 1

    testspace_en = k  
    

    
    return g_max_index, S, P, S_arg, opcount/len(one_mapinfo)
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20

    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
    

    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)
    T = np.empty(n, np.int64)
    

    opcount = 0
    
    

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)

    testspace_en = 1
    

    
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0

    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    T[0] = 0
    c_repeat_weight = 0

    
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        
        T[i] = 0
        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        

        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i            
            c_repeat_weight = coverage_dict[one_mapinfo[i][0]]
            skipcost = oskipcost + c_repeat_weight
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
               
                      


        

        for j in S_arg[:testspace_en][::-1]:
            
              
            
            if(S[j] < (max_scores - one_mapinfo[i][3])):
                
                break
            if(one_mapinfo[i][0] == one_mapinfo[j][0]):
                continue

            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            if((readgap < 0)):
                bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                readgap = 0
                non_overlap_size = one_mapinfo[i][0] - one_mapinfo[j][0]
                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):
                        refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - non_overlap_size)
                    else:
                        refgap = (one_mapinfo[j][1] + one_mapinfo[j][3] - non_overlap_size - one_mapinfo[i][1] - one_mapinfo[i][3])
                else:
                    if(one_mapinfo[j][2] == -1):
                        refgap = one_mapinfo[i][1] + one_mapinfo[j][3] - non_overlap_size - one_mapinfo[j][1]
                    else:
                        refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - non_overlap_size

            else:
                bonus = one_mapinfo[i][3]
                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):
                        refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                    else:
                        refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                else:
                    if(one_mapinfo[j][2] == -1):
                        refgap = one_mapinfo[i][1] - one_mapinfo[j][1]
                    else:
                        refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - one_mapinfo[j][3]
                
            filtered = True



            while(refgap >= 0):
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):

                    break

                gapcost = abs(readgap - refgap)

                if((readgap  > maxgap) or (gapcost > maxdiff)):                 
                    break

                gapcost = gapcost_list[gapcost]

                test_scores = S[j] + bonus - gapcost

                filtered = False

                break
            if(filtered == True):
                gapcost = abs(readgap - refgap)#mark 1
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    test_scores = S[j] - min(50, skipcost) + bonus -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))
                else:    
                    test_scores = S[j] - skipcost + bonus -  min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000))

            if(test_scores > max_scores):

                max_scores = test_scores
                pre_index = j


                            


        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i

      

        
            
        
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    preitem = one_mapinfo[take_index]
    while(True):
        if((P[take_index] == -9999999)):
            break
        take_index = P[take_index]
        nowitem = one_mapinfo[take_index]
        readgap = preitem[0]
        if(preitem[0] >= (nowitem[0] + nowitem[3])):
            path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
        else:
            if(nowitem[2] == 1):
                path.append((nowitem[0], nowitem[1], nowitem[2], preitem[0] - nowitem[0]))
            else:
                path.append((nowitem[0], nowitem[1] + nowitem[3] - preitem[0] + nowitem[0], nowitem[2], preitem[0] - nowitem[0]))
                
        preitem = nowitem

    return g_max_scores, path

@njit
def drop_misplaced_alignment_test1(alignment_list, iloc, debug):
    if((alignment_list[iloc][0][2] == alignment_list[iloc + 1][0][2]) and (alignment_list[iloc][0][2] == alignment_list[iloc + 2][0][2])):


        prealignment_size = alignment_list[iloc][-1][0] + alignment_list[iloc][-1][3] - alignment_list[iloc][0][0]
        midalignment_size = alignment_list[iloc+1][-1][0] + alignment_list[iloc+1][-1][3] - alignment_list[iloc+1][0][0]
        if(midalignment_size > 1000):
            return False
        nowalignment_size = alignment_list[iloc+2][-1][0] + alignment_list[iloc+2][-1][3] - alignment_list[iloc+2][0][0]

        mid_gap = alignment_list[iloc+2][0][0] - (alignment_list[iloc][-1][0] + alignment_list[iloc][-1][3])


        preitem = alignment_list[iloc][-1]
        nowitem = alignment_list[iloc + 1][0]
        readgap = nowitem[0] - preitem[0] - preitem[3]
        if(preitem[2] == 1):
            refgap = nowitem[1] - preitem[1] - preitem[3]

        else:
            refgap = preitem[1]  - nowitem[1] - nowitem[3]

        if(abs(refgap) < 100000):
            DEL = 0
            INS = 0
            if((readgap - refgap) < -30):
                DEL += 1
            elif((readgap - refgap) > 30):
                INS += 1
            else:

                return False
            gap_1 = abs(readgap - refgap)
            preitem = alignment_list[iloc + 1][-1]
            nowitem = alignment_list[iloc + 2][0]

            readgap = nowitem[0] - preitem[0] - preitem[3]
            if(preitem[2] == 1):
                refgap = nowitem[1] - preitem[1] - preitem[3]

            else:
                refgap = preitem[1]  - nowitem[1] - nowitem[3]

            if(abs(refgap) < 100000):
                if((readgap - refgap) < -30):
                    DEL += 1
                elif((readgap - refgap) > 30):
                    INS += 1
                else:

                    return False
                gap_2 = abs(readgap - refgap)
                if(DEL == 1 and INS == 1 and ((min(gap_1, gap_2) / max(gap_1, gap_2)) > 0.7)):

                    alignment_list.pop(iloc+1)
                    return True
    return False
def extend_func(raw_alignment_list, readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, maxdiffratio, debug = False, H = False, nofilter = False, eqx = False):

    TRA_signal = False
    #step 1
    #rebuild chain break
    #remove anchor cause large cost and small alignment
    st = time.time()
    if(H == False):
        alignment_list = rebuild_chain_break(contig2start, raw_alignment_list, large_cost = setting_maxdiff, small_alignment = 50, small_dup = -100)
    else:
        alignment_list = rebuild_chain_break_H(contig2start, raw_alignment_list, large_cost = setting_maxdiff, small_alignment = 30, small_dup = -30)
        nofilter = True
    
        
    tmpiloc = -1
    while((tmpiloc + 1) < len(alignment_list)):
        tmpiloc += 1
        preitem, nowitem = alignment_list[tmpiloc][0], alignment_list[tmpiloc][-1]
        target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, setting_kmersize, contig2seq, contig2start)
        diffratio = edlib.align(query = query, target = target, task = 'distance')['editDistance']/min(len(target), len(query))
        if((diffratio>maxdiffratio)):
            alignment_list.pop(tmpiloc)
            tmpiloc -= 1
            
    #^print_log('step 1: rebuild chain break ', time.time() - st)
    if(debug == True):
        print_log('step 1: rebuild chain break')
        for line in alignment_list:
            tempcontig = pos2contig(line[0][1], contig2start)
            temprefbias = contig2start[tempcontig]
            preitem, nowitem = line[0], line[-1]

            target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, setting_kmersize, contig2seq, contig2start)

            print_log(preitem, nowitem)
            diffratio = edlib.align(query = query, target = target, task = 'distance')['editDistance']/min(len(target), len(query))
            print_log(line[0][0], line[-1][0], line[0][1] - temprefbias, line[-1][1] - temprefbias, diffratio)
            plot_result = np.array(line)
            plt.scatter(plot_result[:,0], plot_result[:,1])
        plt.show()
        
    #step 2
    #extend edge to recover small dup and misplaced alignment
    # and also to merge alignment gaped with deletion or insertion
    st = time.time()
    extend_edge_test(testseq, testseq_len, alignment_list, setting_kmersize, pos2contig, contig2start, contig2seq, san = 1, debug = debug)
    if(debug == True):print_log('step 2: extend edge ', time.time() - st)
    if(debug == True):
        print_log('After extend edge')
        for line in alignment_list:
            tempcontig = pos2contig(line[0][1], contig2start)
            temprefbias = contig2start[tempcontig]
            print_log(line[0][0], line[-1][0], line[0][1] - temprefbias, line[-1][1] - temprefbias)
        
    #step 3
    #remove miss placed alignment which cause del/ins or ins/del in start and end
    st = time.time()
    
    o_alignment_list_len = len(alignment_list)
    filtered = False
    if((len(alignment_list) > 2) and (nofilter == False)):    
        iloc = 0
        while(iloc < (len(alignment_list) - 2)):
            removed = drop_misplaced_alignment_test(alignment_list, iloc, debug = debug)
            if(removed == True):
                continue
            else:
                iloc += 1
                
    
        #extend_edge_drop_test(testseq, testseq_len, alignment_list, setting_kmersize, pos2contig, contig2start, contig2seq, san = 1, debug = debug)
        
    if(debug == True):print_log('step 3: remove miss placed alignment ', time.time() - st)
    if(debug == True):
        print_log('After remove miss placed alignment')
        for line in alignment_list:
            tempcontig = pos2contig(line[0][1], contig2start)
            temprefbias = contig2start[tempcontig]
            print_log(line[0][0], line[-1][0], line[0][1] - temprefbias, line[-1][1] - temprefbias)
            plot_result = np.array(line)
            plt.scatter(plot_result[:,0], plot_result[:,1])
        plt.show()
            

            
    
    if(len(alignment_list) <  o_alignment_list_len):#fill the gap
        filtered = True
        st = time.time()
        extend_edge_test(testseq, testseq_len, alignment_list, setting_kmersize, pos2contig, contig2start, contig2seq, san = 1, debug = debug)
        if(debug == True):print_log('step 4: fill the gap by extend edge ', time.time() - st)
        if(debug == True):
            print_log('After extend edge')
            for line in alignment_list:
                tempcontig = pos2contig(line[0][1], contig2start)
                temprefbias = contig2start[tempcontig]
                print_log(line[0][0], line[-1][0], line[0][1] - temprefbias, line[-1][1] - temprefbias)

        
    #step 4
    #merge conjacent alignment with small readgap or refgap
    st = time.time()
    merge_smallgap = 2000
    too_large_gap = 5000
    if(len(alignment_list) >= 2):
        iloc = 0
        duplist = getdupiloc(alignment_list)
        while((iloc + 1) < len(alignment_list)):
            if(iloc in duplist):
                iloc += 1
                continue
            
            preitem = alignment_list[iloc][-1]
            nowitem = alignment_list[iloc + 1][0]
            if(preitem[2] != nowitem[2] or (pos2contig(preitem[1], contig2start) != pos2contig(nowitem[1], contig2start))):
                iloc += 1
                continue
            presize = alignment_list[iloc][-1][0] + alignment_list[iloc][-1][3] - alignment_list[iloc][0][0]
            nowsize = alignment_list[iloc + 1][-1][0] + alignment_list[iloc + 1][-1][3] - alignment_list[iloc + 1][0][0]
            
            if(alignment_list[iloc][-1][2] == 1):
                ref_presize = alignment_list[iloc][-1][1] + alignment_list[iloc][-1][3] - alignment_list[iloc][0][1]
            else:
                ref_presize = alignment_list[iloc][0][1] + alignment_list[iloc][0][3] - alignment_list[iloc][-1][1]
            if(alignment_list[iloc+1][-1][2] == 1):
                ref_nowsize = alignment_list[iloc+1][-1][1] + alignment_list[iloc+1][-1][3] - alignment_list[iloc+1][0][1]
            else:
                ref_nowsize = alignment_list[iloc+1][0][1] + alignment_list[iloc+1][0][3] - alignment_list[iloc+1][-1][1]
                
            readgap = nowitem[0] - preitem[0] - preitem[3]
            if(preitem[2] == 1):
                refgap = nowitem[1] - preitem[1] - preitem[3]

            else:
                refgap = preitem[1]  - nowitem[1] - nowitem[3]                        
            if(refgap < 0):
                iloc += 1
                continue
            #if((min(readgap, refgap) < merge_smallgap) and (max(readgap, refgap) < too_large_gap)):
            #if((((min(presize, nowsize) > 500) or ((abs(readgap - refgap) / min(presize, nowsize)) < 0.5)) and (max(readgap, refgap) < 20000))):
            if((min(readgap, refgap) < 50) and (abs(readgap - refgap) < 10000)):
                alignment_list[iloc] = List_merge((alignment_list[iloc], alignment_list[iloc + 1]))
                alignment_list.pop(iloc+1)
            else:
                iloc += 1
    if(debug == True):print_log('step 4: merge conjacent alignment ', time.time() - st)
    if(debug == True):
        print_log('After merge conjacent alignment')
        for line in alignment_list:
            tempcontig = pos2contig(line[0][1], contig2start)
            temprefbias = contig2start[tempcontig]
            if(line[0][2] == 1):
                print_log(line[0][0], line[-1][0]+line[-1][3], line[0][1] - temprefbias, line[-1][1]+line[-1][3] - temprefbias)
            else:
                print_log(line[0][0], line[-1][0]+line[-1][3], line[-1][1] - temprefbias, line[0][1]+line[0][3] - temprefbias)
    if(len(alignment_list) > 2):
        iloc = 0
        while(iloc + 2 < len(alignment_list)):
            if(alignment_list[iloc][0][2] == alignment_list[iloc + 2][0][2] and alignment_list[iloc][0][2] != alignment_list[iloc + 1][0][2]):
                if(alignment_list[iloc][0][2] == 1):
                    tempcontig = pos2contig(alignment_list[iloc][0][1], contig2start)
                    temprefbias = contig2start[tempcontig]
                    
                    refen_0 = alignment_list[iloc][-1][1] + alignment_list[iloc][-1][3] - temprefbias
                    readen_0 = alignment_list[iloc][-1][0] + alignment_list[iloc][-1][3]
   
                    
                    refst_1 = alignment_list[iloc + 1][-1][1] - temprefbias
                    readst_1 = alignment_list[iloc + 1][0][0]
                    refen_1 = alignment_list[iloc + 1][0][1] + alignment_list[iloc + 1][0][3] - temprefbias
                    readen_1 = alignment_list[iloc + 1][-1][0] + alignment_list[iloc + 1][-1][3]
                    
                    refst_2 = alignment_list[iloc + 2][0][1] - temprefbias
                    readst_2 = alignment_list[iloc + 2][0][0]
                    #print(readst_1, readen_1)
                    #print(refst_2 - refen_0, refen_1 - refst_1)
                    if(refst_2 - refen_0 == refen_1 - refst_1 and readst_1 - readen_0 + readst_2 - readen_1 == 0):
                        if(refst_1 - refen_0 != 0 and refst_1 - refen_0 + refst_2 - refen_1 == 0):
                            #print( refen_0, refst_1)
                            if(refen_0 > refst_1):
                                tempref = (contig2seq[tempcontig][refen_1: refen_1 + refen_0 - refst_1])
                                tempquery = (testseq[readen_0 - refen_0 + refst_1: readen_0])
                            else:#refen_0 < refst_1
                                tempref = (contig2seq[tempcontig][refen_0: refen_0 - refen_0 + refst_1])
                                tempquery = (testseq[readen_0: readen_0 - refen_0 + refst_1])
                                if(tempref == tempquery):
                                    alignment_list[iloc][-1] = (readen_0 - refen_0 + refst_1, refen_0 - refen_0 + refst_1 + temprefbias, 1, 0)
                                    insertitem = (readen_0 - refen_0 + refst_1, refen_1 + refen_0 - refst_1 + temprefbias, -1, 0)
                                    #print((readen_0 - refen_0 + refst_1, refen_0 - refen_0 + refst_1 + temprefbias, 1, 0))
                                    #print(insertitem)

                                    while(insertitem[0]>=alignment_list[iloc + 1][0][0]):
                                        alignment_list[iloc + 1].pop(0)
                                    alignment_list[iloc + 1].insert(0, insertitem)
                                    #alignment_list[iloc][0] = (readen_0 - refen_0 + refst_1, )
                                    
                            #print(refen_0 )
                            #print(readst_1, readen_1, refst_1 , refen_1 )
                            #print(refst_2)

                                
            iloc += 1
    #step 5
    #split unrelate read
    #note inside alignment there is no way two conjacent anchors
    #has cost large than 30
    #so there is no big gap
    st = time.time()
    new_alignment_list = List()
    cigarlist = []
    for alignment in alignment_list: 

        tmp_alignment_list, tmp_cigarlist = split_alignment_test(alignment, testseq, rc_testseq, testseq_len, kmersize=setting_kmersize , contig2seq = contig2seq, contig2start = contig2start, debug = debug, H = H, eqx = eqx)
        if(debug): print_log(len(tmp_alignment_list), len(tmp_cigarlist))
        iloc = -1
        for alignment in tmp_alignment_list:
            iloc += 1
            new_alignment_list.append(alignment)
            cigarlist.append(tmp_cigarlist[iloc])


    if(debug == True):print_log('step 5: split unrelate read ', time.time() - st)
    if(debug == True):
        print_log('After split unrelate read')
        for line in new_alignment_list:
            tempcontig = pos2contig(line[0][1], contig2start)
            temprefbias = contig2start[tempcontig]
            print_log(line[0][0], line[-1][0], line[0][1] - temprefbias, line[-1][1] - temprefbias)
    
    alignment_list, onemapinfolist = get_onemapinfolist(new_alignment_list, cigarlist, readid, mapq, testseq_len, contig2start, need_reverse)
    return alignment_list, onemapinfolist, TRA_signal, filtered


#testing 20240801
extra = []
for gapcost in range(10000000000):
    extra.append(min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000)))
    if(len(extra) > 1 and extra[-1] == 36):
        break
extra = np.array(extra, dtype = np.float32)

@njit
def get_localmap_multi_all_forDP_inv_guide(raw_alignment_array, testseq, rc_testseq, contig2start, contig2seq, kmersize, skipcost, maxdiff, maxgap, shift = 1):#
    


    def seq2hashtable_multi_test(onelookuptable_s, onelookuptable_m, seq, start, kmersize):
        skiphash = hash('N'*kmersize)
        multi = [0]
        multi.pop(0)
        for iloc in range(0, len(seq) - kmersize + 1, 1):
            hashedkmer = hash(seq[iloc:iloc+kmersize])
            if(skiphash == hashedkmer):
                continue
            if(hashedkmer not in onelookuptable_s):
                onelookuptable_s[hashedkmer] = start + iloc
            else:
                if(hashedkmer in onelookuptable_m):
                    onelookuptable_m[hashedkmer].append(start + iloc)
                else:
                    onelookuptable_m[hashedkmer] = List([onelookuptable_s[hashedkmer], start + iloc])
                    multi.append(hashedkmer)
        for hashedkmer in multi:
            onelookuptable_s.pop(hashedkmer)

    

    readgap = 0
    pre = raw_alignment_array[0]
    for now in raw_alignment_array[1:]:
        if(abs(now[0] - pre[0]) > readgap):
            readgap = abs(now[0] - pre[0])
        pre = now
    readgap += 1000  
    readgap = max(readgap, 5000)
    raw_alignment_array = raw_alignment_array[np.argsort(raw_alignment_array[:, 1])]
    
    startandend = List([(raw_alignment_array[0][1], raw_alignment_array[0][1])])
    for item in raw_alignment_array[1:]:
        if(((item[1] - startandend[-1][1]) < readgap)):
            startandend[-1] = (startandend[-1][0], item[1])
        else:
            if(startandend[-1][0] == startandend[-1][1]):
                startandend.pop(-1)
            startandend.append((item[1], item[1]))
    if(startandend[-1][0] == startandend[-1][1]):
        startandend.pop(-1)
    
    local_lookuptable_s = Dict()
    local_lookuptable_s[0] = 0
    local_lookuptable_s.pop(0)
    
    local_lookuptable_m = Dict()
    local_lookuptable_m[0] = List([0])
    local_lookuptable_m.pop(0)
    
    retry_diffcontig = False
    for item in startandend:
        min_ref, max_ref = item[0], item[1]
        testcontig = pos2contig(min_ref, contig2start)
        if(testcontig != pos2contig(max_ref, contig2start)):
            retry_diffcontig = True
            break
        lookfurther = min(2000, min_ref-contig2start[testcontig])
        min_ref -= lookfurther
        max_ref += 2000
        refseq = contig2seq[testcontig][min_ref-contig2start[testcontig]: max_ref-contig2start[testcontig]]
        seq2hashtable_multi_test(local_lookuptable_s, local_lookuptable_m, refseq, min_ref, kmersize)

    if(retry_diffcontig == True):
        startandend = List([(raw_alignment_array[0][1], raw_alignment_array[0][1])])
        current_contig = pos2contig(raw_alignment_array[0][1], contig2start)
        for item in raw_alignment_array[1:]:
            if(((item[1] - startandend[-1][1]) < readgap) and (current_contig == pos2contig(item[1], contig2start))):
                startandend[-1] = (startandend[-1][0], item[1])
            else:
                if(startandend[-1][0] == startandend[-1][1]):
                    startandend.pop(-1)
                startandend.append((item[1], item[1]))
                current_contig = pos2contig(item[1], contig2start)
        if(startandend[-1][0] == startandend[-1][1]):
            startandend.pop(-1)
            
        local_lookuptable_s = Dict()
        local_lookuptable_s[0] = 0
        local_lookuptable_s.pop(0)

        local_lookuptable_m = Dict()
        local_lookuptable_m[0] = List([0])
        local_lookuptable_m.pop(0)

        retry_diffcontig = False
        for item in startandend:
            min_ref, max_ref = item[0], item[1]
            testcontig = pos2contig(min_ref, contig2start)
            if(testcontig != pos2contig(max_ref, contig2start)):
                retry_diffcontig = True
                break
            lookfurther = min(2000, min_ref-contig2start[testcontig])
            min_ref -= lookfurther
            max_ref += 2000
            refseq = contig2seq[testcontig][min_ref-contig2start[testcontig]: max_ref-contig2start[testcontig]]
            seq2hashtable_multi_test(local_lookuptable_s, local_lookuptable_m, refseq, min_ref, kmersize)

        
    raw_alignment_array = raw_alignment_array[np.argsort(raw_alignment_array[:, 0])]
    readstart = max(0, raw_alignment_array[0][0]-500)
    readend = min(len(testseq)-kmersize+1, raw_alignment_array[-1][0]+500)

    
    one_mapinfo = [(-1, -1, -1, -1)]
    one_mapinfo.pop(0)
    
    iloc = readstart
    iloc -= shift
    

    pointdict = Dict()
    pointdict[21312] = (1, 1, 1, 1)
    pointdict.pop(21312)
    
    pointdict_key = [21312]
    pointdict_key.pop(0)
    
    while(True):
   
        iloc += shift
        if(iloc >= readend):
            break

        if(hash(testseq[iloc: iloc + kmersize]) == hash(rc_testseq[-(iloc + kmersize): -iloc])):
            continue

        biasvalue, biasvalue_1, closest_index, closest_index_1 = findClosest(raw_alignment_array, target = iloc)


        interval = min(biasvalue + biasvalue_1 + 500, 2000)
        upperrefloc =  (raw_alignment_array[closest_index][1] + interval, raw_alignment_array[closest_index_1][1] + interval)
        lowerrefloc =  (raw_alignment_array[closest_index][1] - interval, raw_alignment_array[closest_index_1][1] - interval)
        readgap = abs(iloc - raw_alignment_array[closest_index][0])
   
        hashedkmer = hash(testseq[iloc: iloc + kmersize]) 
        if(hashedkmer in local_lookuptable_s):
            refloc = local_lookuptable_s[hashedkmer]
            refgap = abs(refloc - raw_alignment_array[closest_index][1])
            diff = abs(readgap - refgap)
            if((diff < 500) or (upperrefloc[0] >= refloc and lowerrefloc[0] <= refloc) or (upperrefloc[1] >= refloc and lowerrefloc[1] <= refloc)):
                item = (iloc, refloc, 1)
                point = item[1] - item[0]

                if(point in pointdict):

                    if((pointdict[point][0] + pointdict[point][3]) >= item[0]):
                        bouns = item[0] - (pointdict[point][0] + pointdict[point][3]) + kmersize
                        if(bouns > 0):
                            if((pointdict[point][3] + bouns < 20)):
                                pointdict[point] = (pointdict[point][0], pointdict[point][1], 1, pointdict[point][3] + bouns)
                            else:
                                one_mapinfo.append(pointdict[point])
                                pointdict[point] = (pointdict[point][0] + pointdict[point][3], pointdict[point][1] + pointdict[point][3], 1, bouns)



                    else:    
                        one_mapinfo.append(pointdict[point])
                        pointdict[point] = (item[0], item[1], item[2], kmersize)
                else:
                    pointdict[point] = (item[0], item[1], item[2], kmersize)
                    pointdict_key.append(point)
                    
        elif(hashedkmer in local_lookuptable_m):
            for refloc in local_lookuptable_m[hashedkmer]:
                refgap = abs(refloc - raw_alignment_array[closest_index][1])
                diff = abs(readgap - refgap)
                if((diff < 500) or (upperrefloc[0] >= refloc and lowerrefloc[0] <= refloc) or (upperrefloc[1] >= refloc and lowerrefloc[1] <= refloc)):
                    item = (iloc, refloc, 1)
                    point = item[1] - item[0]

                    if(point in pointdict):

                        if((pointdict[point][0] + pointdict[point][3]) >= item[0]):
                            bouns = item[0] - (pointdict[point][0] + pointdict[point][3]) + kmersize
                            if(bouns > 0):
                                if((pointdict[point][3] + bouns < 20)):
                                    pointdict[point] = (pointdict[point][0], pointdict[point][1], 1, pointdict[point][3] + bouns)
                                else:
                                    one_mapinfo.append(pointdict[point])
                                    pointdict[point] = (pointdict[point][0] + pointdict[point][3], pointdict[point][1] + pointdict[point][3], 1, bouns)



                        else:    
                            one_mapinfo.append(pointdict[point])
                            pointdict[point] = (item[0], item[1], item[2], kmersize)
                    else:
                        pointdict[point] = (item[0], item[1], item[2], kmersize)
                        pointdict_key.append(point)



        hashedkmer = hash(rc_testseq[-(iloc + kmersize): -iloc]) 
        if(hashedkmer in local_lookuptable_s):
            refloc = local_lookuptable_s[hashedkmer]
            refgap = abs(refloc - raw_alignment_array[closest_index][1])
            diff = abs(readgap - refgap)
            if((diff < 500) or (upperrefloc[0] >= refloc and lowerrefloc[0] <= refloc) or (upperrefloc[1] >= refloc and lowerrefloc[1] <= refloc)):
                item = (iloc, refloc, -1)
                point = -(item[1] + item[0])
   
                if(point in pointdict):
                    
                    if((pointdict[point][0] + pointdict[point][3]) >= item[0]):
                        bouns = item[0] - (pointdict[point][0] + pointdict[point][3]) + kmersize
                        if(bouns > 0):
                            if((pointdict[point][3] + bouns < 20)):
                                pointdict[point] = (pointdict[point][0], item[1], -1, pointdict[point][3] + bouns)
                            else:
                                one_mapinfo.append(pointdict[point])
                                pointdict[point] = (pointdict[point][0] + pointdict[point][3], item[1], -1, bouns)
                                    #print_log(pointdict[point][-1])


                    else:    
                        one_mapinfo.append(pointdict[point])
                        pointdict[point] = (item[0], item[1], item[2], kmersize)
                else:
                    pointdict[point] = (item[0], item[1], item[2], kmersize)
                    pointdict_key.append(point)
                    
        elif(hashedkmer in local_lookuptable_m):

            for refloc in local_lookuptable_m[hashedkmer]:
                refgap = abs(refloc - raw_alignment_array[closest_index][1])
                diff = abs(readgap - refgap)
                if((diff < 500) or (upperrefloc[0] >= refloc and lowerrefloc[0] <= refloc) or (upperrefloc[1] >= refloc and lowerrefloc[1] <= refloc)):
                    item = (iloc, refloc, -1)
                    point = -(item[1] + item[0])
   
                    if(point in pointdict):

                        if((pointdict[point][0] + pointdict[point][3]) >= item[0]):
                            bouns = item[0] - (pointdict[point][0] + pointdict[point][3]) + kmersize
                            if(bouns > 0):
                                if((pointdict[point][3] + bouns < 20)):
                                    pointdict[point] = (pointdict[point][0], item[1], -1, pointdict[point][3] + bouns)
                                else:
                                    one_mapinfo.append(pointdict[point])
                                    pointdict[point] = (pointdict[point][0] + pointdict[point][3], item[1], -1, bouns)
                                        #print_log(pointdict[point][-1])


                        else:    
                            one_mapinfo.append(pointdict[point])
                            pointdict[point] = (item[0], item[1], item[2], kmersize)
                    else:
                        pointdict[point] = (item[0], item[1], item[2], kmersize)
                        pointdict_key.append(point)


    for key in pointdict_key:
        one_mapinfo.append(pointdict[key])
    one_mapinfo = np.array(one_mapinfo)
    one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:, 0])]
    #print(one_mapinfo.shape)
  


    return get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo, kmersize = kmersize, skipcost = skipcost, maxdiff = maxdiff, maxgap = maxgap)





@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_fast_all(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500, fast_t = 5):

    extra_size = len(extra) - 1
    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20

    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
            

    

    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    S_i = np.empty(n, np.int64)
    P = np.empty(n, np.int64)

    

    

    
    

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    readlength = one_mapinfo[-1][0] + 1000
    target_arr = np.zeros(n, dtype = np.float64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
        if(one_mapinfo[i][2] == 1):
            target_arr[i] = one_mapinfo[i][1] - one_mapinfo[i][0] + readlength 
        else:
            target_arr[i] = -(one_mapinfo[i][1] + one_mapinfo[i][0] + readlength) 
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)
    testspace_en = 1
    testspace_en_i = 1
    
    
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    S_arg_i = np.empty(n, np.int64)
    S_arg_i[0] = 0
    
    i = 0

    S[i] = one_mapinfo[i][3]
    S_i[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    
    cache_dict = dict()
    cache_dict[0] = 0
    cache_dict.pop(0)
    
    fast = False
    
    
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])

        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        

        
        
        if(prereadloc < one_mapinfo[i][0]):                
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i
                
            k = testspace_en_i
            while(k < i):

                loc_in_sorted_S = smallorequal2target_1d_point_target(S_i, S_i[k], k, S_arg_i, target_arr[k], target_arr) + 1

                S_arg_i[loc_in_sorted_S + 1: k + 1] = S_arg_i[loc_in_sorted_S: k]
                S_arg_i[loc_in_sorted_S] = k
                
                k += 1

            testspace_en_i = k
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
            
            
            fast = False
            if(testspace_en_i > 0):
                temp_endpos = (testspace_en_i - max(testspace_en_i - 1000, 0))

                scorerange = S_i[S_arg_i[testspace_en_i - 1]] - S_i[S_arg_i[max(testspace_en_i - 1000, 0)]] + 1
                
                if(temp_endpos / scorerange > fast_t):
                    fast = True
                    cache_dict = dict()
                    cache_dict[0] = 0
                    cache_dict.pop(0)
            
            
               
                      
        if(fast == True):

            #print_log(i, prereadloc)

            st_loc = testspace_en_i
            en_loc = testspace_en_i

            while(en_loc > 0):
                if(st_loc not in cache_dict):
                    pre_st_loc = st_loc
                    st_loc = smallorequal2target_1d_point(S_i, S_i[S_arg_i[st_loc-1]]-1e-7, en_loc, S_arg_i) + 1
                    cache_dict[pre_st_loc] = st_loc
                else:
                    st_loc = cache_dict[st_loc]

                if(S_i[S_arg_i[st_loc]] < round(max_scores - one_mapinfo[i][3])):

                    break

                S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg_i, one_mapinfo[i][0], one_mapinfo[:,0])




                j = S_arg_i[S_arg_j]
    


                #########
                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                if((readgap < 0)):
                    bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                    readgap = 0
                    non_overlap_size = one_mapinfo[i][0] - one_mapinfo[j][0]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - non_overlap_size)
                        else:
                            refgap = (one_mapinfo[j][1] + one_mapinfo[j][3] - non_overlap_size - one_mapinfo[i][1] - one_mapinfo[i][3])
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] + one_mapinfo[j][3] - non_overlap_size - one_mapinfo[j][1]
                        else:
                            refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - non_overlap_size

                else:
                    bonus = one_mapinfo[i][3]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                        else:
                            refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] - one_mapinfo[j][1]
                        else:
                            refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - one_mapinfo[j][3]

                gapcost = abs(readgap - refgap)
            

                if(one_mapinfo[i][2] == one_mapinfo[j][2] and refgap >= 0 and readgap <= maxgap and gapcost <= maxdiff):


                    test_scores = S[j] + bonus - gapcost_list[gapcost]


                else:
                    if(gapcost > extra_size):
                        gapcost = extra_size
                    test_scores = S[j] - skipcost + bonus -  extra[gapcost]
                if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j
                    


                en_loc = st_loc
        
        else:
            
            for j in S_arg[:testspace_en][::-1]:
            
            

                if(S[j] < (max_scores - one_mapinfo[i][3])):

                    break


                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                if((readgap < 0)):
                    bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                    readgap = 0
                    non_overlap_size = one_mapinfo[i][0] - one_mapinfo[j][0]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - non_overlap_size)
                        else:
                            refgap = (one_mapinfo[j][1] + one_mapinfo[j][3] - non_overlap_size - one_mapinfo[i][1] - one_mapinfo[i][3])
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] + one_mapinfo[j][3] - non_overlap_size - one_mapinfo[j][1]
                        else:
                            refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - non_overlap_size

                else:
                    bonus = one_mapinfo[i][3]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                        else:
                            refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] - one_mapinfo[j][1]
                        else:
                            refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - one_mapinfo[j][3]

                gapcost = abs(readgap - refgap)
            

                if(one_mapinfo[i][2] == one_mapinfo[j][2] and refgap >= 0 and readgap <= maxgap and gapcost <= maxdiff):


                    test_scores = S[j] + bonus - gapcost_list[gapcost]


                else:
                    if(gapcost > extra_size):
                        gapcost = extra_size
                    test_scores = S[j] - skipcost + bonus -  extra[gapcost]
                if(test_scores > max_scores):
                    max_scores = test_scores
                    pre_index = j

        S[i] = max_scores
        S_i[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i

        

    k = testspace_en
    while(k < n):

        loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

        S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
        S_arg[loc_in_sorted_S] = k

        k += 1

    
    return g_max_index, S, P, S_arg

@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_all(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    extra_size = len(extra) - 1
    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20

    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
    



    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)



    

    

    opcount = 0
    
    

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)

    testspace_en = 1
    
    
    test_S = np.empty(0, np.float64)
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0

    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    

    

    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        


        max_scores = one_mapinfo[i][3]

        
        pre_index = -9999999
        

        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            if((opcount/one_mapinfo[i][0] > 1000) and one_mapinfo[i][0] > 500):#mark 1
                return -1, S, P, S_arg, -1.#mark 1
            
            
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i
            
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
               
                      


        
        for j in S_arg[:testspace_en][::-1]:
            
            
            opcount += 1     #mark 1
            
            if(S[j] < (max_scores - one_mapinfo[i][3])):
                
                break


            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            if((readgap < 0)):
                bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                readgap = 0
                non_overlap_size = one_mapinfo[i][0] - one_mapinfo[j][0]
                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):
                        refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - non_overlap_size)
                    else:
                        refgap = (one_mapinfo[j][1] + one_mapinfo[j][3] - non_overlap_size - one_mapinfo[i][1] - one_mapinfo[i][3])
                else:
                    if(one_mapinfo[j][2] == -1):
                        refgap = one_mapinfo[i][1] + one_mapinfo[j][3] - non_overlap_size - one_mapinfo[j][1]
                    else:
                        refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - non_overlap_size

            else:
                bonus = one_mapinfo[i][3]
                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):
                        refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                    else:
                        refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                else:
                    if(one_mapinfo[j][2] == -1):
                        refgap = one_mapinfo[i][1] - one_mapinfo[j][1]
                    else:
                        refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - one_mapinfo[j][3]
                
            gapcost = abs(readgap - refgap)
            

            if(one_mapinfo[i][2] == one_mapinfo[j][2] and refgap >= 0 and readgap <= maxgap and gapcost <= maxdiff):


                test_scores = S[j] + bonus - gapcost_list[gapcost]


            else:
                if(gapcost > extra_size):
                    gapcost = extra_size
                test_scores = S[j] - skipcost + bonus -  extra[gapcost]


            if(test_scores > max_scores):
                max_scores = test_scores
                pre_index = j
                

        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
       


        
    k = testspace_en
    while(k < n):

        loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

        S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
        S_arg[loc_in_sorted_S] = k

        k += 1

    testspace_en = k  
    

    
    return g_max_index, S, P, S_arg, opcount/len(one_mapinfo)
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    extra_size = len(extra) - 1
    
    oskipcost = skipcost
    omaxdiff = maxdiff
    repeat_weight = 20

    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):
        if(gapcost <= 10):
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))
        else:
            gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 2 * np.log2(gapcost))
    

    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)

    

    opcount = 0
    
    

    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
    testspace = np.empty(0, np.int64)

    testspace_en = 1
    

    
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0

    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i

    c_repeat_weight = 0

    
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        

        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        

        
        
        if(prereadloc < one_mapinfo[i][0]):
            
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i            
            c_repeat_weight = coverage_dict[one_mapinfo[i][0]]
            skipcost = oskipcost + c_repeat_weight
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
               
                      


        

        for j in S_arg[:testspace_en][::-1]:
            
              
            
            if(S[j] < (max_scores - one_mapinfo[i][3])):
                
                break


            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            if((readgap < 0)):
                bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                readgap = 0
                non_overlap_size = one_mapinfo[i][0] - one_mapinfo[j][0]
                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):
                        refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - non_overlap_size)
                    else:
                        refgap = (one_mapinfo[j][1] + one_mapinfo[j][3] - non_overlap_size - one_mapinfo[i][1] - one_mapinfo[i][3])
                else:
                    if(one_mapinfo[j][2] == -1):
                        refgap = one_mapinfo[i][1] + one_mapinfo[j][3] - non_overlap_size - one_mapinfo[j][1]
                    else:
                        refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - non_overlap_size

            else:
                bonus = one_mapinfo[i][3]
                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):
                        refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                    else:
                        refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                else:
                    if(one_mapinfo[j][2] == -1):
                        refgap = one_mapinfo[i][1] - one_mapinfo[j][1]
                    else:
                        refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - one_mapinfo[j][1] - one_mapinfo[j][3]
                


            gapcost = abs(readgap - refgap)
            

            if(one_mapinfo[i][2] == one_mapinfo[j][2] and refgap >= 0 and readgap <= maxgap and gapcost <= maxdiff):


                test_scores = S[j] + bonus - gapcost_list[gapcost]


            else:
                if(gapcost > extra_size):
                    gapcost = extra_size
                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    test_scores = S[j] - min(50, skipcost) + bonus -  extra[gapcost]
                else:
                    test_scores = S[j] - skipcost + bonus -  extra[gapcost]

            if(test_scores > max_scores):

                max_scores = test_scores
                pre_index = j


                            


        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i

      

        
            
        
    path = []
    take_index = g_max_index
    path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1] , one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
    preitem = one_mapinfo[take_index]
    while(True):
        if((P[take_index] == -9999999)):
            break
        take_index = P[take_index]
        nowitem = one_mapinfo[take_index]
        readgap = preitem[0]
        if(preitem[0] >= (nowitem[0] + nowitem[3])):
            path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
        else:
            if(nowitem[2] == 1):
                path.append((nowitem[0], nowitem[1], nowitem[2], preitem[0] - nowitem[0]))
            else:
                path.append((nowitem[0], nowitem[1] + nowitem[3] - preitem[0] + nowitem[0], nowitem[2], preitem[0] - nowitem[0]))
                
        preitem = nowitem

    return g_max_scores, path
            
def decode_hit(index_object, index2contig, testseq, testseq_len, contig2start, kmersize, contig2seq, skipcost = (50., 30.), maxdiff = (50, 30), maxgap = 200, check_num = 20, c_bias = 5000, bin_size = 100, overlapprecentage = 0.5, hastra = True, H = False, mid_occ = -1):            
    #redo_ratio = 5#hifi
    redo_ratio = 10#clr, ont
    st = time.time()
    need_reverse, one_mapinfo = get_reversed_chain_numpy_rough(np.array(index_object.map(testseq, check_num = check_num,mid_occ = mid_occ)), testseq_len)
    if(len(one_mapinfo) <= 2):
        return 0, 0., [(0, 0, 0, 0)], -1
    #print_log(need_reverse)
    #print_log('index_object.map', time.time() - st)
    maxgap = 1000
    #print_log('len(one_mapinfo)', len(one_mapinfo))
    st = time.time()
    path_list, primary_index_List, primary_scores_List, all_index_List, mapq, scores_list, factor =  hit2work_1(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage, hastra, H = H)
    if(len(path_list) == 0):
        scores, path = 0., [(0, 0, 0, 0)]
        return 0, scores, path, factor
    #print_log(scores_list)
    #print_log('hit2work_1', time.time() - st)
    #print_log('mapq', mapq)
    st = time.time()
    if(len(primary_index_List) != 0):
        base_iloc = primary_index_List[0]
        
        if(mapq == 0):
            #print_log('len(primary_scores_List)', len(primary_scores_List))
            base_score = primary_scores_List[0][0]
            rc_testseq = str(Seq(testseq).reverse_complement())
            min_diff = 10
            for tmpiloc in range(len(primary_scores_List[0])):
                if(primary_scores_List[0][tmpiloc] / base_score < 0.999):
                    break
                #print_log(all_index_List[0][tmpiloc])
                #print_log(primary_scores_List[0][tmpiloc])
                #print_log(path_list[all_index_List[0][tmpiloc]][-1], ',', path_list[all_index_List[0][tmpiloc]][0])
                #preitem, nowitem = path_list[all_index_List[0][tmpiloc]][-1], path_list[all_index_List[0][tmpiloc]][0]
                preitem, nowitem = return_main_alignment_size(contig2start, np.array(path_list[all_index_List[0][tmpiloc]][::-1]))
                #print_log(preitem, nowitem)
                if(preitem[2] != nowitem[2]):
                    continue
                if(need_reverse == False):
                    target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, kmersize, contig2seq, contig2start)
                else:
                    target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, rc_testseq, testseq, testseq_len, kmersize, contig2seq, contig2start)
                if(min(len(target), len(query)) == 0):
                    print_log(preitem, nowitem)
                    print_log(testseq)
                diff = edlib.align(query = query, target = target, task = 'distance')['editDistance']/min(len(target), len(query))
                #print_log(diff)
                if(diff <= min_diff):
                    min_diff = diff
                    #print_log('change', base_iloc, '->', all_index_List[0][tmpiloc])
                    base_iloc = all_index_List[0][tmpiloc]
                    primary_index_List[0] = base_iloc
 
                    
                #print_log()
        
        scores, path = scores_list[base_iloc], path_list[base_iloc]
    else:
        scores, path = 0., [(0, 0, 0, 0)]
    #print_log('edlib', time.time() - st)
    if(need_reverse == True):
        return mapq, -scores, path, factor
    else:
        return mapq, scores, path, factor
@njit
def hit2work_1(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage = 0.5, hastra = False, H = False):
    def getsecond(x):
        return x[1]
    def getfirst(x):
        return x[0]
    def getlength(x):
        return len(x)

    def get_readloc_set_bin(one_mappos, bin_size):
        return set([i[0] // bin_size for i in one_mappos])
    def get_refloc_set_bin(one_mappos, bin_size):
        return set([i[1] // bin_size for i in one_mappos])
    def get_overlapsize(readloc_set_a, readloc_set_b):
        return len(readloc_set_a&readloc_set_b)/min(len(readloc_set_a), len(readloc_set_b))
    #print_log(len(one_mapinfo))


    
    #one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:,1])] 




    
    hit = False
    minichain_scores = 40
    path_list = List()
    scores_list = []
    scores_list.append(0.)
    scores_list.pop()
    path_list.append([(0, 0, 0, 0)])
    path_list.pop()
    

    max_scores = 0

    
    fast_enable = False
    g_max_index = 0
    factor = -1.



    if(len(one_mapinfo) / testseq_len > 5):
        fast_enable = True
    one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:,0])]
    if(fast_enable == False):
        g_max_index, S, P, S_arg, factor = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_all(one_mapinfo, kmersize = kmersize, skipcost = skipcost[0], maxdiff = maxdiff[0], maxgap = maxgap)


    if(fast_enable == True or g_max_index == -1):
        fast_enable = True
        g_max_index, S, P, S_arg = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_fast_all(one_mapinfo, kmersize = kmersize, skipcost = skipcost[0], maxdiff = maxdiff[0], maxgap = maxgap)

    scores = S[g_max_index]
    #if(scores < minichain_scores):
        #if(hit == False):
            #return 0, 0., [(0, 0, 0, 0)]
        #continue

    ##################
    usedindex = set()

    path = [(0, 0, 0, 0)]
    path.pop()
    take_index = g_max_index
    usedindex.add(take_index)
    score = S[take_index]
    while(True):
        if((P[take_index] == -9999999)):
            path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
            break
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
        take_index = P[take_index]
        usedindex.add(take_index)

    if(score > 40):
        hit = True
        scores_list.append(score)
        path_list.append(path)
    if(scores > max_scores):

        max_scores = scores
    #print_log(score, len(path))

    for take_index in S_arg[::-1]:
        if(take_index in usedindex):
            continue
        path = [(0, 0, 0, 0)]
        path.pop()
        usedindex.add(take_index)
        score = S[take_index]
        while(True):
            if((P[take_index] == -9999999)):
                path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
                break
            path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
            take_index = P[take_index]
            if(take_index in usedindex):
                score = score - S[take_index]
                break
            usedindex.add(take_index)
        if(score > 40):
            scores_list.append(score)
            path_list.append(path)
        if(scores > max_scores):
            max_scores = scores
    ##################
        
        

        

    #print_log(max_scores)
    #print_log()
    if(hit == True and max_scores > 50):#hifi
    #if(hit == True and max_scores > 60):#clr, ont

        order = np.argsort(np.array(scores_list))[::-1]
        #print_log('scores_list[order[0]], len(path_list[order[0]])', scores_list[order[0]], len(path_list[order[0]]))
        primary_rlocset_List = List()
        
        primary_scores_List = List()

        
        primary_index_List = List()

        
        all_index_List = List()

        
        iloc = order[0]
        primary_rlocset_List.append(get_readloc_set_bin(path_list[iloc], bin_size))
        primary_scores_List.append(List([scores_list[iloc]]))
        primary_index_List.append(iloc)
        all_index_List.append(List([iloc]))
        #print_log(scores_list[order[0]], path_list[order[0]][-1][1], path_list[order[0]][0][1])
        for iloc in order[1:]:
            
            #print_log(scores_list[iloc], path_list[iloc][-1][1], path_list[iloc][0][1])
            readloc_set_b = get_readloc_set_bin(path_list[iloc], bin_size)
            maxoverlapsize = 0.
            for p_loc in range(len(primary_rlocset_List)):
                tmp_overlapsize = get_overlapsize(primary_rlocset_List[p_loc], readloc_set_b)
                if(tmp_overlapsize > maxoverlapsize):
                    maxoverlapsize = tmp_overlapsize
                    prefer_p_loc = p_loc
            if(maxoverlapsize < overlapprecentage):
                primary_rlocset_List.append(readloc_set_b)
                primary_scores_List.append(List([scores_list[iloc]]))
                all_index_List.append(List([iloc]))
                primary_index_List.append(iloc)
            else:
                primary_scores_List[prefer_p_loc].append(scores_list[iloc])
                all_index_List[prefer_p_loc].append(iloc)
        
        
        m = len(path_list[order[0]])    
        if(len(primary_scores_List[0]) < 2):
            f1 = primary_scores_List[0][0]
            f2 = 0
        else:
            f1 = primary_scores_List[0][0]
            f2 = primary_scores_List[0][1]
        mapq = min(int(40*(1-f2/f1)*min(1, m/10)*np.log(f1)), 60)
        
        return path_list, primary_index_List, primary_scores_List, all_index_List, mapq, scores_list, factor



    else:
        path_list = List()
        scores_list = []
        scores_list.append(0.)
        scores_list.pop()
        path_list.append([(0, 0, 0, 0)])
        path_list.pop()
        
        primary_scores_List = List()
        primary_scores_List.append(List([0.]))
        primary_scores_List.pop(0)
        
        primary_index_List = List()
        primary_index_List.append(0)
        primary_index_List.pop(0)
        
        all_index_List = List()
        all_index_List.append(List([0]))
        all_index_List.pop(0)
        return path_list, primary_index_List, primary_scores_List, all_index_List, 0, scores_list, factor
