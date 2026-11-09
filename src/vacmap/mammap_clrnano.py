# VACmap - a long-read aligner for structural variation discovery
# Copyright (C) 2023 Hongyu Ding
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
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
from vacmap.output_functions import nm_from_cigar
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
def sort_by_length(x):
    return x[4] - x[3]
def get_overlapsize(readloc_set_a, readloc_set_b):
    return len(readloc_set_a&readloc_set_b)/min(len(readloc_set_a), len(readloc_set_b))


@njit
def get_length(x):
    return len(x)

@njit
def get_readloc_set_bin(one_mappos, bin_size):
    readloc_set = set()
    for item in one_mappos:
        readloc_set.add(item[0]//bin_size)
    return readloc_set


@njit
def seq2hashtable_multi_test(onelookuptable, seq, start, kmersize):

    for iloc in range(0, len(seq) - kmersize + 1, 1):
        hashedkmer = hash(seq[iloc:iloc+kmersize])
        if(hashedkmer in onelookuptable):

            onelookuptable[hashedkmer].append(start + iloc)
        else:
            onelookuptable[hashedkmer] = List([start + iloc])


@njit
def List_merge(List_of_List):
    merged_List = List_of_List[0][:]
    for one_List in List_of_List[1:]:
        for item in one_List:
            merged_List.append(item)
    return merged_List

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
        

##############
#1117mapq
########################
#1212
#print_log(1212, 'coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, 20)')


#print_log(1216)

#print_log(1218)


#print_log(1224)


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
                

#########

@njit
def string2int(x):
    num = 0
    zero = ord('0')
    for c in x:
        num = num * 10 + ord(c) - zero    
    return num
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


########


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


#20240127


#2024312

#20240313

#20240401

#20240415


#20240602
    #f_redo_ratio_list = np.array(f_redo_ratio_list)
    #info = 'median: '+ str(np.median(f_redo_ratio_list)) + ' avg: ' + str(f_redo_ratio_list.mean())
    #logging.info(info)
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


#20240629


#RG

    #f_redo_ratio_list = np.array(f_redo_ratio_list)
    #info = 'median: '+ str(np.median(f_redo_ratio_list)) + ' avg: ' + str(f_redo_ratio_list.mean())
    #logging.info(info)
    
#testing 20240801
extra = []
for gapcost in range(10000000000):
    extra.append(min(36, 30 + 0.5 * np.log(max(gapcost, 1)), min(10, gapcost/100)+min(30, gapcost/1000)))
    if(len(extra) > 1 and extra[-1] == 36):
        break
extra = np.array(extra, dtype = np.float32)


readgapcost_list = np.zeros(100, dtype = np.float32)
for readgapcost in range(1, 100):
    readgapcost_list[readgapcost] = 0.1 * np.log2(readgapcost)

            
@njit
def getdupiloc_numba(alignment_list):

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
def merge_conjacent_alignment(alignment_list, contig2start):
    
    if(len(alignment_list) >= 2):
        merge_smallgap = 2000
        too_large_gap = 5000
        iloc = 0
        duplist = getdupiloc_numba(alignment_list)
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
#0822
        
        
#09032024
@njit
def insertpoint_score_distance(S_i, target_score, k, S_arg_i, target_distance, distance_arr):
    i = 0
    j = k
    if(S_i[S_arg_i[0]] > target_score):
        return 0
    if(S_i[S_arg_i[k - 1]] < target_score):
        return k
    while(i < j):
        
        mid = (i + j) // 2

        now_score = S_i[S_arg_i[mid]]
        if(now_score < target_score):
            i = mid + 1
        elif(now_score > target_score):
            j = mid
        else:

            now_distance = distance_arr[S_arg_i[mid]]

            if(now_distance < target_distance):
                i = mid + 1
            elif(now_distance > target_distance):
                j = mid
            else:
                return mid + 1
    return j
@njit
def closest2targetdistance(target_distance, distance_arr, S_arg_i, st_loc, en_loc):
    i = st_loc
    j = en_loc
    if(distance_arr[S_arg_i[i]] >= target_distance):
        return i
    if(distance_arr[S_arg_i[j - 1]] <= target_distance):
        return j - 1
    while(i < j):
        
        mid = (i + j) // 2

        now_distance = distance_arr[S_arg_i[mid]]
        if(now_distance < target_distance):
            i = mid + 1
        elif(now_distance > target_distance):
            j = mid
        else:
            return mid
    #if(((target_distance - distance_arr[S_arg_i[j - 1]]) < 0) or ((distance_arr[S_arg_i[j]] - target_distance) < 0)):
        #print('error: distance < 0')
    if((target_distance - distance_arr[S_arg_i[j - 1]]) < (distance_arr[S_arg_i[j]] - target_distance)):
        return j - 1
    else:
        return j
    

@njit
def findClosest_1(arr, target):
    n = len(arr)
    if(target <= arr[0]):
        return arr[0] - target, arr[0] - target, 0, 0
    if(target >= arr[n - 1]):
        return target - arr[n - 1], target - arr[n - 1], n-1, n-1

    i = 0
    j = n
    mid = 0
    while(i < j):
        mid = (i + j) // 2

        if(arr[mid] == target):
            return abs(arr[mid] - target), abs(arr[mid] - target), mid, mid

        if(target < arr[mid]):
            j = mid
        else:
            i = mid + 1

    return abs(arr[j - 1] - target), abs(arr[j] - target), j-1, j


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
        elif(info[iloc] == 'S' or info[iloc] == 'H'):
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
        elif(info[iloc] == 'S' or info[iloc] == 'H'):
            continue
        else:
            return '', ''
        preop = info[iloc]
    if(equal_value > 0):
        mdstring += str(equal_value)
    return mdstring, csstring
@njit
def expand_M_ops(info, target, query):
    """Expand 'M' runs into '='/ 'X' runs by base comparison.

    get_MD_CSshort/get_MD_CSlong only understand '=', 'X', 'D' and 'I'.  A
    CIGAR that uses 'M' (VACmap's default unless eqx is enabled) fell into
    their `else: return '', ''` branch, so every alignment shipped an empty
    MD:Z / cs:Z tag.  Expanding here keeps one MD implementation and leaves
    the CIGAR handed to the caller untouched.
    """
    out = []
    refloc = 0
    readloc = 0
    for iloc in range(1, len(info), 2):
        op = info[iloc]
        value = string2int(info[iloc - 1])
        if(op == 'M'):
            run = 0
            runop = ''
            for jloc in range(value):
                cur = '='
                if(target[refloc + jloc] != query[readloc + jloc]):
                    cur = 'X'
                if(cur == runop):
                    run += 1
                else:
                    if(run > 0):
                        out.append(str(run))
                        out.append(runop)
                    run = 1
                    runop = cur
            if(run > 0):
                out.append(str(run))
                out.append(runop)
            refloc += value
            readloc += value
        else:
            out.append(info[iloc - 1])
            out.append(op)
            if(op == '=' or op == 'X' or op == 'D'):
                refloc += value
            if(op == '=' or op == 'X' or op == 'I'):
                readloc += value
    return out


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

    oplist_md = expand_M_ops(oplist, target, query)
    if(shortcs == True):
        mdstring, csstring = get_MD_CSshort(oplist_md, target, query)
    else:
        mdstring, csstring = get_MD_CSlong(oplist_md, target, query)
    return oplist, mdstring, csstring

def mergecigar_md_cs(cigarstring, target, query, shortcs = True):  

    oplist, mdstring, csstring = mergecigar_md_(cigarstring, target, query, shortcs)

    n_cigar = len(oplist)
    return ''.join(oplist), mdstring, csstring, n_cigar
def extend_func(raw_alignment_list, readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, maxdiffratio, option, debug = False, H = False, nofilter = False, eqx = False):

    TRA_signal = False

    alignment_list = rebuild_chain_break(contig2start, raw_alignment_list, large_cost = setting_maxdiff, small_alignment = 50, small_dup = -100)

    
    tmpiloc = -1
    while((tmpiloc + 1) < len(alignment_list)):
        tmpiloc += 1
        preitem, nowitem = alignment_list[tmpiloc][0], alignment_list[tmpiloc][-1]
        target, query, target_st, target_en, query_st, query_en = get_query_target_for_cigar(preitem, nowitem, testseq, rc_testseq, testseq_len, setting_kmersize, contig2seq, contig2start)
        diffratio = edlib.align(query = query, target = target, task = 'distance')['editDistance']/min(len(target), len(query))
        if((diffratio>maxdiffratio)):
            alignment_list.pop(tmpiloc)
            tmpiloc -= 1
            

    extend_edge_test(testseq, testseq_len, alignment_list, setting_kmersize, pos2contig, contig2start, contig2seq, san = 1, debug = debug)


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
                
    
    if(len(alignment_list) <  o_alignment_list_len):#fill the gap
        filtered = True
        st = time.time()
        extend_edge_test(testseq, testseq_len, alignment_list, setting_kmersize, pos2contig, contig2start, contig2seq, san = 1, debug = debug)


    merge_conjacent_alignment(alignment_list, contig2start)

    fix_simple_inv(alignment_list, contig2start, contig2seq, testseq)

    new_alignment_list = List()
    cigarlist = []
    for alignment in alignment_list: 

        tmp_alignment_list, tmp_cigarlist = split_alignment_test(alignment, testseq, rc_testseq, testseq_len, kmersize=setting_kmersize , contig2seq = contig2seq, contig2start = contig2start, debug = debug, H = H, eqx = eqx)

        iloc = -1
        for alignment in tmp_alignment_list:
            iloc += 1
            new_alignment_list.append(alignment)
            cigarlist.append(tmp_cigarlist[iloc])


    alignment_list, onemapinfolist = get_onemapinfolist(new_alignment_list, cigarlist, readid, mapq, testseq_len, contig2start, need_reverse, option['H'])
    return alignment_list, onemapinfolist, TRA_signal, filtered
#RG


#09042024packnew
max_factor = 1000
@njit
def insertpoint_score(S_i, target_score, k, S_arg_i):
    i = 0
    j = k
    if(S_i[S_arg_i[0]] > target_score):
        return 0
    if(S_i[S_arg_i[k - 1]] < target_score):
        return k
    while(i < j):
        
        mid = (i + j) // 2

        now_score = S_i[S_arg_i[mid]]
        if(now_score < target_score):
            i = mid + 1
        elif(now_score > target_score):
            j = mid
        else:
            return mid + 1
    return j

#09082024

#09212024
def P_alignmentstring_comments(infodict, comments):
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
    def checkcomment(comments, infolist, tags):
        if(isinstance(comments, str)):
            comments_info = comments.split('\t')
            for onecomment in comments_info:
                info = onecomment.split(':')
                if((len(info) == 3) and (len(info[0]) == 2) and (info[0] not in tags) and (info[1] in ('A', 'i', 'f', 'Z', 'H', 'B'))):
                    infolist.append(onecomment)
                    tags.add(info[0])
    tags = set(('QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL', 'SA','NM','MD','cs'))
    for key in infodict:
        if(key in name2iloc):
            infolist[name2iloc[key]] = infodict[key]
        else:
            infolist.append(p_other_tag(key, infodict[key]))
            tags.add(key)
    checkcomment(comments, infolist, tags)
    return '\t'.join(infolist)
#09232024
def get_onemapinfolist(new_alignment_list, cigarlist, readid, mapq, testseq_len, contig2start, need_reverse, use_hardclip):
    if(use_hardclip == True):
        clipsyb = 'H'
    else:
        clipsyb = 'S'
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
                    topcigar = str(query_st) + clipsyb
                else:

                    topcigar = ''
                if((testseq_len - query_en) > 0):
                    tailcigar = str(testseq_len - query_en) + clipsyb
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
                    topcigar = str(query_st) + clipsyb
                else:
                    topcigar = ''
                if((testseq_len - query_en) > 0):
                    tailcigar = str(testseq_len - query_en) + clipsyb
                else:
                    tailcigar = ''
                onemapinfolist.append((readid, contig, '-', query_st, query_en, target_st-refbias, target_en-refbias, mapq, topcigar+cigarstring+tailcigar))
        if(use_hardclip == False):
            for line in onemapinfolist:
                if(testseq_len != len(Cigar(line[-1]))):
                    logging.error('testseq_len != Cigar(line[-1])')
                    raise Exception('testseq_len != Cigar(line[-1])')
        else:
            for line in onemapinfolist:
                if((line[4] - line[3]) != len(Cigar(line[-1]))):
                    logging.error('(line[4] - line[3]) != Cigar(line[-1])')
                    raise Exception('(line[4] - line[3]) != Cigar(line[-1])')
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
                    topcigar = str(query_st) + clipsyb
                else:

                    topcigar = ''
                if((testseq_len - query_en) > 0):
                    tailcigar = str(testseq_len - query_en) + clipsyb
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
                    topcigar = str(query_st) + clipsyb
                else:
                    topcigar = ''
                if((testseq_len - query_en) > 0):
                    tailcigar = str(testseq_len - query_en) + clipsyb
                else:
                    tailcigar = ''
                onemapinfolist.append((readid, contig, '+', query_st, query_en, target_st-refbias, target_en-refbias, mapq, topcigar+cigarstring+tailcigar))

        if(use_hardclip == False):
            for line in onemapinfolist:
                if(testseq_len != len(Cigar(line[-1]))):
                    logging.error('testseq_len != Cigar(line[-1])')
                    raise Exception('testseq_len != Cigar(line[-1])')
        else:
            for line in onemapinfolist:
                if((line[4] - line[3]) != len(Cigar(line[-1]))):
                    logging.error('(line[4] - line[3]) != Cigar(line[-1])')
                    raise Exception('(line[4] - line[3]) != Cigar(line[-1])')
        return new_alignment_list, onemapinfolist[::-1]
#20240930

def get_bam_dict_str(mapinfo, query, qual, contig2iloc, contig2seq, md, shortcs, cigar2cg, markunbalancetra, option):
    #'hhk',         ,  '1', '+', 11, 9192, 2767041, 2776138, 60
    #      0            1    2   3    4      5         6      7
    #'18_19897150_+', '18', '+', 0, 4776, 19832244, 19837393, 1]
    if(markunbalancetra == True):
        mapinfo = reassign_mapq(mapinfo)
    else:
        for iloc in range(len(mapinfo)):
            mapinfo[iloc] = list(mapinfo[iloc])

    
    hardclip = option['H']
    
    rc_query = str(Seq(query).reverse_complement())
    mapinfo.sort(key = sort_by_length)
    mapinfo = mapinfo[::-1]
    #mapinfo.sort(key = sortbycontig)
    iloc2nm = dict()
    iloc2md = dict()
    iloc2cs = dict()
    iloc2n_cigar = dict()
    tmpiloc = -1
    fakecigar = option['fakecigar']
    if(fakecigar == True):
        iloc2fakecigar = dict()
        if(hardclip == True):
            clipsyb = 'H'
        else:
            clipsyb = 'S'
    if(md == False):
        for item in mapinfo:
            item[-1], n_cigar = mergecigar_n(item[-1])
            tmpiloc += 1
            if(item[2] == '+'):
                nm = nm_from_cigar(item[8], query, get_refseq(item[1], item[5], item[6], contig2seq))
            else:
                nm = nm_from_cigar(item[8], rc_query, get_refseq(item[1], item[5], item[6], contig2seq))
            iloc2nm[tmpiloc] = nm
            iloc2n_cigar[tmpiloc] = n_cigar
            if(fakecigar == True):
                if(item[3] > 0):
                    onefaketop = str(item[3]) + clipsyb
                else:
                    onefaketop = ''
                if((len(query) - item[4]) > 0):
                    onefaketail = str(len(query) - item[4]) + clipsyb
                else:
                    onefaketail = ''
                tmpdiff = item[4] - item[3] - item[6] + item[5]
                if(tmpdiff > 0):
                    onefakebody = str(item[6] - item[5]) + 'M' + str(tmpdiff) + 'I'
                elif(tmpdiff < 0):
                    onefakebody = str(item[4] - item[3]) + 'M' + str(abs(tmpdiff)) + 'D'
                else:
                    onefakebody = str(item[4] - item[3]) + 'M'
                iloc2fakecigar[tmpiloc] =  ''.join((onefaketop, onefakebody, onefaketail))
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
            nm = nm_from_cigar(cigarstring, query if item[2] == '+' else rc_query, tmp_target)
            item[-1] = cigarstring
            iloc2nm[tmpiloc] = nm
            iloc2md[tmpiloc] = mdstring
            iloc2cs[tmpiloc] = csstring
            iloc2n_cigar[tmpiloc] = n_cigar
            if(fakecigar == True):
                if(item[3] > 0):
                    onefaketop = str(item[3]) + clipsyb
                else:
                    onefaketop = ''
                if((len(query) - item[4]) > 0):
                    onefaketail = str(len(query) - item[4]) + clipsyb
                else:
                    onefaketail = ''
                tmpdiff = item[4] - item[3] - item[6] + item[5]
                if(tmpdiff > 0):
                    onefakebody = str(item[6] - item[5]) + 'M' + str(tmpdiff) + 'I'
                elif(tmpdiff < 0):
                    onefakebody = str(item[4] - item[3]) + 'M' + str(abs(tmpdiff)) + 'D'
                else:
                    onefakebody = str(item[4] - item[3]) + 'M'
                iloc2fakecigar[tmpiloc] =  ''.join((onefaketop, onefakebody, onefaketail))

    if((qual != None) and (len(qual) == len(query))):
        query_qualities = qual
        rc_query_qualities = query_qualities[::-1]
    a_list = []
    primary_iloc = 0
    #if(len(mapinfo) > 1):
        #if(mapinfo[0][7] == 1 and mapinfo[1][7] != 1):
            #primary_iloc = 1

    #QNAME FLAG  RNAME  POS  MAPQ  CIGAR RNEXT  PNEXT  TLEN   SEQ   QUAL
    for iloc in range(len(mapinfo)):
        
        bam_dict = dict()
        if('rg-id' in option):
            bam_dict['RG'] = option['rg-id']
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
                nm = iloc2nm[tmpiloc]
                if(fakecigar == False):
                    salist.append(''.join((item[1], ',', str(item[5]+1), ',', item[2], ',', item[8], ',', str(mq), ',', str(nm)+';')))
                else:
                    salist.append(''.join((item[1], ',', str(item[5]+1), ',', item[2], ',', iloc2fakecigar[tmpiloc], ',', str(mq), ',', str(nm)+';')))


            bam_dict['SA'] = ''.join(salist)
        mq = mapinfo[iloc][7]
        item = primary

        bam_dict['MAPQ'] = str(mq)

        if(item[2] == '+'):
            if(hardclip == False):
                bam_dict['SEQ'] = query
                if((qual != None) and (len(qual) == len(query))):
                    bam_dict['QUAL'] = query_qualities
            else:
                bam_dict['SEQ'] = query[item[3]: item[4]]
                if((qual != None) and (len(qual) == len(query))):
                    bam_dict['QUAL'] = query_qualities[item[3]: item[4]]
                

        else:
            if(hardclip == False):
                bam_dict['SEQ'] = rc_query
                if((qual != None) and (len(qual) == len(query))):
                    bam_dict['QUAL'] = rc_query_qualities
            else:
                bam_dict['SEQ'] = rc_query[item[3]: item[4]]
                if((qual != None) and (len(qual) == len(query))):
                    bam_dict['QUAL'] = rc_query_qualities[item[3]: item[4]]

            
        bam_dict['NM'] = iloc2nm[iloc]
        if(md == True):
            bam_dict['MD'] = iloc2md[iloc]
            bam_dict['cs'] = iloc2cs[iloc]

        a_list.append(P_alignmentstring(bam_dict))
    return a_list
def get_bam_dict_str_comments(mapinfo, query, qual, comments, contig2iloc, contig2seq, md, shortcs, cigar2cg, markunbalancetra, option):
    #'hhk',         ,  '1', '+', 11, 9192, 2767041, 2776138, 60
    #      0            1    2   3    4      5         6      7
    #'18_19897150_+', '18', '+', 0, 4776, 19832244, 19837393, 1]
    if(markunbalancetra == True):
        mapinfo = reassign_mapq(mapinfo)
    else:
        for iloc in range(len(mapinfo)):
            mapinfo[iloc] = list(mapinfo[iloc])

    
    hardclip = option['H']
    
    rc_query = str(Seq(query).reverse_complement())
    mapinfo.sort(key = sort_by_length)
    mapinfo = mapinfo[::-1]
    #mapinfo.sort(key = sortbycontig)
    iloc2nm = dict()
    iloc2md = dict()
    iloc2cs = dict()
    iloc2n_cigar = dict()
    tmpiloc = -1
    fakecigar = option['fakecigar']
    if(fakecigar == True):
        iloc2fakecigar = dict()
        if(hardclip == True):
            clipsyb = 'H'
        else:
            clipsyb = 'S'
    if(md == False):
        for item in mapinfo:
            item[-1], n_cigar = mergecigar_n(item[-1])
            tmpiloc += 1
            if(item[2] == '+'):
                nm = nm_from_cigar(item[8], query, get_refseq(item[1], item[5], item[6], contig2seq))
            else:
                nm = nm_from_cigar(item[8], rc_query, get_refseq(item[1], item[5], item[6], contig2seq))
            iloc2nm[tmpiloc] = nm
            iloc2n_cigar[tmpiloc] = n_cigar
            if(fakecigar == True):
                if(item[3] > 0):
                    onefaketop = str(item[3]) + clipsyb
                else:
                    onefaketop = ''
                if((len(query) - item[4]) > 0):
                    onefaketail = str(len(query) - item[4]) + clipsyb
                else:
                    onefaketail = ''
                tmpdiff = item[4] - item[3] - item[6] + item[5]
                if(tmpdiff > 0):
                    onefakebody = str(item[6] - item[5]) + 'M' + str(tmpdiff) + 'I'
                elif(tmpdiff < 0):
                    onefakebody = str(item[4] - item[3]) + 'M' + str(abs(tmpdiff)) + 'D'
                else:
                    onefakebody = str(item[4] - item[3]) + 'M'
                iloc2fakecigar[tmpiloc] =  ''.join((onefaketop, onefakebody, onefaketail))
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
            nm = nm_from_cigar(cigarstring, query if item[2] == '+' else rc_query, tmp_target)
            item[-1] = cigarstring
            iloc2nm[tmpiloc] = nm
            iloc2md[tmpiloc] = mdstring
            iloc2cs[tmpiloc] = csstring
            iloc2n_cigar[tmpiloc] = n_cigar
            if(fakecigar == True):
                if(item[3] > 0):
                    onefaketop = str(item[3]) + clipsyb
                else:
                    onefaketop = ''
                if((len(query) - item[4]) > 0):
                    onefaketail = str(len(query) - item[4]) + clipsyb
                else:
                    onefaketail = ''
                tmpdiff = item[4] - item[3] - item[6] + item[5]
                if(tmpdiff > 0):
                    onefakebody = str(item[6] - item[5]) + 'M' + str(tmpdiff) + 'I'
                elif(tmpdiff < 0):
                    onefakebody = str(item[4] - item[3]) + 'M' + str(abs(tmpdiff)) + 'D'
                else:
                    onefakebody = str(item[4] - item[3]) + 'M'
                iloc2fakecigar[tmpiloc] =  ''.join((onefaketop, onefakebody, onefaketail))

    if((qual != None) and (len(qual) == len(query))):
        query_qualities = qual
        rc_query_qualities = query_qualities[::-1]
    a_list = []
    primary_iloc = 0
    #if(len(mapinfo) > 1):
        #if(mapinfo[0][7] == 1 and mapinfo[1][7] != 1):
            #primary_iloc = 1

    #QNAME FLAG  RNAME  POS  MAPQ  CIGAR RNEXT  PNEXT  TLEN   SEQ   QUAL
    for iloc in range(len(mapinfo)):
        
        bam_dict = dict()
        if('rg-id' in option):
            bam_dict['RG'] = option['rg-id']
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
                nm = iloc2nm[tmpiloc]
                if(fakecigar == False):
                    salist.append(''.join((item[1], ',', str(item[5]+1), ',', item[2], ',', item[8], ',', str(mq), ',', str(nm)+';')))
                else:
                    salist.append(''.join((item[1], ',', str(item[5]+1), ',', item[2], ',', iloc2fakecigar[tmpiloc], ',', str(mq), ',', str(nm)+';')))


            bam_dict['SA'] = ''.join(salist)
        mq = mapinfo[iloc][7]
        item = primary

        bam_dict['MAPQ'] = str(mq)


        if(item[2] == '+'):
            if(hardclip == False):
                bam_dict['SEQ'] = query
                if((qual != None) and (len(qual) == len(query))):
                    bam_dict['QUAL'] = query_qualities
            else:
                bam_dict['SEQ'] = query[item[3]: item[4]]
                if((qual != None) and (len(qual) == len(query))):
                    bam_dict['QUAL'] = query_qualities[item[3]: item[4]]
                

        else:
            if(hardclip == False):
                bam_dict['SEQ'] = rc_query
                if((qual != None) and (len(qual) == len(query))):
                    bam_dict['QUAL'] = rc_query_qualities
            else:
                bam_dict['SEQ'] = rc_query[item[3]: item[4]]
                if((qual != None) and (len(qual) == len(query))):
                    bam_dict['QUAL'] = rc_query_qualities[item[3]: item[4]]
        bam_dict['NM'] = iloc2nm[iloc]
        if(md == True):
            bam_dict['MD'] = iloc2md[iloc]
            bam_dict['cs'] = iloc2cs[iloc]

        a_list.append(P_alignmentstring_comments(bam_dict, comments))
    return a_list
def get_reversed_chain_numpy_rough(raw_alignment_array, testseq_len):
    if(len(raw_alignment_array) < 3):
        return False, raw_alignment_array
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

                
                cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 2, gap_open_2 = 24, gap_extend_2 = 1, bw=-1, zdropvalue=-1, eqx = eqx)
                
                new_alignment[-1].append(nowitem)
                cigarlist[-1].append(cigarstring)
                    #print_log(cigarstring)#
            else:
                logging.error(str(preitem) + ' '+ str(nowitem))
                logging.error("ERROR: Failed to compute CIGAR ")
                raise Exception("ERROR: Failed to compute CIGAR ")

            preitem = nowitem
            iloc += 1
        if(cigarlist[-1] == []):
            logging.error('cigarlist[-1] == []')
            logging.error("ERROR: Failed to compute CIGAR ")
            raise Exception("ERROR: Failed to compute CIGAR ")


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
                cigarstring, zdropedcode, q_e, t_e, tmpdelcount, tmpinscount = mp.k_cigar(target, query, match = 2, mismatch = -4, gap_open_1 = 4, gap_extend_1 = 2, gap_open_2 = 24, gap_extend_2 = 1, bw=-1, zdropvalue=-1, eqx = eqx)
                
                new_alignment[-1].append(nowitem)
                cigarlist[-1].append(cigarstring)
                    #print_log(cigarstring)#
            else:
                logging.error(str(preitem) + ' '+ str(nowitem))
                logging.error("ERROR: Failed to compute CIGAR ")
                raise Exception("ERROR: Failed to compute CIGAR ")


            preitem = nowitem
            iloc += 1
        if(cigarlist[-1] == []):
            logging.error('cigarlist[-1] == []')
            logging.error("ERROR: Failed to compute CIGAR ")
            raise Exception("ERROR: Failed to compute CIGAR ")

        return new_alignment, cigarlist
    
#scar
#print_log('SCAR')
#more chain
@njit
def get_localmap_multi_all_forDP_inv_guide_1(oreadstart, oreadend, one_mapinfo, raw_alignment_array, testseq, rc_testseq, contig2start, contig2seq, kmersize, skipcost, maxdiff, maxgap, shift = 1):#


    def seq2hashtable_multi_test(multi, onelookuptable_s, onelookuptable_m, seq, start, kmersize):
        skiphash = hash('N'*kmersize)
        
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
        

    if(raw_alignment_array[0][0] < raw_alignment_array[-1][0]):
        look_span = max(raw_alignment_array[0][0], len(testseq) - raw_alignment_array[-1][0])+500
    else:
        look_span = max(raw_alignment_array[-1][0], len(testseq) - raw_alignment_array[0][0])+500
    look_span = 7000
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
    
    multi = [0]
    multi.pop(0)

    retry_diffcontig = False
    for item in startandend:
        min_ref, max_ref = item[0], item[1]
        testcontig = pos2contig(min_ref, contig2start)
        if(testcontig != pos2contig(max_ref, contig2start)):
            retry_diffcontig = True
            break
        lookfurther = min(look_span, min_ref-contig2start[testcontig])
        min_ref -= lookfurther
        max_ref += look_span
        refseq = contig2seq[testcontig][min_ref-contig2start[testcontig]: max_ref-contig2start[testcontig]]
        seq2hashtable_multi_test(multi, local_lookuptable_s, local_lookuptable_m, refseq, min_ref, kmersize)
    for hashedkmer in multi:
        local_lookuptable_s.pop(hashedkmer)

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
        
        multi = [0]
        multi.pop(0)

        retry_diffcontig = False
        for item in startandend:
            min_ref, max_ref = item[0], item[1]
            testcontig = pos2contig(min_ref, contig2start)
            if(testcontig != pos2contig(max_ref, contig2start)):
                retry_diffcontig = True
                break
            lookfurther = min(look_span, min_ref-contig2start[testcontig])
            min_ref -= lookfurther
            max_ref += look_span
            refseq = contig2seq[testcontig][min_ref-contig2start[testcontig]: max_ref-contig2start[testcontig]]
            seq2hashtable_multi_test(multi, local_lookuptable_s, local_lookuptable_m, refseq, min_ref, kmersize)
        for hashedkmer in multi:
            local_lookuptable_s.pop(hashedkmer)


    raw_alignment_array = raw_alignment_array[np.argsort(raw_alignment_array[:, 0])]
    if(oreadend < 0):
        readstart = max(0, raw_alignment_array[0][0]-500)
        readend = min(len(testseq)-kmersize+1, raw_alignment_array[-1][0]+500)
    else:
        readstart, readend = oreadstart, oreadend

    readstart = max(0, raw_alignment_array[0][0]-look_span)
    readend = min(len(testseq)-kmersize+1, raw_alignment_array[-1][0]+look_span)


    iloc = readstart
    iloc -= shift


    pointdict = Dict()
    pointdict[21312] = (1, 1, 1, 1)
    pointdict.pop(21312)

    pointdict_key = [21312]
    pointdict_key.pop(0)
    
    readposarr = raw_alignment_array[:, 0].astype(np.int32)

    while(True):

        iloc += shift
        if(iloc >= readend):
            break
        forwardhashedkmer, reversehashedkmer = hash(testseq[iloc: iloc + kmersize]), hash(rc_testseq[-(iloc + kmersize): -iloc])
        if(forwardhashedkmer == reversehashedkmer):
            continue

        biasvalue, biasvalue_1, closest_index, closest_index_1 = findClosest_1(readposarr, target = iloc)


        interval = min(biasvalue + biasvalue_1 + 500, 2000)
        closest_refpos_1 = raw_alignment_array[closest_index][1]
        closest_refpos_2 = raw_alignment_array[closest_index_1][1]
        upperrefloc =  (closest_refpos_1 + interval, closest_refpos_2 + interval)
        lowerrefloc =  (closest_refpos_1 - interval, closest_refpos_2 - interval)
        readgap = abs(iloc - raw_alignment_array[closest_index][0])

        hashedkmer = forwardhashedkmer
        if(hashedkmer in local_lookuptable_s):
            refloc = local_lookuptable_s[hashedkmer]
            refgap = abs(refloc - closest_refpos_1)
            diff = abs(readgap - refgap)
            if((diff < 500) or (upperrefloc[0] >= refloc and lowerrefloc[0] <= refloc) or (upperrefloc[1] >= refloc and lowerrefloc[1] <= refloc)):
                item = (iloc, refloc, 1)
                point = item[1] - item[0]

                if(point in pointdict):
                    cache_item_0, cache_item_1, cache_item_2, cache_item_3  = pointdict[point]
                    if((cache_item_0 + cache_item_3) >= item[0]):
                        bouns = item[0] - (cache_item_0 + cache_item_3) + kmersize
                        if(bouns > 0):
                            if((cache_item_3 + bouns < 20)):
                                pointdict[point] = (cache_item_0, cache_item_1, 1, cache_item_3 + bouns)
                            else:
                                one_mapinfo.append((cache_item_0, cache_item_1, cache_item_2, cache_item_3))
                                pointdict[point] = (cache_item_0 + cache_item_3, cache_item_1 + cache_item_3, 1, bouns)

                                
                    else:    
                        one_mapinfo.append((cache_item_0, cache_item_1, cache_item_2, cache_item_3))
                        pointdict[point] = (item[0], item[1], item[2], kmersize)
                else:
                    pointdict[point] = (item[0], item[1], item[2], kmersize)
                    pointdict_key.append(point)

        elif(hashedkmer in local_lookuptable_m):
            for refloc in local_lookuptable_m[hashedkmer]:
                refgap = abs(refloc - closest_refpos_1)
                diff = abs(readgap - refgap)
                if((diff < 500) or (upperrefloc[0] >= refloc and lowerrefloc[0] <= refloc) or (upperrefloc[1] >= refloc and lowerrefloc[1] <= refloc)):
                    item = (iloc, refloc, 1)
                    point = item[1] - item[0]

                    if(point in pointdict):
                        cache_item_0, cache_item_1, cache_item_2, cache_item_3  = pointdict[point]
                        if((cache_item_0 + cache_item_3) >= item[0]):
                            bouns = item[0] - (cache_item_0 + cache_item_3) + kmersize
                            if(bouns > 0):
                                if((cache_item_3 + bouns < 20)):
                                    pointdict[point] = (cache_item_0, cache_item_1, 1, cache_item_3 + bouns)
                                else:
                                    one_mapinfo.append((cache_item_0, cache_item_1, cache_item_2, cache_item_3))
                                    pointdict[point] = (cache_item_0 + cache_item_3, cache_item_1 + cache_item_3, 1, bouns)


                        else:    
                            one_mapinfo.append((cache_item_0, cache_item_1, cache_item_2, cache_item_3))
                            pointdict[point] = (item[0], item[1], item[2], kmersize)
                    else:
                        pointdict[point] = (item[0], item[1], item[2], kmersize)
                        pointdict_key.append(point)


        hashedkmer = reversehashedkmer
        if(hashedkmer in local_lookuptable_s):
            refloc = local_lookuptable_s[hashedkmer]
            refgap = abs(refloc - closest_refpos_1)
            diff = abs(readgap - refgap)
            if((diff < 500) or (upperrefloc[0] >= refloc and lowerrefloc[0] <= refloc) or (upperrefloc[1] >= refloc and lowerrefloc[1] <= refloc)):
                item = (iloc, refloc, -1)
                point = -(item[1] + item[0])

                if(point in pointdict):
                    cache_item_0, cache_item_1, cache_item_2, cache_item_3  = pointdict[point]
                    if((cache_item_0 + cache_item_3) >= item[0]):
                        bouns = item[0] - (cache_item_0 + cache_item_3) + kmersize
                        if(bouns > 0):
                            if((cache_item_3 + bouns < 20)):
                                pointdict[point] = (cache_item_0, item[1], -1, cache_item_3 + bouns)
                            else:
                                one_mapinfo.append((cache_item_0, cache_item_1, cache_item_2, cache_item_3))
                                pointdict[point] = (cache_item_0 + cache_item_3, item[1], -1, bouns)
                                    #print_log(pointdict[point][-1])


                    else:    
                        one_mapinfo.append((cache_item_0, cache_item_1, cache_item_2, cache_item_3))
                        pointdict[point] = (item[0], item[1], item[2], kmersize)
                else:
                    pointdict[point] = (item[0], item[1], item[2], kmersize)
                    pointdict_key.append(point)

        elif(hashedkmer in local_lookuptable_m):

            for refloc in local_lookuptable_m[hashedkmer]:
                refgap = abs(refloc - closest_refpos_1)
                diff = abs(readgap - refgap)
                if((diff < 500) or (upperrefloc[0] >= refloc and lowerrefloc[0] <= refloc) or (upperrefloc[1] >= refloc and lowerrefloc[1] <= refloc)):
                    item = (iloc, refloc, -1)
                    point = -(item[1] + item[0])

                    if(point in pointdict):
                        cache_item_0, cache_item_1, cache_item_2, cache_item_3  = pointdict[point]
                        if((cache_item_0 + cache_item_3) >= item[0]):
                            bouns = item[0] - (cache_item_0 + cache_item_3) + kmersize
                            if(bouns > 0):
                                if((cache_item_3 + bouns < 20)):
                                    pointdict[point] = (cache_item_0, item[1], -1, cache_item_3 + bouns)
                                else:
                                    one_mapinfo.append((cache_item_0, cache_item_1, cache_item_2, cache_item_3))
                                    pointdict[point] = (cache_item_0 + cache_item_3, item[1], -1, bouns)
                                        #print_log(pointdict[point][-1])


                        else:    
                            one_mapinfo.append((cache_item_0, cache_item_1, cache_item_2, cache_item_3))
                            pointdict[point] = (item[0], item[1], item[2], kmersize)
                    else:
                        pointdict[point] = (item[0], item[1], item[2], kmersize)
                        pointdict_key.append(point)
                        
    for key in pointdict_key:
        one_mapinfo.append(pointdict[key])
    return readstart, readend
    

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


            if((abs(readgap - refgap) <= large_cost) and (refgap >= (-20)) and (readgap < 100)):
                if(pos2contig(preitem[1], contig2start) == pos2contig(nowitem[1], contig2start)):

                    if(refgap >= 0):
                        alignment_list[-1].append(nowitem)
                        preitem = nowitem
                        continue
                    else:
                        if(readgap <= 20):
                            continue
                        else:
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
    def select_secondary_alignment(order, path_list, scores_list, S_arr, testseq_len):
        secondary_path_list = List()
        secondary_path_list.append([(0, 0, 0, 0)])
        secondary_path_list.pop()
        if(len(path_list) > 1):
            en_loc = testseq_len
            loc2score = np.zeros(testseq_len)
            for anchor, score in zip(path_list[0], S_arr):
                st_loc = anchor[0]
                loc2score[st_loc: en_loc] = score
                en_loc = st_loc
            for iloc in order[1:]:
                one_path, f2 = path_list[iloc], scores_list[iloc]
                en_loc, st_loc = one_path[0][0], one_path[-1][0]
                if(en_loc - st_loc < 50):
                    continue
                f1 = max(loc2score[en_loc] - loc2score[st_loc], 1)
                if(f2 / f1 > 0.9 or abs(f1 - f2) < 40):
                    skip = False
                    for pri_path in secondary_path_list:
                        pri_en_loc, pri_st_loc = pri_path[0][0], pri_path[-1][0]
                        overlap_size = max(min(en_loc, pri_en_loc) - max(pri_st_loc, st_loc), 0)
                        if((overlap_size / (en_loc - st_loc)) > 0.5):
                            skip = True
                            break
                    if(skip == False):
                        secondary_path_list.append(one_path)
  
            
            return secondary_path_list
                    
                
        else:
            return secondary_path_list
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
    S_arr = [0.]
    S_arr.pop(0)
    take_index = g_max_index
    usedindex.add(take_index)
    score = S[take_index]
    while(True):
        if((P[take_index] == -9999999)):
            path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
            S_arr.append(S[take_index])
            break
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
        S_arr.append(S[take_index])
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
    #if(hit == True and max_scores > 40):#hifi
    if(hit == True and max_scores > 60):#clr, ont

        order = np.argsort(np.array(scores_list))[::-1]
        #print_log('scores_list[order[0]], len(path_list[order[0]])', scores_list[order[0]], len(path_list[order[0]]))
        primary_rlocset_List = List()
        
        primary_scores_List = List()

        
        primary_index_List = List()

        
        all_index_List = List()
        
        if(order[0] != 0):
            for iloc in range(len(order)):
                if(order[iloc] == 0):
                    order[iloc] = order[0]
                    order[0] = 0
                    break

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


        return path_list, primary_index_List, primary_scores_List, all_index_List, mapq, scores_list, factor, select_secondary_alignment(order, path_list, scores_list, S_arr, testseq_len)


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
        
        secondary_path_list = List()
        secondary_path_list.append([(0, 0, 0, 0)])
        secondary_path_list.pop()
        return path_list, primary_index_List, primary_scores_List, all_index_List, 0, scores_list, factor, secondary_path_list
    
@njit
def hit2work_1_64(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage = 0.5, hastra = False, H = False):
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
    def select_secondary_alignment(order, path_list, scores_list, S_arr, testseq_len):
        secondary_path_list = List()
        secondary_path_list.append([(0, 0, 0, 0)])
        secondary_path_list.pop()
        if(len(path_list) > 1):
            en_loc = testseq_len
            loc2score = np.zeros(testseq_len)
            for anchor, score in zip(path_list[0], S_arr):
                st_loc = anchor[0]
                loc2score[st_loc: en_loc] = score
                en_loc = st_loc
            for iloc in order[1:]:
                one_path, f2 = path_list[iloc], scores_list[iloc]
                en_loc, st_loc = one_path[0][0], one_path[-1][0]
                if(en_loc - st_loc < 50):
                    continue
                f1 = max(loc2score[en_loc] - loc2score[st_loc], 1)
                if(f2 / f1 > 0.9 or abs(f1 - f2) < 40):
                    skip = False
                    for pri_path in secondary_path_list:
                        pri_en_loc, pri_st_loc = pri_path[0][0], pri_path[-1][0]
                        overlap_size = max(min(en_loc, pri_en_loc) - max(pri_st_loc, st_loc), 0)
                        if((overlap_size / (en_loc - st_loc)) > 0.5):
                            skip = True
                            break
                    if(skip == False):
                        secondary_path_list.append(one_path)
            
            return secondary_path_list
                    
                
        else:
            return secondary_path_list
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
        g_max_index, S, P, S_arg, factor = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_all_64(one_mapinfo, kmersize = kmersize, skipcost = skipcost[0], maxdiff = maxdiff[0], maxgap = maxgap)


    if(fast_enable == True or g_max_index == -1):
        fast_enable = True
        g_max_index, S, P, S_arg = get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_fast_all_64(one_mapinfo, kmersize = kmersize, skipcost = skipcost[0], maxdiff = maxdiff[0], maxgap = maxgap)

    scores = S[g_max_index]
    #if(scores < minichain_scores):
        #if(hit == False):
            #return 0, 0., [(0, 0, 0, 0)]
        #continue

    ##################
    usedindex = set()

    path = [(0, 0, 0, 0)]
    path.pop()
    S_arr = [0.]
    S_arr.pop(0)
    take_index = g_max_index
    usedindex.add(take_index)
    score = S[take_index]
    while(True):
        if((P[take_index] == -9999999)):
            path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
            S_arr.append(S[take_index])
            break
        path.append((one_mapinfo[take_index][0], one_mapinfo[take_index][1], one_mapinfo[take_index][2], one_mapinfo[take_index][3]))
        S_arr.append(S[take_index])
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
    #if(hit == True and max_scores > 40):#hifi
    if(hit == True and max_scores > 60):#clr, ont

        order = np.argsort(np.array(scores_list))[::-1]
        #print_log('scores_list[order[0]], len(path_list[order[0]])', scores_list[order[0]], len(path_list[order[0]]))
        primary_rlocset_List = List()
        
        primary_scores_List = List()

        
        primary_index_List = List()

        
        all_index_List = List()

        if(order[0] != 0):
            for iloc in range(len(order)):
                if(order[iloc] == 0):
                    order[iloc] = order[0]
                    order[0] = 0
                    break
        
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


        return path_list, primary_index_List, primary_scores_List, all_index_List, mapq, scores_list, factor, select_secondary_alignment(order, path_list, scores_list, S_arr, testseq_len)


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
        
        secondary_path_list = List()
        secondary_path_list.append([(0, 0, 0, 0)])
        secondary_path_list.pop()
        return path_list, primary_index_List, primary_scores_List, all_index_List, 0, scores_list, factor, secondary_path_list
    
def decode_hit(index_object, index2contig, testseq, testseq_len, contig2start, kmersize, contig2seq, skipcost = (50., 30.), maxdiff = (50, 30), maxgap = 200, check_num = 20, c_bias = 5000, bin_size = 100, overlapprecentage = 0.5, hastra = True, H = False, mid_occ = -1):            
    #redo_ratio = 5#hifi
    redo_ratio = 10#clr, ont
    st = time.time()
    need_reverse, one_mapinfo = get_reversed_chain_numpy_rough(np.array(index_object.map(testseq, check_num = check_num,mid_occ = mid_occ)), testseq_len)
    if(len(one_mapinfo) <= 2):
        scores, path = 0., []
        return 0, scores, path, 0, []
    #print_log(need_reverse)
    #print_log('index_object.map', time.time() - st)
    maxgap = 1000
    #print_log('len(one_mapinfo)', len(one_mapinfo))
    st = time.time()
    if(len(one_mapinfo) < 2147483638):
        path_list, primary_index_List, primary_scores_List, all_index_List, mapq, scores_list, factor, secondary_path_list =  hit2work_1(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage, hastra, H = H)
    else:
        path_list, primary_index_List, primary_scores_List, all_index_List, mapq, scores_list, factor, secondary_path_list =  hit2work_1_64(one_mapinfo, index2contig, contig2start, testseq_len, skipcost, maxdiff, maxgap, check_num, c_bias, bin_size, kmersize, overlapprecentage, hastra, H = H)
    if(len(path_list) == 0):
        scores, path = 0., []
        return 0, scores, path, factor, []
    #print_log(scores_list)
    #print_log('hit2work_1', time.time() - st)
    #print_log('mapq', mapq)
    st = time.time()
    if(len(primary_index_List) != 0):
        base_iloc = primary_index_List[0]
        
        
        scores, path = scores_list[base_iloc], path_list[base_iloc]
        return_path_list = [path]
        for one_path in secondary_path_list:
            return_path_list.append(one_path)
    else:
        return_path_list = []
        scores, path = 0., [(0, 0, 0, 0)]
    #print_log('edlib', time.time() - st)
    if(need_reverse == True):
        return mapq, -scores, path, factor, return_path_list
    else:
        return mapq, scores, path, factor, return_path_list


def get_readmap_DP_test(readid, testseq, contig2start, contig2seq, index_object, index2contig, option, kmersize = 15, local_kmersize = 9, local_maxdiff = 30, refine = True, local_skipcost = 59., golbal_skipcost = (45., 45.),  golbal_maxdiff = (30, 30), check_num = 20, bin_size = 100, hastra = False, debug = False, H = False, mid_occ = -1, redo_ratio = 5, eqx = False):
    kmersize = index_object.k
    golbal_skipcost = (option['golbal_skipcost'], option['golbal_skipcost'])
    golbal_maxdiff = (option['golbal_maxdiff'], option['golbal_maxdiff'])
    local_skipcost = option['local_skipcost']
    local_maxdiff = option['local_maxdiff']

    setting_kmersize = kmersize
    setting_maxdiff = golbal_maxdiff[1]

    
    rc_testseq = str(Seq(testseq).reverse_complement())
    testseq_len = len(testseq)
                           
    mapq, scores, raw_alignment_list, redo_flag, path_list = decode_hit(index_object, index2contig, testseq, testseq_len, contig2start, kmersize, contig2seq, skipcost = golbal_skipcost, maxdiff = golbal_maxdiff, maxgap = 200, check_num = check_num, c_bias = 5000, bin_size = bin_size, overlapprecentage = 0.5, hastra = hastra, H = H, mid_occ = mid_occ)


    if(scores == 0.):
        return [], ([], []), 0, redo_flag
    
    
    if((scores < 0.)):
        need_reverse = True
    else:
        need_reverse = False

    setting_maxdiff = local_maxdiff
    setting_kmersize = local_kmersize


    new_path_list = List()
    for one_path in path_list:
        new_path_list.append(np.array(one_path))

    if(need_reverse == False):
        scores, raw_alignment_list = get_localmap_multi_all_forDP_inv_guide_list(new_path_list, testseq, rc_testseq, contig2start, contig2seq, kmersize = setting_kmersize, skipcost = local_skipcost, maxdiff = setting_maxdiff, maxgap = 99, shift = 1)

    else:
        scores, raw_alignment_list = get_localmap_multi_all_forDP_inv_guide_list(new_path_list, rc_testseq, testseq, contig2start, contig2seq, kmersize = setting_kmersize, skipcost = local_skipcost, maxdiff = setting_maxdiff, maxgap = 99, shift = 1)
        testseq, rc_testseq = rc_testseq, testseq

    if(len(raw_alignment_list) <= 1):
        return [], ([], []), 0, redo_flag


    alignment_list, onemapinfolist, TRA_signal, filtered = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, maxdiffratio = option['maxdivergence'], option = option, debug = debug, H = H, nofilter = option['nodiscard'], eqx = eqx)
    if(len(onemapinfolist) == 0):
        return [], ([], []), 0, redo_flag

    
    if(option['nodiscard'] == False and filtered == True and pairedindel(List([line[-1] for line in onemapinfolist]), indelsize = 30) == True):
        alignment_list, onemapinfolist, TRA_signal, filtered = extend_func(List(raw_alignment_list[::-1]), readid, mapq, testseq, rc_testseq, testseq_len, setting_kmersize, pos2contig, contig2start, contig2seq, setting_maxdiff, need_reverse, maxdiffratio = option['maxdivergence'], option = option, debug = debug, H = H, nofilter = True, eqx = eqx)
        
    
    return onemapinfolist, (alignment_list, raw_alignment_list), TRA_signal, redo_flag
#output cache
def get_list_of_readmap_stdout(raw_queue, cooked_queue, minimap, contig2seq, hastra, H, header, option):
    cache_size = 100
    cache_size_chr = 2000000
    a_list = []
    a_chars = 0
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

    
    rt_list = []


    n_in = n_nomap = n_dperr = n_outerr = 0
    while(True):
        readidandseq = raw_queue.get()
        if(type(readidandseq) == int):
            break


        n_in += 1
        try:
            onemapinfolist, (alignment_list,raw_alignment_list), TRA_signal, f_redo_ratio = get_readmap_DP_test(readidandseq[0], readidandseq[1], contig2start, contig2seq, minimap, index2contig, option, hastra = False, redo_ratio = redo_ratio, eqx = option['eqx'], check_num = option['c'])
        except Exception as error:
            n_dperr += 1

            if(option['debug'] == True):
                logging.error(error)
                logging.info(readidandseq[0])
                logging.info(readidandseq[1])

            continue


        if(len(onemapinfolist) == 0):
            n_nomap += 1
        if(len(onemapinfolist) != 0):
            try:
                tmp_a_list = get_bam_dict_str(onemapinfolist, readidandseq[1], readidandseq[2], contig2iloc, contig2seq, option['md'], option['shortcs'], option['cigar2cg'], option['markunbalancetra'], option)
            except Exception as error:
                n_outerr += 1
                if(option['debug'] == True):
                    logging.info(error)
                    logging.info(readidandseq[0])
                    logging.info(readidandseq[1])
                continue
            if((tmp_a_list) == None):
                n_outerr += 1
                continue
            else:
                for a_line in tmp_a_list:
                    a_chars += len(a_line)
                a_list += tmp_a_list
                if(a_chars > cache_size_chr):
                    cooked_queue.put(a_list)
                    a_list = []
                    a_chars = 0
                else:
                    continue
    if(n_in > 0 and (n_nomap > 0 or n_dperr > 0 or n_outerr > 0)):
        logging.info("[AUDIT] mode-H output worker: input_reads=%d no_alignment=%d dp_error=%d output_error=%d emitted=%d" % (n_in, n_nomap, n_dperr, n_outerr, n_in - n_nomap - n_dperr - n_outerr))
    if(len(a_list) > 0):
        cooked_queue.put(a_list)
        
def get_list_of_readmap_stdout_comments(raw_queue, cooked_queue, minimap, contig2seq, hastra, H, header, option):
    cache_size = 100
    cache_size_chr = 2000000
    a_list = []
    a_chars = 0
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

    
    rt_list = []


    while(True):
        readidandseq = raw_queue.get()
        if(type(readidandseq) == int):
            break


        try:
            onemapinfolist, (alignment_list,raw_alignment_list), TRA_signal, f_redo_ratio = get_readmap_DP_test(readidandseq[0], readidandseq[1], contig2start, contig2seq, minimap, index2contig, option, hastra = False, redo_ratio = redo_ratio, eqx = option['eqx'], check_num = option['c'])
        except Exception as error:

            if(option['debug'] == True):
                logging.error(error)
                logging.info(readidandseq[0])
                logging.info(readidandseq[1])

            continue


        if(len(onemapinfolist) != 0):
            try:
                tmp_a_list = get_bam_dict_str_comments(onemapinfolist, readidandseq[1], readidandseq[2], readidandseq[3], contig2iloc, contig2seq, option['md'], option['shortcs'], option['cigar2cg'], option['markunbalancetra'], option)
            except Exception as error:
                if(option['debug'] == True):
                    logging.info(error)
                    logging.info(readidandseq[0])
                    logging.info(readidandseq[1])
                continue
            if((tmp_a_list) == None):
                continue
            else:
                for a_line in tmp_a_list:
                    a_chars += len(a_line)
                a_list += tmp_a_list
                if(a_chars > cache_size_chr):
                    cooked_queue.put(a_list)
                    a_list = []
                    a_chars = 0
                else:
                    continue
    if(len(a_list) > 0):
        cooked_queue.put(a_list)
#@njit
def fix_simple_inv(alignment_list, contig2start, contig2seq, testseq):
    def get_rc(s):


        rc_s = ''
        for c in s[::-1]:
            if(c == 'A'):
                rc_s += 'T'
            elif(c == 'T'):
                rc_s +=  'A'
            elif(c == 'G'):
                rc_s +=  'C'
            elif(c == 'C'):
                rc_s +=  'G'
            else:
                rc_s +=  'N'

        return rc_s
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
                    '''print('readen_0, refen_0', readen_0, refen_0)
                    print('readst_1, refst_1', readst_1, refst_1)
                    print()
                    print('readen_1, refen_1', readen_1, refen_1)
                    print('readst_2, refst_2', readst_2, refst_2)
                    print()'''
     
                    if(refst_2 - refen_0 == refen_1 - refst_1 and readst_1 - readen_0 + readst_2 - readen_1 == 0):
                        if(refst_1 - refen_0 != 0 and refst_1 - refen_0 + refst_2 - refen_1 == 0):
                            #print('refen_0, refst_1', refen_0, refst_1)
                            if(refen_0 > refst_1):
                                tempref = get_rc((contig2seq[tempcontig][refen_1: refen_1 + refen_0 - refst_1]))
                                tempquery = (testseq[readen_0 - refen_0 + refst_1: readen_0])
                                if(tempref == tempquery):
                                    #print('tempref == tempquery')
                                    bias = refen_0 - refst_1

                                    alignment_list[iloc+2][0] = (readst_2 - bias, refst_2 - bias + temprefbias, 1, 0)
                                    insertitem = (readst_2 - bias, refen_0 + temprefbias, -1, 0)
               

                                    while(insertitem[0] <= (alignment_list[iloc + 1][-1][0] + alignment_list[iloc + 1][-1][3])):
                                        alignment_list[iloc + 1].pop(-1)
                                    assert len(alignment_list[iloc + 1]) > 0
                                    alignment_list[iloc + 1].append(insertitem)
                                    #print(alignment_list[iloc + 1][-5:])
                            else:#refen_0 < refst_1
                                tempref = (contig2seq[tempcontig][refen_0: refen_0 - refen_0 + refst_1])
                                tempquery = (testseq[readen_0: readen_0 - refen_0 + refst_1])
                                if(tempref == tempquery):
                                    #print('tempref == tempquery')
                                    alignment_list[iloc][-1] = (readen_0 - refen_0 + refst_1, refen_0 - refen_0 + refst_1 + temprefbias, 1, 0)
                                    insertitem = (readen_0 - refen_0 + refst_1, refen_1 + refen_0 - refst_1 + temprefbias, -1, 0)
                                    #print((readen_0 - refen_0 + refst_1, refen_0 - refen_0 + refst_1 + temprefbias, 1, 0))
                                    #print(insertitem)

                                    while(insertitem[0]>=alignment_list[iloc + 1][0][0]):
                                        alignment_list[iloc + 1].pop(0)
                                    assert len(alignment_list[iloc + 1]) > 0
                                    alignment_list[iloc + 1].insert(0, insertitem)
                                        #alignment_list[iloc][0] = (readen_0 - refen_0 + refst_1, )
                            '''print(tempref, tempquery)        
                            print(refen_0 )
                            print(readst_1, readen_1, refst_1 , refen_1 )
                            print(refst_2)
                    print('end')'''

                                
            iloc += 1
#new eqution
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_all_64(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    extra_size = len(extra) - 1
    
    oskipcost = skipcost#mark_a
    omaxdiff = maxdiff#mark_a
    repeat_weight = 20#mark_a

    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):

        gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))

    
    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)


    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)

            
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)

    
    testspace_en = 1
    
    
    S_arg = np.empty(n, np.int64)
    S_arg[0] = 0
    
    i = 0

    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    

    opcount = 0
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        

        max_scores = one_mapinfo[i][3]

        
        pre_index = -9999999
        

        if(prereadloc < one_mapinfo[i][0]):
            
            if((opcount/i) > max_factor):#mark 1
                return -1, S, P, S_arg, -1.#mark 1
            
            
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = insertpoint_score(S, S[k], k, S_arg)

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
               
        for j in S_arg[:testspace_en][::-1]:
            
            
            if(S[j] > (max_scores - one_mapinfo[i][3])):
                opcount += 1


                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                if((readgap < 0)):
                    bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                    readgap = 0
                    overlap_size = one_mapinfo[j][0] +one_mapinfo[j][3] - one_mapinfo[i][0]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = one_mapinfo[i][1] + overlap_size - (one_mapinfo[j][1] + one_mapinfo[j][3])
                        else:
                            refgap = one_mapinfo[j][1] - (one_mapinfo[i][1] + bonus)
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] + overlap_size - one_mapinfo[j][1] + 1
                        else:
                            refgap = one_mapinfo[i][1] + bonus - 1 - (one_mapinfo[j][1] + one_mapinfo[j][3])

                else:
                    bonus = one_mapinfo[i][3]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                        else:
                            refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] - one_mapinfo[j][1] + 1
                        else:
                            refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - 1 - one_mapinfo[j][1] - one_mapinfo[j][3]


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
            else:
                break
                

        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
       

    k = testspace_en
    while(k < n):

        loc_in_sorted_S = insertpoint_score(S, S[k], k, S_arg)

        S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
        S_arg[loc_in_sorted_S] = k

        k += 1

    testspace_en = k  
    

    return g_max_index, S, P, S_arg, 0
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_fast_all_64(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500, fast_t = 5):

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

        gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))

            
    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    S_i = np.empty(n, np.int64)
    P = np.empty(n, np.int64)

    
    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    readlength = one_mapinfo[-1][0] + 1000
    target_arr = np.zeros(n, dtype = np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
        if(one_mapinfo[i][2] == 1):
            target_arr[i] = one_mapinfo[i][1] - one_mapinfo[i][0] + readlength 
        else:
            target_arr[i] = -(one_mapinfo[i][1] + one_mapinfo[i][0] + readlength) 
            

    prereadloc = one_mapinfo[0][0]

    
    testspace_en = 1
    testspace_en_i = 1
    
    
    S_arg_i = np.empty(n, np.int64)
    S_arg_i[0] = 0
    
    i = 0

    S[i] = one_mapinfo[i][3]
    S_i[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    

    fast = False
    
    S_i_count = np.zeros(one_mapinfo[-1][0] + 50, dtype = np.int64)
    S_i_count[one_mapinfo[i][3]] = 1
    
    max_score_i = 0
    

    readpos_arr = one_mapinfo[:,0]
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])

        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        

        if(prereadloc < one_mapinfo[i][0]):                

                
            k = testspace_en_i
            while(k < i):
                
                S_i_count[S_i[k]] += 1
                if(S_i[k] > max_score_i):
                    max_score_i = S_i[k]
                
                loc_in_sorted_S = insertpoint_score_distance(S_i, S_i[k], k, S_arg_i, target_arr[k], target_arr)

                S_arg_i[loc_in_sorted_S + 1: k + 1] = S_arg_i[loc_in_sorted_S: k]
                S_arg_i[loc_in_sorted_S] = k
                
                k += 1

            testspace_en_i = k
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)

            
            prereadloc = one_mapinfo[i][0]
            
            
        c_score_i = max_score_i     
        st_loc = testspace_en_i
        en_loc = testspace_en_i
        f_kmersize = one_mapinfo[i][3] + 1
        while(c_score_i > (max_scores - f_kmersize)):
            now_count = S_i_count[c_score_i]
            if(now_count == 0):
                c_score_i -= 1
                continue
            st_loc = en_loc - now_count
            '''tmp_set = set()
            tmp_set.add(-1)
            if(st_loc < 0 or en_loc <= 0):
                print('st_loc, en_loc', st_loc, en_loc)
            for tmp_item in S[S_arg[st_loc: en_loc]]:
                tmp_set.add(tmp_item//1)
            if(len(tmp_set) != 2):
                print(S[S_arg[st_loc: en_loc]])'''
            if(now_count > fast_t):
                

                #S_arg_j = closest2targetdistance(target_arr[i], target_arr, S_arg_i, st_loc, en_loc)
                #S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg_i, one_mapinfo[i][0], readpos_arr)
                j = S_arg_i[closest2targetdistance(target_arr[i], target_arr, S_arg_i, st_loc, en_loc)]
    

                #########
                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                if((readgap < 0)):
                    bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                    readgap = 0
                    overlap_size = one_mapinfo[j][0] +one_mapinfo[j][3] - one_mapinfo[i][0]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = one_mapinfo[i][1] + overlap_size - (one_mapinfo[j][1] + one_mapinfo[j][3])
                        else:
                            refgap = one_mapinfo[j][1] - (one_mapinfo[i][1] + bonus)
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] + overlap_size - one_mapinfo[j][1] + 1
                        else:
                            refgap = one_mapinfo[i][1] + bonus - 1 - (one_mapinfo[j][1] + one_mapinfo[j][3])

                else:
                    bonus = one_mapinfo[i][3]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                        else:
                            refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] - one_mapinfo[j][1] + 1
                        else:
                            refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - 1 - one_mapinfo[j][1] - one_mapinfo[j][3]


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
                    

            else:
                tmp_j = en_loc - 1
                while(tmp_j >= st_loc):
                    

                    j = S_arg_i[tmp_j]
                    tmp_j -= 1
                #for j in S_arg[st_loc: en_loc][::-1]:
            
                    #S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg_i, one_mapinfo[i][0], readpos_arr)

                    readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                    if((readgap < 0)):
                        bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                        readgap = 0
                        overlap_size = one_mapinfo[j][0] +one_mapinfo[j][3] - one_mapinfo[i][0]
                        if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                            if(one_mapinfo[i][2] == 1):
                                refgap = one_mapinfo[i][1] + overlap_size - (one_mapinfo[j][1] + one_mapinfo[j][3])
                            else:
                                refgap = one_mapinfo[j][1] - (one_mapinfo[i][1] + bonus)
                        else:
                            if(one_mapinfo[j][2] == -1):
                                refgap = one_mapinfo[i][1] + overlap_size - one_mapinfo[j][1] + 1
                            else:
                                refgap = one_mapinfo[i][1] + bonus - 1 - (one_mapinfo[j][1] + one_mapinfo[j][3])

                    else:
                        bonus = one_mapinfo[i][3]
                        if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                            if(one_mapinfo[i][2] == 1):
                                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                            else:
                                refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                        else:
                            if(one_mapinfo[j][2] == -1):
                                refgap = one_mapinfo[i][1] - one_mapinfo[j][1] + 1
                            else:
                                refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - 1 - one_mapinfo[j][1] - one_mapinfo[j][3]


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
            
            c_score_i -= 1

            
        S[i] = max_scores
        S_i[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i

        
    k = testspace_en_i
    while(k < n):

        S_i_count[S_i[k]] += 1
        if(S_i[k] > max_score_i):
            max_score_i = S_i[k]

        loc_in_sorted_S = insertpoint_score_distance(S_i, S_i[k], k, S_arg_i, target_arr[k], target_arr)

        S_arg_i[loc_in_sorted_S + 1: k + 1] = S_arg_i[loc_in_sorted_S: k]
        S_arg_i[loc_in_sorted_S] = k

        k += 1

    
    return g_max_index, S, P, S_arg_i
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_d_all(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500):

    extra_size = len(extra) - 1
    
    oskipcost = skipcost#mark_a
    omaxdiff = maxdiff#mark_a
    repeat_weight = 20#mark_a

    
    gap_arr = np.empty(one_mapinfo[-1][0])#mark 1
    gap_arr[0] = 200#mark 1
    
    g_max_scores = 0.
    g_max_index = -1
     
    gapcost_list = np.zeros(maxdiff + 1, dtype = np.float64)
    for gapcost in range(1, maxdiff + 1):

        gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))

    
    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int32)


    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)

            
    prereadloc = one_mapinfo[0][0]
    skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)

    
    testspace_en = 1
    
    
    S_arg = np.empty(n, np.int32)
    S_arg[0] = 0
    
    i = 0

    S[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    

    opcount = 0

    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])
        

        max_scores = one_mapinfo[i][3]

        
        pre_index = -9999999
        

        if(prereadloc < one_mapinfo[i][0]):
            
            if((opcount/i) > max_factor):#mark 1
                return -1, S, P, S_arg, -1.#mark 1
            
            
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = insertpoint_score(S, S[k], k, S_arg)

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0]
            
               
        for j in S_arg[:testspace_en][::-1]:
            
            
            if(S[j] > (max_scores - one_mapinfo[i][3])):
                
                opcount += 1

                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                if((readgap < 0)):
                    bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                    readgap = 0
                    overlap_size = one_mapinfo[j][0] +one_mapinfo[j][3] - one_mapinfo[i][0]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = one_mapinfo[i][1] + overlap_size - (one_mapinfo[j][1] + one_mapinfo[j][3])
                        else:
                            refgap = one_mapinfo[j][1] - (one_mapinfo[i][1] + bonus)
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] + overlap_size - one_mapinfo[j][1] + 1
                        else:
                            refgap = one_mapinfo[i][1] + bonus - 1 - (one_mapinfo[j][1] + one_mapinfo[j][3])

                else:
                    bonus = one_mapinfo[i][3]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                        else:
                            refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] - one_mapinfo[j][1] + 1
                        else:
                            refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - 1 - one_mapinfo[j][1] - one_mapinfo[j][3]


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
            else:
                break
                

        S[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i
       

    k = testspace_en
    while(k < n):

        loc_in_sorted_S = insertpoint_score(S, S[k], k, S_arg)

        S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
        S_arg[loc_in_sorted_S] = k

        k += 1

    testspace_en = k  
    

    return g_max_index, S, P, S_arg, 0
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

        gapcost_list[gapcost] = (0.01 * kmersize * gapcost + 0.5 * np.log2(gapcost))

            
    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    S_i = np.empty(n, np.int64)
    P = np.empty(n, np.int32)

    
    coverage_dict = np.zeros(one_mapinfo[-1][0] + 1, np.int64)
    readlength = one_mapinfo[-1][0] + 1000
    target_arr = np.zeros(n, dtype = np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
        if(one_mapinfo[i][2] == 1):
            target_arr[i] = one_mapinfo[i][1] - one_mapinfo[i][0] + readlength 
        else:
            target_arr[i] = -(one_mapinfo[i][1] + one_mapinfo[i][0] + readlength) 
            

    prereadloc = one_mapinfo[0][0]

    
    testspace_en = 1
    testspace_en_i = 1
    
    
    S_arg_i = np.empty(n, np.int32)
    S_arg_i[0] = 0
    
    i = 0

    S[i] = one_mapinfo[i][3]
    S_i[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    

    fast = False
    
    S_i_count = np.zeros(one_mapinfo[-1][0] + 50, dtype = np.int64)
    S_i_count[one_mapinfo[i][3]] = 1
    
    max_score_i = 0
    

    readpos_arr = one_mapinfo[:,0]
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])

        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        

        if(prereadloc < one_mapinfo[i][0]):                

                
            k = testspace_en_i
            while(k < i):
                
                S_i_count[S_i[k]] += 1
                if(S_i[k] > max_score_i):
                    max_score_i = S_i[k]
                
                loc_in_sorted_S = insertpoint_score_distance(S_i, S_i[k], k, S_arg_i, target_arr[k], target_arr)

                S_arg_i[loc_in_sorted_S + 1: k + 1] = S_arg_i[loc_in_sorted_S: k]
                S_arg_i[loc_in_sorted_S] = k
                
                k += 1

            testspace_en_i = k
            
            
            gap_arr[prereadloc] = maxgap + 2*(one_mapinfo[i][0] - prereadloc)#mark 1
            
            skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)

            
            prereadloc = one_mapinfo[i][0]
            
            
        c_score_i = max_score_i     
        st_loc = testspace_en_i
        en_loc = testspace_en_i
        f_kmersize = one_mapinfo[i][3] + 1
        while(c_score_i > (max_scores - f_kmersize)):
            now_count = S_i_count[c_score_i]
            if(now_count == 0):
                c_score_i -= 1
                continue
            st_loc = en_loc - now_count
            '''tmp_set = set()
            tmp_set.add(-1)
            if(st_loc < 0 or en_loc <= 0):
                print('st_loc, en_loc', st_loc, en_loc)
            for tmp_item in S[S_arg[st_loc: en_loc]]:
                tmp_set.add(tmp_item//1)
            if(len(tmp_set) != 2):
                print(S[S_arg[st_loc: en_loc]])'''
            if(now_count > fast_t):
                

                #S_arg_j = closest2targetdistance(target_arr[i], target_arr, S_arg_i, st_loc, en_loc)
                #S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg_i, one_mapinfo[i][0], readpos_arr)
                j = S_arg_i[closest2targetdistance(target_arr[i], target_arr, S_arg_i, st_loc, en_loc)]
    

                #########
                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                if((readgap < 0)):
                    bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                    readgap = 0
                    overlap_size = one_mapinfo[j][0] +one_mapinfo[j][3] - one_mapinfo[i][0]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = one_mapinfo[i][1] + overlap_size - (one_mapinfo[j][1] + one_mapinfo[j][3])
                        else:
                            refgap = one_mapinfo[j][1] - (one_mapinfo[i][1] + bonus)
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] + overlap_size - one_mapinfo[j][1] + 1
                        else:
                            refgap = one_mapinfo[i][1] + bonus - 1 - (one_mapinfo[j][1] + one_mapinfo[j][3])

                else:
                    bonus = one_mapinfo[i][3]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                        else:
                            refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] - one_mapinfo[j][1] + 1
                        else:
                            refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - 1 - one_mapinfo[j][1] - one_mapinfo[j][3]


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
                    

            else:
                tmp_j = en_loc - 1
                while(tmp_j >= st_loc):
                    

                    j = S_arg_i[tmp_j]
                    tmp_j -= 1
                #for j in S_arg[st_loc: en_loc][::-1]:
            
                    #S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg_i, one_mapinfo[i][0], readpos_arr)

                    readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                    if((readgap < 0)):
                        bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                        readgap = 0
                        overlap_size = one_mapinfo[j][0] +one_mapinfo[j][3] - one_mapinfo[i][0]
                        if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                            if(one_mapinfo[i][2] == 1):
                                refgap = one_mapinfo[i][1] + overlap_size - (one_mapinfo[j][1] + one_mapinfo[j][3])
                            else:
                                refgap = one_mapinfo[j][1] - (one_mapinfo[i][1] + bonus)
                        else:
                            if(one_mapinfo[j][2] == -1):
                                refgap = one_mapinfo[i][1] + overlap_size - one_mapinfo[j][1] + 1
                            else:
                                refgap = one_mapinfo[i][1] + bonus - 1 - (one_mapinfo[j][1] + one_mapinfo[j][3])

                    else:
                        bonus = one_mapinfo[i][3]
                        if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                            if(one_mapinfo[i][2] == 1):
                                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                            else:
                                refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                        else:
                            if(one_mapinfo[j][2] == -1):
                                refgap = one_mapinfo[i][1] - one_mapinfo[j][1] + 1
                            else:
                                refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - 1 - one_mapinfo[j][1] - one_mapinfo[j][3]


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
            
            c_score_i -= 1

            
        S[i] = max_scores
        S_i[i] = max_scores
        P[i] = pre_index

        if(max_scores > g_max_scores):

            g_max_scores = max_scores
            g_max_index = i

        
    k = testspace_en_i
    while(k < n):

        S_i_count[S_i[k]] += 1
        if(S_i[k] > max_score_i):
            max_score_i = S_i[k]

        loc_in_sorted_S = insertpoint_score_distance(S_i, S_i[k], k, S_arg_i, target_arr[k], target_arr)

        S_arg_i[loc_in_sorted_S + 1: k + 1] = S_arg_i[loc_in_sorted_S: k]
        S_arg_i[loc_in_sorted_S] = k

        k += 1

    
    return g_max_index, S, P, S_arg_i
#fast in second chain


###############

readgapcost_list = np.zeros(100, dtype = np.float32)
for readgapcost in range(1, 100):
    readgapcost_list[readgapcost] = 0.1 * np.log2(readgapcost + 1)
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_fast_64(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500, fast_t = 5):

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
    S_i = np.empty(n, np.int64)
    P = np.empty(n, np.int64)

    
    coverage_dict = np.zeros(one_mapinfo[-1][0] + 5000, np.int64)
    readlength = one_mapinfo[-1][0] + 1000
    target_arr = np.zeros(n, dtype = np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
        if(one_mapinfo[i][2] == 1):
            target_arr[i] = one_mapinfo[i][1] - one_mapinfo[i][0] + readlength 
        else:
            target_arr[i] = -(one_mapinfo[i][1] + one_mapinfo[i][0] + readlength) 
            

    prereadloc = one_mapinfo[0][0] + one_mapinfo[0][3]

    
    testspace_en = 1
    testspace_en_i = 1
    
    
    S_arg_i = np.empty(n, np.int64)
    S_arg_i[0] = 0
    
    i = 0

    S[i] = one_mapinfo[i][3]
    S_i[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    

    fast = False
    
    S_i_count = np.zeros(one_mapinfo[-1][0] + 50, dtype = np.int64)
    S_i_count[one_mapinfo[i][3]] = 1
    
    max_score_i = 0
    

    readpos_arr = one_mapinfo[:,0]
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])

        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        

        if(prereadloc < one_mapinfo[i][0] + one_mapinfo[i][3]):                

                
            k = testspace_en_i
            while(k < i):
                
                S_i_count[S_i[k]] += 1
                if(S_i[k] > max_score_i):
                    max_score_i = S_i[k]
                
                loc_in_sorted_S = insertpoint_score_distance(S_i, S_i[k], k, S_arg_i, target_arr[k], target_arr)

                S_arg_i[loc_in_sorted_S + 1: k + 1] = S_arg_i[loc_in_sorted_S: k]
                S_arg_i[loc_in_sorted_S] = k
                
                k += 1

            testspace_en_i = k
            
            
            #skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            #maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)

            
            prereadloc = one_mapinfo[i][0] + one_mapinfo[i][3]
            
            
        c_score_i = max_score_i     
        st_loc = testspace_en_i
        en_loc = testspace_en_i
        f_kmersize = one_mapinfo[i][3] + 1
        while(c_score_i > (max_scores - f_kmersize)):
            now_count = S_i_count[c_score_i]
            if(now_count == 0):
                c_score_i -= 1
                continue
            st_loc = en_loc - now_count
            '''tmp_set = set()
            tmp_set.add(-1)
            if(st_loc < 0 or en_loc <= 0):
                print('st_loc, en_loc', st_loc, en_loc)
            for tmp_item in S[S_arg[st_loc: en_loc]]:
                tmp_set.add(tmp_item//1)
            if(len(tmp_set) != 2):
                print(S[S_arg[st_loc: en_loc]])'''
            if(now_count > fast_t):
                

                #S_arg_j = closest2targetdistance(target_arr[i], target_arr, S_arg_i, st_loc, en_loc)
                #S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg_i, one_mapinfo[i][0], readpos_arr)
                j = S_arg_i[closest2targetdistance(target_arr[i], target_arr, S_arg_i, st_loc, en_loc)]
    

                #########
                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                ####
                if((readgap < 0)):
                    bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                    if(bonus <= 0):
                        continue
                    readgap = 0
                    overlap_size = one_mapinfo[j][0] +one_mapinfo[j][3] - one_mapinfo[i][0]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = one_mapinfo[i][1] + overlap_size - (one_mapinfo[j][1] + one_mapinfo[j][3])
                        else:
                            refgap = one_mapinfo[j][1] - (one_mapinfo[i][1] + bonus)
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] + overlap_size - one_mapinfo[j][1] + 1
                        else:
                            refgap = one_mapinfo[i][1] + bonus - 1 - (one_mapinfo[j][1] + one_mapinfo[j][3])

                else:
                    bonus = one_mapinfo[i][3]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                        else:
                            refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] - one_mapinfo[j][1] + 1
                        else:
                            refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - 1 - one_mapinfo[j][1] - one_mapinfo[j][3]
                ####

                gapcost = abs(readgap - refgap)
            

############################################################################
                if(one_mapinfo[i][2] == one_mapinfo[j][2] and refgap >= 0 and readgap <= maxgap and gapcost <= maxdiff):


                    test_scores = S[j] + bonus - gapcost_list[gapcost] - readgapcost_list[readgap]

                    if(test_scores > max_scores):

                        max_scores = test_scores
                        pre_index = j
                        

                else:
                    if(gapcost > extra_size):
                        gapcost = extra_size


                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                        tmp_penalty = min(50, skipcost) + extra[gapcost]
                        test_scores = S[j] + bonus - tmp_penalty
                    else:
                        tmp_penalty = skipcost + extra[gapcost]
                        test_scores = S[j] + bonus - tmp_penalty


                    if(test_scores > max_scores):

                        max_scores = test_scores
                        pre_index = j
                        
############################################################################
                    

            else:
                tmp_j = en_loc - 1
                while(tmp_j >= st_loc):
                    

                    j = S_arg_i[tmp_j]
                    tmp_j -= 1
                #for j in S_arg[st_loc: en_loc][::-1]:
            
                    #S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg_i, one_mapinfo[i][0], readpos_arr)

                    readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                    ####
                    if((readgap < 0)):
                        bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                        if(bonus <= 0):
                            continue
                        readgap = 0
                        overlap_size = one_mapinfo[j][0] +one_mapinfo[j][3] - one_mapinfo[i][0]
                        if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                            if(one_mapinfo[i][2] == 1):
                                refgap = one_mapinfo[i][1] + overlap_size - (one_mapinfo[j][1] + one_mapinfo[j][3])
                            else:
                                refgap = one_mapinfo[j][1] - (one_mapinfo[i][1] + bonus)
                        else:
                            if(one_mapinfo[j][2] == -1):
                                refgap = one_mapinfo[i][1] + overlap_size - one_mapinfo[j][1] + 1
                            else:
                                refgap = one_mapinfo[i][1] + bonus - 1 - (one_mapinfo[j][1] + one_mapinfo[j][3])

                    else:
                        bonus = one_mapinfo[i][3]
                        if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                            if(one_mapinfo[i][2] == 1):
                                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                            else:
                                refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                        else:
                            if(one_mapinfo[j][2] == -1):
                                refgap = one_mapinfo[i][1] - one_mapinfo[j][1] + 1
                            else:
                                refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - 1 - one_mapinfo[j][1] - one_mapinfo[j][3]
                    ####

                    gapcost = abs(readgap - refgap)


############################################################################
                    if(one_mapinfo[i][2] == one_mapinfo[j][2] and refgap >= 0 and readgap <= maxgap and gapcost <= maxdiff):


                        test_scores = S[j] + bonus - gapcost_list[gapcost] - readgapcost_list[readgap]

                        if(test_scores > max_scores):

                            max_scores = test_scores
                            pre_index = j


                    else:
                        if(gapcost > extra_size):
                            gapcost = extra_size

            
                        if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                            tmp_penalty = min(50, skipcost) + extra[gapcost]
                            test_scores = S[j] + bonus - tmp_penalty
                        else:
                            tmp_penalty = skipcost + extra[gapcost]
                            test_scores = S[j] + bonus - tmp_penalty


                        if(test_scores > max_scores):

                            max_scores = test_scores
                            pre_index = j

############################################################################
                        
                        
            en_loc = st_loc
            
            c_score_i -= 1

            
        S[i] = max_scores
        S_i[i] = max_scores
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
        if(preitem[0] < (nowitem[0] + nowitem[3])):
            overlap_size = nowitem[0] + nowitem[3] - preitem[0]
            if(preitem[2] == 1):
                path[-1] = ((preitem[0] + overlap_size, preitem[1] + overlap_size, preitem[2], preitem[3] - overlap_size))
            else:
                path[-1] = ((preitem[0] + overlap_size, preitem[1], preitem[2], preitem[3] - overlap_size))
        path.append((nowitem[0], nowitem[1], nowitem[2], nowitem[3]))

                
        preitem = nowitem
    #print('fast', g_max_scores)
    return g_max_scores, path
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_fast(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500, fast_t = 5):

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
    S_i = np.empty(n, np.int64)
    P = np.empty(n, np.int32)

    
    coverage_dict = np.zeros(one_mapinfo[-1][0] + 5000, np.int64)
    readlength = one_mapinfo[-1][0] + 1000
    target_arr = np.zeros(n, dtype = np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
        if(one_mapinfo[i][2] == 1):
            target_arr[i] = one_mapinfo[i][1] - one_mapinfo[i][0] + readlength 
        else:
            target_arr[i] = -(one_mapinfo[i][1] + one_mapinfo[i][0] + readlength) 
            

    prereadloc = one_mapinfo[0][0] + one_mapinfo[0][3]

    
    testspace_en = 1
    testspace_en_i = 1
    
    
    S_arg_i = np.empty(n, np.int32)
    S_arg_i[0] = 0
    
    i = 0

    S[i] = one_mapinfo[i][3]
    S_i[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    

    fast = False
    
    S_i_count = np.zeros(one_mapinfo[-1][0] + 50, dtype = np.int64)
    S_i_count[one_mapinfo[i][3]] = 1
    
    max_score_i = 0
    

    readpos_arr = one_mapinfo[:,0]
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])

        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        

        if(prereadloc < one_mapinfo[i][0] + one_mapinfo[i][3]):                

                
            k = testspace_en_i
            while(k < i):
                
                S_i_count[S_i[k]] += 1
                if(S_i[k] > max_score_i):
                    max_score_i = S_i[k]
                
                loc_in_sorted_S = insertpoint_score_distance(S_i, S_i[k], k, S_arg_i, target_arr[k], target_arr)

                S_arg_i[loc_in_sorted_S + 1: k + 1] = S_arg_i[loc_in_sorted_S: k]
                S_arg_i[loc_in_sorted_S] = k
                
                k += 1

            testspace_en_i = k
            
            
            #skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            #maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)

            
            prereadloc = one_mapinfo[i][0] + one_mapinfo[i][3]
            
            
        c_score_i = max_score_i     
        st_loc = testspace_en_i
        en_loc = testspace_en_i
        f_kmersize = one_mapinfo[i][3] + 1
        while(c_score_i > (max_scores - f_kmersize)):
            now_count = S_i_count[c_score_i]
            if(now_count == 0):
                c_score_i -= 1
                continue
            st_loc = en_loc - now_count
            '''tmp_set = set()
            tmp_set.add(-1)
            if(st_loc < 0 or en_loc <= 0):
                print('st_loc, en_loc', st_loc, en_loc)
            for tmp_item in S[S_arg[st_loc: en_loc]]:
                tmp_set.add(tmp_item//1)
            if(len(tmp_set) != 2):
                print(S[S_arg[st_loc: en_loc]])'''
            if(now_count > fast_t):
                

                #S_arg_j = closest2targetdistance(target_arr[i], target_arr, S_arg_i, st_loc, en_loc)
                #S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg_i, one_mapinfo[i][0], readpos_arr)
                j = S_arg_i[closest2targetdistance(target_arr[i], target_arr, S_arg_i, st_loc, en_loc)]
    

                #########
                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                ####
                if((readgap < 0)):
                    bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                    if(bonus <= 0):
                        continue
                    readgap = 0
                    overlap_size = one_mapinfo[j][0] +one_mapinfo[j][3] - one_mapinfo[i][0]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = one_mapinfo[i][1] + overlap_size - (one_mapinfo[j][1] + one_mapinfo[j][3])
                        else:
                            refgap = one_mapinfo[j][1] - (one_mapinfo[i][1] + bonus)
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] + overlap_size - one_mapinfo[j][1] + 1
                        else:
                            refgap = one_mapinfo[i][1] + bonus - 1 - (one_mapinfo[j][1] + one_mapinfo[j][3])

                else:
                    bonus = one_mapinfo[i][3]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                        else:
                            refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] - one_mapinfo[j][1] + 1
                        else:
                            refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - 1 - one_mapinfo[j][1] - one_mapinfo[j][3]
                ####

                gapcost = abs(readgap - refgap)
            

############################################################################
                if(one_mapinfo[i][2] == one_mapinfo[j][2] and refgap >= 0 and readgap <= maxgap and gapcost <= maxdiff):


                    test_scores = S[j] + bonus - gapcost_list[gapcost] - readgapcost_list[readgap]

                    if(test_scores > max_scores):

                        max_scores = test_scores
                        pre_index = j
                        

                else:
                    if(gapcost > extra_size):
                        gapcost = extra_size


                    if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                        tmp_penalty = min(50, skipcost) + extra[gapcost]
                        test_scores = S[j] + bonus - tmp_penalty
                    else:
                        tmp_penalty = skipcost + extra[gapcost]
                        test_scores = S[j] + bonus - tmp_penalty


                    if(test_scores > max_scores):

                        max_scores = test_scores
                        pre_index = j
                        
############################################################################
                    

            else:
                tmp_j = en_loc - 1
                while(tmp_j >= st_loc):
                    

                    j = S_arg_i[tmp_j]
                    tmp_j -= 1
                #for j in S_arg[st_loc: en_loc][::-1]:
            
                    #S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg_i, one_mapinfo[i][0], readpos_arr)

                    readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                    ####
                    if((readgap < 0)):
                        bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                        if(bonus <= 0):
                            continue
                        readgap = 0
                        overlap_size = one_mapinfo[j][0] +one_mapinfo[j][3] - one_mapinfo[i][0]
                        if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                            if(one_mapinfo[i][2] == 1):
                                refgap = one_mapinfo[i][1] + overlap_size - (one_mapinfo[j][1] + one_mapinfo[j][3])
                            else:
                                refgap = one_mapinfo[j][1] - (one_mapinfo[i][1] + bonus)
                        else:
                            if(one_mapinfo[j][2] == -1):
                                refgap = one_mapinfo[i][1] + overlap_size - one_mapinfo[j][1] + 1
                            else:
                                refgap = one_mapinfo[i][1] + bonus - 1 - (one_mapinfo[j][1] + one_mapinfo[j][3])

                    else:
                        bonus = one_mapinfo[i][3]
                        if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                            if(one_mapinfo[i][2] == 1):
                                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                            else:
                                refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                        else:
                            if(one_mapinfo[j][2] == -1):
                                refgap = one_mapinfo[i][1] - one_mapinfo[j][1] + 1
                            else:
                                refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - 1 - one_mapinfo[j][1] - one_mapinfo[j][3]
                    ####

                    gapcost = abs(readgap - refgap)


############################################################################
                    if(one_mapinfo[i][2] == one_mapinfo[j][2] and refgap >= 0 and readgap <= maxgap and gapcost <= maxdiff):


                        test_scores = S[j] + bonus - gapcost_list[gapcost] - readgapcost_list[readgap]

                        if(test_scores > max_scores):

                            max_scores = test_scores
                            pre_index = j


                    else:
                        if(gapcost > extra_size):
                            gapcost = extra_size

            
                        if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                            tmp_penalty = min(50, skipcost) + extra[gapcost]
                            test_scores = S[j] + bonus - tmp_penalty
                        else:
                            tmp_penalty = skipcost + extra[gapcost]
                            test_scores = S[j] + bonus - tmp_penalty


                        if(test_scores > max_scores):

                            max_scores = test_scores
                            pre_index = j

############################################################################
                        
                        
            en_loc = st_loc
            
            c_score_i -= 1

            
        S[i] = max_scores
        S_i[i] = max_scores
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
        if(preitem[0] < (nowitem[0] + nowitem[3])):
            overlap_size = nowitem[0] + nowitem[3] - preitem[0]
            if(preitem[2] == 1):
                path[-1] = ((preitem[0] + overlap_size, preitem[1] + overlap_size, preitem[2], preitem[3] - overlap_size))
            else:
                path[-1] = ((preitem[0] + overlap_size, preitem[1], preitem[2], preitem[3] - overlap_size))
        path.append((nowitem[0], nowitem[1], nowitem[2], nowitem[3]))

                
        preitem = nowitem
    #print('fast', g_max_scores)
    return g_max_scores, path
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
    
    
    coverage_dict = np.zeros(one_mapinfo[-1][0] + 5000, np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    prereadloc = one_mapinfo[0][0] + one_mapinfo[0][3]
    #skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    #maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
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
        

        if(prereadloc < (one_mapinfo[i][0] + one_mapinfo[i][3])):
            if(opcount > 100000 and opcount / prereadloc > 1000):
                if(n < np.iinfo(np.int32).max-10):
                    return get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_fast(one_mapinfo, kmersize = kmersize, skipcost = skipcost, maxdiff = maxdiff, maxgap = maxgap, fast_t = 5)
                else:
                    return get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_fast_64(one_mapinfo, kmersize = kmersize, skipcost = skipcost, maxdiff = maxdiff, maxgap = maxgap, fast_t = 5)
                    
            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i
            #c_repeat_weight = coverage_dict[one_mapinfo[i][0]]
            #skipcost = oskipcost + c_repeat_weight
            #maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0] + one_mapinfo[i][3]
            
               
        for j in S_arg[:testspace_en][::-1]:
            opcount += 1
              
            
            if(S[j] < (max_scores - one_mapinfo[i][3])):
                
                break


            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            ####
            if((readgap < 0)):
                bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                if(bonus <= 0):
                    continue
                readgap = 0
                overlap_size = one_mapinfo[j][0] +one_mapinfo[j][3] - one_mapinfo[i][0]
                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):
                        refgap = one_mapinfo[i][1] + overlap_size - (one_mapinfo[j][1] + one_mapinfo[j][3])
                    else:
                        refgap = one_mapinfo[j][1] - (one_mapinfo[i][1] + bonus)
                else:
                    if(one_mapinfo[j][2] == -1):
                        refgap = one_mapinfo[i][1] + overlap_size - one_mapinfo[j][1] + 1
                    else:
                        refgap = one_mapinfo[i][1] + bonus - 1 - (one_mapinfo[j][1] + one_mapinfo[j][3])

            else:
                bonus = one_mapinfo[i][3]
                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):
                        refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                    else:
                        refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                else:
                    if(one_mapinfo[j][2] == -1):
                        refgap = one_mapinfo[i][1] - one_mapinfo[j][1] + 1
                    else:
                        refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - 1 - one_mapinfo[j][1] - one_mapinfo[j][3]

            ####

                
            gapcost = abs(readgap - refgap)
            

            if(one_mapinfo[i][2] == one_mapinfo[j][2] and refgap >= 0 and readgap <= maxgap and gapcost <= maxdiff):


                test_scores = S[j] + bonus - gapcost_list[gapcost] - readgapcost_list[readgap]
                
                if(test_scores > max_scores):

                    max_scores = test_scores
                    pre_index = j
                    

            else:
                if(gapcost > extra_size):
                    gapcost = extra_size

                if(one_mapinfo[i][2] != one_mapinfo[j][2]):
                    tmp_penalty = min(50, skipcost) + extra[gapcost]
                    test_scores = S[j] + bonus - tmp_penalty
                else:
                    tmp_penalty = skipcost + extra[gapcost]
                    test_scores = S[j] + bonus - tmp_penalty


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
        if(preitem[0] < (nowitem[0] + nowitem[3])):
            overlap_size = nowitem[0] + nowitem[3] - preitem[0]
            if(preitem[2] == 1):
                path[-1] = ((preitem[0] + overlap_size, preitem[1] + overlap_size, preitem[2], preitem[3] - overlap_size))
            else:
                path[-1] = ((preitem[0] + overlap_size, preitem[1], preitem[2], preitem[3] - overlap_size))
        path.append((nowitem[0], nowitem[1], nowitem[2], nowitem[3]))

                
        preitem = nowitem
    #print('normal', g_max_scores)
    return g_max_scores, path
#testreadgap
log2cache = np.array([0.5 * np.log2((gapcost + 1)) for gapcost in range(100000)])
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_mismatch_fast_64(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500, large_readgap = 30, fast_t = 5):

    extra_size = len(extra) - 1
    log2cache_size = len(log2cache) - 1
    

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


    large_readgapcost_list = np.zeros(maxgap + 1, dtype=np.float32)
    for readgap in range(1, maxgap + 1):
        if(large_readgap <= readgap):
            large_readgapcost_list[readgap] = 0.5 * readgap
        else:
            large_readgapcost_list[readgap] = 0.1 * np.log2(readgap + 1)

    
    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    S_i = np.empty(n, np.int64)
    P = np.empty(n, np.int64)

    
    coverage_dict = np.zeros(one_mapinfo[-1][0] + 5000, np.int64)
    readlength = one_mapinfo[-1][0] + 1000
    target_arr = np.zeros(n, dtype = np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
        if(one_mapinfo[i][2] == 1):
            target_arr[i] = one_mapinfo[i][1] - one_mapinfo[i][0] + readlength 
        else:
            target_arr[i] = -(one_mapinfo[i][1] + one_mapinfo[i][0] + readlength) 
            

    prereadloc = one_mapinfo[0][0] + one_mapinfo[0][3]

    
    testspace_en = 1
    testspace_en_i = 1
    
    
    S_arg_i = np.empty(n, np.int64)
    S_arg_i[0] = 0
    
    i = 0

    S[i] = one_mapinfo[i][3]
    S_i[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    

    fast = False
    
    S_i_count = np.zeros(one_mapinfo[-1][0] + 50, dtype = np.int64)
    S_i_count[one_mapinfo[i][3]] = 1
    
    max_score_i = 0
    

    readpos_arr = one_mapinfo[:,0]
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])

        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        

        if(prereadloc < one_mapinfo[i][0] + one_mapinfo[i][3]):                

                
            k = testspace_en_i
            while(k < i):
                
                S_i_count[S_i[k]] += 1
                if(S_i[k] > max_score_i):
                    max_score_i = S_i[k]
                
                loc_in_sorted_S = insertpoint_score_distance(S_i, S_i[k], k, S_arg_i, target_arr[k], target_arr)

                S_arg_i[loc_in_sorted_S + 1: k + 1] = S_arg_i[loc_in_sorted_S: k]
                S_arg_i[loc_in_sorted_S] = k
                
                k += 1

            testspace_en_i = k
            
            
            #skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            #maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)

            
            prereadloc = one_mapinfo[i][0] + one_mapinfo[i][3]
            
            
        c_score_i = max_score_i     
        st_loc = testspace_en_i
        en_loc = testspace_en_i
        f_kmersize = one_mapinfo[i][3] + 1
        while(c_score_i > (max_scores - f_kmersize)):
            now_count = S_i_count[c_score_i]
            if(now_count == 0):
                c_score_i -= 1
                continue
            st_loc = en_loc - now_count
            '''tmp_set = set()
            tmp_set.add(-1)
            if(st_loc < 0 or en_loc <= 0):
                print('st_loc, en_loc', st_loc, en_loc)
            for tmp_item in S[S_arg[st_loc: en_loc]]:
                tmp_set.add(tmp_item//1)
            if(len(tmp_set) != 2):
                print(S[S_arg[st_loc: en_loc]])'''
            if(now_count > fast_t):
                

                #S_arg_j = closest2targetdistance(target_arr[i], target_arr, S_arg_i, st_loc, en_loc)
                #S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg_i, one_mapinfo[i][0], readpos_arr)
                j = S_arg_i[closest2targetdistance(target_arr[i], target_arr, S_arg_i, st_loc, en_loc)]
    

                #########
                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                ####
                if((readgap < 0)):
                    bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                    if(bonus <= 0):
                        continue
                    readgap = 0
                    overlap_size = one_mapinfo[j][0] +one_mapinfo[j][3] - one_mapinfo[i][0]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = one_mapinfo[i][1] + overlap_size - (one_mapinfo[j][1] + one_mapinfo[j][3])
                        else:
                            refgap = one_mapinfo[j][1] - (one_mapinfo[i][1] + bonus)
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] + overlap_size - one_mapinfo[j][1] + 1
                        else:
                            refgap = one_mapinfo[i][1] + bonus - 1 - (one_mapinfo[j][1] + one_mapinfo[j][3])

                else:
                    bonus = one_mapinfo[i][3]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                        else:
                            refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] - one_mapinfo[j][1] + 1
                        else:
                            refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - 1 - one_mapinfo[j][1] - one_mapinfo[j][3]
                ####

                gapcost = abs(readgap - refgap)
            

############################################################################
                if(one_mapinfo[i][2] == one_mapinfo[j][2] and refgap >= 0 and readgap <= maxgap and gapcost <= maxdiff):


                    test_scores = S[j] + bonus - gapcost_list[gapcost] - large_readgapcost_list[readgap]

                    if(test_scores > max_scores):

                        max_scores = test_scores
                        pre_index = j


                else:

                    tmp_penalty = skipcost + log2cache[min(log2cache_size, gapcost)]#+ extra[gapcost]
                    test_scores = S[j] + bonus - tmp_penalty


                    if(test_scores > max_scores):

                        max_scores = test_scores
                        pre_index = j
                        
############################################################################
                    

            else:
                tmp_j = en_loc - 1
                while(tmp_j >= st_loc):
                    

                    j = S_arg_i[tmp_j]
                    tmp_j -= 1
                #for j in S_arg[st_loc: en_loc][::-1]:
            
                    #S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg_i, one_mapinfo[i][0], readpos_arr)

                    readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                    ####
                    if((readgap < 0)):
                        bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                        if(bonus <= 0):
                            continue
                        readgap = 0
                        overlap_size = one_mapinfo[j][0] +one_mapinfo[j][3] - one_mapinfo[i][0]
                        if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                            if(one_mapinfo[i][2] == 1):
                                refgap = one_mapinfo[i][1] + overlap_size - (one_mapinfo[j][1] + one_mapinfo[j][3])
                            else:
                                refgap = one_mapinfo[j][1] - (one_mapinfo[i][1] + bonus)
                        else:
                            if(one_mapinfo[j][2] == -1):
                                refgap = one_mapinfo[i][1] + overlap_size - one_mapinfo[j][1] + 1
                            else:
                                refgap = one_mapinfo[i][1] + bonus - 1 - (one_mapinfo[j][1] + one_mapinfo[j][3])

                    else:
                        bonus = one_mapinfo[i][3]
                        if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                            if(one_mapinfo[i][2] == 1):
                                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                            else:
                                refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                        else:
                            if(one_mapinfo[j][2] == -1):
                                refgap = one_mapinfo[i][1] - one_mapinfo[j][1] + 1
                            else:
                                refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - 1 - one_mapinfo[j][1] - one_mapinfo[j][3]
                    ####

                    gapcost = abs(readgap - refgap)


############################################################################
                    if(one_mapinfo[i][2] == one_mapinfo[j][2] and refgap >= 0 and readgap <= maxgap and gapcost <= maxdiff):


                        test_scores = S[j] + bonus - gapcost_list[gapcost] - large_readgapcost_list[readgap]

                        if(test_scores > max_scores):

                            max_scores = test_scores
                            pre_index = j


                    else:

                        tmp_penalty = skipcost + log2cache[min(log2cache_size, gapcost)]#+ extra[gapcost]
                        test_scores = S[j] + bonus - tmp_penalty


                        if(test_scores > max_scores):

                            max_scores = test_scores
                            pre_index = j

############################################################################
                        
                        
            en_loc = st_loc
            
            c_score_i -= 1

            
        S[i] = max_scores
        S_i[i] = max_scores
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
        if(preitem[0] < (nowitem[0] + nowitem[3])):
            overlap_size = nowitem[0] + nowitem[3] - preitem[0]
            if(preitem[2] == 1):
                path[-1] = ((preitem[0] + overlap_size, preitem[1] + overlap_size, preitem[2], preitem[3] - overlap_size))
            else:
                path[-1] = ((preitem[0] + overlap_size, preitem[1], preitem[2], preitem[3] - overlap_size))
        path.append((nowitem[0], nowitem[1], nowitem[2], nowitem[3]))

                
        preitem = nowitem
    #print('fast', g_max_scores)
    return g_max_scores, path
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_mismatch_fast(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500, large_readgap = 30, fast_t = 5):

    extra_size = len(extra) - 1
    log2cache_size = len(log2cache) - 1
    

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

    large_readgapcost_list = np.zeros(maxgap + 1, dtype=np.float32)
    for readgap in range(1, maxgap + 1):
        if(large_readgap <= readgap):
            large_readgapcost_list[readgap] = 0.5 * readgap
        else:
            large_readgapcost_list[readgap] = 0.1 * np.log2(readgap + 1)
            

    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    S_i = np.empty(n, np.int64)
    P = np.empty(n, np.int32)

    
    coverage_dict = np.zeros(one_mapinfo[-1][0] + 5000, np.int64)
    readlength = one_mapinfo[-1][0] + 1000
    target_arr = np.zeros(n, dtype = np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
        if(one_mapinfo[i][2] == 1):
            target_arr[i] = one_mapinfo[i][1] - one_mapinfo[i][0] + readlength 
        else:
            target_arr[i] = -(one_mapinfo[i][1] + one_mapinfo[i][0] + readlength) 
            

    prereadloc = one_mapinfo[0][0] + one_mapinfo[0][3]

    
    testspace_en = 1
    testspace_en_i = 1
    
    
    S_arg_i = np.empty(n, np.int32)
    S_arg_i[0] = 0
    
    i = 0

    S[i] = one_mapinfo[i][3]
    S_i[i] = one_mapinfo[i][3]
    P[i] = -9999999
    g_max_scores = one_mapinfo[i][3]
    g_max_index = i
    

    fast = False
    
    S_i_count = np.zeros(one_mapinfo[-1][0] + 50, dtype = np.int64)
    S_i_count[one_mapinfo[i][3]] = 1
    
    max_score_i = 0
    

    readpos_arr = one_mapinfo[:,0]
    
    for i in range(1, n):
        #print_log('start: S[0:'+str(i)+']', S[0:i])

        P[i] = -9999999
        max_scores = one_mapinfo[i][3]
        tmp_target_score = max_scores
        
        pre_index = -9999999
        

        if(prereadloc < one_mapinfo[i][0] + one_mapinfo[i][3]):                

                
            k = testspace_en_i
            while(k < i):
                
                S_i_count[S_i[k]] += 1
                if(S_i[k] > max_score_i):
                    max_score_i = S_i[k]
                
                loc_in_sorted_S = insertpoint_score_distance(S_i, S_i[k], k, S_arg_i, target_arr[k], target_arr)

                S_arg_i[loc_in_sorted_S + 1: k + 1] = S_arg_i[loc_in_sorted_S: k]
                S_arg_i[loc_in_sorted_S] = k
                
                k += 1

            testspace_en_i = k
            
            
            #skipcost = oskipcost + coverage_dict[one_mapinfo[i][0]]
            #maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)

            
            prereadloc = one_mapinfo[i][0] + one_mapinfo[i][3]
            
            
        c_score_i = max_score_i     
        st_loc = testspace_en_i
        en_loc = testspace_en_i
        f_kmersize = one_mapinfo[i][3] + 1
        while(c_score_i > (max_scores - f_kmersize)):
            now_count = S_i_count[c_score_i]
            if(now_count == 0):
                c_score_i -= 1
                continue
            st_loc = en_loc - now_count
            '''tmp_set = set()
            tmp_set.add(-1)
            if(st_loc < 0 or en_loc <= 0):
                print('st_loc, en_loc', st_loc, en_loc)
            for tmp_item in S[S_arg[st_loc: en_loc]]:
                tmp_set.add(tmp_item//1)
            if(len(tmp_set) != 2):
                print(S[S_arg[st_loc: en_loc]])'''
            if(now_count > fast_t):
                

                #S_arg_j = closest2targetdistance(target_arr[i], target_arr, S_arg_i, st_loc, en_loc)
                #S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg_i, one_mapinfo[i][0], readpos_arr)
                j = S_arg_i[closest2targetdistance(target_arr[i], target_arr, S_arg_i, st_loc, en_loc)]
    

                #########
                readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                ####
                if((readgap < 0)):
                    bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                    if(bonus <= 0):
                        continue
                    readgap = 0
                    overlap_size = one_mapinfo[j][0] +one_mapinfo[j][3] - one_mapinfo[i][0]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = one_mapinfo[i][1] + overlap_size - (one_mapinfo[j][1] + one_mapinfo[j][3])
                        else:
                            refgap = one_mapinfo[j][1] - (one_mapinfo[i][1] + bonus)
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] + overlap_size - one_mapinfo[j][1] + 1
                        else:
                            refgap = one_mapinfo[i][1] + bonus - 1 - (one_mapinfo[j][1] + one_mapinfo[j][3])

                else:
                    bonus = one_mapinfo[i][3]
                    if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                        if(one_mapinfo[i][2] == 1):
                            refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                        else:
                            refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                    else:
                        if(one_mapinfo[j][2] == -1):
                            refgap = one_mapinfo[i][1] - one_mapinfo[j][1] + 1
                        else:
                            refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - 1 - one_mapinfo[j][1] - one_mapinfo[j][3]
                ####

                gapcost = abs(readgap - refgap)
            

############################################################################
                if(one_mapinfo[i][2] == one_mapinfo[j][2] and refgap >= 0 and readgap <= maxgap and gapcost <= maxdiff):


                    test_scores = S[j] + bonus - gapcost_list[gapcost] - large_readgapcost_list[readgap]

                    if(test_scores > max_scores):

                        max_scores = test_scores
                        pre_index = j


                else:

                    tmp_penalty = skipcost + log2cache[min(log2cache_size, gapcost)]#+ extra[gapcost]
                    test_scores = S[j] + bonus - tmp_penalty


                    if(test_scores > max_scores):

                        max_scores = test_scores
                        pre_index = j
                        
############################################################################
                    

            else:
                tmp_j = en_loc - 1
                while(tmp_j >= st_loc):
                    

                    j = S_arg_i[tmp_j]
                    tmp_j -= 1
                #for j in S_arg[st_loc: en_loc][::-1]:
            
                    #S_arg_j = closest2target_1d_point_pos(target_arr, target_arr[i], st_loc, en_loc, S_arg_i, one_mapinfo[i][0], readpos_arr)

                    readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

                    ####
                    if((readgap < 0)):
                        bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                        if(bonus <= 0):
                            continue
                        readgap = 0
                        overlap_size = one_mapinfo[j][0] +one_mapinfo[j][3] - one_mapinfo[i][0]
                        if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                            if(one_mapinfo[i][2] == 1):
                                refgap = one_mapinfo[i][1] + overlap_size - (one_mapinfo[j][1] + one_mapinfo[j][3])
                            else:
                                refgap = one_mapinfo[j][1] - (one_mapinfo[i][1] + bonus)
                        else:
                            if(one_mapinfo[j][2] == -1):
                                refgap = one_mapinfo[i][1] + overlap_size - one_mapinfo[j][1] + 1
                            else:
                                refgap = one_mapinfo[i][1] + bonus - 1 - (one_mapinfo[j][1] + one_mapinfo[j][3])

                    else:
                        bonus = one_mapinfo[i][3]
                        if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                            if(one_mapinfo[i][2] == 1):
                                refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                            else:
                                refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                        else:
                            if(one_mapinfo[j][2] == -1):
                                refgap = one_mapinfo[i][1] - one_mapinfo[j][1] + 1
                            else:
                                refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - 1 - one_mapinfo[j][1] - one_mapinfo[j][3]
                    ####

                    gapcost = abs(readgap - refgap)


############################################################################
                    if(one_mapinfo[i][2] == one_mapinfo[j][2] and refgap >= 0 and readgap <= maxgap and gapcost <= maxdiff):


                        test_scores = S[j] + bonus - gapcost_list[gapcost] - large_readgapcost_list[readgap]

                        if(test_scores > max_scores):

                            max_scores = test_scores
                            pre_index = j


                    else:

                        tmp_penalty = skipcost + log2cache[min(log2cache_size, gapcost)]#+ extra[gapcost]
                        test_scores = S[j] + bonus - tmp_penalty


                        if(test_scores > max_scores):

                            max_scores = test_scores
                            pre_index = j

############################################################################
                        
                        
            en_loc = st_loc
            
            c_score_i -= 1

            
        S[i] = max_scores
        S_i[i] = max_scores
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
        if(preitem[0] < (nowitem[0] + nowitem[3])):
            overlap_size = nowitem[0] + nowitem[3] - preitem[0]
            if(preitem[2] == 1):
                path[-1] = ((preitem[0] + overlap_size, preitem[1] + overlap_size, preitem[2], preitem[3] - overlap_size))
            else:
                path[-1] = ((preitem[0] + overlap_size, preitem[1], preitem[2], preitem[3] - overlap_size))
        path.append((nowitem[0], nowitem[1], nowitem[2], nowitem[3]))

                
        preitem = nowitem
    #print('fast', g_max_scores)
    return g_max_scores, path
@njit            #one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500
def get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_mismatch(one_mapinfo, kmersize = 15, skipcost = 50., maxdiff = 30, maxgap = 500, large_readgap = 30):

    extra_size = len(extra) - 1
    log2cache_size = len(log2cache) - 1
    
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

    large_readgapcost_list = np.zeros(maxgap + 1, dtype=np.float32)
    for readgap in range(1, maxgap + 1):
        if(large_readgap <= readgap):
            large_readgapcost_list[readgap] = 0.5 * readgap
        else:
            large_readgapcost_list[readgap] = 0.1 * np.log2(readgap + 1)
    

    n = len(one_mapinfo)
    S = np.empty(n, np.float64)
    P = np.empty(n, np.int64)

  
    opcount = 0
    
    
    coverage_dict = np.zeros(one_mapinfo[-1][0] + 5000, np.int64)
    for i in range(n):

        coverage_dict[one_mapinfo[i][0]] = min(coverage_dict[one_mapinfo[i][0]]+1, repeat_weight)
            

    prereadloc = one_mapinfo[0][0] + one_mapinfo[0][3]
    #skipcost = oskipcost + coverage_dict[one_mapinfo[0][0]]
    #maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[0][0]], 10)
    
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
        

        if(prereadloc < (one_mapinfo[i][0] + one_mapinfo[i][3])):
            if(opcount > 100000 and opcount / prereadloc > 1000):
                if(n < np.iinfo(np.int32).max-10):
                    return get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_mismatch_fast(one_mapinfo, kmersize = kmersize, skipcost = skipcost, maxdiff = maxdiff, maxgap = maxgap, large_readgap = large_readgap, fast_t = 5)
                else:
                    return get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_mismatch_fast_64(one_mapinfo, kmersize = kmersize, skipcost = skipcost, maxdiff = maxdiff, maxgap = maxgap, large_readgap = large_readgap, fast_t = 5)
                

            k = testspace_en
            while(k < i):
                
                loc_in_sorted_S = smallorequal2target_1d_point(S, S[k], k, S_arg) + 1

                S_arg[loc_in_sorted_S + 1: k + 1] = S_arg[loc_in_sorted_S: k]
                S_arg[loc_in_sorted_S] = k
                
                k += 1
            
            testspace_en = i
            #c_repeat_weight = coverage_dict[one_mapinfo[i][0]]
            #skipcost = oskipcost + c_repeat_weight
            #maxdiff = max(omaxdiff - coverage_dict[one_mapinfo[i][0]], 10)
            
            prereadloc = one_mapinfo[i][0] + one_mapinfo[i][3]
            
               
        for j in S_arg[:testspace_en][::-1]:
            opcount += 1


            if(S[j] < (max_scores - one_mapinfo[i][3])):

                break


            readgap = (one_mapinfo[i][0] - one_mapinfo[j][0] - one_mapinfo[j][3])

            ####
            if((readgap < 0)):
                bonus = one_mapinfo[i][0] + one_mapinfo[i][3] - one_mapinfo[j][0] - one_mapinfo[j][3]
                if(bonus <= 0):
                    continue
                readgap = 0
                overlap_size = one_mapinfo[j][0] +one_mapinfo[j][3] - one_mapinfo[i][0]
                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):
                        refgap = one_mapinfo[i][1] + overlap_size - (one_mapinfo[j][1] + one_mapinfo[j][3])
                    else:
                        refgap = one_mapinfo[j][1] - (one_mapinfo[i][1] + bonus)
                else:
                    if(one_mapinfo[j][2] == -1):
                        refgap = one_mapinfo[i][1] + overlap_size - one_mapinfo[j][1] + 1
                    else:
                        refgap = one_mapinfo[i][1] + bonus - 1 - (one_mapinfo[j][1] + one_mapinfo[j][3])

            else:
                bonus = one_mapinfo[i][3]
                if(one_mapinfo[i][2] == one_mapinfo[j][2]):
                    if(one_mapinfo[i][2] == 1):
                        refgap = (one_mapinfo[i][1] - one_mapinfo[j][1] - one_mapinfo[j][3])
                    else:
                        refgap = (one_mapinfo[j][1] - one_mapinfo[i][1] - one_mapinfo[i][3])
                else:
                    if(one_mapinfo[j][2] == -1):
                        refgap = one_mapinfo[i][1] - one_mapinfo[j][1] + 1
                    else:
                        refgap = one_mapinfo[i][1] + one_mapinfo[i][3] - 1 - one_mapinfo[j][1] - one_mapinfo[j][3]

            ####


            gapcost = abs(readgap - refgap)


            if(one_mapinfo[i][2] == one_mapinfo[j][2] and refgap >= 0 and readgap <= maxgap and gapcost <= maxdiff):


                test_scores = S[j] + bonus - gapcost_list[gapcost] - large_readgapcost_list[readgap]

                if(test_scores > max_scores):

                    max_scores = test_scores
                    pre_index = j


            else:

                tmp_penalty = skipcost + log2cache[min(log2cache_size, gapcost)]#+ extra[gapcost]
                test_scores = S[j] + bonus - tmp_penalty


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
        if(preitem[0] < (nowitem[0] + nowitem[3])):
            overlap_size = nowitem[0] + nowitem[3] - preitem[0]
            if(preitem[2] == 1):
                path[-1] = ((preitem[0] + overlap_size, preitem[1] + overlap_size, preitem[2], preitem[3] - overlap_size))
            else:
                path[-1] = ((preitem[0] + overlap_size, preitem[1], preitem[2], preitem[3] - overlap_size))
        path.append((nowitem[0], nowitem[1], nowitem[2], nowitem[3]))

                
        preitem = nowitem
    #print('normal', g_max_scores)
    return g_max_scores, path

@njit
def get_localmap_multi_all_forDP_inv_guide_list(raw_alignment_array_list, testseq, rc_testseq, contig2start, contig2seq, kmersize, skipcost, maxdiff, maxgap, shift = 1):#
    def getcounts(x):
        return 1/len(x)
    def drop_somechains(return_path_list):
        iloclist = np.zeros(len(return_path_list) - 1, dtype = np.int64)

        distance = np.ones(len(return_path_list) - 1, dtype = np.int64) * np.iinfo(np.int64).max
        strandcount = np.zeros((len(return_path_list) - 1, 2))
        chains_strandcount = np.zeros((len(return_path_list) - 1, 2))
        for item in return_path_list[0]:
            for chainiloc in range(len(return_path_list) - 1):
                chain = return_path_list[chainiloc + 1]
                if(item[0] >= chain[-1][0] and item[0] <= chain[0][0]):
                    if(item[2] == 1):
                        strandcount[chainiloc][0] += 1
                    else:
                        strandcount[chainiloc][1] += 1
                while(chain[iloclist[chainiloc]][0] > item[0]):
                    if(iloclist[chainiloc] < (len(chain) - 1)):
                        iloclist[chainiloc] += 1
                    else:
                        break

                testitem = chain[iloclist[chainiloc]]
                if(abs(item[1] - testitem[1]) < distance[chainiloc]):
                    distance[chainiloc] = abs(item[1] - testitem[1])
                #print(item, testitem)
            #print()
        for chainiloc in range(len(return_path_list) - 1):
            chain = return_path_list[chainiloc + 1]
            for item in chain:
                if(item[2] == 1):
                    chains_strandcount[chainiloc][0] += 1
                else:
                    chains_strandcount[chainiloc][1] += 1
        new_return_path_list = List([return_path_list[0]])
        iloc = 0
        for dis in distance: 
            if((strandcount[iloc][0] > strandcount[iloc][1]) and (chains_strandcount[iloc][0] > chains_strandcount[iloc][1])):
                keep = True
            elif((strandcount[iloc][0] < strandcount[iloc][1]) and (chains_strandcount[iloc][0] < chains_strandcount[iloc][1])):
                keep = True
            else:
                keep = False
            iloc += 1
            if((keep == False and (dis < 500)) or ((return_path_list[iloc][0][0] - return_path_list[iloc][-1][0]) < 100)):
                continue
            else:
                new_return_path_list.append(return_path_list[iloc])
        return new_return_path_list
    def merge_chain(ochains):
        def get_start(r_chain):
            return r_chain[-1][0]
        def get_length(r_chain):
            return len(r_chain)
        #print(len(ochains))
        chains = ochains[1:]
        chains.sort(key = get_start)
        iloc = 0
        while(iloc < len(chains) - 1):
            jloc = iloc + 1
            while(jloc < len(chains)):
                if(chains[iloc][0][0] + chains[iloc][0][3] <= chains[jloc][-1][0]):
                    if(chains[iloc][0][2] == chains[jloc][-1][2]):
                        readgap = chains[jloc][-1][0] - chains[iloc][0][0] - chains[iloc][0][3]
                        
                        if(chains[iloc][0][2] == 1):
                            refgap = chains[jloc][-1][1] - chains[iloc][0][1] - chains[iloc][0][3]
                        else:
                            refgap = chains[iloc][0][1] - chains[jloc][-1][1] - chains[jloc][-1][3]
                        if(abs(readgap - refgap) < 500):
                            #print(readgap, refgap)
                            #print(chains[iloc][0])
                            #print(chains[jloc][-1])
                            #print()

                            chains[iloc] = np.concatenate((chains[jloc], chains[iloc]))
                            #for item in chains[iloc]:
                                #print(item)
                            #print(chains[iloc])
                            chains.pop(jloc)
                            continue
                        
                jloc += 1
            iloc += 1
        chains.sort(key = get_length)
        return_chains = List([ochains[0]])
        for chain in chains:
            return_chains.append(chain)
        #print(len(return_chains))
        return return_chains
    raw_alignment_array_list = merge_chain(raw_alignment_array_list)
    one_mapinfo = [(-1, -1, -1, -1)]
    one_mapinfo.pop(0)
    raw_alignment_array_list = drop_somechains(raw_alignment_array_list)
    raw_alignment_array_list.sort(key = getcounts)
    readstart, readend = get_localmap_multi_all_forDP_inv_guide_1(-1, -1, one_mapinfo, raw_alignment_array_list[0], testseq, rc_testseq, contig2start, contig2seq, kmersize, skipcost, maxdiff, maxgap, shift = 1)
    count = 2
    if(len(raw_alignment_array_list) > 1):
        for raw_alignment_array in raw_alignment_array_list[1:]:
            tmp_readstart, tmp_readend = get_localmap_multi_all_forDP_inv_guide_1(readstart, readend, one_mapinfo, raw_alignment_array, testseq, rc_testseq, contig2start, contig2seq, kmersize, skipcost, maxdiff, maxgap, shift = 1)
            count += 1
            if(count > 5):
                break

    one_mapinfo = np.array(one_mapinfo)
    one_mapinfo = one_mapinfo[np.argsort(one_mapinfo[:, 0] + one_mapinfo[:, 3])]
    if(len(raw_alignment_array_list) > 1):
        return get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list_mismatch( one_mapinfo, kmersize = kmersize, skipcost = skipcost, maxdiff = maxdiff, maxgap = maxgap, large_readgap = 30)
    else:
        return get_optimal_chain_sortbyreadpos_forSV_inv_test_merged_fine_list(one_mapinfo, kmersize = kmersize, skipcost = skipcost, maxdiff = maxdiff, maxgap = maxgap)
