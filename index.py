import mappy as mp
import sys
import os
pstring = sys.argv
if(len(pstring) != 3):
    print('Usage')
    print('python index.py reference_genome_path output_index_path')
elif(os.path.isfile(pstring[1]) == False):
    print('reference genome not find')
else:
    mp.Aligner(fn_idx_in=pstring[1], k=15, w=10, fn_idx_out=pstring[2])