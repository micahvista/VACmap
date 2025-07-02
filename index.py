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
