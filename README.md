# VACmap
VACmap-—an advanced long-read mapping method explicitly tailored for accurately mapping complex structural variations.


Installation
------------

    #Install from github (requires Python 3.6.* or newer)
    git clone https://github.com/micahvista/VACmap.git
    cd VACmap
    conda create -n VACmap heapdict edlib pandas pysam numba biopython matplotlib cigar psutil samtools
    conda activate VACmap
    pip install setuptools==57
    cd index
    python setup.py install

Input
-----

MAMnet takes sorted and indexed alignment files in BAM format as inputs. And MAMnet has been successfully tested on PacBio CLR, PacBio HiFi (CCS) and Oxford Nanopore data and alignment files produced by the read aligners `minimap2 <https://github.com/lh3/minimap2>`_, `pbmm2 <https://github.com/PacificBiosciences/pbmm2/>`_ , `NGMLR <https://github.com/philres/ngmlr>`_, and BWA-MEM.

Output
------

MAMnet produces SV calls in the Variant Call Format (VCF).

Usage
----------------------
    Index refernece
    ----------------------
    python index.py reference_genome_path output_index_path # VACmap uses the implementation in minimap2 to build an index of the reference sequence. 
    
    Mapping long reads
    ----------------------
    python mammap.py -ref output_index_path -read /read.fasta -outputdir /outputdir -mode H or L -maxworker number_of_threads
    #-read The path of long read 
    #-outputdir The output path of 
    #-outputpath the output path of called vcf file
    #-threads the number of threads to use. (default: all available thread)
    #-step data shift size [1-200]. (default: 50)
    #-includecontig the list of contig to preform detection. (default: [], all contig are used)


Changelog
---------


Contact
-------

If you experience any problems or have suggestions please create an issue or a pull request.

Citation
---------


License
-------

The project is licensed under the GNU General Public License.
