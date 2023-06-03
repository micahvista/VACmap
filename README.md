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

VACmap takes long read and reference sequence as inputs. And VACmap has been successfully tested on PacBio CLR, PacBio HiFi (CCS) and Oxford Nanopore data. Unlike other aligners use different parameter setting for different data type, VACmap uses same parameter setting for different data type. But depend on applications, we prodive two fine tune setting for downstream structural variation analyse. -mode H is fine tune for complex structural variation or mixed(both simple and complex structural variation) analyse. -mode L is fine tune for simple structural variation analyse. Comparing to -mode H, -mode L use more stringent settings to discard suspicious alignments. However, these alignments could be correct or misaligned.

Output
------

VACmap output alignment file in compressed BAM format.

Usage
----------------------
    Index refernece
    ----------------------
    conda activate VACmap
    python index.py reference_genome_path output_index_path # VACmap uses the implementation in minimap2 to build an index of the reference sequence. 
    
    Mapping long reads
    ----------------------
    conda activate VACmap
    python mammap.py -ref output_index_path -read /read.fasta -outputdir /outputdir -mode H or L -maxworker number_of_threads
    #-ref The index file produced by minimap2. The index file can be reuse for different data type. 
    #-read The path of long read. 
    #-outputdir The path to store output alignments. 
    #-maxworker The number of threads to use. 
    #-mode H best for complex structual variation or mixed (simple and complex structual variation) analyse. L best for simple structual variation analyse.
    
    Optinal
    ----------------------
    samtools merge -@64 merged.bam /outputdir/*bam
    samtools sort -@4 merged.bam > merged.sorted.bam
    samtools index -@64 merged.sorted.bam


# VACsim
VACsim-—a tools to simulate and implant simple and complex structual variation to reference genome.

Input
-----

VACmap takes long read and reference sequence as inputs. And VACmap has been successfully tested on PacBio CLR, PacBio HiFi (CCS) and Oxford Nanopore data. Unlike other aligners use different parameter setting for different data type, VACmap uses same parameter setting for different data type. But depend on applications, we prodive two fine tune setting for downstream structural variation analyse. -mode H is fine tune for complex structural variation or mixed(both simple and complex structural variation) analyse. -mode L is fine tune for simple structural variation analyse. Comparing to -mode H, -mode L use more stringent settings to discard suspicious alignments. However, these alignments could be correct or misaligned.

Output
------

VACmap output alignment file in compressed BAM format.

Usage
----------------------
    Index refernece
    ----------------------
    conda activate VACmap
    python index.py reference_genome_path output_index_path # VACmap uses the implementation in minimap2 to build an index of the reference sequence. 
    
    Mapping long reads
    ----------------------
    conda activate VACmap
    python mammap.py -ref output_index_path -read /read.fasta -outputdir /outputdir -mode H or L -maxworker number_of_threads
    #-ref The index file produced by minimap2. The index file can be reuse for different data type. 
    #-read The path of long read. 
    #-outputdir The path to store output alignments. 
    #-maxworker The number of threads to use. 
    #-mode H best for complex structual variation or mixed (simple and complex structual variation) analyse. L best for simple structual variation analyse.
    
    Optinal
    ----------------------
    samtools merge -@64 merged.bam /outputdir/*bam
    samtools sort -@4 merged.bam > merged.sorted.bam
    samtools index -@64 merged.sorted.bam



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
