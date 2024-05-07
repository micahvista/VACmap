# VACmap
VACmap-—a long-read aligner specifically designed for complex structural variation discovery.


Installation
------------

    git clone https://github.com/micahvista/VACmap.git
    cd VACmap
    conda env create --name vacmap_env --file VACmap_environment.yml
    conda activate vacmap_env
    python setup.py install

Usage
----------------------    
    
    Map long reads
    ----------------------
    conda activate vacmap_env
    vacmap -ref /ref.fasta -read /read.fasta -mode H or L -t number_of_threads | samtools sort -@4 > alignments.sorted.bam
    samtools index -@4 alignments.sorted.bam
    
    #-ref reference sequence. 
    #-read The path of long reads. 
    #-t The number of threads to use. 
    #-mode H For aligning high error rate long read (Pacbio CLR, ONT). 
    #-mode L For aligning low error rate long read (Pacbio HiFi). 
    #--eqx Output =/X CIGAR operators for sequence match/mismatch.
    #--MD Output the MD tag.
    #--cs[=short|long] Output the cs tag. (deflaut: short cs).
    






Changelog
---------


Contact
-------

If you experience any problems or have suggestions please create an issue or a pull request.

Citation
---------

https://www.biorxiv.org/content/10.1101/2023.08.03.551566v2

License
-------

The project is licensed under the GNU General Public License.
