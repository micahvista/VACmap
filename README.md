# VACmap
VACmap-—a long-read aligner specifically designed for complex structural variation discovery.

Important
---------
If you are interested in detecting complex variants, we recommend using the '-mode S' option to increase sensitivity to small rearrangements. However, existing structural variant callers may fail to extract variant signals from the alignments produced by this mode.

# Demo
This is a demonstration of the non-linear alignment algorithm used in VACmap: http://64.64.240.35:8050/. I hope this tool proves helpful for your research!


Installation
------------

    git clone https://github.com/micahvista/VACmap.git
    cd VACmap
    conda env create --name vacmap_env --file VACmap_environment.yml
    conda activate vacmap_env
    python setup.py install

Usage
----------------------    
    
    conda activate vacmap_env
    vacmap -ref /ref.fasta -read /read.fasta -mode S -t number_of_threads | samtools sort -@4 > alignments.sorted.bam
    samtools index -@4 alignments.sorted.bam
    
    #-ref The path of reference sequence. 
    #-read The path of long reads. 
    #-t The number of threads to use. 
    #-mode S For discovering complex variants (Pacbio CLR, ONT, HiFi). 
    #-mode H For aligning high error rate long read (Pacbio CLR, ONT). 
    #-mode L For aligning low error rate long read (Pacbio HiFi). 
    #--eqx Output =/X CIGAR operators for sequence match/mismatch.
    #--MD Output the MD tag.
    #--cs[=short|long] Output the cs tag. (deflaut: short cs).
    #-k k-mer size (no larger than 28, deflaut: 15) # set -k 19 -w 10 for HiFi data 
        to reduce run time (2X faster) but there is very small decrease in accuracy.
    #-w minimizer window size. (deflaut: 10)
    






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
