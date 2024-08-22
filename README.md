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
    
    conda activate vacmap_env
    vacmap -ref /ref.fasta -read /read.fasta -mode S -t number_of_threads | samtools sort -@4 > alignments.sorted.bam
    samtools index -@4 alignments.sorted.bam
    
    -ref The path of reference sequence. 
    -read The path of long reads. 
    -t The number of threads to use. 
    

    Mapping:
        -mode H For aligning high error rate long read (Pacbio CLR, ONT). 
        -mode L For aligning low error rate long read (Pacbio HiFi). 
        -mode S Increase the sensitivity for small variants. (<100bp). (Pacbio CLR, ONT, HiFi). 
        -mode R Use a fixed value for the variation penalty. (Pacbio CLR, ONT, HiFi). 
        
        -k k-mer size (no larger than 28, deflaut: 15) # set -k 19 -w 10 for HiFi data 
            to reduce run time (2X faster) but there is very small decrease in accuracy.
        -w minimizer window size. (deflaut: 10)

    
    Output: 

        --eqx Output =/X CIGAR operators for sequence match/mismatch.
        
        --MD Output the MD tag.
        
        --cs[=short|long] Output the cs tag. (deflaut: short cs).
        
        Copy FASTA/Q comments to output. [--copycomments]
    
        --rg-id <string> Adds RG:Z:<string> to all alignments
        --rg-sm <string> RG header: Sample 
        --rg-lb <string> RG header: Library 
        --rg-pl <string> RG header: Platform
        --rg-ds <string> RG header: Description
        --rg-dt <string> RG header: Date (format: YYYY-MM-DD)
        --rg-pu <string> RG header: Platform unit 
        --rg-pi <string> RG header: Median insert size
        --rg-pg <string> RG header: Programs 
        --rg-cn <string> RG header: sequencing center
        --rg-fo <string> RG header: Flow order 
        --rg-ks <string> RG header: Key sequence 
        --rg-pm <string> Platform model. Free-form text providing further details of the platform/technology used.
        --rg-bc <string> Barcode sequence identifying the sample or library.
    






Changelog
---------

Copy FASTA/Q comments to output. Aug 22, 2024

Improved speed, 40% faster. Aug 9, 2024

Fix bug. Aug 5, 2024

Enabled overlapping anchors in non-linear chaining. Jul 4, 2024

Contact
-------

If you experience any problems or have suggestions please create an issue or a pull request.

Demo
-------
This is a demonstration of the non-linear alignment algorithm used in VACmap: http://64.64.240.35:8050/. I hope this tool proves helpful for your research!

Citation
---------

https://www.biorxiv.org/content/10.1101/2023.08.03.551566v2

License
-------

The project is licensed under the GNU General Public License.
