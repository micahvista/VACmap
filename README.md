# VACmap
VACmap-—a long-read aligner specifically designed for complex structural variation discovery.


Requirements
------------

- Linux OS (tested on Ubuntu)
- Python

Installation
------------

    git clone https://github.com/micahvista/VACmap.git
    cd VACmap
    conda env create --name vacmap_env --file VACmap_environment.yml
    conda activate vacmap_env
    python setup.py install

Usage
----------------------  

    Map long genomic reads:
    ---------------------- 
    
    vacmap -ref /ref.fasta -read /read.fasta -mode H -t 8 | samtools sort -@4 > alignments.sorted.bam
    samtools index -@4 alignments.sorted.bam
    
    *Memory usage: <20GB
    *Speed: Processing HG002 50X ONT data (170GB) took 5.8 hours using 40 threads.
    
    Map assembly:
    ---------------------- 

    vacmap -ref /ref.fasta -read /read.fasta -mode asm -t 8  -workdir /home/usr/workdir/ --H --fakecigar | samtools sort -@4 > alignments.sorted.bam
    samtools index -@4 alignments.sorted.bam
    
    or (lower memory usage)
    
    vacmap -ref /ref.fasta -read /read.fasta -mode asm -t 8  -workdir /home/usr/workdir/ --H --fakecigar > alignments.sam
    samtools sort -@4 alignments.sam > alignments.sorted.bam
    samtools index -@4 alignments.sorted.bam
    
    *Memory usage: 30GB-70GB Using more threads requires more memory.
    *Speed: Processing HG002 v1.0 assembly(48 contigs) took 2.3 hours using 40 threads.

Parameter
----------------------  
    Mandatory parameter:
    
        -ref The path of reference sequence. 
        -read The path of long reads. 
        -t The number of threads to use. (Note: In asm mode, using fewer threads results in lower memory usage.)
        -mode S|L|H|R|asm
        -workdir Only for asm mode, store temporary data. If the folder not existexist, VACmap will create it automatically.

    Input:
        Supported read input types: FASTA, FASTA.GZ, FASTQ, FASTQ.GZ, BAM.
            [Warning: Duplicate reads (based on read names) will be processed only once.]
        Multiple read input options:
            -read read_1.fa read_2.fa ...
            -read read_1.fa -read read_2.fa ...

    Mapping:
        -mode H For aligning high error rate long read (Pacbio CLR, ONT). 
        -mode L For aligning low error rate long read (Pacbio HiFi). 
        -mode S Increase the sensitivity for small variants. (<100bp).
        -mode R Use a fixed value for the variation penalty, more sensitive to translocation events, 
            such as gene conversion. (Pacbio HiFi). 
        -mode asm For full genome alignment.   
        
        -k k-mer size (no larger than 28, deflaut: 15) # set -k 19 -w 10 for HiFi data 
            to reduce run time (2X faster) but there is very small decrease in accuracy.
        -w minimizer window size. (deflaut: 10)

    
    Output: 

        --eqx Output =/X CIGAR operators for sequence match/mismatch.
        
        --MD Output the MD tag.
        
        --cs[=short|long] Output the cs tag. (deflaut: short cs).
        
        --copycomments Copy FASTA/Q comments to output. 

        --H Use hard-clipping instead of soft-clipping for clipped sequences.
            [Warning: SV callers may fail to detect split-read events.]

        --Q Ignore base quality in the input file.

        --fakecigar Use approximate CIGAR in the SA tag.
    
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

Validation Installation
---------------------- 

To verify a successful installation, use the sample test data provided in the testdata folder. If the installation is successful, VACmap will generate three alignments.

Limitations
----------------------  

    Single Reference Sequence Size Limit: 
        2,147,483,647 bases.
    Single Query Sequence Size Limit:
        2,147,483,647 bases for Modes L, H, S, R.
        Unlimited for Mode asm.
        

Demo
-------
This is a demonstration of the non-linear alignment algorithm used in VACmap: [http://154.12.33.179:8050](http://154.12.33.179:8050/)/. I hope this tool proves helpful for your research!


Changelog
---------

Full genome alignment. Oct 3, 2024

Reduce runtime (45% faster for aligning HiFi data) and include options to decrease output size. Seq 3, 2024

Reduce runtime (45% faster for aligning HiFi data) and include options to decrease output size. Seq 3, 2024

Copy FASTA/Q comments to output. Aug 22, 2024

Improved speed, 40% faster. Aug 9, 2024

Fix bug. Aug 5, 2024

Enabled overlapping anchors in non-linear chaining. Jul 4, 2024


TODO
---------

    Reimplement VACmap in C:
        Seed: Completed
        Non-linear Chaining: Completed
        Extend: Completed
        Testing

Contact
-------

If you experience any problems or have suggestions please create an issue or a pull request.


Citation
---------

[https://www.biorxiv.org/content/10.1101/2023.08.03.551566v3](https://www.biorxiv.org/content/10.1101/2023.08.03.551566v3)

License
-------

The project is licensed under the GNU General Public License.
