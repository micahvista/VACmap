# VACmap
VACmap-—a long-read aligner specifically designed for complex structural variation discovery.


Installation
------------

    git clone https://github.com/micahvista/VACmap.git
    cd VACmap
    conda env create --name vacmap_env --file VACmap_environment.yml
    conda activate vacmap_env
    python setup.py install

Input
-----

VACmap takes long reads and reference sequences as inputs. 

Output
------

VACmap outputs alignment files in compressed BAM format.

Usage
----------------------
    Index reference
    ----------------------

    conda create -n minimap minimap2
    conda activate minimap
    minimap2 -d output_index_path reference_genome_path  -w 10 -k 15
    
    
    Map long reads
    ----------------------
    conda activate vacmap_env
    vacmap -ref output_index_path -read /read.fasta -outputdir /outputdir -mode H or L -maxworker number_of_threads
    #-ref The index file produced by minimap2. 
    #-read The path of long reads. 
    #-outputdir The path to store output alignments. 
    #-maxworker The number of threads to use. 
    #-mode H For aligning high error rate long read (Pacbio CLR, ONT). 
    #-mode L For aligning low error rate long read (Pacbio CCS). 
    
    Merge and sort BAM file
    ----------------------
    samtools merge -@8 merged.bam /outputdir/*bam
    samtools sort -@4 merged.bam > merged.sorted.bam
    samtools index -@8 merged.sorted.bam


# VACsim
VACsim-—a tool to simulate and implant simple and complex structural variation to reference genome.

Input
-----

VACsim takes a configuration file used to generate simulated structural variation and a reference sequence as inputs. And VACsim allows user to specify or randomly create the type of complex structural variation.

Output
------

VACsim output altered reference sequence.

Usage
----------------------
    Configure the parameterfile
    ----------------------
    #NML: Normal sequence, DEL: Deletion event, INS: Insertion event, INV: Inversion event, DUP: Duplication event, TRA: Translocation event 
    #Define the composition of complex structural variations. For "DUP:100:200:0:3", DUP refers to duplication events, 100: 200 refers to the lower bound (include) and upper bound (exclude) to the size of the duplication sequence, 0 refers to this duplication sequence is not inverted (1 for inverted), 3 refer to the duplication sequence is repeated three times.    
    #number: The number of simulated structural variations.
    Specified{DEL:100:200,INS:100:1000,INV:100:200,DUP:100:200:0:3,TRA:200:400:1;number=2}
    Specified{INV:100:200,INV:100:200,TRA:200:400:1;number=2}
    Specified{INV:100:200,NML:100:200,TRA:200:400:1;number=2}
    Specified{INV:100:200;number=2}
    #eventcount=[1(include),5(exclude)] VACsim will randomly choose a number as the event count of simulated structural variation.
    #eventset VACsim will randomly choose event count times to determine the composition of structural variations. 
    Random{eventset=["DEL:100:200","INS:100:1000","INV:100:200","DUP:100:200","TRA:200:400"];eventcount=[1,5];number=2}
    Random{eventset=["DEL:100:200,INV:100:200","INS:100:1000,NML:100:200","NML:100:200,INV:100:200","DUP:100:200","TRA:200:400"];eventcount=[4,20];number=2}
        
    Simulating SV
    ----------------------
    conda activate VACmap
    cd VACmap
    python mamsim.py -parameterfilepath parameterfile -inputgenomepath reference_genome_path -altedgenomepath alted_genome_path -outputvcfpath ouput.vcf
    #-parameterfilepath The configured parameter file path. 
    #-inputgenomepath The path of input genome sequence path. 
    #-altedgenomepath The output path of altered genome sequence path. 
    #-outputvcfpath The output path of VCF file. 

    
    Add SNPs and simulate long read
    ----------------------
    git clone https://github.com/fritzsedlazeck/SURVIVOR.git
    cd SURVIVOR/Debug
    make
    cd ..
    unzip HG002_Pac_error_profile_bwa.txt.zip
    
    ~/SURVIVOR/Debug/SURVIVOR simSV alted_genome_path parameter_file 0.001 0 alted_genome_path.snp
    
    ~/SURVIVOR/Debug/SURVIVOR simreads alted_genome_path.snp.fasta ~/SURVIVOR/HG002_Pac_error_profile_bwa.txt 40 PAC40x.fasta




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
