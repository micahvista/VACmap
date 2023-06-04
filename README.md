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
    
    OR
    
    conda create -n minimap minimap2
    conda activate minimap
    minimap2 -d output_index_path reference_genome_path  -w 10 -k 15
    
    
    Map long reads
    ----------------------
    conda activate VACmap
    python mammap.py -ref output_index_path -read /read.fasta -outputdir /outputdir -mode H or L -maxworker number_of_threads
    #-ref The index file produced by minimap2. The index file can be reuse for different data type. 
    #-read The path of long read. 
    #-outputdir The path to store output alignments. 
    #-maxworker The number of threads to use. 
    #-mode H best for complex structual variation or mixed (simple and complex structual variation) analyse. L best for simple structual variation analyse.
    
    Merge and sort BAM file
    ----------------------
    samtools merge -@64 merged.bam /outputdir/*bam
    samtools sort -@4 merged.bam > merged.sorted.bam
    samtools index -@64 merged.sorted.bam


# VACsim
VACsim-—a tools to simulate and implant simple and complex structual variation to reference genome.

Input
-----

VACsim takes a configuration file used to generate simulated structual variation and a reference sequence as inputs. And VACsim allow user to specify or randomly create the type of complex structural variation.

Output
------

VACsim output altered reference sequence.

Usage
----------------------
    Configure the parameterfile
    ----------------------
    #NML: Normal sequence, DEL: Deletion event, INS: Insertion event, INV: Inversion event, DUP: Duplication event, TRA: Traslocation event 
    #Define the composition of complex structural variations. For "DUP:100:200:0:3", DUP refer to duplication events, 100: 200 refer to lower bound (include) and upper bound (exclude) to the size of duplication sequence, 0 refer to this duplication sequence is not inverted (1 for inverted), 3 refer to the duplication sequence is repeated three times.    
    #number: The number of simulated structural variations.
    Specified{DEL:100:200,INS:100:1000,INV:100:200,DUP:100:200:0:3,TRA:200:400:1;number=2}
    Specified{INV:100:200,INV:100:200,TRA:200:400:1;number=2}
    Specified{INV:100:200,NML:100:200,TRA:200:400:1;number=2}
    Specified{INV:100:200;number=2}
    #eventcount=[1(include),5(exclude)] VACsim will randomly choose a number as the event count of simluated structural variation.
    #eventset VACsim will randomly choose eventcount times to determine the composition of structural variations. Important: If you provide a nest event (e.g. "DEL:100:200,INV:100:200") the final eventcount may exceed defined eventcount. User can identify the true eventcount in output VCF using bp annotation.
    Random{eventset=["DEL:100:200","INS:100:1000","INV:100:200","DUP:100:200","TRA:200:400"];eventcount=[1,5];number=2}
    Random{eventset=["DEL:100:200,INV:100:200","INS:100:1000,NML:100:200","NML:100:200,INV:100:200","DUP:100:200","TRA:200:400"];eventcount=[4,20];number=2}
        
    Simulating SV
    ----------------------
    conda activate VACmap
    cd VACmap
    python mamsim.py -parameterfilepath parameterfile -inputgenomepath reference_genome_path -altedgenomepath alted_genome_path -outputvcfpath ouput.vcf
    #-parameterfilepath The configured the parameter file path. 
    #-inputgenomepath The path of input genome sequence path. 
    #-altedgenomepath The output path of alted genome sequence path. 
    #-outputvcfpath The output path of VCF file. 

    
    Add SNPs and simulating long read
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
