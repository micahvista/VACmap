# VACsim
VACsim-—a tool to simulate and implant simple and complex structural variation to reference genome.
The synthetic data used in VACmap paper can be download from: https://drive.google.com/drive/folders/1pJ-6DNNRuFKX-tgcRIKkfIsSMchzJTkX?usp=sharing.

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
    conda activate vacmap_env
    cd VACmap
    python vacsim.py -parameterfilepath parameterfile -inputgenomepath reference_genome_path -altedgenomepath alted_genome_path -outputvcfpath ouput.vcf
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
    
    ~/SURVIVOR/Debug/SURVIVOR simreads alted_genome_path.snp.fasta ~/SURVIVOR/HG002_Pac_error_profile_bwa.txt 20 PAC40x.fasta
