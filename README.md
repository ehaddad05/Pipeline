# Pipeline


# Bugs that are still being worked on:
  -Featurecounts output says that there are no reads for the accession I proccessed (Unassigned_NoFeatures	43855675; Unassigned_Unmapped	1065304). This is for sure a problem with the pipeline, because I tested the pipeline commands with the same accession's sorted bam file, one that was made prior to the creation of this pipeline, and it yielded correct results. The bowtie2 alignment output is equivalent to the alignment done in the past. Howver, looking at the file sizes, the pipeline produced a sorted bam file that was about 100mb larger than the one produced in the past. I will look into the cause.

# Features that are untested:
The combining of multiple pieces of reads into one read
Analyzing with descriptive accession name input (due to the featurecount bug)
  
# Features that were tested and work:
Accession processing (however the bug with featurecounts is making the output empty)

Cache System

Analyzing with .fcount input (I used .fcount outputs from our prior work to test this)
  
# Current Path of Pipeline:
Processing: Accession -> Trim Accession -> BuildIndex using input genbank -> Align accession to genome -> make SAM to BAM -> Sort BAM -> Run FeatureCounts -> CacheReferences

Analyzing:
Take input -> Prepare KEGG table -> Make Dataframe of reads -> (If accession files are in pieces, sum the reads) -> Normalize reads -> Calculate TE of each condition -> Merge with KEGG

# Plans for the next update:
Implement function to get RBS 

Implement function to get CAI

Implement functions that create figures displaying relationships with log_te

Merge the -p and -a arguments




# Process an accession file

-p : This argument needs to be filled to tell the pipeline you are processing. In the future, I am going to merge -p with -a.

-a : This argument takes the accession [If you already have the accession file downloaded, make sure to put it in SRA_INPUT_FILES so that the pipeline does not download it again!]

-n : This argument allows you to give the accession a descriptive name

-g : This argument takes the genbank file for that accession

python3 program.py -p y -a ERR2736130 -n RNA_Cyano_Condition1 -g ~/Cyanobacteria.gbk

# Creating a new project out of certain accession files
-u: Takes the name of the project

-q: Takes a list of descriptive names for the input files (this option is only used if you are directly providing the fcount files)

-l: Takes a list of either descriptive names of accession files OR a list of .fcount references

-k: Takes the kegg of the organism

*If using -q to provide a list of descriptive names for the input .fcount files, make sure you order the lists respective to one another.
*Make sure that the descriptive name contains RNA (to denote rna-seq) or RPF (to denote ribo-seq) in the front of the name!

Using fcount references:

python3 program.py -u ProjectName -q RNA_Cyano_Condition1 RPF_Cyano_Condition1 -l ~/Pipeline/OUTPUT/ERR2736130.fcount ~/Pipeline/OUTPUT/ERR2736134.fcount -k ~/Pipeline/CyanoKEGG.txt

Using descriptive accession references:

python3 program.py -u ProjectName -l RNA_Cyano_Condition1 RPF_Cyano_Condition1 -k ~/Pipeline/CyanoKEGG.txt

*If one of the samples have multiple reads for one condition, make sure to append "Piece#" at the end of the descriptive name so that the pipeline can sum up their reads into one dataframe column

Example:
python3 program.py -u ProjectName -l RNA_Cyano_Condition1_Piece1 RNA_Cyano_Condition1_Piece2 RPF_Cyano_Condition1_Piece1 RPF_Cyano_Condition1_Piece2 -k ~/Pipeline/CyanoKEGG.txt

# Quick tools
-t: Takes the name of what you want to get. Current choices are: GBK amino acid fasta (gbk_fa), GBK nucleotide fasta (gbk_na)

-i: Takes the reference of the genbank file

Example:

python3 program.py -t gbk_fa -i ~/Cyanobacteria.gbk

python3 program.py -t gbk_na -i ~/Cyanobacteria.gbk
