### smMIPS-Detection-Pipeline

## Introduction
smMIPS-Detection-Pipeline is a bioinformatics analysis pipeline to identify somatic mutations from Illumina Next Generation Sequencing data for libraries prepared using smMIPS.
The pipeline is built using Nextflow (version 20.10.0).
The pipeline trims adapters, creates paired end assembly, maps the reads to the human genome, uses GATK Best Practices to create BAM files. CoverView is used to report read depth coverage and quality metrics of the data. A variety of variant callers are used to call mutations followed by a machine learning based algorithm in SomaticSeq to create a consensus based VCF. Platypus and Freebayes are used to call INDELs. FLT3-ITD variants are called using Get-ITD. The variants are then annotated with ANNOVAR and CAVA and formatted using custom scripts.


## Pipeline Summary

1. Adaptor Trimming (fastq-mcf)
2. Merge paired-end reads (PEAR)
3. Alignment (bwa mem)
4. SAMtools conversion
5. Generating Final BAM files based on GATK Best Practices
    i. RealignerTargetCreator
    ii. BaseRecalibrator
    iii. PrintReads
    iv. SAMtools sort and index on aligned recaliberated BAM files 
6. Variant Calling using-
    a. Mutect2
    b. Freebayes
    c. Platypus
    d. VarDict
    e. VarScan
    f. Strelka
    g. LoFreq
    h. SomaticSeq with inputs from-
            i. Mutect2
            ii. VarDict
            iii. VarScan
            iv. Lofreq
            v. Strelka

7. Freebayes and Platypus VCF files are combined (GATK CombineVariants)
8. Annotation of SomaticSeq and Combined(Freebayes+Platypus) VCF files. (Using ANNOVAR)
9. getITD
10. Coverview
11. CAVA
12. ANNOVAR Annotated files are formatted using custom python scripts
13. Excel Sheet generated for each sample in Final_Output Directory 
14. Temporary files for each sample are deleted



## Software Dependencies
fastq-mcf = 1.05
PEAR = v0.9.10
Picard = 2.17.1
BWA = 0.7.12
Samtools = 1.7
GATK = 3.8
Picard = 2.17.1
Freebayes = v1.3.2
Lofreq =  2.1.4
VarDictJava = 1.8
Varscan = v2.3.9
CAVA = 1.2.3
Platypus

## Other Dependencies
Python 3.6
Python 2.7
pandas 0.20.4


## Usage
### Prerequisites-
1. Install dependencies listed above.
2. Create a sample sheet in .CSV format, with a list of sample IDs.
    (For example: "20NGS1234" for sample 20NGS1234_S4_L001_R1_001.fastq.gz  20NGS1234_S4_L001_R2_001.fastq.gz)
    
3. Create BEDfiles as follows-
    a. Create a BED4 coordinate file
        (chr    start-coordinate    stop-coordinate    gene/exon-name
        Example: chr9	133738125	133738246	ABL1Exon4_1)
        
        Update path to this bedfile in nextflow.config for the "bedfile" parameter.
        
    b. Generate platypus regions file as follows-
        
        awk 'BEGIN{OFS=""}{print $1,":",$2,"-",$3}' original_file.bed > original_file_regions.txt
      Add path to this bedfile in nextflow.config file for the "regions" parameter.
      
    c. Sort, compress and index bed file
        
        sort -k 1,1 -k 2,2n -k 3,3n original_bedfile.bed | bgzip -c > original_bedfile.bed.gz
        tabix -p bed original_bedfile.bed.gz
    
    (Note: Make sure the 4 files below have the same name: example bedfile.bed , bedfile_regions.txt, bedfile.bed.gz , bedfile.bed.gz.tbi)
    
4. Upload all fastq.gz files to a directory. Add the path to the directory in the sequences parameter in nextflow.config or provide the path on command line with the --sequences option.
5. Update config file (nextflow.config) to add the paths to the input sample sheet, directory containing sequences, BEDfile and to the tools required by the pipeline.

### To run-
    
    cd /path/to/directory/containing/main.nf/
    nextflow run main.nf [ options ] -entry [workflow]

## Optional Arguments (not required if updated in nextflow.config file)
    --sequences             path to directory containing fastq files.
    --input                 path to csv file containing list of sample IDs.
    --genome                path to reference genome
    --adaptors              path to FASTA file of adaptors for trimming
    --bedfile               path to BED4 coordinate BEDfile (without the extension)
## Workflow
    --entry MIPS            to run the MIPS pipeline
    --entry ALL             to run ALL pipeline
    --entry CLL             to run CLL pipeline
    --entry AMPLICON        for amplicon based pipelines