# smMIPS-Detection-Pipeline

Introduction
smMIPS-Detection-Pipeline is a bioinformatics analysis pipeline to identify somatic mutations from Illumina Next Generation Sequencing data.
The pipeline is built using Nextflow (version 20.04.1).

Pipeline Summary

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
12. Qualimap
13. Multianno files from ANNOVAR formatted using custom python scripts
14. Excel Sheet generated for each sample in Final_Output Directory 
15. Temporary files for each sample are deleted

Software Versions on which pipeline was developed-
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
Platypus

Usage
Prerequesites-
1. Install dependencies listed above.
2. Create a sample sheet in .CSV format, with a list of sample IDs.
3. Create BEDfiles as follows-
    a. Create a BED4 coordinate file
        (chr    start-coordinate    stop-coordinate    gene/exon-name)
        Example: chr9	133738125	133738246	ABL1Exon4_1
        Update path to this bedfile in nextflow.config for the "bedfile" parameter.
        
    b. Create a BED6 coordinate file for qualimap by adding two extra columns with '0' and '.' for each region.
        Example: chr9	133738125	133738246	ABL1Exon4_1   0   .
        Update path to this bedfile in nextflow.config for the "qualimap_bedfile" parameter.
        
    c. Generate platypus regions file as follows-
        
        awk 'BEGIN{OFS=""}{print $1,":",$2,"-",$3}' original_file.bed > original_file_regions.txt
      Add path to this bedfile in nextflow.config file for the "regions" parameter.
      
    d. Sort, compress and index bed file
        
        sort -k 1,1 -k 2,2n -k 3,3n original_bedfile.bed | bgzip -c > original_bedfile.bed.gz
        tabix -p bed original_bedfile.bed.gz

4. Update config file (nextflow.config) to update the paths to the input sample sheet, directory containing sequences, BEDfiles and to the tools required by the pipeline.

To run-
    
    cd /path/to/directory/containing/main.nf/
    nextflow run main.nf [ options ]

Optional Arguments (not required if updated in nextflow.config file)
    --sequences             path to directory containing fastq files.
    --input                 path to csv file containing list of sample IDs.
    --genome                path to reference genome
    --adaptors              path to FASTA file of adaptors for trimming
    --bedfile               path to BED4 coordinate BEDfile
    --qualimap_bedfile      path to BED6 qualimap BEDfile
    --regions               path to coordinate regions in .txt format
    

