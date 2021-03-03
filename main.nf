#!/usr/bin/env nextflow


Channel
    .fromPath(params.input)
    .splitCsv(header:false)
    .map{ it }
    .set { samples_ch }

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=

Sample list: ${params.input}
BED file: ${params.bedfile}
Sequences in:${params.sequences}

"""

process trimming {
	maxForks 15
	publishDir "${PWD}/${Sample[0]}/processed_reads/", mode: 'copy'
	input:
		val Sample from samples_ch	
	output:
		tuple(Sample, file("*.fastq")) into trimmed_ch
	script:
	"""
	${params.ea_utils_path}/fastq-mcf -o ${Sample[0]}.R1.trimmed.fastq -o ${Sample[0]}.R2.trimmed.fastq -l 53 -k 0 -q 0 ${params.adaptors} ${params.sequences}/${Sample[0]}*_R1_*.fastq.gz ${params.sequences}/${Sample[0]}*_R2_*.fastq.gz
	"""
}

process gzip{
	executor="local"
	publishDir "${PWD}/${Sample[0]}/processed_reads/", mode: 'copy'
	input:
		set Sample, file(trimmedFiles) from trimmed_ch
	output:
		val Sample into gzipped_ch
	script:
	"""
	mkdir "$PWD/${Sample[0]}/Annovar_Modified/" 
	gzip -f ${PWD}/${Sample[0]}/processed_reads/${trimmedFiles[0]}
	gzip -f ${PWD}/${Sample[0]}/processed_reads/${trimmedFiles[1]}
	"""
}

process pair_assembly{
	publishDir "${PWD}/${Sample[0]}/assembled_reads/", mode: 'copy'
	input:
		val Sample from gzipped_ch
	output:
		tuple(Sample, file('*')) into paired_ch
	script:
	"""
	${params.pear_path} -f ${PWD}/${Sample[0]}/processed_reads/${Sample[0]}.R1.trimmed.fastq.gz -r ${PWD}/${Sample[0]}/processed_reads/${Sample[0]}.R2.trimmed.fastq.gz -o ${Sample[0]} -n 53 -j 15
	"""
}

process mapping_reads{
	publishDir "${PWD}/${Sample[0]}/mapped_reads/", mode: 'copy'
	input:
		set Sample, file(pairAssembled) from paired_ch
	output:
		tuple (Sample, file ("*.sam")) into bwa_ch
	script:
	"""
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample[0]}\\tPI:200" -M -t 20 ${params.genome} ${pairAssembled[0]} > ${Sample[0]}.sam
	"""
}

process sam_conversion{
	publishDir "$PWD/${Sample[0]}/mapped_reads/", mode: 'copy', pattern: '*.fxd_sorted.bam'
	publishDir "$PWD/${Sample[0]}/mapped_reads/", mode: 'copy', pattern: '*.fxd_sorted.bam.bai'

	input:
		set (Sample, file(samFile)) from bwa_ch
	output:
		tuple (Sample, file ("*.fxd_sorted.bam"), file ("*.fxd_sorted.bam.bai")) into realign_ch
		tuple (Sample, file ("*.fxd_sorted.bam"), file ("*.fxd_sorted.bam.bai")) into indelRealign_ch
	
	script:
	"""
	java -Xmx8G -jar ${params.picard_path} FixMateInformation I= ${samFile} O= ${Sample[0]}.fxd.sam VALIDATION_STRINGENCY=SILENT
	${params.samtools} view -bT ${params.genome} ${Sample[0]}.fxd.sam > ${Sample[0]}.fxd.bam
	${params.samtools} sort ${Sample[0]}.fxd.bam > ${Sample[0]}.fxd_sorted.bam
	${params.samtools} index ${Sample[0]}.fxd_sorted.bam > ${Sample[0]}.fxd_sorted.bam.bai
	"""
}

process RealignerTargetCreator {
	publishDir "${PWD}/${Sample[0]}/gatk38_processing/", mode: 'copy', pattern: '*.intervals'
	
	input:
		set (Sample, file (bamFiles)) from realign_ch
	output:
		tuple (Sample, file ("*.intervals")) into realignTargetCreater_ch
	script:
	"""
	java -Xmx8G -jar ${params.GATK38_path} -T RealignerTargetCreator -R ${params.genome} -nt 10 -I ${bamFiles[0]} --known ${params.site1} -o ${Sample[0]}.intervals
	"""
}

process IndelRealigner{
	publishDir "${PWD}/${Sample[0]}/gatk38_processing/", mode: 'copy', pattern: '*.realigned.bam'
	input:
		set (Sample, file (targetIntervals), file(bamFiles)) from realignTargetCreater_ch.join(indelRealign_ch)
	output:
		tuple (Sample, file ("*.realigned.bam")) into IndelRealigner_ch
		tuple (Sample, file ("*.realigned.bam")) into printreads_ch
	script:
	"""
	java -Xmx8G -jar ${params.GATK38_path} -T IndelRealigner -R ${params.genome} -I ${bamFiles[0]} -known ${params.site1} --targetIntervals ${targetIntervals} -o ${Sample[0]}.realigned.bam
	"""
}

process BaseRecalibrator{
	publishDir "${PWD}/${Sample[0]}/gatk38_processing/", mode: 'copy', pattern: '*.recal_data.table'
	input:
		set (Sample, file (bamFile)) from IndelRealigner_ch
	output:
		tuple (Sample, file ("*.recal_data.table")) into BaseRecalibrator_ch
	script:
	"""
	java -Xmx8G -jar ${params.GATK38_path} -T BaseRecalibrator -R ${params.genome} -I ${bamFile} -knownSites ${params.site2} -knownSites ${params.site3} -o ${Sample[0]}.recal_data.table
	"""
}

process PrintReads{
	publishDir "${PWD}/${Sample[0]}/gatk38_processing/", mode: 'copy', pattern: '*.aligned.recalibrated.bam'
	input:
		set (Sample, file (realignedBam), file (recal_dataTable)) from printreads_ch.join(BaseRecalibrator_ch)
	output:
		tuple (Sample, file ("*.aligned.recalibrated.bam")) into printedReads_ch
	script:
	"""
	java -Xmx8G -jar ${params.GATK38_path} -T PrintReads -R ${params.genome} -I ${realignedBam} --BQSR ${recal_dataTable} -o ${Sample[0]}.aligned.recalibrated.bam
	"""
}

process generatefinalbam{
	publishDir "$PWD/${Sample[0]}/gatk38_processing/", mode: 'copy', pattern: '*.final.bam'
	publishDir "$PWD/Final_Output/${Sample[0]}/", mode: 'copy', pattern: '*.final.bam'
	publishDir "$PWD/${Sample[0]}/gatk38_processing/", mode: 'copy', pattern: '*.final.bam.bai'
	publishDir "$PWD/Final_Output/${Sample[0]}/", mode: 'copy', pattern: '*.final.bam.bai'
	
	input:
		set (Sample , file(alignedRecalibratedBam)) from printedReads_ch
	output:
		tuple (Sample, file ("*.final.bam"),  file ("*.final.bam.bai")) into bamGenerated_ch
		tuple (Sample, file ("*")) into runMutect_ch
		tuple (Sample, file ("*")) into runFreebayes_ch
		tuple (Sample, file ("*.*")) into runPlatypus_ch
		tuple (Sample, file ("*.*")) into runVardict_ch
		tuple (Sample, file ("*.*")) into runVarscan_ch
		tuple (Sample, file ("*.*")) into runStrelka_ch
		tuple (Sample, file ("*.*")) into runLofreq_ch
		tuple (Sample, file ("*.*")) into rungetITD_ch
		tuple (Sample, file ("*.*")) into runCoverage_ch
		tuple (Sample, file ("*.*")) into runCoverview_ch
		tuple (Sample) into runSomaticSeq_ch
		
	script:
	"""
	${params.samtools} sort ${alignedRecalibratedBam} > ${Sample[0]}.final.bam
	${params.samtools} index ${Sample[0]}.final.bam > ${Sample[0]}.final.bam.bai
	"""
}

process mutect2_run{
	publishDir "$PWD/${Sample[0]}/variants/", mode: 'copy', pattern: '*.mutect2.vcf'
	
	input:
		set (Sample , file(finalBams)) from runMutect_ch
	output:
		tuple (Sample, file ("*.mutect2.vcf")) into mutect2Done_ch
		val Sample into mutect2somaticseq_ch

	script:
	"""
	java -Xmx10G -jar ${params.GATK38_path} -T MuTect2 -R ${params.genome} -I:tumor ${finalBams[0]} -o ${Sample[0]}.mutect2.vcf --dbsnp ${params.site2} -L ${params.bedfile} -nct 30 -contamination 0.02 -mbq 30
	"""
}

process freebayes_run{
	publishDir "$PWD/${Sample[0]}/variants/", mode: 'copy', pattern: '*.freebayes.vcf'
	
	input:
		set (Sample , file(finalBams)) from runFreebayes_ch
	output:
		file ("*.freebayes.vcf")
		tuple (Sample, file ("*.freebayes.vcf")) into freebayesDone_ch

	script:
	"""
	${params.freebayes_path} -f ${params.genome} -b ${finalBams[0]} -t ${params.bedfile} > ${Sample[0]}.freebayes.vcf 	
	"""
}

process platypus_run{
	executor="local"
	publishDir "$PWD/${Sample[0]}/variants/", mode: 'copy', pattern: '*.platypus.vcf'
	input:
		set (Sample , file(finalBams)) from runPlatypus_ch
	output:
		tuple (Sample, file ("*.platypus.vcf")) into platypusDone_ch
	script:
	"""
	python ${params.platypus_path} callVariants --bamFiles=${finalBams[0]} --refFile=${params.genome} --output=${Sample[0]}.platypus.vcf --nCPU=15 --minFlank=0 --filterDuplicates=0 --maxVariants=6 --minReads=6 --regions=${params.regions}
	"""
}

process vardict_run{
	publishDir "$PWD/${Sample[0]}/variants/", mode: 'copy', pattern: '*.vardict.vcf'
	
	input:
		set (Sample , file(finalBams)) from runVardict_ch
	output:
		tuple (Sample, file ("*.vardict.vcf")) into vardictDone_ch
		val Sample into vardict2somaticseq_ch
	script:
	"""
	VarDict -G ${params.genome} -f 0.03 -N ${Sample[0]} -b ${finalBams[0]} -c 1 -S 2 -E 3 -g 4 ${params.bedfile} | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample[0]} -E -f 0.03 > ${Sample[0]}.vardict.vcf
	"""
}

process varscan_run{
	publishDir "$PWD/${Sample[0]}/variants/", mode: 'copy', pattern: '*.varscan_snp.vcf'
	publishDir "$PWD/${Sample[0]}/variants/", mode: 'copy', pattern: '*.varscan_indel.vcf'
	publishDir "$PWD/${Sample[0]}/variants/", mode: 'copy', pattern: '*.varscan_snp.vcf.gz'
	publishDir "$PWD/${Sample[0]}/variants/", mode: 'copy', pattern: '*.varscan_indel.vcf.gz'
	publishDir "$PWD/${Sample[0]}/variants/", mode: 'copy', pattern: '*.varscan.vcf'
	
	input:
		set (Sample , file(finalBams)) from runVarscan_ch
	output:
		tuple (Sample, file ("*.varscan_snp.vcf.gz"), file ("*.varscan_indel.vcf.gz")) into varscanDone_ch
		tuple (file ("*.varscan_snp.vcf"),  file ("*.varscan_indel.vcf"), file("*.varscan.vcf"))
		val Sample into varscan2somaticseq_ch
		
	script:
	"""
	${params.samtools} mpileup -f ${params.genome} ${finalBams[0]} > ${Sample[0]}.mpileup
	java -jar ${params.varscan_path} mpileup2snp ${Sample[0]}.mpileup --min-coverage 10 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.003 --p-value 1e-4 --output-vcf 1 > ${Sample[0]}.varscan_snp.vcf
	java -jar ${params.varscan_path} mpileup2indel ${Sample[0]}.mpileup --min-coverage 10 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.003 --p-value 1e-4 --output-vcf 1 > ${Sample[0]}.varscan_indel.vcf
	bgzip -c ${Sample[0]}.varscan_snp.vcf > ${Sample[0]}.varscan_snp.vcf.gz
	bgzip -c ${Sample[0]}.varscan_indel.vcf > ${Sample[0]}.varscan_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample[0]}.varscan_snp.vcf.gz
	${params.bcftools_path} index -t ${Sample[0]}.varscan_indel.vcf.gz
	${params.bcftools_path} concat -a ${Sample[0]}.varscan_snp.vcf.gz ${Sample[0]}.varscan_indel.vcf.gz -o ${Sample[0]}.varscan.vcf
	"""
}

process strelka_run{
	publishDir "$PWD/${Sample[0]}/variants/strelka", mode: 'copy', pattern: '*.strelka.vcf'
	
	input:
		set (Sample , file(finalBams)) from runStrelka_ch
	output:
		val (Sample) into strelkaDone_ch
		val Sample into strelka2somaticseq_ch
	script:
	"""
	${params.strelka_path}/configureStrelkaGermlineWorkflow.py --bam ${finalBams[0]} --referenceFasta ${params.genome} --callRegions  ${params.bedfile}.gz --targeted --runDir ${PWD}/${Sample[0]}/variants/strelka/
	${PWD}/${Sample[0]}/variants/strelka/runWorkflow.py -m local -j 20
	gunzip -f ${PWD}/${Sample[0]}/variants/strelka/results/variants/variants.vcf.gz
	mv ${PWD}/${Sample[0]}/variants/strelka/results/variants/variants.vcf $PWD/${Sample[0]}/variants/${Sample[0]}.strelka.vcf
	
	${params.strelka_path}/configureStrelkaSomaticWorkflow.py --normalBam ${params.NA12878_bam}  --tumorBam ${finalBams[0]} --referenceFasta ${params.genome} --callRegions ${params.bedfile}.gz --targeted --runDir ${PWD}/${Sample[0]}/variants/strelka-somatic/
	${PWD}/${Sample[0]}/variants/strelka-somatic/runWorkflow.py -m local -j 20
	
	${params.bcftools_path} concat -a ${PWD}/${Sample[0]}/variants/strelka-somatic/results/variants/somatic.indels.vcf.gz ${PWD}/${Sample[0]}/variants/strelka-somatic/results/variants/somatic.snvs.vcf.gz -o ${PWD}/${Sample[0]}/variants/${Sample[0]}.strelka-somatic.vcf
	"""
}

process lofreq_run{
	publishDir "$PWD/${Sample[0]}/variants/", mode: 'copy', pattern: '*.lofreq.filtered.vcf'
	
	input:
		set (Sample , file(finalBams)) from runLofreq_ch
	output:
		tuple (Sample, file ("*.lofreq.filtered.vcf")) into lofreqDone_ch
		val Sample into lofreq2somaticseq_ch
	script:
	"""
	${params.lofreq_path} viterbi -f ${params.genome} -o ${Sample[0]}.lofreq.pre.bam ${finalBams[0]}
	${params.samtools} sort ${Sample[0]}.lofreq.pre.bam > ${Sample[0]}.lofreq.bam
	${params.lofreq_path} call -b dynamic -C 50 -a 0.00005 -q 30 -Q 30 -m 50 -f ${params.genome} -l ${params.bedfile} -o ${Sample[0]}.lofreq.vcf ${Sample[0]}.lofreq.bam
	${params.lofreq_path} filter -a 0.01 -i ${Sample[0]}.lofreq.vcf -o ${Sample[0]}.lofreq.filtered.vcf
	"""
}

process somaticSeq_run {
	publishDir "$PWD/${Sample_id[0]}/variants/", mode: 'copy', pattern: '*.somaticseq.vcf'
	publishDir "$PWD/${Sample_id[0]}/variants/", mode: 'copy', pattern: '*.avinput'
	publishDir "$PWD/${Sample_id[0]}/ANNOVAR/", mode: 'copy', pattern: '*.hg19_multianno.csv'
	input:
		each (Sample_id) from runSomaticSeq_ch.toList()
		val (Sample) from mutect2somaticseq_ch.toList()
		val (Sample) from vardict2somaticseq_ch.toList()
		val (Sample) from varscan2somaticseq_ch.toList()
		val (Sample) from lofreq2somaticseq_ch.toList()
		val (Sample) from strelka2somaticseq_ch.toList()
	output:
		tuple (Sample_id, file ("*.somaticseq.vcf")) into somaticseqDone_ch
		tuple (Sample_id, file ("*.avinput"), file ("*.hg19_multianno.csv"))
		tuple (Sample_id, file ("*.hg19_multianno.csv")) into formatSomaticseq_ch
	script:
	"""
	python3 /scratch/hematopath/programs/somaticseq/somaticseq/somaticseq_parallel.py --output-directory ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.somaticseq --genome-reference ${params.genome} --inclusion-region ${params.bedfile} --threads 25 --algorithm xgboost  --dbsnp-vcf  ~/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file ${PWD}/${Sample_id[0]}/gatk38_processing/${Sample_id[0]}.final.bam --mutect2-vcf ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.mutect2.vcf --vardict-vcf ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.vardict.vcf --varscan-vcf ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.varscan.vcf --lofreq-vcf ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.lofreq.vcf	--strelka-vcf ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.strelka.vcf  --sample-name ${Sample_id[0]}
	
	grep "^#" ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.somaticseq/Consensus.sSNV.vcf > ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.somaticseq/somaticseq_snv.vcf
	grep -v "^#" ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.somaticseq/Consensus.sSNV.vcf | sort -k1,1V -k2,2g >> ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.somaticseq/somaticseq_snv.vcf
	bgzip -c ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.somaticseq/somaticseq_snv.vcf > ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.somaticseq/somaticseq_snv.vcf.gz
	${params.bcftools_path} index -t ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.somaticseq/somaticseq_snv.vcf.gz
	
	grep "^#" ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.somaticseq/Consensus.sINDEL.vcf > ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.somaticseq/somaticseq_indel.vcf
	grep -v "^#" ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.somaticseq/Consensus.sINDEL.vcf | sort -k1,1V -k2,2g >> ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.somaticseq/somaticseq_indel.vcf
	bgzip -c ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.somaticseq/somaticseq_indel.vcf > ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.somaticseq/somaticseq_indel.vcf.gz
	${params.bcftools_path} index -t ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.somaticseq/somaticseq_indel.vcf.gz
	
	${params.bcftools_path} concat -a ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.somaticseq/somaticseq_snv.vcf.gz ${PWD}/${Sample_id[0]}/variants/${Sample_id[0]}.somaticseq/somaticseq_indel.vcf.gz -o ${Sample_id[0]}.somaticseq.vcf
	
	perl ${params.annovar2_path}/convert2annovar.pl -format vcf4 ${Sample_id[0]}.somaticseq.vcf  --outfile ${Sample_id[0]}.somaticseq.avinput --withzyg --includeinfo
	
	perl ${params.annovar2_path}/table_annovar.pl ${Sample_id[0]}.somaticseq.avinput --out ${Sample_id[0]}.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all --operation g,r,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovar2_path}/humandb/ --xreffile ${params.annovar2_path}/example/gene_fullxref.txt
	"""
}

process coverage {
	publishDir "$PWD/${Sample[0]}/coverage/", mode: 'copy'
	input:
		set (Sample , file(finalBams)) from runCoverage_ch
	output:
		file ("*")
		
	script:
	"""
	${params.bedtools} bamtobed -i ${finalBams[0]} > ${Sample[0]}.bed
	${params.bedtools} coverage -counts -a ${params.bedfile} -b ${Sample[0]}.bed > ${Sample[0]}.counts.bed
	
	"""
}

process getITD_run{
	publishDir "$PWD/${Sample[0]}/getITD/", mode: 'copy'
	input:
		set (Sample , file(finalBams)) from rungetITD_ch
	output:
		file ("*")
		val Sample into getITD_Done_ch
	script:
	"""
	gunzip -c ${params.sequences}/${Sample[0]}*_R1_*.fastq.gz > ${Sample[0]}_R1_.fastq
	gunzip -c ${params.sequences}/${Sample[0]}*_R2_*.fastq.gz > ${Sample[0]}_R2_.fastq
	
	python3 ${params.getitd_path}/getitd.py ${Sample[0]} ${Sample[0]}_R1_.fastq ${Sample[0]}_R2_.fastq -require_indel_free_primers False -anno ${params.getitd_path}/anno/amplicon_kayser.tsv -reference ${params.getitd_path}/anno/amplicon.txt -forward_primer ${params.forward_primer} -reverse_primer ${params.reverse_primer} -forward_adapter  ATACGAGATCCGTAATCGGGAAGCTGAAG -reverse_adapter ACACGCACGATCCGACGGTAGTGT -nkern 25
	"""
}

process combine_variants{
	publishDir "$PWD/${Sample[0]}/variants/", mode: 'copy'
	publishDir "$PWD/${Sample[0]}/variants/", mode: 'copy', pattern: '*.avinput'
	publishDir "$PWD/${Sample[0]}/ANNOVAR/", mode: 'copy', pattern: '*.hg19_multianno.csv'
	
	input:
		set (Sample, file(freebayesVcf), file(platypusVcf)) from freebayesDone_ch.join(platypusDone_ch)
	output:
		file ("*")
		tuple (Sample, file ("*.combined.vcf")) into combineVariantsDone_ch
		tuple (Sample, file ("*.hg19_multianno.csv")) into formatCombined_ch
	script:
	"""
	grep "^#" ${PWD}/${Sample[0]}/variants/${Sample[0]}.freebayes.vcf > ${Sample[0]}.freebayes.sorted.vcf
	grep -v "^#" ${PWD}/${Sample[0]}/variants/${Sample[0]}.freebayes.vcf | sort -k1,1V -k2,2g >> ${Sample[0]}.freebayes.sorted.vcf
	
	grep "^#" ${PWD}/${Sample[0]}/variants/${Sample[0]}.platypus.vcf > ${Sample[0]}.platypus.sorted.vcf
	grep -v "^#" ${PWD}/${Sample[0]}/variants/${Sample[0]}.platypus.vcf | sort -k1,1V -k2,2g >> ${Sample[0]}.platypus.sorted.vcf
	
	java -jar ${params.GATK38_path} -T CombineVariants -R ${params.genome} --variant ${Sample[0]}.freebayes.sorted.vcf --variant ${Sample[0]}.platypus.sorted.vcf -o ${Sample[0]}.combined.vcf -genotypeMergeOptions UNIQUIFY
	
	perl ${params.annovar2_path}/convert2annovar.pl -format vcf4 ${Sample[0]}.combined.vcf  --outfile ${Sample[0]}.combined.avinput --withzyg --includeinfo
	
	perl ${params.annovar2_path}/table_annovar.pl ${Sample[0]}.combined.avinput --out ${Sample[0]}.combined --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all --operation g,r,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovar2_path}/humandb/ --xreffile ${params.annovar2_path}/example/gene_fullxref.txt
	"""
}

process coverview_run {
	executor="local"
	publishDir "$PWD/${Sample[0]}/Coverview/", mode: 'copy'
	input:
		set (Sample , file(finalBam)) from runCoverview_ch
	output:
		val (Sample) into coverviewDone_ch
	script:
	"""
	mkdir $PWD/${Sample[0]}/Coverview/
	${params.coverview_path}/coverview -i ${finalBam[0]} -b ${params.bedfile} -c ${params.coverview_path}/config/config.txt -o ${PWD}/${Sample[0]}/Coverview/${Sample[0]}.coverview
	python3 ${params.coverview_script_path} ${PWD}/${Sample[0]}/Coverview/${Sample[0]}.coverview_regions.txt ${PWD}/${Sample[0]}/Coverview/${Sample[0]}.coverview_regions.csv
	cp ${PWD}/${Sample[0]}/Coverview/${Sample[0]}.coverview_regions.csv ${PWD}/Coverview/${Sample[0]}.coverview_regions.csv
	"""
}

process coverview_report {
	executor="local"
	input:
		set (Sample) from coverviewDone_ch.toList()
	output:
		val Sample into coverviewreportDone_ch
	script:
	"""
	python3.6 ${params.coverview_report_path} ${PWD}/Coverview/ ${PWD}/Final_Output/
	"""
}

process cava {
	publishDir "$PWD/${Sample[0]}/CAVA/", mode: 'copy'
	
	input:
		tuple (Sample, file (somaticVcf), file(combinedVcf)) from somaticseqDone_ch.join(combineVariantsDone_ch)
	
	output:
		file ("*")
		val Sample into cavaDone_ch
	script:
	"""
	echo true
	${params.cava_path}/cava -c ${params.cava_path}/config.txt -t 10 -i $PWD/${Sample[0]}/variants/${Sample[0]}.somaticseq.vcf -o ${Sample[0]}.somaticseq
	${params.cava_path}/cava -c ${params.cava_path}/config.txt -t 10 -i $PWD/${Sample[0]}/variants/${Sample[0]}.combined.vcf -o ${Sample[0]}.combined
	python3.6 ${params.cava_script_path} ${Sample[0]}.somaticseq.txt ${Sample[0]}.combined.txt ${Sample[0]}.cava.csv
	"""
}

process Qualimap{
	executor="local"
	errorStrategy 'retry'
	input:
		set (Sample, file (finalBam)) from bamGenerated_ch
	output:
		val Sample into qualimapDone_ch
	shell:
	"""
	echo ${PWD}
	${params.qualimap_path}/qualimap bamqc -bam ${finalBam[0]} -gff ${params.qualimap_bedfile} -c -outdir ${PWD}/${Sample[0]}/Qualimap/ -outfile ${Sample[0]}.qc_report.pdf
	"""
}

process format_somaticseq {
	executor="local"
	input:
		set (Sample , file(hgMultiannoFile)) from formatSomaticseq_ch
	output:
		val Sample into somaticseqMerge_ch
	script:
	"""
	python3.6 ${params.format_somaticseq_script} ${PWD}/${Sample[0]}/ANNOVAR/${hgMultiannoFile} ${PWD}/${Sample[0]}/Annovar_Modified/${Sample[0]}.somaticseq.csv
	"""
}

process format_combined {
	executor="local"
	input:
		set (Sample , file(hgMultiannoFile)) from formatCombined_ch
	output:
		val Sample into combinedMerge_ch
	script:
	"""
	python3.6 ${params.format_combined_script} ${PWD}/${Sample[0]}/ANNOVAR/${hgMultiannoFile} ${PWD}/${Sample[0]}/Annovar_Modified/${Sample[0]}.combined.csv
	"""
}

process merge_csv {
	executor="local"
	input:
		each (Sample) from somaticseqMerge_ch.toList()
		set (Sample2) from combinedMerge_ch.toList()
		set (Sample3) from cavaDone_ch.toList()
	output:
		val Sample into mergeDone_ch
	script:
	"""
	echo true
	python3.6 ${params.merge_csvs_script} ${Sample[0]} ${PWD}/${Sample[0]}/Annovar_Modified/ ${PWD}/Final_Output/${Sample[0]}/${Sample[0]}.xlsx ${PWD}/${Sample[0]}/CAVA/ ${PWD}/${Sample[0]}/Coverview/${Sample[0]}.coverview_regions.csv
	"""
}

process Final_Output {
	executor="local"
	input:
		val Sample from qualimapDone_ch.join(getITD_Done_ch)
	output:
		val Sample into finalOutput_ch
	script:
	"""
	python3.6 /scratch/hematopath/hemoseq_v2/mutation_detector_nextflow/scripts/coverageplot.py ${Sample} $PWD/${Sample}/coverage/${Sample}.counts.bed $PWD/${Sample}/coverage/
	cp ${PWD}/${Sample}/Qualimap/${Sample}.qc_report.pdf ${PWD}/Final_Output/${Sample}/
	cp -r ${PWD}/${Sample}/getITD/${Sample}_getitd/ ${PWD}/Final_Output/${Sample}/
	cp ${PWD}/${Sample}/coverage/${Sample}.Low_Coverage.png ${PWD}/Final_Output/${Sample}/
	"""
}

process remove_files{
	input:
		each Sample from finalOutput_ch.toList()
		val Sample1 from coverviewreportDone_ch.toList()
		val Sample2 from mergeDone_ch.toList()
	script:
	"""
	rm -rf ${PWD}/${Sample[0]}/
	"""
}



