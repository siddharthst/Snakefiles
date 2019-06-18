import os
shell.executable("/bin/bash")

# Set working directory 
workdir: ''

# Set path to paired end reads
pathToPEReads = ''

# Set path to snakemake index
pathToHISAT2index = ''

# Set path to temp
pathToTemp = 'temp/'

# Set path for BAM output
pathToBAM = ''

# This step will generate list of files present in directory 
# containging sequencing reads. 
# Change the suffix accordingly.
readsPE, = glob_wildcards(pathToPEReads + '{sample}_1.fastq.gz')

# Number of threads for samtools
samToolsThreads = 10

# Define cleaning parameters for trim_galore
quality = 20
minLength = 80
stringency = 2 
RightClip_Read1 = 11
RightClip_Read2 = 11

rule start:
    input:
        expand(pathToBAM + '{sample}.bam', sample=readsPE)

# Since this makefile deals with alignment using raw data
# we need to clean the resulting alignment files.
# Modify it if using anything other than standard illumina
# TruSeq
rule trimgalorePE:
    threads:2
	input:
		fastq1 = pathToPEReads + '{sample}_1.fastq.gz',
		fastq2 = pathToPEReads + '{sample}_2.fastq.gz'
	output:
		output1 = pathToTemp + '{sample}_1_val_1.fq.gz',
		output2 = pathToTemp + '{sample}_2_val_2.fq.gz'
	shell:'''
        trim_galore --quality {quality} --paired --length {minLength} --illumina --stringency {stringency} \
        --clip_R1 {RightClip_Read1} --clip_R2 {RightClip_Read2} -o {pathToTemp} \ 
        --paired {input.fastq1} {input.fastq2}
        '''

# Do remember to configure the HISAT2 parameters. 
# Read the manual here : https://ccb.jhu.edu/software/hisat2/manual.shtml
rule align_hisat2:
    input:
        fastq1 = pathToTemp + '{sample}_1_val_1.fq.gz',
        fastq2 = pathToTemp + '{sample}_2_val_2.fq.gz'
    output: 
        output = 'processed_data/BAM/{sample}.bam'
    threads: 10
    shell:'''
        hisat2 -q -p {threads} --no-discordant --no-mixed --rna-strandness RF --dta -x {pathToHISAT2index} \ 
        -1 {input.fastq1} -2 {input.fastq2} | samtools sort -@ {samToolsThreads} -o {output}
        '''

rule clean:
    shell:
        '''
        rm temp/*
        '''
