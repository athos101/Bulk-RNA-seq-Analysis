import subprocess
import os
import sys

output_directory = sys.argv[3]
genome_directory = sys.argv[1]
fastq_directory = sys.argv[2]

subprocess.call(f"STAR --runThreadN 64 --runMode genomeGenerate --genomeDir {genome_directory} --genomeFastaFiles {genome_directory}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --sjdbGTFfile {genome_directory}/Homo_sapiens.GRCh38.113.gtf", shell=True)

os.mkdir(output_directory)
for fastq in os.listdir(fastq_directory):
    if fastq.endswith('.fastq.gz'):
        prefix=fastq.strip(".fastq.gz") + "_output"
        os.mkdir(output_directory + prefix)
        print ("Currently mapping: " + fastq)
        subprocess.call("STAR --runThreadN 64 --genomeDir /" + genome_directory + " --readFilesCommand zcat --outFilterType BySJout --outFilterMismatchNoverLmax 0.04 --outFilterMismatchNmax 999 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --outFilterMultimapNmax 20 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesIn /"+ fastq_directory + fastq + " --clip3pAdapterSeq GATCGGAAGAGCACACGTCTGAACTCCAGTCAC --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix " + output_directory + prefix + "/", shell=True)
