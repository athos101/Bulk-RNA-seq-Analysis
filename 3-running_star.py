import subprocess
import os

output_directory = "STAR_output/"
genome_directory = "/data/genomes/h38/STAR/"
fastq_directory = "/data/analysis/hypoxia/fastq/"

os.mkdir(output_directory)
for fastq in os.listdir(fastq_directory):
    if fastq.endswith('.fastq.gz'):
        prefix=fastq.strip(".fastq.gz") + "_output"
        os.mkdir(output_directory + prefix)
        print ("Currently mapping: " + fastq)
        subprocess.call("STAR --runThreadN 64 --genomeDir " + genome_directory + " --readFilesCommand zcat --outFilterType BySJout --outFilterMismatchNoverLmax 0.04 --outFilterMismatchNmax 999 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --outFilterMultimapNmax 20 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesIn "+ fastq_directory + fastq + " --clip3pAdapterSeq GATCGGAAGAGCACACGTCTGAACTCCAGTCAC --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix " + output_directory + prefix + "/", shell=True)