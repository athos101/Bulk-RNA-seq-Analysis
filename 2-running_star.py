import subprocess
import os
import sys

output_directory = "bulk_data/control/counts"
genome_directory = "bulk_data/Genome"
fastq_directory = "bulk_data/control/fastq"

subprocess.call(f"STAR --runThreadN 16 --runMode genomeGenerate --genomeDir {genome_directory} --genomeFastaFiles {genome_directory}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --sjdbGTFfile {genome_directory}/Homo_sapiens.GRCh38.gtf", shell=True)
if (not os.path.isdir(output_directory)):
    os.mkdir(output_directory)

f=[]

for fastq in os.listdir(fastq_directory):
    for file in os.listdir(fastq_directory+"/"+fastq):
        name=fastq_directory+"/"+fastq+"/"+file
        f.append(name)
        print(name)
    print ("Currently mapping: " + fastq)
    # run STAR on the current fastq file
    subprocess.call("STAR --runThreadN 16 --genomeDir " + genome_directory + " --outFilterType BySJout --outFilterMismatchNoverLmax 0.04 --outFilterMismatchNmax 999 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --outFilterMultimapNmax 20 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesIn " + f[0] +" "+ f[1] + " --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix " + output_directory +"/"+ fastq + "/", shell=True)
