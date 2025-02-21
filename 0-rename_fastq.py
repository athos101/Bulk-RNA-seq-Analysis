import os
import sys
import re

sample_path = sys.argv[1]

def rename_fastq(name):
    name=re.sub(r".*\[(.*?)\]",r"\1",name)
    if("Adj" in name):
    	pattern=r"Pul0?(\d+)-Adj_S\d+_L0{2}(\d+)_R(\d+)_\d+"
    	replacement=r"Pul\1_Adj_L\2_R\3"
    elif("Tumor" in name):
    	pattern=r"Pul0?(\d+)-Tumor_S\d+_L0{2}(\d+)_R(\d+)_\d+"
    	replacement=r"Pul\1_Tumor_L\2_R\3"
    
    return re.sub(pattern, replacement, name)

for f in os.listdir(sample_path):
    old_path=os.path.join(sample_path,f)
    new_path=os.path.join(sample_path, rename_fastq(f))
    os.rename(old_path, new_path)
    print(f"Renamed: {f} --> {rename_fastq(f)}")
