import os
import sys

sample_path = sys.argv[1]

for f in os.listdir(sample_path):
    print(f)