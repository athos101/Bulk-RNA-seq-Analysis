import os
import pandas as pd

# Diretório onde as pastas estão localizadas
base_directory = "/caminho/do/diretorio"  # Substitua pelo caminho correto

# Coluna desejada para extrair as contagens de genes
# 0 = GeneID, 1 = Unstranded, 2 = Forward, 3 = Reverse
desired_column = 1  # Exemplo: Unstranded
samples_counts = {}

for folder_name in os.listdir(base_directory):
    if folder_name.endswith("_output"):
        folder_path = os.path.join(base_directory, folder_name)
        file_path = os.path.join(folder_path, "ReadsPerGene.out.tab")
        if os.path.exists(file_path):
            df = pd.read_csv(file_path, sep="\t", header=None)
            gene_ids = df[0]
            counts = df[desired_column]
            sample_name = folder_name.replace("_output", "")
            samples_counts[sample_name] = counts

result_df = pd.DataFrame(samples_counts, index=gene_ids)
result_df = result_df[~result_df.index.str.startswith("N_")]
output_file = os.path.join(base_directory, "raw_counts.csv")
result_df.to_csv(output_file, sep=",")

print(f"Arquivo 'raw_counts.csv' salvo em: {output_file}")