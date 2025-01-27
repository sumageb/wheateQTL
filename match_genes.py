#!/usr/bin/env python3

# Set input file paths
file1 = "modified_geneversion.txt"
file2 = "hq_homeologous_genes.txt"
output_file = "matched_genes.txt"

# Read all genes from file1 into a set for fast lookup
gene_set = set()
with open(file1, 'r') as f1:
    for line in f1:
        gene_id = line.strip()
        gene_set.add(gene_id)

# Debugging: Print the size of the gene set
print(f"Total unique genes in file1: {len(gene_set)}")

# Open the output file
with open(output_file, 'w') as out_f:
    # Read each line from file2 and check if all genes are present in the set
    with open(file2, 'r') as f2:
        header = f2.readline().strip()  # Read the header
        out_f.write(header + '\n')  # Write the header to the output file

        for line in f2:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                genes1 = parts[0].strip('"').split(",")
                genes2 = parts[1].strip('"').split(",")
                genes3 = parts[2].strip('"').split(",")

                # Check if all genes in each part are in the gene_set
                all_genes = genes1 + genes2 + genes3
                if all(gene in gene_set for gene in all_genes if gene):
                    # Debugging: Print the matched line
                    print(f"Matched line: {line.strip()}")
                    out_f.write(line.strip() + '\n')

print(f"Matching complete. Results are saved in {output_file}")

