dataset = {'A': {'ALA': 0, 'ARG': 15, 'ASN': 0, 'ASP': 22, 'CYS': 0, 'GLN': 0, 'GLY': 0, 'HIS': 0, 'ILE': 0, 'LEU': 2, 'LYS': 8, 'MET': 12, 'PHE': 58, 'PRO': 0, 'SER': 18, 'THR': 8, 'TRP': 0, 'TYR': 0, 'VAL': 8, 'GLU': 0}, 'U': {'ALA': 17, 'ARG': 26, 'ASN': 16, 'ASP': 12, 'CYS': 0, 'GLN': 0, 'GLY': 0, 'HIS': 0, 'ILE': 2, 'LEU': 21, 'LYS': 44, 'MET': 7, 'PHE': 0, 'PRO': 0, 'SER': 0, 'THR': 0, 'TRP': 0, 'TYR': 0, 'VAL': 8, 'GLU': 0}, 'C': {'ALA': 7, 'ARG': 1, 'ASN': 0, 'ASP': 0, 'CYS': 0, 'GLN': 13, 'GLY': 0, 'HIS': 0, 'ILE': 0, 'LEU': 0, 'LYS': 69, 'MET': 0, 'PHE': 11, 'PRO': 0, 'SER': 0, 'THR': 4, 'TRP': 0, 'TYR': 39, 'VAL': 0, 'GLU': 0}, 'G': {'ALA': 0, 'ARG': 17, 'ASN': 11, 'ASP': 37, 'CYS': 0, 'GLN': 42, 'GLY': 6, 'HIS': 0, 'ILE': 1, 'LEU': 1, 'LYS': 16, 'MET': 0, 'PHE': 0, 'PRO': 0, 'SER': 6, 'THR': 0, 'TRP': 0, 'TYR': 4, 'VAL': 4, 'GLU': 0}} 
total_interactions = 0
for rna_base in dataset.keys():
    for aa in dataset[rna_base].keys():
        total_interactions += dataset[rna_base][aa]

print("Total number of interactions observed (sumNIpq):", total_interactions)

for aa in dataset['A'].keys():
    total_interactions = 0
    for rna_base in dataset.keys():
        if aa in dataset[rna_base]:
            total_interactions += dataset[rna_base][aa]
    print("Total number of interactions of", aa, "with any RNA molecule (Nsp):", total_interactions)

