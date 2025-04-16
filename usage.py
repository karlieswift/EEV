"""
@Env: /anaconda3/python3.11
@Time: 2025/4/16-10:51
@Auth: karlieswift
@File: usage.py
@Desc: 
"""
from eev.eev import EnergyEntropy

# DNA
dna_seq = 'AGCTGAGTAAGATCA'
energy_values = 2
eev = EnergyEntropy(data_type='dna', energy_values=energy_values)
vector = eev.seq2vector(dna_seq)
print("DNA Sequence Vector:")
print(vector)

# Protein
protein_seq = 'PGGSLTLSCAASEPVFEANTMGWYRQAPGKQRELVASIRNVGGTNYADSVKGRF'
eev = EnergyEntropy(data_type='protein', energy_values=energy_values)
vector = eev.seq2vector(protein_seq)
print("\nProtein Sequence Vector:")
print(vector)

# RNA
rna_seq = 'ACGUAAAUCG'
eev = EnergyEntropy(data_type='rna', energy_values=energy_values)
vector = eev.seq2vector(rna_seq)
print("\nRNA Sequence Vector:")
print(vector)
