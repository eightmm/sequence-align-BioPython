import os
import argparse
import numpy as np
import pandas as pd

from tqdm import tqdm

from Bio.PDB import PDBParser
from Bio.Align import substitution_matrices, PairwiseAligner

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

def clean_sequence(seq):
    return "".join([aa for aa in seq if aa in VALID_AA])

def get_sequence_from_pdb(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    sequence = []

    aa_three_to_one = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'CYX': 'C',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'HIE': 'H', 'HIP': 'H', 'HID': 'H',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ' and residue.resname in aa_three_to_one:
                    sequence.append(residue.resname)

    return clean_sequence(''.join([aa_three_to_one[res] for res in sequence]))

def extract_sequences_from_fasta(file_path):
    sequences, names = [], []
    with open(file_path, 'r') as f:
        seq, name = "", ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    seq = clean_sequence(seq)
                    if seq:
                        sequences.append(seq)
                        names.append(name)
                name = line[1:].split("|")[0]
                seq = ""
            else:
                seq += line
        if seq:
            seq = clean_sequence(seq)
            if seq:
                sequences.append(seq)
                names.append(name)
    return sequences, names

def extract_sequences_from_text(file_path):
    sequences, names = [], []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if ',' in line:
                name, seq = line.split(',', 1)
                seq = clean_sequence(seq)
                if seq:
                    names.append(name.strip())
                    sequences.append(seq.strip())
            elif line:
                seq = clean_sequence(line)
                if seq:
                    sequences.append(seq)
                    names.append(f"seq{len(names) + 1}")
    return sequences, names

def extract_sequences(file_path):
    file_ext = os.path.splitext(file_path)[1].lower()
    if file_ext == '.pdb':
        return [get_sequence_from_pdb(file_path)], [os.path.basename(file_path).split('.')[0]]
    elif file_ext in ['.fa', '.fasta']:
        return extract_sequences_from_fasta(file_path)
    else:
        return extract_sequences_from_text(file_path)

def calculate_identity(seq1, seq2):
    substitution_matrices.load()
    matrix = substitution_matrices.load("BLOSUM62")
    aligner = PairwiseAligner()
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(seq1, seq2)[0]
    return (alignment.counts().identities / alignment.length) * 100

def compare_sequences(seq_list1, seq_list2):
    output = np.zeros((len(seq_list1), len(seq_list2)))
    for i, seq1 in tqdm( enumerate(seq_list1) ):
        for j, seq2 in enumerate(seq_list2):
            output[i, j] = round(calculate_identity(seq1, seq2), 2)
    return output

def main(args):
    temp_sequences, temp_names = extract_sequences(args.temp)
    query_sequences, query_names = extract_sequences(args.query)
    identity_matrix = compare_sequences(temp_sequences, query_sequences)
    df = pd.DataFrame(identity_matrix, index=temp_names, columns=query_names)
    df.to_csv(args.output, sep='\t', index=True)
    print(f"Identity matrix saved to {args.output}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate sequence identity from a temp and query file.')
    parser.add_argument('--temp', type=str, required=True, help='Input temp (sequence, FASTA, or PDB file).')
    parser.add_argument('--query', type=str, required=True, help='Input query (sequence, FASTA, or PDB file).')
    parser.add_argument('--output', type=str, required=True, help='Output tsv file.')
    args = parser.parse_args()
    main(args)

