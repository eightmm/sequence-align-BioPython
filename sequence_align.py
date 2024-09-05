import os
import argparse
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.Align import substitution_matrices, PairwiseAligner

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

    return ''.join([aa_three_to_one[res] for res in sequence])

def calculate_identity(seq1, seq2):
    substitution_matrices.load()
    matrix = substitution_matrices.load("BLOSUM62")

    aligner = PairwiseAligner()
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(seq1, seq2)
    best_alignment = alignment[0]
    aligned_length = best_alignment.length
    identity_score = best_alignment.counts().identities
    return (identity_score / aligned_length) * 100

def compare_sequences(seq_list1, seq_list2):
    output = np.zeros((len(seq_list1), len(seq_list2)))

    for i, seq1 in enumerate(seq_list1):
        for j, seq2 in enumerate(seq_list2):
            output[i, j] = round(calculate_identity(seq1, seq2), 2)  # 소숫점 두 자리까지 반올림
    return output

def extract_sequences(input_path):
    if os.path.splitext(input_path)[1].lower() == '.pdb':
        return [get_sequence_from_pdb(input_path)]
    else:
        with open(input_path, 'r') as f:
            sequences = [line.strip() for line in f if line.strip()]
        return sequences


def extract_sequences(input_path):
    if os.path.splitext(input_path)[1].lower() == '.pdb':
        return [get_sequence_from_pdb(input_path)], [os.path.basename(input_path).split('.')[0]]
    else:
        with open(input_path, 'r') as f:
            sequences = [line.strip() for line in f if line.strip()]
        return sequences, [f"{i}" for i in range(len(sequences))]


def main(args):
    temp_sequences, temp_names = extract_sequences(args.temp)
    query_sequences, query_names = extract_sequences(args.query)

    if len(temp_names) > 1:
        temp_names = [f"seq{i+1}" for i in range(len(temp_sequences))]
    
    if len(query_names) > 1:
        query_names = [f"seq{i+1}" for i in range(len(query_sequences))]

    identity_matrix = compare_sequences(temp_sequences, query_sequences)

    df = pd.DataFrame(identity_matrix, index=temp_names, columns=query_names)
    df.to_csv(args.output, sep='\t', index=True)

    print(f"Identity matrix saved to {args.output}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate sequence identity from a temp and query file.')
    parser.add_argument('--temp', type=str, required=True, help='Input temp (sequence or PDB file).')
    parser.add_argument('--query', type=str, required=True, help='Input query (sequence or PDB file).')
    parser.add_argument('--output', type=str, required=True, help='Output tsv file.')

    args = parser.parse_args()
    main(args)

