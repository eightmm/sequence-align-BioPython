
# sequence-identity-calculation-BioPython

This repository contains a Python script to calculate sequence identity between sequences extracted from PDB files or sequence files using BioPython. The identity matrix is saved as a `.tsv` file.

## Project Structure

```
.
├── LICENSE
├── README.md
├── sequence_align.py
├── example
│   ├── 3BEX.pdb
│   ├── seq.txt
│   ├── output.tsv
│   └── run.sh
```

- `sequence_align.py`: Main Python script to calculate sequence identity.
- `example/3BEX.pdb`: Example PDB file to extract protein sequences.
- `example/seq.txt`: Example file containing multiple protein sequences.
- `example/output.tsv`: Output file with the identity matrix.
- `example/run.sh`: Example bash script to run the alignment process.

## Installation

Before running the script, ensure that you have Python installed and the required dependencies.

1. Install Python (version 3.x recommended).
2. Install the required Python packages using the following command:

```bash
pip install biopython pandas numpy
```

## Usage

You can calculate sequence identity between two files (PDB or sequence text files) by running the script ` sequence-identity.py`.

### Example Usage

In the `example` folder, an example of how to run the script using a PDB file (`3BEX.pdb`) and a sequence file (`seq.txt`) is provided.

1. **Run the Example Script**

You can use the provided `run.sh` script to run the example:

```bash
bash example/run.sh
```

This will execute the `sequence_align.py` script using the example files.

### Command Line Usage

You can also run the script directly from the command line as follows:

```bash
python sequence_align.py --temp example/3BEX.pdb --query example/seq.txt --output example/output.tsv
```

- `--temp`: Path to the first input file (either a PDB file or a sequence text file or a fasta file).
- `--query`: Path to the second input file (either a PDB file or a sequence text file or a fasta file).
- `--output`: Path to save the output `.tsv` file.

### Output

The identity matrix will be saved in the specified `.tsv` file. For example, after running the provided example, you can view the output in `example/output.tsv`.

### Example Output

```
    seq1   seq2   seq3
seq1  100.00  73.44   58.62
seq2  73.44   100.00  71.43
seq3  58.62   71.43   100.00
```

This matrix represents the percentage identity between the sequences extracted from the input files.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

