
python ../sequence_align.py --temp seq.txt --query seq.txt --out out0.tsv
python ../sequence_align.py --temp 3BEX.pdb --query seq.txt --out out1.tsv
python ../sequence_align.py --temp seq.txt --query 3BEX.pdb --out out2.tsv
python ../sequence_align.py --temp 3BEX.pdb --query 3BEX.pdb --out out3.tsv
