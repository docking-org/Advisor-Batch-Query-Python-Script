# Advisor Batch Query

A script that calls arthor.docking.org's similarity search API on the aggregator dataset.

## Instructions

1. Install rdkit
   - `pip install rdkit`
1. Enter a list of smiles separated by new lines and put them in "smiles.txt" file.
1. Run the script
   - `python advisor_batch.py`
1. Choose an affinity from strong, standard, or weak.
1. The script should create a file called results.txt.
