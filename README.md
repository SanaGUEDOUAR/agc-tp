# OTU Clustering Tool

This tool performs OTU (Operational Taxonomic Unit) clustering from amplicon sequences using a greedy algorithm based on sequence abundance and identity. It allows for the dereplication of sequences and the identification of unique OTUs.

## Table of Contents
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Output](#output)
- [Quality Assessment](#quality-assessment)
- [License](#license)
- [Contact](#contact)

## Features
- Reads compressed FASTA files (.fasta.gz).
- Dereplicates sequences based on length and count.
- Identifies OTUs using a greedy clustering approach.
- Outputs results in a FASTA format.
- Allows for quality assessment using `vsearch` against reference databases.

## Requirements
- Python 3.x
- Required Python packages:
  - `argparse`
  - `gzip`
  - `textwrap`
  - `pathlib`
  - `collections`
  - `nwalign3`

### Create a virtual environment
It is recommended to create a virtual environment to install the required packages. You can create a virtual environment using `conda`:
```bash
conda create -n agc python=3.9 pip numpy
conda activate agc
pip install nwalign3 pytest pylint pytest-cov
```

## Installation
1. Clone the repository:
```bash
git clone git@github.com:SanaGUEDOUAR/agc-tp.git
cd agc-tp
```

2. Install any additional dependencies if needed.

## Usage
```bash
python agc/agc.py -i <amplicon_file> -s <minseqlen> -m <mincount> -o <output_file>
```
### Arguments
- -i, --amplicon_file: Path to the compressed FASTA file containing the amplicon sequences (required).
- -s, --minseqlen: Minimum sequence length for dereplication (default: 400).
- -m, --mincount: Minimum count for dereplication (default: 10).
- -o, --output_file: Path to the output file (default: OTU.fasta).

### Example
```bash
python agc/agc.py -i amplicons.fasta.gz -s 400 -m 10 -o results/OTU.fasta
```

## Output
The program generates a FASTA file containing the identified OTUs. Each OTU is labeled with an occurrence count.

## Quality Assessment
To assess the quality of the resulting OTUs, you can use the vsearch tool. Install vsearch using:
```bash
conda install vsearch
```
Then run the following command to align the OTUs against a reference database:
```bash
vsearch --usearch_global OTU.fasta --db data/mock_16S.fasta --id 0.8 --blast6out results/resultat.tsv
```
This command produces a resultat.tsv file containing the alignment results.

## License
This project is licensed under the GNU General Public License v3.0. See the [LICENSE](https://www.gnu.org/licenses/gpl-3.0.html) file for more details.

## Contact
For questions or feedback, please contact:

Author: Sana Guedouar
Email: [sana.guedouar@etu.u-paris.fr](mailto:sana.guedouar@etu.u-paris.fr)
University: Université Paris Cité