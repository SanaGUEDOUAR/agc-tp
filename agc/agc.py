#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import gzip
import textwrap
from pathlib import Path
from collections import Counter
from typing import Iterator, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Sana GUEDOUAR"
__copyright__ = "Universite Paris Cité"
__credits__ = ["Sana GUEEDOUAR"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Sana GUEDOUAR"
__email__ = "sana.guedouar@etu.u-paris.fr"
__status__ = "Developpement"

#==============================================================
def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=f"{sys.argv[0]} -h")

    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=Path,
                        default=Path("OTU.fasta"), help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """
    header = ""
    sequence = ""
    for line in gzip.open(amplicon_file, 'rt'):
        if line.startswith(">"):
            if len(sequence) >= minseqlen:
                yield sequence
            header = line.strip()
            sequence = ""
        else:
            sequence += line.strip()
    if len(sequence) >= minseqlen:  # Yield the last sequence
        yield sequence


def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequences and return unique sequences sorted by count.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum count for the sequence to be considered.
    :return: A generator object that provides a [sequence, count] list with sequences
             having a count >= mincount, sorted by occurrence in descending order.
    """
    # Step 1: Count the occurrences of each sequence that meets the length requirement
    sequence_counts = Counter(read_fasta(amplicon_file, minseqlen))

    # Step 2: Sort the sequences by count in descending order
    sorted_sequences = sorted(sequence_counts.items(), key=lambda x: x[1], reverse=True)

    # Step 3: Yield only the sequences with counts >= mincount
    for sequence, count in sorted_sequences:
        if count >= mincount:
            yield [sequence, count]


def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences
                            in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The percentage of identity between the two sequences.
    """
    seq1, seq2 = alignment_list
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != '-' and b != '-')  # Ignore gaps
    identity = matches / len(seq1) * 100
    return identity



def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int,
                                chunk_size: int, kmer_size: int) -> List:
    """Compute an abundance greedy clustering
        regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    otu_list = []
    for sequence, count in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        is_new_otu = True
        for otu in otu_list:
            alignment = nw.global_align(sequence, otu[0], gap_open=-1, gap_extend=-1,
                                        matrix=str(Path(__file__).parent / "MATCH"))
            identity = get_identity(alignment)
            if identity >= 97:  # Vérifier l'identité supérieure à 97%
                is_new_otu = False
                break
        if is_new_otu:
            otu_list.append([sequence, count])  # Ajout de la nouvelle OTU
    return otu_list


def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, "w", encoding="utf-8") as f:
        for i, (sequence, count) in enumerate(OTU_list):
            f.write(f">OTU_{i+1} occurrence:{count}\n")
            f.write(f"{textwrap.fill(sequence, width=80)}\n")

#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Step 1: Perform OTU clustering
    print("Performing OTU clustering...")
    OTU_list = abundance_greedy_clustering(
        args.amplicon_file, args.minseqlen, args.mincount, chunk_size=100, kmer_size=8
    )
    # Step 2: Write OTU results to the output file
    print(f"Writing OTU results to {args.output_file}...")
    write_OTU(OTU_list, args.output_file)
    print("Done.")


if __name__ == '__main__':
    main()
