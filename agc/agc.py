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
import os
import gzip
import statistics
import textwrap
import numpy as np
from pathlib import Path
from collections import Counter
from typing import Iterator, Dict, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw
np.int = int

__author__ = "Nicolas Salaun"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Nicolas Salaun"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Nicolas Salaun"
__email__ = "nicolas.salaun@orange.fr"
__status__ = "Developpement"



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
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
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
    with gzip.open(amplicon_file, "rt") as f:
        seq = ""
        for line in f:
            if line.startswith(">"):
                if len(seq) >= minseqlen:
                    yield seq
                seq = ""
            else:
                seq += line.strip()
        yield seq if len(seq) >= minseqlen else ""



def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of sequence with a count >= mincount and a length >= minseqlen.
    """
    seq_count = Counter()
    for sequence in read_fasta(amplicon_file, minseqlen):
        seq_count[sequence] += 1

    for seq, count in seq_count.items():
        if count >= mincount:
            yield [seq, count]

    
    
    
def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    return sum(1 for a, b in zip(*alignment_list) if a == b) / len(alignment_list[0]) * 100
    

def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int, chunk_size: int, kmer_size: int) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    seqs_gen = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    OTU_bank = [next(seqs_gen)]
    
    for i, seq_info in enumerate(seqs_gen):
        if i % 100 == 0:
            print(i)
        if all(get_identity(nw.global_align(seq_info[0], otu_seq_info[0], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))) <= 97 for otu_seq_info in OTU_bank):
            OTU_bank.append([seq_info[0],seq_info[1]])

    return OTU_bank


def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, "w") as f:
        for i, otu in enumerate(OTU_list):
            f.write(f">OTU_{i+1} occurrence:{otu[1]}\n")
            f.write(f"{otu[0]}\n")



#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    args = get_arguments()
    otu = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, 0, 0)
    write_OTU(otu, args.output_file)


if __name__ == '__main__':
    main()
