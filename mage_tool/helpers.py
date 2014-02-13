# -*- coding: utf-8 -*-

try:
    from string import maketrans
except ImportError:
    maketrans = str.maketrans

import re
from math import log


"""
Constants
"""

amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

degenerate_nucleotides = (
    ("A", frozenset(["A"])),
    ("T", frozenset(["T"])),
    ("G", frozenset(["G"])),
    ("C", frozenset(["C"])),
    ("R", frozenset(["A", "G"])),
    ("M", frozenset(["A", "C"])),
    ("Y", frozenset(["C", "T"])),
    ("K", frozenset(["G", "T"])),
    ("S", frozenset(["C", "G"])),
    ("W", frozenset(["A", "T"])),
    ("B", frozenset(["C", "G", "T"])),
    ("D", frozenset(["A", "G", "T"])),
    ("H", frozenset(["A", "C", "T"])),
    ("V", frozenset(["A", "C", "G"])),
    ("N", frozenset(["A", "C", "G", "T"]))
)

nts_to_dgn = {nts: dgn for dgn, nts in degenerate_nucleotides}
dgn_to_nts = {dgn: nts for dgn, nts in degenerate_nucleotides}

default_codon_table = {
    'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N', 'ACA': 'T', 'ACC': 'T',
    'ACG': 'T', 'ACT': 'T', 'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
    'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I', 'CAA': 'Q', 'CAC': 'H',
    'CAG': 'Q', 'CAT': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'CTA': 'L', 'CTC': 'L',
    'CTG': 'L', 'CTT': 'L', 'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GGA': 'G', 'GGC': 'G',
    'GGG': 'G', 'GGT': 'G', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'TAA': '$', 'TAC': 'Y', 'TAG': '$', 'TAT': 'Y', 'TCA': 'S', 'TCC': 'S',
    'TCG': 'S', 'TCT': 'S', 'TGA': '$', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
    'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}

default_dgn_table = {'A': 'GCN', 'C': 'TGY', 'E': 'GAR', '$': 'TRR',
                     'G': 'GGN', 'F': 'TTY', 'I': 'ATH', 'H': 'CAY',
                     'K': 'AAR', 'M': 'ATG', 'L': 'YTN', 'N': 'AAY',
                     'Q': 'CAR', 'P': 'CCN', 'S': 'WSN', 'R': 'MGN',
                     'T': 'ACN', 'W': 'TGG', 'V': 'GTN', 'Y': 'TAY',
                     'D': 'GAY'}

#E. coli default codon usage
default_codon_usage = {'CTT': ['L', 0.12], 'ATG': ['M', 1.0],
                       'ACA': ['T', 0.13], 'ACG': ['T', 0.24],
                       'ATC': ['I', 0.35], 'ATA': ['I', 0.07],
                       'AGG': ['R', 0.03], 'CCT': ['P', 0.17],
                       'AGC': ['S', 0.33], 'AGA': ['R', 0.02],
                       'ATT': ['I', 0.58], 'CTG': ['L', 0.46],
                       'CTA': ['L', 0.05], 'ACT': ['T', 0.16],
                       'CCG': ['P', 0.55], 'AGT': ['S', 0.14],
                       'CCA': ['P', 0.14], 'CCC': ['P', 0.13],
                       'TAT': ['Y', 0.53], 'GGT': ['G', 0.29],
                       'CGA': ['R', 0.07], 'CGC': ['R', 0.44],
                       'CGG': ['R', 0.07], 'GGG': ['G', 0.12],
                       'TAG': ['$', 0.0], 'GGA': ['G', 0.13],
                       'TAA': ['$', 0.64], 'GGC': ['G', 0.46],
                       'TAC': ['Y', 0.47], 'CGT': ['R', 0.36],
                       'GTA': ['V', 0.17], 'GTC': ['V', 0.18],
                       'GTG': ['V', 0.4], 'GAG': ['E', 0.3],
                       'GTT': ['V', 0.25], 'GAC': ['D', 0.35],
                       'GAA': ['E', 0.7], 'AAG': ['K', 0.27],
                       'AAA': ['K', 0.73], 'AAC': ['N', 0.53],
                       'CTC': ['L', 0.1], 'CAT': ['H', 0.55],
                       'AAT': ['N', 0.47], 'CAC': ['H', 0.45],
                       'CAA': ['Q', 0.3], 'CAG': ['Q', 0.7],
                       'TGT': ['C', 0.42], 'TCT': ['S', 0.11],
                       'GAT': ['D', 0.65], 'TTT': ['F', 0.57],
                       'TGC': ['C', 0.58], 'TGA': ['$', 0.36],
                       'TGG': ['W', 1.0], 'TTC': ['F', 0.43],
                       'TCG': ['S', 0.16], 'TTA': ['L', 0.15],
                       'TTG': ['L', 0.12], 'TCC': ['S', 0.11],
                       'ACC': ['T', 0.47], 'TCA': ['S', 0.15],
                       'GCA': ['A', 0.21], 'GCC': ['A', 0.31],
                       'GCG': ['A', 0.38], 'GCT': ['A', 0.11]}

"""
DNA string tools
"""

def reverse_complement(x):
    """Reverse complement.

    Checks for an internal reverse_complement() method of the supplied object,
    otherwise performs a string translate. String translate retains case.

        >>> s = "AaTtGgCc"
        >>> reverse_complement(s)
        'gGcCaAtT'

    Internal method:

        >>> class Whatever:
        ...     def reverse_complement(self):
        ...         return "reverse complement"
        >>> w = Whatever()
        >>> reverse_complement(w)
        'reverse complement'

    """
    try:
        return x.reverse_complement()
    except AttributeError:
        #return x.translate(maketrans("ATGCatgc", "TACGtacg"))[::-1]
        return str(x).translate(maketrans("ATGCRYMKBVDHatgcrymkbvdh",
                                          "TACGYRKMVBHDtacgyrkmvbhd"))[::-1]


def reverse_complement_dgn(x):
    """Reverse complement a degenerate DNA sequence."""
    print("DeprecationWarning: use normal reverse_complement()")
    return reverse_complement(x)


def valid_na(seq):
    """Check for valid DNA/RNA.

        >>> valid_na("ATCGacTCUuUGACT")
        True

    """
    return bool(re.match(r"^[ATGCU]+$", str(seq), re.IGNORECASE))


def valid_rna(seq):
    """Check for valid RNA.

        >>> valid_rna("ATAGCAT")
        False
        >>> valid_rna("AUuagU")
        True

    """
    return bool(re.match(r"^[AGCU]+$", seq, re.IGNORECASE))


def valid_dna(seq):
    """Check for valid DNA.

        >>> valid_dna("AtagCAT")
        True
        >>> valid_dna("AUuagU")
        False

    """
    return bool(re.match(r"^[ATGC]+$", seq, re.IGNORECASE))


def valid_dgn(seq):
    """Check for valid degenerate DNA.

        >>> valid_dgn("NnYDAtagCAT")
        True
        >>> valid_dgn("AUuagU")
        False

    """
    return bool(re.match(r"^[ATGCRYMKSWBDHVN]+$", seq, re.IGNORECASE))


def valid_pept(seq):
    """Check for valid amino acids"""
    return bool(re.match(r"^[ACDEFGHIKLMNPQRSTVWY*]+$", seq))


def translate(seq, codon_table):
    seq = str(seq).upper().replace("U", "T")
    if not valid_dna(seq):
        raise ValueError("Not valid DNA/RNA: {}".format(seq))

    pept = list()
    for i in range(0, len(seq), 3):
        pept.append(codon_table[seq[i:i+3]])

    return "".join(pept)


def reverse_translate(seq, dgn_table):
    """Reverse translate a peptide sequence into a degenerate DNA sequence."""
    if not valid_pept(seq):
        raise ValueError("Not valid protein: {}".format(seq))

    dgn = list()
    for AA in seq:
        dgn.append(dgn_table[AA])

    return "".join(dgn)


def cds_to_wobble(seq, codon_table, dgn_table):
    """Convert a DNA sequence to a wobble (degenerate) DNA sequence."""
    seq = str(seq).upper().replace("U", "T")
    if not valid_dna(seq):
        raise ValueError("Not valid DNA/RNA: {}".format(seq))

    dgn = list()
    for i in range(0, len(seq), 3):
        AA = codon_table[seq[i:i+3]]
        dgn.append(dgn_table[AA])

    return "".join(dgn)


def seqs_to_degenerate(seqs):
    """Turn a list of sequences into one degenerate DNA sequence"""

    if not seqs:
        raise ValueError("No sequences supplied to create degenerate sequences.")

    tested = list()
    length = len(seqs[0])
    for seq in seqs:
        if len(seq) != length:
            raise ValueError("All sequences must be same length, expected {}, found {} in {}".format(length, len(seq), seq))

        seq = str(seq).upper().replace("U", "T")
        if not valid_dna(seq):
            raise ValueError("Not valid DNA/RNA: {}".format(seq))

        tested.append(seq)

    dgn_seq = list()
    for pool in zip(*seqs):
        dgn_seq.append(nts_to_dgn[frozenset(pool)])

    return "".join(dgn_seq)


def is_inside(db_start, db_end, q_start, q_end):
    """Checks whether q_start:q_end is inside db_start:db_end

    Checks for all possibilities:

        >>> is_inside(10, 15, 8, 12)
        True
        >>> is_inside(10, 15, 12, 18)
        True
        >>> is_inside(10, 15, 11, 14)
        True
        >>> is_inside(10, 15, 10, 15)
        True
        >>> is_inside(10, 15, 10, 18)
        True
        >>> is_inside(10, 15, 8, 17)
        True

    Return False if slices do not cross:

        >>> is_inside(10, 15, 5, 8)
        False
        >>> is_inside(10, 15, 16, 18)
        False
        >>> is_inside(10, 15, 5, 10)
        False
        >>> is_inside(10, 15, 15, 17)
        False

    If start > end, it means the string is circular:

        >>> is_inside(10, 15, 12, 2)
        True
        >>> is_inside(10, 15, 8, 2)
        True
        >>> is_inside(10, 15, 17, 2)
        False
        >>> is_inside(10, 2, 12, 18)
        True
        >>> is_inside(10, 2, 12, 1)
        True
        >>> is_inside(10, 4, 6, 8)
        False

    Same goes for negative values:

        >>> is_inside(10, 15, -6, 8)
        False

    """
    if db_start > db_end:
        db_end += max(db_start, q_start, q_end)
    if q_start > q_end:
        q_end += db_end

    if q_start < db_start < q_end:
        return True
    if q_start < db_end < q_end:
        return True
    if q_start >= db_start and q_end <= db_end:
        return True

    return False


def extract_circular(parent, start, end):
    """Extract start:end from a circular sequence.

    TODO: Doctest

    """
    if start >= 0:
        extr = parent[start:end]
    else:
        #Leader extends backwards beyond 0
        extr = parent[start:] + parent[:end]

    #Leader extends beyond len(genome)
    if end > len(parent):
        extr += parent[:end-len(parent)]

    return extr


default_primer_constants = {
    #      kcal/mol   cal/k*mol
    "AA": {"H": -7.9, "S": -22.2},
    "AT": {"H": -7.2, "S": -20.4},
    "TA": {"H": -7.2, "S": -21.3},
    "CA": {"H": -8.5, "S": -22.7},
    "GT": {"H": -8.4, "S": -22.4},
    "CT": {"H": -7.8, "S": -21.0},
    "GA": {"H": -8.2, "S": -22.2},
    "CG": {"H": -10.6,"S": -27.2},
    "GC": {"H": -9.8, "S": -24.4},
    "GG": {"H": -8.0, "S": -19.9},
    #Initiation w GC
    "G": {"H": 0.1, "S": -2.8},
    #Initiation w AT
    "A": {"H": 2.3, "S":  4.1},
    #Symmetry penalty
    "sym": {"H": 0,   "S": -1.4},
    "dS_salt": 0.368,
}


def primer_tm(primer, salt=0.05, c=5e-8, constants=None):
    """Calculate primer Tm according to SantaLucia 1998.

    salt is the concentration of monovalent cations in mol/L, and c is the
    primer concentration in mol/L. constants are given as a dict of neighbors.
    See default_primer_constants for an example.

    Full citation:
        SantaLucia, J. A unified view of polymer, dumbbell, and oligonucleotide
        DNA nearest-neighbor thermodynamics. Proceedings of the National
        Academy of Sciences of the United States of America 95, 1460–5 (1998).

    """
    prm = str(primer).upper()
    if not valid_dna(prm):
        raise ValueError("A valid DNA sequence (no degenerate NTs) must be used "
                         "with primer_tm. Got: '{}'.".format(primer))

    #Load constants or default constants
    cnst = constants if constants else default_primer_constants

    #First termial pairing penalty
    st = prm[0] if prm[0] in cnst else reverse_complement(prm[0])
    if st in cnst:
        dH = cnst[st]["H"]
        dS = cnst[st]["S"]

    #Symmetry penalty
    half = int(len(prm)/2)
    if prm[:half] == prm[-half:][::-1]:
        dH += cnst["sym"]["H"]
        dS += cnst["sym"]["S"]

    #Nearest neighbors.
    #Last terminal penalty is included.
    for i in range(len(prm)):
        pair = prm[i:i+2]
        pair = pair if pair in cnst else reverse_complement(pair)
        if pair not in cnst: continue
        dH += cnst[pair]["H"]
        dS += cnst[pair]["S"]

    NN = len(prm) - 1
    dS_salt = dS + cnst["dS_salt"] * NN * log(salt)

    R = 1.987 #cal/molK
    return dH*1000/(dS_salt + R*log(c/4)) - 273.15


def make_primer(ref_seq, temp, start, end, salt_c=0.05, primer_c=5e-8, max_len=200):
    """Extract a primer from ref_seq with Tm of temp."""
    for i in range(max_len):
        pr = extract_circular(ref_seq, start, end)
        Tm = primer_tm(pr, salt_c, primer_c)
        if Tm > temp:
            break
        #Keep going
        start -= 1
    else:
        raise Exception("Target temp ({}) not reached for forward primer. "
                        "(Reached: {})".format(temp, Tm))
    return pr, Tm


def make_rev_primer(ref_seq, temp, start, end, salt_c=0.05, primer_c=5e-8, max_len=200):
    """Extract a primer from ref_seq with Tm of temp."""
    for i in range(max_len):
        pr = extract_circular(ref_seq, start, end)
        Tm = primer_tm(pr, salt_c, primer_c)
        if Tm > temp:
            break
        #Keep going
        end += 1
    else:
        raise Exception("Target temp ({}) not reached for forward primer. "
                        "(Reached: {})".format(temp, Tm))
    return pr, Tm


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    import doctest
    doctest.testmod()
