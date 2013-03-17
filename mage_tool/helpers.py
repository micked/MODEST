from string import maketrans
import re

"""
Constants
"""

amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

degenerate_nucleotides = {
    "R": ["A", "G"],
    "M": ["A", "C"],
    "Y": ["C", "T"],
    "K": ["G", "T"],
    "S": ["C", "G"],
    "W": ["A", "T"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"]
}

"""
DNA string tools
"""
    
def complement(x):
    """Complement a DNA string"""
    try:
        return x.complement()
    except AttributeError:
        return x.translate(maketrans("ATGCatgc", "TACGtacg"))


def reverse_complement(x):
    """Reverse complement."""
    try:
        return x.reverse_complement()
    except AttributeError:
        return x.translate(maketrans("ATGCatgc", "TACGtacg"))[::-1]


def reverse_complement_dgn(x):
    """Reverse complement a degenerate DNA sequence."""
    return str(x).translate(maketrans("ATGCRYMKBVDHatgcrymkbvdh",
                                      "TACGYRKMVBHDtacgyrkmvbhd"))[::-1]


def valid_na(seq):
    """Check for valid DNA/RNA"""
    return bool(re.match(r"^[ATGCU]+$", seq, re.IGNORECASE))


def valid_rna(seq):
    """Check for valid RNA"""
    return bool(re.match(r"^[AGCU]+$", seq, re.IGNORECASE))


def valid_dna(seq):
    """Check for valid DNA"""
    return bool(re.match(r"^[ATGC]+$", seq, re.IGNORECASE))


def valid_dgn(seq):
    """Check for valid degenerate DNA"""
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


def cds_to_wobble(seq, options):
    """Convert a DNA sequence to a wobble (degenerate) DNA sequence."""
    seq = str(seq).upper().replace("U", "T")
    if not valid_dna(seq):
        raise ValueError("Not valid DNA/RNA: {}".format(seq))

    codon_table = options["codon_table"]
    dgn_table = options["dgn_table"]

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

    dgn_lookup = {
        frozenset(["A"]): "A",
        frozenset(["T"]): "T",
        frozenset(["G"]): "G",
        frozenset(["C"]): "C",
        frozenset(["A", "G"]): "R",
        frozenset(["A", "C"]): "M",
        frozenset(["C", "T"]): "Y",
        frozenset(["G", "T"]): "K",
        frozenset(["C", "G"]): "S",
        frozenset(["A", "T"]): "W",
        frozenset(["C", "G", "T"]): "B",
        frozenset(["A", "G", "T"]): "D",
        frozenset(["A", "C", "T"]): "H",
        frozenset(["A", "C", "G"]): "V",
        frozenset(["A", "C", "G", "T"]): "N"
    }

    dgn_seq = list()
    for pool in zip(*seqs):
        dgn_seq.append(dgn_lookup[frozenset(pool)])

    return "".join(dgn_seq)