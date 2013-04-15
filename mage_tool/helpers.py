try:
    from string import maketrans
except ImportError:
    maketrans = str.maketrans

import re

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


if __name__ == "__main__":
    import doctest
    doctest.testmod()