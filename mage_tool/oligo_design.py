#!/usr/bin/env python

"""
Module for designing DNA oligos from mutation objects
"""

from Bio.Seq import Seq

from helpers import reverse_complement


def mut_to_oligo(mut, ref_genome, oligo_len=90):
    
    post_seq_len = (oligo_len - len(mut.after) +0)/2
    pre_seq_len = oligo_len - post_seq_len - len(mut.after)
    
    offset = 1
    #if not len(mut.before): offset = 0
    pre_seq = ref_genome[ mut.position-pre_seq_len-1 : mut.position-offset ]
    post_seq = ref_genome[ mut.position+len(mut.before)-offset : mut.position+len(mut.before)+post_seq_len-1 ]
    
    out = pre_seq + mut.after.lower() + post_seq
    
    return out


def barcoding(oligo, library_name, barcode_library):
    return oligo


def target_lagging_strand(oligo, mut_pos, ori, ter):
    """TODO"""
    #ori = 3886229
    #ter = 1493223
    if 0 < mut_pos < ter or ori < mut_pos:
        #reverse complemented to target lagging strand
        #swapcase to id the ones that are reverse complemented
        out = reverse_complement(oligo)
        try:
            out = out.swapcase()
        except AttributeError:
            out = Seq(str(out).swapcase())
        return out
    elif ter < mut_pos < ori:
        #do nothing, the oligo is targetting lagging strand already
        return oligo
    else:
        raise Exception("Error locating lagging stand")
