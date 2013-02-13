from Bio.Seq import Seq

from helpers import reverse_complement

def mut_to_oligo(mut, ref_genome, oligo_len=90):
    genome_slice_len = oligo_len - len(mut.after)
    
    out = "AAA"
    
    out = target_lagging_strand(out, mut.position, ori=3886229, ter=1493223)
    
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
