from helpers import reversecomplement

def mut_to_oligo(mut):
    pass


def barcoding(oligo, library_name, barcode_library):
    pass


def target_lagging_strand(oligo, pos, ori, ter):
    """TODO"""
	#ori = 3886229
	#ter = 1493223
	if 0 < mut_pos < ter or ori < mut_pos:
		#reverse complemented to target lagging strand
		return reversecomplement(oligo).swapcase() #swapcase to id the ones that are reverse complemented
	elif ter < mut_pos < ori:
		#do nothing, the oligo is targetting lagging strand already
		return oligo_mut 
	else:
		raise Exception("Error locating lagging stand")
