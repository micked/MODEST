#!/usr/bin/env python

"""
DOC
"""

import argparse
import re

from mage_tool.RBS.RBS_MC_Design import Monte_Carlo_Design
from mage_tool.RBS.RBS_Calculator import NoRBSError


if __name__ == "__main__":
    #Read command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("pre_seq", help="Pre/leader sequence")
    parser.add_argument("cds", help="Coding sequence")
    parser.add_argument("target", help="Target Initialisation rate", type=float)
    args = parser.parse_args()

    #Check Sequence and TIR Inputs
    valid_start_codons = ('ATG','AUG','GTG','GUG')

    TIR = args.target
    pre_seq = args.pre_seq
    post_seq = args.cds

    if TIR <= 0: raise ValueError("Invalid Target Initialisation Rate.")

    nuc_regex = re.compile(r"^[ATGCU]+$", re.IGNORECASE)
    if not nuc_regex.match(pre_seq):
        raise ValueError("Invalid pre sequence")
    if not nuc_regex.match(post_seq) or post_seq[0:3].upper() not in valid_start_codons:
        raise ValueError("Invalid CDS/post sequence")

    total_iterations = 0
    max_iter = 10000
    iterations = max_iter
    while iterations == max_iter:
        try:
            TIR_out, RBS_out, estimator, iterations = Monte_Carlo_Design(pre_seq, post_seq, RBS_init=None, TIR_target=TIR, dG_target=None, MaxIter=max_iter, verbose=False)
            total_iterations += iterations
        except NoRBSError:
            print "caught noRBS"
            continue
        print total_iterations, "iterations"

    print "Program Executed" + "\n" + RBS_out + "\n" + str(TIR_out)