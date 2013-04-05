#!/usr/bin/env python

"""
DOC
"""

import argparse

from mage_tool.translation import RBS_predict
from mage_tool.translation import dG_to_AU

if __name__ == "__main__":
    #Read command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("pre", help="")
    parser.add_argument("cds", help="")
    #parser.add_argument("--old", "--start", help="", type=int, default=0)
    args = parser.parse_args()

    dG = RBS_predict(args.pre, args.cds, verbose=True)

    print
    print "Expression:", dG_to_AU(dG)