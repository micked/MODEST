#!/usr/bin/env python

"""
DOC
"""

import argparse

from mage_tool.RBS.RBS_Calculator import RBS_Calculator

if __name__ == "__main__":
    #Read command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("seq", help="")
    parser.add_argument("-p", "--start", help="", type=int, default=0)
    args = parser.parse_args()

    seq = args.seq
    start = args.start

    if not start:
        start_range = [0, len(seq)]
    else:
        start_range = [int(start), int(start)+1]

    name = "no name"

    #Create instance of RBS Calculator
    calcObj = RBS_Calculator(seq, start_range, name)
    calcObj.calc_dG()

    dG_total_list = calcObj.dG_total_list[:]
    start_pos_list = calcObj.start_pos_list[:]
    kinetic_score_list = calcObj.kinetic_score_list[:]

    expr_list = []
    for dG in dG_total_list:
        expr_list.append(calcObj.calc_expression_level(dG))

    print len(expr_list)
    for (expr,start_pos,ks) in zip(expr_list,start_pos_list,kinetic_score_list):
        print start_pos, expr, ks
    