#!/usr/bin/env python2

"""
Commandline interface to <>

Yadda yadda yadda
"""
from __future__ import print_function
import os.path
import argparse

from mage_tool.IO import OligoLibraryReport


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("csvfile", help="CSV project report")
    parser.add_argument("-o", help="PDF output", default=None)
    parser.add_argument("--HTML", help="Generate HTML output instead", action='store_true')
    args = parser.parse_args()

    if args.o:
        output_pdf = args.o
    else:
        output_pdf = os.path.splitext(args.csvfile)[0]

    if not args.HTML:
        output_pdf += '.pdf'

    print("Writing report to {}..".format(output_pdf))


    project = os.path.basename(os.path.splitext(args.csvfile)[0])
    report = OligoLibraryReport(project)
    report.parse_and_generate(args.csvfile, csv_file=True)
    if args.HTML:
        report.write_html(output_pdf)
    else:
        report.write_pdf(output_pdf)
