#!/usr/bin/env python2

"""
General interface to <>
"""
from __future__ import print_function
import signal
import logging
from multiprocessing import Pool, Value, Lock

import mage_tool.run_control as rc
from mage_tool.operations import OPERATIONS
from mage_tool.IO import ParserError

#Define a log
log = logging.getLogger("MODEST")
log.addHandler(logging.NullHandler())

"""
Constants
"""
counter = Value('i', 0)
lock = Lock()

def parse_adjustments(adjlist, genes, config, barcoding_lib):
    """Parse an adjlist to operation objects.

    adjlist is a list of dictionaries.

    """
    parsed_operations = list()
    error_list = list()

    for adj in adjlist:
        #Check existence of gene and operation
        gene_str = adj["gene"]
        op_str = adj["op"]
        i = adj["line_id"]
        old_error_len = len(error_list)

        try:
            gene = genes[gene_str]
            op = OPERATIONS[op_str]
        except KeyError:
            if op_str not in OPERATIONS:
                error_list.append("Operation '{}' not found in line {}.".format(op_str, i))
            if gene_str not in genes:
                error_list.append("Gene '{}' not found in line {}.".format(gene_str, i))
            continue

        current_operation = op(i, gene, adj["options"], config)

        #Validate existance of barcodes
        for bc in current_operation.barcodes:
            if bc not in barcoding_lib:
                error_list.append("Barcode '{}' not found in barcode lib".format(bc))

        #Operation has an error
        if not current_operation:
            for e in current_operation.errorlist:
                error_list.append("{} error: {}".format(current_operation, e))
            continue

        #No new errors
        if len(error_list) == old_error_len:
            parsed_operations.append(current_operation)

    if error_list:
        raise ParserError("Error making operations", error_list)

    return parsed_operations


def run_adjustments(oplist, genome, project, barcoding_lib, threaded=True):
    """Run a list of operation objects.

    Returns a list of Oligo objects.

    """
    if threaded:
        pool = Pool(rc.CONF["processes"], process_initializer)

    results = list()
    for op in oplist:
        #increase oligo counter
        with lock:
            counter.value += 1
        if threaded:
            kwargs = {"genome": genome, "project": project, "barcoding_lib": barcoding_lib, "count": counter.value}
            results.append(pool.apply_async(create_oligos_decorator, (op, kwargs,)))
        else:
            try:
                results.append(op.create_oligos(genome, project, barcoding_lib, counter.value))
            except KeyboardInterrupt:
                log.error("Computation manually stopped")
                break

    oligos = list()
    for r in results:
        try:
            if threaded:
                res = r.get(999999)
            else:
                res = r

            if res is None:
                log.error("Caught exception doing operation. "
                          "Exception printed to stdout.")
            else:
                oligos.extend(res)
        except KeyboardInterrupt:
            print()
            log.error("Computation manually stopped.")
            return oligos

    return oligos


def process_initializer():
    """Ignores KeyboardInterrupt. Useful for subprocesses."""
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def create_oligos_decorator(op, kwargs):
    """Catches all exceptions from create_oligo and prints to screen."""
    try:
        return op.create_oligos(**kwargs)
    except Exception:
        import traceback
        print()
        traceback.print_exc()
        return None


"""
Unittests
~~~~~~~~~
"""

import unittest
import oligo_design

class InterfaceTests(unittest.TestCase):
    def setUp(self):
        gene_cds = ("ATGTCGTGTGAAGAACTGGAAATTGTCTGGAACAATATTAAAGCCGAAGCCAGAACGCTG"
                    "GCGGACTGTGAGCCAATGCTGGCCAGTTTTTACCACGCGACGCTACTCAAGCACGAAAAC"
                    "CTTGGCAGTGCACTGAGCTACATGCTGGCGAACAAGCTGTCATCGCCAATTATGCCTGCT"
                    "ATTGCTATCCGTGAAGTGGTGGAAGAAGCCTACGCCGCTGACCCGGAAATGATCGCCTCT"
                    "GCGGCCTGTGATATTCAGGCGGTGCGTACCCGCGACCCGGCAGTCGATAAATACTCAACC"
                    "CCGTTGTTATACCTGAAGGGTTTTCATGCCTTGCAGGCCTATCGCATCGGTCACTGGTTG"
                    "TGGAATCAGGGGCGTCGCGCACTGGCAATCTTTCTGCAAAACCAGGTTTCTGTGACGTTC"
                    "CAGGTCGATATTCACCCGGCAGCAAAAATTGGTCGCGGTATCATGCTTGACCACGCGACA"
                    "GGCATCGTCGTTGGTGAAACGGCGGTGATTGAAAACGACGTATCGATTCTGCAATCTGTG"
                    "ACGCTTGGCGGTACGGGTAAATCTGGTGGTGACCGTCACCCGAAAATTCGTGAAGGTGTG"
                    "ATGATTGGCGCGGGCGCGAAAATCCTCGGCAATATTGAAGTTGGGCGCGGCGCGAAGATT"
                    "GGCGCAGGTTCCGTGGTGCTGCAACCGGTGCCGCCGCATACCACCGCCGCTGGCGTTCCG"
                    "GCTCGTATTGTCGGTAAACCAGACAGCGATAAGCCATCAATGGATATGGACCAGCATTTC"
                    "AACGGTATTAACCATACATTTGAGTATGGGGATGGGATCTAA")
        self.gene = oligo_design.Gene("cysE", 3780585, -1, gene_cds)
        self.genome = oligo_design.Gene("genome", 0, 1, "A")

    def test_bool(self):
        op = OPERATIONS["residue_mutation"]
        config = {}
        opF1 = op(0, self.gene, "", config)
        opF2 = op(1, self.gene, "mut=A", config)
        opT1 = op(2, self.gene, "mut=E166G", config)
        opT2 = op(3, self.gene, "", config, mut="E166G")
        self.assertFalse(opF1)
        self.assertFalse(opF2)
        self.assertTrue(opT1)
        self.assertTrue(opT2)


if __name__ == "__main__":
    unittest.main()
