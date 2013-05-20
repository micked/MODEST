#!/usr/bin/env python2

"""
General interface to <>
"""

import re
import signal
import logging
from multiprocessing import Pool, Value, Lock

from oligo_design import Oligo
import translation
import manual

#Define a log
log = logging.getLogger("MODEST")
log.addHandler(logging.NullHandler())

"""
Config and constants
"""
NUM_PROCESSES = None

counter = Value('i', 0)
lock = Lock()

#Operations are gradually added to this dict
OPERATIONS = dict()

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

        try:
            gene = genes[gene_str]
            op = OPERATIONS[op_str]
        except KeyError:
            if op_str not in OPERATIONS:
                error_list.append("Operation {} not found in line {}."
                                  "".format(op_str, i))
            if gene_str not in genes:
                error_list.append("Gene {} not found in line {}."
                                  "".format(gene_str, i))
            continue

        #Validate existance of barcodes
        barcodes = adj["barcodes"]
        for bc in barcodes:
            if bc not in barcoding_lib:
                error_list.append("Barcode {} not found in barcode lib."
                                  "".format(bc))

        current_operation = op(i, gene, adj["options"], config, barcodes)

        if not current_operation:
            for e in current_operation.errorlist:
                error_list.append("{} error: {}".format(op, e))
            continue

        if not error_list:
            parsed_operations.append(current_operation)

    return parsed_operations, error_list


def run_adjustments(oplist, genome, project, barcoding_lib, threaded=True):
    """Run a list of operation objects.

    Returns a list of Oligo objects.

    """
    if threaded:
        pool = Pool(NUM_PROCESSES, process_initializer)

    results = list()
    for op in oplist:
        if threaded:
            kwargs = {"genome": genome, "project": project, "barcoding_lib": barcoding_lib}
            results.append(pool.apply_async(create_oligos_decorator, (op, kwargs,)))
        else:
            try:
                results.append(op.create_oligos(genome, project, barcoding_lib))
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
            print
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
        print
        traceback.print_exc()
        return None


class BaseOperation(object):
    """
    Sets the basis for an operation. Does nothing by itself, but can be
    subclassed to create new operations which overload .run().
    """

    #Default options and type
    # default_options = {"option": (type, "default_value")}
    default_options = {}

    #Iterable of required options
    required = ()

    #Set this to false to disable operation on whole genome
    genome_allowed = True

    op_str = "base_operation"

    def __init__(self, line_id, gene, opt_str, config, barcodes=[], **options):
        self.ok = True
        self.errorlist = list()
        self.opt_str = "NotEvaluated"

        self.barcodes = barcodes
        self.line_id = line_id
        self.gene = gene
        self.config = config

        if opt_str:
            options.update(self.parse_options(opt_str))

        self.validate_options(options)

        if not self.genome_allowed and str(gene) == "genome":
            self.error("Cannot use on genome")

        self.create_opt_str()
        if self:
            self.post_init()

    """
    Overload these
    """

    def post_init(self):
        """Called after a successful init."""
        pass

    def run(self):
        return None

    """
    Status
    """

    def __nonzero__(self):
        return self.ok

    def __str__(self):
        return "[{}/{}] line {} {}".format(self.op_str, self.gene, self.line_id, self.opt_str)

    def __repr__(self):
        return self.__str__()

    """
    Option parsing
    """

    def parse_options(self, opt_str):
        """Parse options string."""
        if not opt_str:
            options = []
        else:
            options = opt_str.split(",")

        #Parse supplied keywords and values
        parsed_options = dict()
        for option in options:
            tmp = option.split("=", 1)
            if len(tmp) < 2:
                self.error("Invalid option without value <{}>".format(option))
                self.ok = False
            else:
                parsed_options[tmp[0].lower()] = tmp[1]

        return parsed_options

    def validate_options(self, options):
        """Validate a dictionary of options.

        Fill out remaining options with default options.

        """
        self.options = dict()

        #Iterate default options
        for o in self.default_options:
            #Option supplied, parse
            if o in options:
                try:
                    tp = self.default_options[o][0]
                    self.options[o] = tp(options[o])
                except ValueError:
                    self.error("Invalid value '{}' for {}".format(options[o], self.default_options[o][0]))
            elif o in self.required:
                self.error("Required option '{}' not found.".format(o))
            #Use default value
            else:
                self.options[o] = self.default_options[o][1]

    def create_opt_str(self):
        self.opt_str = ""
        for o, v in self.options.items():
            if isinstance(v, list) or isinstance(v, tuple):
                v = ";".join(v)
            self.opt_str += "," + o + "=" + str(v)

        self.opt_str = self.opt_str.strip(",")

    def error(self, e):
        """Add an error."""
        self.errorlist.append(e)
        self.ok = False

    def create_oligos(self, genome, project, barcoding_lib):
        """Run an operation and create oligos"""
        muts = self.run()
        gene = self.gene
        config = self.config
        barcodes = self.barcodes

        if muts is None:
            log.warn(str(self) + " did not make any mutations.")
            return []

        oligos = list()
        for mut, code, operation, values in muts:
            #reset barcode counter
            j = 1
            #increase oligo counter
            with lock:
                counter.value += 1
            number = "{:0>4}".format(counter.value) # + barcoding
            oligo = Oligo(mut, gene, project, number, oligo_len=90)
            oligo.set_oligo(genome, optimise=True, threshold=-20.0)
            oligo.target_lagging_strand(config["replication"]["ori"], config["replication"]["ter"])
            #Back tracing
            oligo.code = code
            oligo.operation = operation
            oligo.operation_values = values
            #Add to log
            log.info(" ".join([operation, ">>", oligo.short_id()]))
            log.info(oligo.id())
            #Add barcodes
            for barcode_ids in barcodes:
                temp_oligo = oligo.copy()
                temp_oligo.number = "{}.{}".format(number, j) # + barcoding
                temp_oligo.add_barcodes(barcode_ids, barcoding_lib)
                #Add
                oligos.append(temp_oligo)
                j += 1

        return oligos


"""
Custom types
~~~~~~~~~~~~
"""

def truefalse(t):
    """Parsed bool"""
    if t.lower() in ["true", "1", "yes"]:
        return True
    elif t.lower() in ["false", "0", "no"]:
        return False
    else:
        raise ValueError("Could not parse: {} as true/false value".format(t))


def rbs_method(s):
    """RBS library method."""
    s = s.lower()
    if s in ["exp", "fuzzy"]:
        return s
    else:
        raise ValueError("Unknown RBS_library method: {}".format(s))


def residue_mutlist(s):
    """Several residue mutations."""
    muts = s.split(";")
    for mut in muts:
        m = re.match("^([A-Z*$])(\d+)([a-z]?)([A-Z*$])$", mut)
        if not m:
            raise ValueError("Invalid mutation: {}".format(mut))
    return muts


def dna_mut(s):
    """Mutation to DNA."""
    t1 = re.match(r"^\[(\w*)=(\w*)\]\.(-?\d*)$", s)
    t2 = re.match(r"^(\w*)\[(\w*)=(\w*)\](\w*)$", s)
    if not t1 and t2:
        raise ValueError("Invalid DNA mutation: {}".format(s))
    return s


def fold_list(s):
    """A list of improvement folds."""
    s = s.split(";")
    remax = re.compile(r"^max\d*$")
    remin = re.compile(r"^min\d*$")
    for i,f in enumerate(s):
        f = f.lower()
        if remax.match(f) or remin.match(f):
            s[i] = f
        else:
            s[i] = float(f)
    return s


"""
Custom mutation operations
~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

class DNAMutation(BaseOperation):

    default_options = {"mut": (dna_mut, None)}
    required = ("mut",)
    genome_allowed = True
    op_str = "dna_mutation"

    def run(self):
        """Allows for a desired mutation.

        Options:
        - ``mut=upstream[before=after]downstream``

        I.e. ``mut=TATCAACGCC[GCTC=A]GCTTTCATGACT`` changes
        TATCAACGCC\ **GCTCG**\ CTTTCATGACT to TATCAACGCC\ **A**\ GCTTTCATGACT.

        """
        mut = self.options["mut"]
        mut = manual.gene_mutation(self.gene, mut)
        if not mut:
            log.error(str(self) + " Not mutating, did not find mutation box")
            return None

        code = "DNAMut"
        return [(mut, code, str(self), [])]

OPERATIONS[DNAMutation.op_str] = DNAMutation


class ResidueMutation(BaseOperation):

    default_options = {"mut": (residue_mutlist, None)}
    required = ("mut",)
    genome_allowed = False
    op_str = "residue_mutation"

    def run(self):
        """``residue_mutation``: Mutating a residue.

        Options:
        - ``mut=``

        I.e. ``mut=`` changes

        """
        cdn_tbl = self.config["codon_table"]
        dgn_tbl = self.config["dgn_table"]
        cdn_usage = self.config["codon_usage"]
        mut = self.options["mut"]
        mut = manual.residue_mutation(self.gene, mut, cdn_tbl, dgn_tbl, cdn_usage)
        if not mut:
            log.error(str(self) + " Not mutating.")
            return None

        code = "ResMut"
        return [(mut, code, str(self), [])]

OPERATIONS[ResidueMutation.op_str] = ResidueMutation


"""
Translation modifications
~~~~~~~~~~~~~~~~~~~~~~~~~
"""

class StartCodonOptimal(BaseOperation):

    default_options = {}
    required = ()
    genome_allowed = False
    op_str = "start_codon_optimal"

    def run(self):
        """Mutates a start codon to the optimal start codon.

        Tries to mutate the start codon to the optimal start codon defined in the
        strain config file (Usually ``ATG``\ ).

        Options and default values:
        - None

        """
        mut = translation.replace_start_codon(self.gene, self.config["start_codons"][0])
        if not mut:
            log.debug(str(self) + " Not mutating, optimal start codon found.")
            return []

        code = "OptStart{}".format(len(mut.before))
        return [(mut, code, str(self), [])]

OPERATIONS[StartCodonOptimal.op_str] = StartCodonOptimal


class TranslationalKnockout(BaseOperation):

    default_options = {"ko_frame": (int, 10),
                       "ko_mutations": (int, 3)}
    required = ()
    genome_allowed = False
    op_str = "translational_knockout"

    def run(self):
        """Gene knock-out by premature stop-codon.

        Tries to knock out a gene by introducing a number of early stop-codons in
        the CDS with the least amount of mutations possible.

        Options and default values:
        - ``ko_frame`` Number of codons that are applicable to be mutated.
            E.g. a value of 10 means the operation will try to mutate stop codons
            into the CDS within 10 codons of the start codon. Default is within one
            half of the length of the CDS.
        - ``ko_mutations`` number of stop codons to introduce. Default (and
            minimum) is the number of different stop codons available in the genome
            configuration file (normally 3).

        """
        stop_codons = self.config["stop_codons"]
        ko_mutations = self.options["ko_mutations"]
        ko_frame = self.options["ko_frame"]
        mut = translation.translational_KO(self.gene, stop_codons, ko_mutations, ko_frame)
        code = "TransKO{}".format(mut._codon_offset)
        return [(mut, code, str(self), [0, mut._codon_offset])]

OPERATIONS[TranslationalKnockout.op_str] = TranslationalKnockout


class RBSLibrary(BaseOperation):

    default_options = {"target": (float, 5000000.),
                       "n": (int, 10),
                       "max_mutations": (int, 10),
                       "method": (rbs_method, "exp"),
                       "m": (float, 0)}
    required = ()
    genome_allowed = False
    op_str = "RBS_library"

    def post_init(self):
        #Correct very low target values
        if self.options["target"] < 0.1:
            val = self.options["target"]
            self.options["target"] = 0.1
            self.create_opt_str()
            log.debug("{} adjusting very low target value '{}' to 0.1".format(self, val))

    def run(self):
        """Create a library of different RBS expression levels.

        Options and default values:
        - ``target=5000000`` Target expression level to reach for in AU. If
            target is reached, computation is stopped, and library will be created.
            if target is not reached within the specified number of mutations, a
            library of expression levels closest to target as possible will be
            created.
        - ``n=10`` Number of library sequences to create.
        - ``max_mutations=10`` Maximum number of mutations to attempt.
        - ``method=exp`` How to create the library. Two methods are available:
            ``exp`` and ``fuzzy``. ``exp`` creates a library where each new
            sequence is an ``m``-fold improvement over the last. ``m`` can either
            be supplied via the ``m``-parameter, or calculated automatically to
            create an evenly spaced library from wt level to target. The ``exp``
            method runs multiple Monte Carlo simulations to reach each target,
            however, it uses information from previous runs to more quickly reach
            subsequent targets. ``fuzzy`` tries to replicate the ``exp`` library,
            only the Monte Carlo simulation is only run once, and inbetween
            are collected along the way. This method yields a less precise library,
            but is quicker. Additionally, ``fuzzy`` enables picking out the best
            sequences below a certain mutation count by using the ``m`` parameter.
            Fx. using ``m=6``, ``fuzzy`` will collect the best possible sequences
            with a maximum of 1, 2, .. 6 mutations. It will then try to fill out
            the rest of the library with evenly spaced sequences.
        - ``m=0`` see ``method`` for explanation on ``m``.

        This operation will run an Monte-Carlo simulation in an attempt to reach
        the specified target value within a number of mutations. Lower numbers of
        mutations are tried first and are always prioritised over similar
        expression levels with a higher mutation count.

        """
        method = self.options["method"]
        target = self.options["target"]
        n = self.options["n"]
        m = self.options["m"]
        max_mutations = self.options["max_mutations"]
        if method == "exp":
            muts = translation.RBS_library(self.gene, target, n, max_mutations, m)
        elif method == "fuzzy":
            muts = translation.RBS_library_fuzzy(self.gene, target, n, max_mutations, m)

        if not muts:
            return None

        muts_out = list()
        for i, m in enumerate(muts):
            code = "RBSlib{}_{:.1f}/{:.1f}({})".format(i, m._AU, m._orgAU, m._n)
            l_op = str(self) + " {:.3f} (wt: {:.2f})".format(m._AU, m._orgAU)
            muts_out.append((m, code, l_op, [m._orgAU, m._AU]))

        return muts_out

OPERATIONS[RBSLibrary.op_str] = RBSLibrary


"""
Transcription modifications
~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

class PromoterLibrary(BaseOperation):

    default_options = {"targets": (fold_list, ["max"]),
                       "max_mutations": (int, 10),
                       "matrix": (str, "sigma70")}
    required = ()
    genome_allowed = False
    op_str = "promoter_library"

    def run(self):
        """Create a library of different promoter expression levels.

        """
        targets = self.options["targets"]
        max_mutations = self.options["max_mutations"]
        matrix = self.options["matrix"]
        muts = promoter.promoter_library(self.gene, targets, max_mutations, matrix)

        if not muts:
            return None

        muts_out = list()
        for i, m in enumerate(muts):
            code = "promoterlib{}_{:.1f}({})".format(i, m._fold, m._n)
            l_op = "{} {:.3f})".format(self, m._fold)
            muts_out.append((m, code, l_op, [m._fold]))

        return muts_out

OPERATIONS[PromoterLibrary.op_str] = PromoterLibrary


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
