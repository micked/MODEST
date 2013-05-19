#!/usr/bin/env python

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
Config
"""
NUM_PROCESSES = None

counter = Value('i', 0)
lock = Lock()


def parse_adjustments(adjfilehandle, genes, barcode_lib=None):
    """Pass"""
    error_list = list()
    parsed_operations = list()
    #Iterate lines:
    for i, line in enumerate(adjfilehandle, 1):
        if not line.strip() or line[0] == "#":
            continue
        line = line.split()
        if len(line) < 3:
            #Append error and quit line
            error_list.append("Too few arguments in line {} in adjustmentlist."
                              "".format(i))
            continue
        elif len(line) == 3:
            options = ""
        else:
            options = line[3]

        gene_str, op_str, barcodes = line[0:3]
        barcodes = barcodes.split(",")

        #Check existence of gene and operation
        try:
            gene = genes[gene_str]
            op, op_kwargs = OPERATIONS[op_str]
        except KeyError:
            if op_str not in OPERATIONS:
                error_list.append("Operation {} not found in line {}."
                                  "".format(op_str, i))
            if gene_str not in genes:
                error_list.append("Gene {} not found in line {}."
                                  "".format(gene_str, i))
            continue

        #Validate existance of barcodes
        for bc in barcodes:
            if barcode_lib and bc not in barcode_lib:
                error_list.append("Barcode {} not found in barcode lib."
                                  "".format(bc))

        current_operation = {"op": op, "gene": gene, "barcodes": barcodes}
        options, options_str = parse_options(options, op_kwargs)

        if options is None:
            error_list.append("Options parser error line {}: {}"
                              "".format(i, options_str))
            continue

        op_str = "[{}/{}] line {} {}".format(op_str, gene, i, options_str)

        current_operation.update({"options": options, "op_str": op_str})
        parsed_operations.append(current_operation)

    return parsed_operations, error_list


def run_adjustments(adjfilehandle, genes, genome, config, project,
                    barcoding_lib, threaded=True):
    """Pass"""
    parser_args = (adjfilehandle, genes, barcoding_lib)
    parsed_operations, errors = parse_adjustments(*parser_args)
    if errors:
        return [], errors

    oligo_kwargs = {"genome": genome, "config": config, "project": project,
                    "barcoding_lib": barcoding_lib}

    if threaded:
        pool = Pool(NUM_PROCESSES, process_initializer)

    results = list()
    for run_op in parsed_operations:
        op_kwargs = oligo_kwargs.copy()
        op_kwargs.update(run_op)
        if threaded:
            results.append(pool.apply_async(create_oligos_decorator, (op_kwargs,)))
        else:
            try:
                results.append(create_oligos(**op_kwargs))
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
            return oligos, []

    return oligos, []


def process_initializer():
    """Ignores all exceptions. Useful for subprocesses"""
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def create_oligos_decorator(kwargs):
    """TODO"""
    try:
        return create_oligos(**kwargs)
    except Exception:
        import traceback
        print
        traceback.print_exc()
        return None


def create_oligos(genome, op, gene, config, options, op_str, project, barcodes, barcoding_lib):
    """Run an operation and create oligos"""
    oligos = list()
    muts = op(gene, op_str, config, **options)

    if muts is None:
        log.warn(op_str + " did not make any mutations.")
        return []

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
==================
Mutation functions
==================

Each function described here must take a gene object, an op string (for logging
purposes), a strain config dict, and additional kwargs described in global
OPERATIONS dict.

Each function must output a list of tuples, each tuple being a mutation in the
form: (mut, code, op, valuelist)

Where:
    mut       = Mutation object
    code      = Compressed code string for oligo ID (length ~ 10 chars)
    op        = String with operation and output status
    valuelist = Status/output values in form [wt, altered, *other]

List of status/output values can be empty, but wt value must always be at pos
0, and altered at 1. If wt value is not applicable, return 0 or "" as wt value.

If an error occured within the function, the function should return None and
log a message stating the error.


Custom Mutations
~~~~~~~~~~~~~~~~
"""

def dna_mutation(gene, op, config, mut):
    """Allows for a desired mutation.

    Options:
      - ``mut=upstream[before=after]downstream``

    I.e. ``mut=TATCAACGCC[GCTC=A]GCTTTCATGACT`` changes
    TATCAACGCC\ **GCTCG**\ CTTTCATGACT to TATCAACGCC\ **A**\ GCTTTCATGACT.

    """
    if not mut:
        log.error(op + " No mutation found.")
        return None

    mut = manual.gene_mutation(gene, mut)
    if not mut:
        log.error(op + " Not mutating, did not find mutation box")
        return None

    code = "DNAMut"
    return [(mut, code, op, [])]


def residue_mutation(gene, op, config, mut):
    """Allows for a mutating a desired residue

    Options:
      - ``mut=``

    I.e. ``mut=`` changes

    """
    if str(gene) == "genome":
        log.error(op + " Cannot use residue_mutation on genome")
        return None

    if not mut:
        log.error(op + " No mutation specified!")
        return None

    mut = manual.residue_mutation(gene, mut, config["codon_table"],
                                  config["dgn_table"], config["codon_usage"])
    if not mut:
        log.error(op + " Not mutating, did not find mutation box")
        return None

    code = "ResMut"
    return [(mut, code, op, [])]


"""
Translation modifications
~~~~~~~~~~~~~~~~~~~~~~~~~
"""

def start_codon_optimal(gene, op, config):
    """Mutates a start codon to the optimal start codon.

    Tries to mutate the start codon to the optimal start codon defined in the
    strain config file (Usually ``ATG``\ ).

    Options and default values:
      - None

    """
    if str(gene) == "genome":
        log.error(op + " Cannot use start_codon_optimal on genome")
        return None

    mut = translation.replace_start_codon(gene, config["start_codons"][0])
    if not mut:
        log.debug(op + " Not mutating, optimal start codon found.")
        return []

    code = "OptStart{}".format(len(mut.before))
    return [(mut, code, op, [])]


def translational_knockout(gene, op, config, ko_mutations, ko_frame):
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
    if str(gene) == "genome":
        log.error(op + " Cannot use translational_knockout on genome")
        return None

    stop_codons = config["stop_codons"]
    mut = translation.translational_KO(gene, stop_codons, ko_mutations, ko_frame)
    code = "TransKO{}".format(mut._codon_offset)
    return [(mut, code, op, [0, mut._codon_offset])]


def RBS_library(gene, op, config, target, n, max_mutations, method, m):
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
    if str(gene) == "genome":
        log.error(op + " Cannot use RBS_library on genome")
        return None

    if method == "exp":
        muts = translation.RBS_library(gene, target, n, max_mutations, m)
    elif method == "fuzzy":
        muts = translation.RBS_library_fuzzy(gene, target, n, max_mutations, m)

    if not muts:
        return None

    muts_out = list()
    for i, m in enumerate(muts):
        code = "RBSlib{}_{:.1f}/{:.1f}({})".format(i, m._AU, m._orgAU, m._n)
        l_op = op + " {:.3f} (wt: {:.2f})".format(m._AU, m._orgAU)
        muts_out.append((m, code, l_op, [m._orgAU, m._AU]))

    return muts_out


"""
Transcription modifications
~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

def promoter_library(gene, op, config, targets, max_mutations, matrix):
    """Create a library of different promoter expression levels.

    """
    if str(gene) == "genome":
        log.error(op + " Cannot use promoter_library on genome")
        return None

    muts = promoter.promoter_library(gene, targets, max_mutations, matrix)


    if not muts:
        return None

    muts_out = list()
    for i, m in enumerate(muts):
        code = "promoterlib{}_{:.1f}({})".format(i, m._fold, m._n)
        l_op = op + " {:.3f})".format(m._fold)
        muts_out.append((m, code, l_op, [m._fold]))

    return muts_out


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
            s[ı] = float(f)
    return s


"""
Dict to map all allowed operations and associated options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

OPERATIONS = {
    "start_codon_optimal":    (start_codon_optimal, {}),
    "translational_knockout": (translational_knockout, {"ko_frame": (int, 10),
                                                        "ko_mutations": (int, 3)}),
    "RBS_library":            (RBS_library, {"target": (float, 5000000.),
                                             "n": (int, 10),
                                             "max_mutations": (int, 10),
                                             "method": (rbs_method, "exp"),
                                             "m": (float, 0)}),
    "dna_mutation":           (dna_mutation, {"mut": (dna_mut, None)}),
    "residue_mutation":       (residue_mutation, {"mut": (residue_mutlist, None)}),
    "promoter_library":       (promoter_library, {"targets": (fold_list, ["max"]),
                                                  "max_mutations": (int, 10),
                                                  "matrix": (str, "sigma70")})
    }

"""
General options parser
"""

def parse_options(options, kwds):
    """Parse options string according to keywords and types, return option dict

    kwds is suppled as a dict of "keyword: [type, default_value]"

    options string is comma_seperated, no spaces: opt=4,opt2=something
    Possible values:
        int: 3
        float: 3.14
        string: something
        truefalse: True|yes|1 or False|no|0
        int_list: 1;2;3
        float_list: 3.14;5.13;5.0
        string_list: some;thing
        codon_list: ATG;ATC
        mut: [GTG=ATC].10
        ...?
    """
    if not options:
        options = []
    else:
        options = options.split(",")

    #Parse supplied keywords and values
    supp_opt = dict()
    for option in options:
        tmp = option.split("=", 1)
        if len(tmp) < 2:
            return None, "Invalid option without value <{}>".format(option)
        supp_opt[tmp[0].lower()] = "".join(tmp[1:])

    opts_out = dict()

    #Iterate requested keywords
    for o in kwds:
        #Option supplied, parse
        if o.lower() in supp_opt:
            try:
                opts_out[o] = kwds[o][0](supp_opt[o])
            except ValueError:
                return None, "Invalid value \"{}\" for {}".format(supp_opt[o], kwds[o][0])
        #Use default value
        else:
            opts_out[o] = kwds[o][1]

    opt_str = ""
    for o, v in opts_out.items():
        if isinstance(v, list) or isinstance(v, tuple):
            v = ";".join(v)
        opt_str += "," + o + "=" + str(v)

    opt_str = opt_str.strip(",")
    #opt_str = ",".join([o + "=" + str(v) for o,v in opts_out.items()])

    return opts_out, opt_str
