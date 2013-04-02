#!/usr/bin/env python

"""
General interface to <>
"""

import logging
from multiprocessing import Pool, Value, Lock
import signal

from oligo_design import Oligo
import translation
import custom

#Define a log
log = logging.getLogger("MODEST")
log.addHandler(logging.NullHandler())


"""
Config
"""
NUM_PROCESSES = None

counter = Value('i', 0)
lock = Lock()


def parse_adjustments(adjfilehandle, genes):
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
            error_list.append("Too few arguments in line {} in adjustmentlist.".format(i))
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
                error_list.append("Operation {} not found in line {}.".format(op_str, i))
            if gene_str not in genes:
                error_list.append("Gene {} not found in line {}.".format(gene, i))
            continue

        current_operation = {"op": op, "gene": gene, "barcodes": barcodes}
        options, options_str = parse_options(options, op_kwargs)

        if options is None:
            error_list.append("Options parser error line {}: {}".format(i, options_str))
            continue

        op_str = "[{}/{}] line {} {}".format(op_str, gene, i, options_str)

        current_operation.update({"options": options, "op_str": op_str})
        parsed_operations.append(current_operation)

    return parsed_operations, error_list


def run_adjustments(adjfilehandle, genes, genome, config, project, barcoding_lib):
    """Pass"""
    parsed_operations, errors = parse_adjustments(adjfilehandle, genes)
    if errors:
        return [], errors

    oligo_kwargs = {"genome": genome, "config": config, "project": project, "barcoding_lib": barcoding_lib}

    pool = Pool(NUM_PROCESSES)
    results = list()

    for run_op in parsed_operations:
        op_kwargs = oligo_kwargs.copy()
        op_kwargs.update(run_op)
        results.append(pool.apply_async(create_oligos_decorator, (op_kwargs,)))

    oligos = list()
    for r in results:
        try:
            res = r.get(999999)
            if res is None:
                log.error("Caught exception doing operation. Exception printed to stdout.")
            else:
                oligos.extend(res)
        except KeyboardInterrupt:
            print
            log.error("Computation manually stopped.")
            return oligos, []
        except Exception:
            raise 

    return oligos, []


def run_adjustments_unthreaded(adjfilehandle, genes, genome, config, project, barcoding_lib):
    """Pass"""
    parsed_operations, errors = parse_adjustments(adjfilehandle, genes)
    if errors:
        return [], errors

    oligo_kwargs = {"genome": genome, "config": config, "project": project, "barcoding_lib": barcoding_lib}

    results = list()

    for run_op in parsed_operations:
        op_kwargs = oligo_kwargs.copy()
        op_kwargs.update(run_op)
        results.append(create_oligos(**op_kwargs))

    oligos = list()
    for r in results:
        oligos.extend(r)

    return oligos, []


def create_oligos_decorator(kwargs):
    """TODO"""
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    try:
        return create_oligos(**kwargs)
    except KeyboardInterrupt:
        raise
    except Exception as exception:
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
        log.info(" ".join([operation, str(mut), ">>", oligo.id()]))
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

Each function must output a list of tuples, each tuple buing a mutation and
containing (Mutation object, small code string, operation string + status, 
            list of status/output values: [wt value, altered value, *other])

List of status/output values can be empty, but wt value must always be at 0, and
altered at 1. If wt valueis not applicable, return 0 or "" as wt value.


Custom Mutations
~~~~~~~~~~~~~~~~
"""

def custom_mutation(gene, op, config, mut):
    """Allows for a desired mutation.

    Options:
      - ``mut=upstream[before=after]downstream``

    I.e. ``mut=TATCAACGCC[GCTC=A]GCTTTCATGACT`` changes
    TATCAACGCC\ **GCTCG**\ CTTTCATGACT to TATCAACGCC\ **A**\ GCTTTCATGACT.

    """
    mut = custom.custom_mutation(gene, mut)
    if not mut:
        log.debug(op + " Not mutating, did not find mutation box")
        return []
        
    code = "CustomMut"
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
    mut = translation.replace_start_codon(gene, config["start_codons"][0])
    if not mut:
        log.debug(op + " Not mutating, optimal start codon found.")
        return []

    code = "OptStart{}".format(len(mut.before))
    return [(mut, code, op, [])]


def translational_KO(gene, op, config, ko_frame):
    """Gene knock-out by premature stop-codon.

    Tries to knock out a gene by introducing an early stop-codon in the CDS with
    the least amount of mutations possible.

    Options and default values:
      - ``ko_frame=10`` Number of codons that are applicable to be mutated.
        E.g. a value of 10 means the operation will try to mutate a stop codon
        into the CDS within 10 codons of the start codon.

    Mutations are prioritised as:
        - Single mutation, NNN -> XNN
        - Double concurrent mutation, NNN -> XXN
        - Double gapped mutation, NNN -> XNX
        - Triple mutation, NNN -> XXX

    """
    stop_codons = config["stop_codons"]
    mut = translation.translational_KO(gene, stop_codons, ko_frame)
    code = "TransKO{}".format(mut._codon_offset)
    return [(mut, code, op, [0, mut._codon_offset])]


def RBSopt_single(gene, op, config, insert, delete, top, maximise):
    """Depreceated"""
    muts = translation.RBS_single_mutation(gene, insert, delete, top, maximise)
    if not muts:
        return None

    muts_out = list()
    for m in muts:
        code = "RBSoptSingle{:.2f}".format(m._adjustment)
        l_op = op + " {:.4f}X wt".format(m._adjustment)
        muts_out.append((m, code, l_op, [m._adjustment]))

    return muts_out


def RBS_library(gene, op, config, target, n, max_mutations, passes):
    """Create a library of different RBS expression levels.

    Options and default values:
      - ``target=5000000`` Target expression level to reach for, in AU. If
        target is reached, computation is stopped, and library will be created.
        if target is not reached within the specified number of mutations, a
        library of expression levels closest to target as possible will be
        created.
      - ``n=10`` Number of mutations in library to create.
      - ``max_mutations=10`` Maximum number of mutations to attempt.
      - ``passes=1`` Number of computations. Specifying more passes will result
        in the computation restarting. Doing more passes creates more diversity
        and will result in a library with a higher resolution.
      - ``id=-`` Specify an ID for the library, which will be included in the
        code/ID for the resulting oligo. Usefull for tracing important libraries.
        Leave at ``-`` to have an empty id.

    This operation will run an Monte-Carlo simulation in an attempt to reach the
    specified target value within a number of mutations. Lower numbers of
    mutations are tried first and are always prioritised over similar expression
    levels with a higher mutation count.

    If target is quickly reached (usually happens when lowering expression
    levels), then RBS_library may not result in a full library. This can helped
    by doing more passes.

    """
    muts = translation.generate_RBS_library(gene, target=target, n=n,
                                            max_mutations=max_mutations,
                                            passes=passes)

    if not muts:
        return None

    muts_out = list()
    for i, m in enumerate(muts):
        code = "RBSlib{}_{:.1f}/{:.1f}({})".format(i, m._expr, m._org_expr, m._n_muts)
        l_op = op + " {:.3f} (wt: {:.2f})".format(m._expr, m._org_expr)
        muts_out.append((m, code, l_op, [m._org_expr, m._expr]))

    return muts_out


"""
Custom types
"""

def int_list(lst):
    raise Exception("TODO")

def truefalse(t):
    """Parsed bool"""
    if t.lower() in ["true", "1", "yes"]:
        return True
    elif t.lower() in ["false", "0", "no"]:
        return False
    else:
        raise ValueError("Could not parse: {} as true/false value".format(t))


"""
Dict to map all allowed operations and associated options.
"""

OPERATIONS = {
    "start_codon_optimal":  (start_codon_optimal, {}),
    "translational_KO":     (translational_KO, {"ko_frame": (int, 10)}),
    "RBSopt_single":        (RBSopt_single, {"insert": (truefalse, False),
                                             "delete": (truefalse, False),
                                             "top": (int, 3),
                                             "maximise": (truefalse, True)}),
    "RBS_library":          (RBS_library, {"target": (float, 5000000.),
                                           "n": (int, 10),
                                           "max_mutations": (int, 10),
                                           "passes": (int, 1)}),
    "custom_mutation":      (custom_mutation, {"mut": (str, "[=]")})
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

    opt_str = ",".join([o + "=" + str(v) for o,v in opts_out.items()])

    return opts_out, opt_str