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

NUM_PROCESSES = None #Use max available processes.

counter = Value('i', 0)
lock = Lock()
    

def initializer(*args):
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def create_oligos(genome, op, gene, config, options, op_str, project, multi_oligo_barcodes, barcoding_lib):
    """Run an operation and create oligos"""
    oligos = list()

    muts = op(gene, config, options, op_str)

    if muts is None:
        log.warn(op_str + " did not make any mutations.")

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
        for barcode_ids in multi_oligo_barcodes:
            temp_oligo = oligo.copy()
            temp_oligo.number = "{}.{}".format(number, j) # + barcoding
            temp_oligo.add_barcodes(barcode_ids, barcoding_lib)
            #Add
            oligos.append(temp_oligo)
            j += 1

    return oligos


def interface(adjustments, genes, genome, config, barcoding_lib, project=None):
    """General interface to <>

    adjustments is a parsed list
    genes is a dictionary with Gene objects
    genome is a Biopython Seq object or a string

    """
    p_pool = Pool(NUM_PROCESSES, initializer, (counter, lock))
    oligos_list = list()

    for a in adjustments:
        op = a["operation"]
        gene = a["gene"]
        op_str = "[{}/{}] line {}".format(op, gene, a["line"])
        
        try:
            #Collect gene
            gene = genes[gene]
            #Collect operation
            op = operations[op]
        except KeyError:
            if op not in operations:
                log.error("Operation {} not found. Doing nothing.".format(op_str))
            if gene not in genes:
                log.error("Gene {} not found. Doing nothing.".format(op_str))
        else:
            #Barcode id
            multi_oligo_barcodes = a["barcode_id"].split(',')

            #Do operation
            args = (genome, op, gene, config, a["options"], op_str, project, multi_oligo_barcodes, barcoding_lib)
            olis = p_pool.apply_async(create_oligos, args)
            oligos_list.append(olis)

    oligos = list()
    for oligos_out in oligos_list:
        try:
            oligos.extend(oligos_out.get(999999))
        except KeyboardInterrupt:
            #Kills program .. No output.
            print
            log.error("Computation manually stopped.")
            return oligos

    return oligos


"""
Custom Mutations
"""

def custom_mutation(gene, config, options, op):
    """Allows for a desired mutation.

    Options:
      - ``mut=upstream[before=after]downstream``

    I.e. ``mut=TATCAACGCC[GCTC=A]GCTTTCATGACT`` changes
    TATCAACGCC\ **GCTCG**\ CTTTCATGACT to TATCAACGCC\ **A**\ GCTTTCATGACT.

    """

    options, opt_str = parse_options(options, {"mut": (str, "[=]")}, op)
    op += " " + opt_str
    mut = custom.custom_mutation(gene, options["mut"])
    if not mut:
        log.debug(op + " Not mutating, did not find mutation box")
        return []
        
    code = "CustomMut"
    return [(mut, code, op, [])]


"""
Translation modifications
"""

def start_codon_optimal(gene, config, options, op):
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


def translational_KO(gene, config, options, op):
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
    options, opt_str = parse_options(options, {"ko_frame": (int, 10)}, op)
    op += " " + opt_str
    if not options:
        return None

    KO_frame = options["ko_frame"]
    stop_codons = config["stop_codons"]
    try:
        mut = translation.translational_KO(gene, stop_codons, KO_frame)
        code = "TransKO{}".format(mut._codon_offset)
        return [(mut, code, op, [mut._codon_offset])]
    except Exception as e:
        log.error(op + " Caught exception: " + str(e))
        return None


def RBSopt_single(gene, config, options, op):
    """Depreceated"""
    kwds = {"insert": (truefalse, False),
            "delete": (truefalse, False),
            "top": (int, 3),
            "maximise": (truefalse, True)
    }
    options, opt_str = parse_options(options, kwds, op)
    if not options:
        return None
    muts = translation.RBS_single_mutation(gene, **options)
    if not muts:
        return None

    op += " " + opt_str

    muts_out = list()
    for m in muts:
        code = "RBSoptSingle{:.2f}".format(m._adjustment)
        l_op = op + " {:.4f}X wt".format(m._adjustment)
        muts_out.append((m, code, l_op, [m._adjustment]))

    return muts_out


def RBS_library(gene, config, options, op):
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
    kwds = {"target": (float, 5000000.),
            "n": (int, 10),
            "max_mutations": (int, 10),
            "passes": (int, 1),
            "id": (str, "-")
    }
    options, opt_str = parse_options(options, kwds, op)
    if not options:
        return None

    muts = translation.generate_RBS_library(gene, target=options["target"],
                                            n=options["n"],
                                            max_mutations=options["max_mutations"],
                                            passes=options["passes"])

    op += " " + opt_str
    ID = options["id"]
    if ID == "-":
        ID = ""

    if not muts:
        return None

    muts_out = list()
    for i, m in enumerate(muts):
        code = "RBSlib{}{}_{:.1f}/{:.1f}({})".format(ID, i, m._expr, m._org_expr, m._n_muts)
        l_op = op + " {:.3f} (wt: {:.2f})".format(m._expr, m._org_expr)
        muts_out.append((m, code, l_op, [m._expr, m._org_expr]))

    return muts_out


"""
General options parser
"""

def parse_options(options, kwds, func):
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
    if not kwds:
        return None

    if not options:
        options = []
    else:
        options = options.split(",")

    #Parse supplied keywords and values
    supp_opt = dict()
    for option in options:
        tmp = option.split("=", 1)
        if len(tmp) < 2:
            log.error("{} invalid option without value <{}>".format(func, option))
            return None, ""
        supp_opt[tmp[0].lower()] = "".join(tmp[1:])

    opts_out = dict()

    #Iterate requested keywords
    for o in kwds:
        #Option supplied, parse
        if o.lower() in supp_opt:
            try:
                opts_out[o] = kwds[o][0](supp_opt[o])
            except ValueError:
                log.error("{} invalid value \"{}\" for {}".format(func, supp_opt[o], kwds[o][0]))
                return None, ""
        #Use default value
        else:
            opts_out[o] = kwds[o][1]

    opt_str = ",".join([o + "=" + str(v) for o,v in opts_out.items()])

    return opts_out, opt_str


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
Dict to map all allowed operations
"""

operations = {
    "start_codon_optimal":  start_codon_optimal,
    "translational_KO":     translational_KO,
    "RBSopt_single":        RBSopt_single,
    "RBS_library":          RBS_library,
    "custom_mutation":      custom_mutation,
    }