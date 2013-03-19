#!/usr/bin/env python

"""
General interface to <>
"""

import logging
from multiprocessing import Pool

from oligo_design import Oligo
from translation import replace_start_codon
from translation import translational_KO
from translation import RBS_single_mutation
from translation import generate_RBS_library

#Define a log
log = logging.getLogger("MODEST")
log.addHandler(logging.NullHandler())


def create_oligos(genome, op, gene, config, options, op_str, project, barcode_id, barcoding_lib, i):
    """Run an operation and create oligos"""
    muts = op(gene, config, options, op_str)

    oligos = list()

    if muts is None:
        log.warn(op_str + " did not make any mutations.")

    j = 0
    for mut, code, operation, values in muts:
        number = "{:0>4}.{}".format(i, j) # + barcoding
        oligo = Oligo(mut, gene, project, number, oligo_len=90)
        oligo.set_oligo(genome, optimise=True, threshold=-20.0)
        oligo.target_lagging_strand(config["replication"]["ori"], config["replication"]["ter"])
        #Add barcodes
        oligo.add_barcode(barcode_id, barcoding_lib)
        #Back tracing
        oligo.code = code
        oligo.operation = operation
        oligo.operation_values = values
        #Add to log
        log.info(" ".join([operation, str(mut), ">>", oligo.id()]))
        #Add
        oligos.append(oligo)
        j += 1

    return oligos


def interface(adjustments, genes, genome, config, barcoding_lib, project=None):
    """General interface to <>

    adjustments is a parsed list
    genes is a dictionary with Gene objects
    genome is a Biopython Seq object or a string

    """
    p_pool = Pool(processes=4)
    oligos_list = list()

    i = 0
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
            barcode_id = a["barcode_id"]

            #Do operation
            args = (genome, op, gene, config, a["options"], op_str, project, barcode_id, barcoding_lib, i)
            oligos_list.append(p_pool.apply_async(create_oligos, args))
            i += 1

    #This does not work:
    #DO SOMETHINGNGNGNGN!
    try:
        p_pool.close()
    except KeyboardInterrupt:
        log.warn("Computation manually stopped.")

    oligos = list()
    for oligos_out in oligos_list:
        oligos.extend(oligos_out.get())

    return oligos


"""
Translation modifications
"""

def MAGE_start_codon_optimal(gene, config, options, op):
    mut = replace_start_codon(gene, config["start_codons"][0])
    if not mut:
        log.debug(op + " Not mutating, optimal start codon found.")
        return []

    code = "OptStart{}".format(len(mut.before))
    return [(mut, code, op, [])]

def MAGE_translational_KO(gene, config, options, op):
    options, opt_str = parse_options(options, {"ko_frame": (int, 10)}, op)
    op += " " + opt_str
    if not options:
        return None

    KO_frame = options["ko_frame"]
    stop_codons = config["stop_codons"]
    try:
        mut = translational_KO(gene, stop_codons, KO_frame)
        code = "TransKO{}".format(mut._codon_offset)
        return [(mut, code, op, [mut._codon_offset])]
    except Exception as e:
        log.error(op + " Caught exception: " + str(e))
        return None


def MAGE_RBSopt_single(gene, config, options, op):
    kwds = {"insert": (truefalse, False),
            "delete": (truefalse, False),
            "top": (int, 3),
            "maximise": (truefalse, True)
    }
    options, opt_str = parse_options(options, kwds, op)
    if not options:
        return None
    muts = RBS_single_mutation(gene, **options)
    if not muts:
        return None

    op += " " + opt_str

    muts_out = list()
    for m in muts:
        code = "RBSoptSingle{:.2f}".format(m._adjustment)
        l_op = op + " {:.4f}X wt".format(m._adjustment)
        muts_out.append((m, code, l_op, [m._adjustment]))

    return muts_out


def MAGE_RBS_library(gene, config, options, op):
    kwds = {"target": (float, 5000000.),
            "n": (int, 10),
            "max_mutations": (int, 10),
            "passes": (int, 1),
            "id": (str, "-")
    }
    options, opt_str = parse_options(options, kwds, op)
    if not options:
        return None

    muts = generate_RBS_library(gene, target=options["target"],
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
        tmp = option.split("=")
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
    "start_codon_optimal":  MAGE_start_codon_optimal,
    "translational_KO":     MAGE_translational_KO,
    "RBSopt_single":        MAGE_RBSopt_single,
    "RBS_library":          MAGE_RBS_library,
    }
