#!/usr/bin/env python

"""
General interface to <>
"""

import logging

from oligo_design import Oligo
from translation import replace_start_codon
from translation import translational_KO
from translation import RBS_single_mutation

#Define a log
log = logging.getLogger("MODEST")
log.addHandler(logging.NullHandler())


def interface(adjustments, genes, genome, config, project=None):
    """General interface to <>

    adjustments is a parsed list
    genes is a dictionary with Gene objects
    genome is a Biopython Seq object or a string
    """

    oligos = list()

    i = 0
    for a in adjustments:
        #Collect gene
        gene = genes[a["gene"]]
        op = a["operation"]
        op_str = "[{}/{}] line {}".format(op, gene, a["line"])
        #Try to find operation
        try:
            op = operations[op]
        except KeyError:
            log.error("{} not found. Doing nothing.".format(op_str))
        else:
            #Do operation
            muts = op(gene, config, a["options"], op_str)

            #Functions may return several oligos
            for mut, code, operation in muts:
                oligo = Oligo(mut, gene, project, i, oligo_len=90)
                oligo.set_oligo(genome, optimise=True, threshold=-20.0)
                oligo.target_lagging_strand(config["replication"]["ori"], config["replication"]["ter"])
                #Back tracing
                oligo.code = code
                oligo.operation = operation
                #Add to log
                log.info(" ".join([operation, str(mut), ">>", oligo.id()]))
                #Add
                oligos.append(oligo)
                i += 1

    return oligos


"""
Translation modifications
"""

def MAGE_start_codon_optimal(gene, config, options, op):
    mut = replace_start_codon(gene, config["start codons"][0])
    if not mut:
        log.debug(op + " Not mutating, optimal start codon found.")
        return []

    code = "OptStart{}".format(len(mut.before))
    return [(mut, code, op)]

def MAGE_translational_KO(gene, config, options, op):
    options, opt_str = parse_options(options, {"ko_frame": (int, 10)}, op)
    op += " " + opt_str
    if not options:
        return []

    KO_frame = options["ko_frame"]
    stop_codons = config["stop codons"]
    try:
        mut = translational_KO(gene, stop_codons, KO_frame)
        code = "TransKO{}".format(mut._codon_offset)
        return [(mut, code, op)]
    except Exception as e:
        log.error(op + " Caught exception: " + str(e))
        return []


def MAGE_RBSopt_single(gene, config, options, op):
    kwds = {"insert": (truefalse, False),
            "delete": (truefalse, False),
            "top": (int, 3),
            "maximise": (truefalse, True)
    }
    options, opt_str = parse_options(options, kwds, op)
    if not options:
        return []
    muts = RBS_single_mutation(gene, **options)
    if not muts:
        return []

    op += " " + opt_str

    muts_out = list()
    for m in muts:
        code = "RBSoptSingle{:.2f}".format(m._adjustment)
        l_op = op + " {:.4f}X wt".format(m._adjustment)
        muts_out.append((m, code, l_op))

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
        mut: [GTG=ATC]
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
    "RBSopt_single":        MAGE_RBSopt_single
    }