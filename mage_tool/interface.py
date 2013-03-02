#!/usr/bin/env python

"""
General interface to <>
"""

import logging

from oligo_design import Oligo
from translation import replace_start_codon
from translation import translational_KO

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
        #Try to find operation
        try:
            op = operations[a["operation"]]
        except KeyError:
            log.error("Operation [{}] not found. Doing nothing.".format(a["operation"]))

        #Collect gene and do operation
        gene = genes[a["gene"]]
        muts = op(gene, config, a["options"], a["line"])
        
        #Functions may return several oligos
        for mut, code, operation in muts:
            oligo = Oligo(mut, gene, project, i, oligo_len=90)
            oligo.make_oligo(genome)
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

def MAGE_start_codon_optimal(gene, config, options, line):
    op = "[start_codon_optimal/{}] line {}".format(gene, line)

    mut = replace_start_codon(gene, "ATG")
    if not mut:
        log.debug(op + " Not mutating, optimal start codon found.")
        return []

    code = "OptStart{}".format(len(mut.before))
    return [(mut, code, op)]


def MAGE_translational_KO(gene, config, options, line):
    op = "[translational_KO/{}] line {}".format(gene, line)

    options, opt_str = parse_options(options, {"ko_frame": (int, 10)}, op)
    op += " " + opt_str
    if not options:
        return []

    KO_frame = options["ko_frame"]
    stop_codons = config["stop codons"]
    try:
        mut = translational_KO(gene, stop_codons, KO_frame)
        code = "TransKO{}".format(len(mut.before))
        return [(mut, code, op)]
    except Exception as e:
        log.error(op + " Caught exception: " + str(e))

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
            log.error("{}invalid option without value <{}>".format(func, option))
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
                log.error("{}invalid value {} for {}".format(func, supp_opt[o], kwds[o][0]))
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


"""
Dict to map all allowed operations
"""

operations = {
    "start_codon_optimal":  MAGE_start_codon_optimal,
    "translational_KO":     MAGE_translational_KO,
    }