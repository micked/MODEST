
"""
mage_tool operations

stateful module, which holds all operations
"""

from __future__ import absolute_import

import logging
import os.path
import pkgutil

from mage_tool.oligo_design import Oligo
import mage_tool.run_control as rc

#Define a log
log = logging.getLogger("MODEST")
log.addHandler(logging.NullHandler())


#Operations are gradually added to this dict
OPERATIONS = dict()


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

    #If computation is longer than a minute, set this to True
    heavy = False

    op_str = "base_operation"

    def __init__(self, line_id, gene, opt_str, config, **options):
        self.ok = True
        self.errorlist = list()
        self.opt_str = "NotEvaluated"

        self.line_id = line_id
        self.gene = gene
        self.config = config

        if opt_str:
            options.update(self.parse_options(opt_str))

        #Set custom ID
        self.custom_id = None
        if "id" in options:
            self.custom_id = options["id"]

        #Set barcodes
        self.barcodes = []
        if "barcodes" in options:
            self.barcodes = options["barcodes"]
            if isinstance(self.barcodes, basestring):
                self.barcodes = self.barcodes.split(";")

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
                    self.error("Invalid value '{}' for '{}'".format(options[o], self.default_options[o][0].__name__))
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

        if self.barcodes:
            self.opt_str += ',barcodes=' + ';'.join(self.barcodes)

    def error(self, e):
        """Add an error."""
        self.errorlist.append(e)
        self.ok = False

    def create_oligos(self, genome, project, barcoding_lib, count=None):
        """Run an operation and create oligos"""
        muts = self.run()
        gene = self.gene
        config = self.config
        barcodes = self.barcodes

        if muts is None:
            log.warn(str(self) + " did not make any mutations.")
            return []

        if count is None:
            count = self.line_id
        oligos = list()
        for i, (mut, code, operation, values) in enumerate(muts, 1):
            number = "{:0>4}.{}".format(count, i)
            oligo = Oligo(mut, gene, project, number, oligo_len=rc.CONF["oligo_length"])
            oligo.set_oligo(genome, optimise=True, threshold=-20.0)
            oligo.target_lagging_strand(config["replication"]["ori"], config["replication"]["ter"])

            #Back tracing
            oligo.code = code
            oligo.operation = operation
            oligo.operation_values = values
            #Set custom ID
            if self.custom_id:
                oligo.set_custom_id(self.custom_id)

            #Add to log
            log.info(" ".join([operation, ">>", oligo.short_id()]))
            log.info(oligo.id("full"))

            #reset barcode counter
            j = 1
            #Add barcodes
            if not barcodes:
                oligos.append(oligo)
            else:
                for barcode_ids in barcodes:
                    temp_oligo = oligo.copy()
                    temp_oligo.number = "{}.{}".format(number, j)
                    temp_oligo.add_barcodes(barcode_ids, barcoding_lib)
                    #Add
                    oligos.append(temp_oligo)
                    j += 1
        return oligos


def register_module(module):
    """Imports a module and registers all operations.

    An operation is defined as a class that inherits from BaseOperation.

    """
    new_ops = __import__(module, globals(), locals(), ["OPERATIONS"])
    OPERATIONS.update(new_ops.OPERATIONS)


#Import all internal operations
pkgpath = os.path.dirname(__file__)
for _, name, _ in pkgutil.iter_modules([pkgpath]):
    register_module("mage_tool.operations." + name)


#Import user operations
for module in rc.CONF["operations"]:
    register_module(module)
