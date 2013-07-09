
"""
Classes dealing with input/output and format changes
"""

from __future__ import absolute_import

import re
import os
import csv
import codecs
import logging

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

from mage_tool.oligo_design import Gene
from mage_tool.oligo_design import Sequence
from mage_tool.oligo_design import WobbleError
from mage_tool.helpers import is_inside
from mage_tool.helpers import extract_circular
from mage_tool.helpers import cds_to_wobble
from mage_tool.helpers import seqs_to_degenerate
from mage_tool.helpers import reverse_complement
from mage_tool.helpers import reverse_complement_dgn
import mage_tool.run_control as rc

log = logging.getLogger("MODEST")
log.addHandler(logging.NullHandler())


class ParserError(Exception):
    def __init__(self, message, error_list=[]):
        Exception.__init__(self, message)
        self.error_list = error_list


"""
Adjustment file parsers
~~~~~~~~~~~~~~~~~~~~~~~
"""

def raw_adjlist_to_adjlist(adjfilehandle):
    """Convert an adjustmentlist in raw text to adjustment list dict."""
    error_list = list()
    adjlist = list()

    for i, line in enumerate(adjfilehandle, 1):
        if not line.strip() or line[0] == "#":
            continue

        line = line.split()
        if len(line) < 3:
            #Append error and quit line
            error_list.append("Too few arguments in line {}.".format(i))
            continue

        elif len(line) == 3:
            options = ""
        else:
            options = line[3]

        gene, op, barcodes = line[0:3]
        barcodes = barcodes.split(",")
        adjlist.append({"options": options, "gene": gene, "op": op, "barcodes": barcodes, "line_id": i})

    #Raise ParserError with error_list
    if error_list:
        raise ParserError("Error parsing list of adjustments.", error_list)

    return adjlist



"""
Genelist parsers
~~~~~~~~~~~~~~~~
"""


#generic_cfg = {'Definition': 'Generic',
    #'start_codons': ['ATG', 'GTG', 'TTG', 'ATT', 'CTG'],
    #'replication': {'ter extended': [598894, 2375400],
                    #'ter': [1339769, 1682272], 'ori': [3923767, 3923998]},
    #'stop_codons': ['TAG', 'TAA', 'TGA'],
    #'dgn_table': {'A': 'GCN', 'C': 'TGY', 'E': 'GAR', '$': 'TRR', 'G': 'GGN',
                  #'F': 'TTY', 'I': 'ATH', 'H': 'CAY', 'K': 'AAR', 'M': 'ATG',
                  #'L': 'YTN', 'N': 'AAY', 'Q': 'CAR', 'P': 'CCN', 'S': 'WSN',
                  #'R': 'MGN', 'T': 'ACN', 'W': 'TGG', 'V': 'GTN', 'Y': 'TAY',
                  #'D': 'GAY'},
    #'codon_table': {'CTT': 'L', 'ATG': 'M', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I',
                    #'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'CTC': 'L',
                    #'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N',
                    #'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'ACT': 'T', 'CAC': 'H',
                    #'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q',
                    #'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R',
                    #'CAG': 'Q', 'CGC': 'R', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F',
                    #'TGC': 'C', 'GGG': 'G', 'TAG': '$', 'GGA': 'G', 'TAA': '$',
                    #'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S', 'TTA': 'L',
                    #'GAC': 'D', 'CGT': 'R', 'GAA': 'E', 'TGG': 'W', 'GCA': 'A',
                    #'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V',
                    #'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TGA': '$',
                    #'TTG': 'L', 'TCC': 'S', 'TCA': 'S', 'TCT': 'S'}}


def seqIO_to_genelist(genome, config, include_genes=None, include_genome=False,
                      exclude_genes=None, leader_len=rc.CONF["leader_length"],
                      promoter_len=rc.CONF["promoter_length"]):
    """TODO"""

    genes = dict()

    for g in genome.features:
        if g.type == "CDS":
            name = gene_name = g.qualifiers["gene"][0]
            locus_tag = g.qualifiers["locus_tag"][0]

            if include_genes and name not in include_genes:
                if locus_tag in include_genes:
                    name = locus_tag
                else:
                    continue
            elif exclude_genes:
                if name in exclude_genes or locus_tag in exclude_genes:
                    continue

            #Skip pseudo genes unless they are specifically requested
            if not include_genes and "pseudo" in g.qualifiers:
                continue

            #Little bit of warning:
            #Bio converts positions to 0-index internally
            #This is mostly a good thing
            strand = g.location.strand
            cds = g.extract(genome).seq

            start = g.location.start
            end = g.location.end

            if strand not in [1, -1]:
                raise Exception("Invalid strand: {}. Use 1, or -1".format(strand))

            #Get leader and leader position
            leader, l_start, l_end = get_leader(genome.seq, start, end, strand, leader_len)
            try:
                leader = find_wobble_seq(genome, leader, l_start, l_end, config["codon_table"], config["dgn_table"])
            except WobbleError as e:
                log.error("Wobble error in gene.leader: {}; {}".format(name, e))

            #Detect operons
            in_operon = None
            if "operons" in config:
                for operon in config["operons"]:
                    if locus_tag in config["operons"][operon]["genes"]:
                        if in_operon is not None:
                            raise Exception("Gene {}/{} found in more than one operon. "
                                            "Check your config file.".format(gene_name, locus_tag))
                        in_operon = config["operons"][operon]

            if in_operon is None:
                operon_strand = strand
                operon_start = start
                operon_end = end
            else:
                operon_strand = in_operon["strand"]
                operon_start = in_operon["start"] - 1
                operon_end = in_operon["end"]

            #Extract promoter sequences
            promoter, p_start, p_end = get_leader(genome.seq, operon_start, operon_end, operon_strand, promoter_len)
            try:
                promoter = find_wobble_seq(genome, promoter, p_start, p_end, config["codon_table"], config["dgn_table"])
            except WobbleError as e:
                log.error("Wobble error in gene.promoter: {}; {}".format(name, e))

            if strand == 1:
                pos = int(start)
                promoter_pos = p_start - pos
            else:
                pos = int(end)
                promoter_pos = pos - p_end
                leader = leader.reverse_complement()
                promoter = promoter.reverse_complement()

            #If gene is already found and is specifically asked for.
            if name in genes and include_genes and name in include_genes:
                #log.warn("Gene {} found more than once.".format(name))
                raise ParserError("Gene {} found more than once. Use locus_tag "
                                  "[{}] instead.".format(name, locus_tag))
            else:
                gene = Gene(name, pos, strand, cds, leader, promoter, promoter_pos)
                gene.gene_name = gene_name
                gene.locus_tag = locus_tag
                gene.in_operon = in_operon is not None

                genes[name] = gene
                #Mixing genes and locus_tags are ok
                #If no include_genes are supplied, duplicate all genes
                if not include_genes or locus_tag in include_genes:
                    genes[locus_tag] = gene

    if include_genome:
        #Circumvent the default Sequence object to save memory
        genes["genome"] = Gene("genome", 0, 1, "A")
        genes["genome"].cds = genome.seq

    return genes


def get_leader(parent_seq, start, end, strand, length):
    """Extract a leader.

    Take strand into account, and wrap-around sequences that are close to
    either start or end of genome.

    """
    if strand not in [1, -1]:
        raise Exception("Invalid strand: {}. Use eiter 1 or -1.".format(strand))

    #Gene start, and leader start/end
    if strand == 1:
        leader_end = start
        leader_start = leader_end - length
    else:
        leader_start = end
        leader_end = leader_start + length

    leader = extract_circular(parent_seq, leader_start, leader_end)

    #Leader in genome context (+1)
    return Sequence(leader), leader_start, leader_end


def find_wobble_seq(genome, leader, l_start, l_end, codon_table, dgn_table):
    """Find a wobble sequence between l_start and l_end.

    Returns None if no wobble sequence if found

    """
    #Look for CDS in leader
    for c in genome.features:
        if c.type == "CDS":
            #Skip all pseudo genes
            if "pseudo" in c.qualifiers:
                continue
            if is_inside(l_start, l_end, c.location.start, c.location.end):
                w_cds = c.extract(genome).seq
                try:
                    w_seq = cds_to_wobble(w_cds, codon_table, dgn_table)
                except KeyError as e:
                    c_name = c.qualifiers["locus_tag"][0]
                    raise WobbleError("Invalid translation for [{}]: {}".format(c_name, e))

                #Wobble sequence in genome context (+1)
                if c.location.strand == -1:
                    w_seq = reverse_complement(w_seq)

                #TODO:
                #If there is disjonted range, this will fail.
                #I.e. join(3948583..3949566,3949565..3950227)
                #Calculate start offset
                start_offset = c.location.start - l_start
                leader.add_wobble(w_seq, start_offset)

    return leader

"""
Barcoding parsers
~~~~~~~~~~~~~~~~~
"""

def parse_barcode_library(barcode_filehandle):
    """Parse an open barcoding library file.

    Barcoding library file has the following syntax:

        >>> barcode_file = '''#Comments
        ...
        ... #Empty lines are allowed
        ... #First header is for primers:
        ... >PRIMERS
        ... #ID SEQUENCE
        ... F1 CTACCTTGCA #Comment on primer F1
        ... F2 GGAATTGAGA#Comment on this primer
        ... R1 CCGTCCGTTA
        ... R2 ATTTCCCTTG
        ...
        ... #Second header is for the library
        ... >LIBRARY
        ... #ID FWD REV
        ... ID1 F1 R1#Comment on ID1
        ... ID2 F2 R2 #Comment in ID2
        ... #Primers can be mixed
        ... ID3 F1 R2'''.splitlines()


    barcode_file must be open file-like iterable:

        >>> barcodes = parse_barcode_library(barcode_file)


    Barcodes are returned as a dictionary:

        >>> for id in sorted(barcodes):
        ...     (id,
        ...     barcodes[id]["forward"],
        ...     barcodes[id]["reverse"])
        ...
        ('ID1', 'CTACCTTGCA', 'TAACGGACGG')
        ('ID2', 'GGAATTGAGA', 'CAAGGGAAAT')
        ('ID3', 'CTACCTTGCA', 'CAAGGGAAAT')

    """
    primers = dict()
    barcodes = dict()

    primer_flag = False
    barcodes_flag = False

    error_list = list()

    for line in barcode_filehandle:
        if not line.strip() or line[0] == "#":
            continue

        #Filter comments
        line = line.split("#")[0].strip()

        #Raise parser flags
        if line == ">PRIMERS":
            primer_flag = True
            continue
        elif line == ">LIBRARY":
            barcodes_flag = True
            continue

        line = line.split()
        if barcodes_flag:
            #Get forward
            try:
                forward = primers[line[1]]
            except KeyError:
                error_list.append("Primer '{}' not found in barcode file".format(line[1]))
                continue
            #Get reverse
            try:
                reverse = reverse_complement(primers[line[2]])
            except KeyError:
                error_list.append("Primer '{}' not found in barcode file".format(line[2]))
                continue

            #Add to library
            barcodes[line[0]] = {"forward": forward, "reverse": reverse}

        #Collect primers
        elif primer_flag:
            primers[line[0]] = line[1]
        else:
            raise ParserError("Malformed barcode file")

    if error_list:
        raise ParserError("Malformed barcode file", error_list)

    return barcodes

"""
configuration parsers
~~~~~~~~~~~~~~~~~~~~~
"""

def create_config_tables(config, cfg_basedir="./"):
    """Create additional lookup tables from the parsed config yaml file."""
    #Translation table
    AAtable = dict() #{AA: [] for AA in amino_acids}
    config["codon_table"] = dict()
    config["stop_codons"] = list()
    for triplet, (AA, usage) in config["codon_usage"].items():
        config["codon_table"][triplet] = AA
        if AA in AAtable:
            AAtable[AA].append(triplet)
        else:
            AAtable[AA] = [triplet]
        if AA == "$":
            config["stop_codons"].append(triplet)

    #Reverse translation table
    config["dgn_table"] = dict()
    for AA in AAtable:
        config["dgn_table"][AA] = seqs_to_degenerate(AAtable[AA])

    if "DOOR_operon_source" in config:
        opr = os.path.join(cfg_basedir, config["DOOR_operon_source"])
        if not os.path.isfile(opr):
            raise Exception("Could not find opr file: {}".format(opr))
        with open(opr) as oprfile:
            operons = parse_DOOR_oprfile(oprfile)
        if "operons" not in config:
            config["operons"] = dict()
        config["operons"].update(operons)

    if "promoters" in config:
        for prom in config["promoters"]:
            for i, l in enumerate(config["promoters"][prom]):
                ntdict = dict()
                ntdict["A"],ntdict["C"],ntdict["G"],ntdict["T"] = float(l[0]), float(l[1]), float(l[2]), float(l[3])
                config["promoters"][prom][i] = ntdict

    return config


def parse_DOOR_oprfile(oprfilehandle):
    """Parse an operonfile in DOOR format."""
    header = None
    operons = dict()
    for line in oprfilehandle:
        if not line.strip(): continue
        if line[0] == "#": continue
        if not header:
            header = [h.lower() for h in line.split()]
            continue
        line = dict(zip(header, line.split()))
        if not line["operonid"] in operons:
            operons[line["operonid"]] = {"start": float("inf"), "end": float("-inf"),
                                         "strand": 1, "genes": []}
        operon = operons[line["operonid"]]
        #Convert strand
        operon["strand"] = 1 if line["strand"] == "+" else -1
        operon["start"] = min(operon["start"], int(line["start"]))
        operon["end"] = max(operon["end"], int(line["end"]))
        operon["genes"].append(line["synonym"])

    return operons


"""
Oligolist output
~~~~~~~~~~~~~~~~
"""

def oligolist_to_tabfile(oligolist, output):
    """Writeout oligolist to a tab separated file.

    Supply string filename or open file

    """
    cls = False
    if not hasattr("write", output):
        output = open(output, "w")
        cls = True

    for o in sorted(oligolist, key=lambda x:x.number):
        output.write(o.id() + "\t" + o.output() + "\n")

    if cls:
        output.close()


def oligolist_to_csv(oligolist, output=None):
    """Print report from oligolist in CSV format."""
    csvoutlist = list()

    op_reg = re.compile(r"\[(\w+)/(\w+)\] (line \d+)\s?(\S*)")

    for o in oligolist:
        m = re.match(op_reg, o.operation)
        try:
            op_cmd, gene, line, op_options = m.groups()
        except AttributeError:
            print("Operation not recognized: {}".format(o.operation))
            return
        csvoutlist.append([o.short_id(), op_cmd, gene, line, op_options,
                           str(o.mut), "+".join(o.barcode_ids), o.output(),
                           o.dG_fold] + o.operation_values)

    csvoutlist.sort(key=lambda x: (x[1], x[0]))

    if output:
        cls = False
        if not hasattr("write", output):
            output = open(output, "wb")
            cls = True

        #For Libreoffice Calc
        output.write(codecs.BOM_UTF8)
        csv_w = csv.writer(output)

        #TODO add additional headers
        headers = ["id", "operation", "gene", "line", "options", "mutation",
                   "barcodes", "oligo", "mfe", "wt", "altered"]
        csv_w.writerow(headers)

        for n in csvoutlist:
            csv_w.writerow(n)

        if cls:
            output.close()

    return [headers] + csvoutlist


def oligolist_to_mascfile(oligolist, masc_kwargs, mascfile=None):
    """Generate MASC PCR primers.

    Write primers to mascfile in plain text format if mascfile is supplied.

    """
    masc_primers = dict()
    #unique oligos
    un_oligos = list()
    for oligo in oligolist:
        ID = oligo.short_id().split(".")[0]
        if ID not in un_oligos:
            try:
                masc_primers[ID] = oligo.mut.MASC_primers(**masc_kwargs)
                un_oligos.append(ID)
            except Exception as ex:
                log.warning("Could not create MASC primers for {}/{}: {}".format(ID, oligo.mut, ex))
                continue
            masc_primers[ID]["mut"] = oligo.mut

    if mascfile:
        cls = False
        if not hasattr(mascfile, "write"):
            mascfile = open(mascfile, "wb")
            cls = True

        #For Libreoffice Calc
        mascfile.write(codecs.BOM_UTF8)
        csv_w = csv.writer(mascfile)

        headers = ["id", "mut", "forward(wt)", "forward(mut)"]
        for l in sorted(masc_kwargs["lengths"]):
            headers.append(str(l))
        csv_w.writerow(headers)

        for ID in un_oligos:
            prm = masc_primers[ID]
            line = [ID, prm["mut"].__str__(idx=1)]
            line.append(str(prm["fpwt"][0]))
            line.append(str(prm["fpmut"][0]))
            for l in sorted(masc_kwargs["lengths"]):
                line.append(str(prm[l][0]))
            csv_w.writerow(line)

        if cls:
            mascfile.close()

    return masc_primers


class OligoLibraryReport:
    """Print a report from a csv report.

    This is meant to be general, and new file format writers can easily be
    plugged into the class

    """
    def __init__(self, project):
        self.project = project
        self.sections = dict()
        self.parsed = False

    def parse_and_generate(self, csv_in, csv_file=False):
        """Wrapper to parse csv and generate report"""
        if csv_file:
            self.parse_csvfile(csv_in)
        else:
            self.parse_csvlist(csv_in)

        self.generate_report()

    def parse_csvlist(self, csvlist):
        """Parse input csvlist"""
        self.oplib = dict()
        header = False
        for line in csvlist:
            if not header:
                header = line
                header[0] = header[0].strip(codecs.BOM_UTF8)
                continue
            entry = {h: line[i] for i,h in enumerate(header) if len(line) > i}

            op = entry["operation"]
            gene = entry["gene"]
            if op not in self.oplib:
                self.oplib[op] = dict()
            if gene not in self.oplib[op]:
                self.oplib[op][gene] = list()

            self.oplib[op][gene].append(entry)

        self.parsed = True

    def parse_csvfile(self, csvfile):
        """Parse input csv file"""
        with open(csvfile) as csvf:
            reader = csv.reader(csvf)
            self.parse_csvlist(reader)

    def generate_report(self):
        """TODO"""
        if not self.parsed:
            raise Exception("generate_report called with parsed data")

        self.generate_rbs_lib_section()

    """
    Specific sections
    """

    def generate_rbs_lib_section(self, prefix="01", plot_log=True):
        """RBS_library stats generation"""
        from matplotlib import pyplot as plt
        import math

        if "RBS_library" not in self.oplib:
            return

        elements = [("p", "Ribosome binding sites libraries.")]

        bar_width = 0.3
        #sort by altered
        sba = lambda gene: max([float(oligo["altered"]) for oligo in self.oplib["RBS_library"][gene]])
        genes = sorted(self.oplib["RBS_library"].keys(), key=sba, reverse=True)

        #Genes per plot
        gpp = min(range(6,11), key=lambda i: len(genes) % i)
        if len(genes) % gpp != 0:
            gpp = max(range(6,11), key=lambda i: len(genes) % i)
        fig_w = 8.
        fig_h = 3.
        gw = fig_w/gpp

        for i in range(0, len(genes), gpp):
            curr_genes = genes[i:i+gpp]
            x = [0]
            y1 = list()
            y2 = list()
            labels = list()
            mi = 1
            for gene in curr_genes:
                wt = float(self.oplib["RBS_library"][gene][0]["wt"])
                expr = max([float(oligo["altered"]) for oligo in self.oplib["RBS_library"][gene]])
                min_expr = min([float(oligo["altered"]) for oligo in self.oplib["RBS_library"][gene]])
                x.append(x[-1]+1)
                y1.append(wt)
                y2.append(expr)
                labels.append(gene)
                mi = min([mi, min_expr])

            del(x[0])

            fig = plt.figure(figsize=(fig_w, fig_h))
            ax = fig.add_subplot(111)
            bar_wt = plt.bar(x, y1, bar_width, color="#9FF33D", linewidth=0, bottom=0, log=plot_log)
            bar_lib = plt.bar([x1+bar_width for x1 in x], y2, bar_width, color="#007633", bottom=0, linewidth=0, log=plot_log)

            ax.set_xticklabels(labels)
            ax.set_xticks([x1+bar_width for x1 in x])



            if plot_log:
                ma = max(y1 + y2)*1.1#round(max(y2), len(str(max(y2))))
                ax.set_ylim(mi, ma*2)

            ax.set_xlim(1, len(curr_genes)+1)

            for gene, j in zip(curr_genes, x):
                step = 0.5 / (len(self.oplib["RBS_library"][gene]))
                clr = 0.5
                for oligo in self.oplib["RBS_library"][gene]:
                    color = "{:.1f}".format(clr)
                    y3 = float(oligo["altered"])
                    ax.annotate("", xy=(j+1.5*bar_width, y3), xytext=(j+0.75+bar_width*0.5, y3),
                                arrowprops={"arrowstyle": "wedge", "linewidth": 0,
                                            "color": color})
                    clr -= step

            plt.subplots_adjust(top=0.85, left=0.08, right=0.95*len(curr_genes)/gpp)

            plt.ylabel("log(expression [AU])" if plot_log else "expression [AU]")
            t = fig.text(0.08, 0.90, "RBS library expression levels",
                         horizontalalignment='left', fontsize=13)
            #plt.title("RBS library expression levels")

            prx_arr = plt.Rectangle((0, 0), 1, 1, fc="0.5", linewidth=0)

            if 0.99 < len(curr_genes)/gpp < 1.01:
                plt.legend([bar_wt, bar_lib, prx_arr], ["wt", "lib", "oligos"],
                    bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
                    ncol=3, borderaxespad=0., frameon=False)
            else:
                plt.legend([bar_wt, bar_lib, prx_arr], ["wt", "lib", "oligos"],
                    bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=False)

            imgdata_png = StringIO()
            imgdata_pdf = StringIO()
            fig.savefig(imgdata_png, format="PNG", dpi=300)
            fig.savefig(imgdata_pdf, format="PDF")
            imgdata_png.seek(0)
            imgdata_pdf.seek(0)

            elements.append(("png/pdf", imgdata_png, imgdata_pdf, "page", fig_h/fig_w))

        self.sections[prefix + "RBS_library"] = elements

    """
    Output formats
    """

    def write_pdf(self, report_out):
        """Write report directly to PDF"""
        #Local imports since reports are not strictly necessary
        from reportlab import platypus as RL
        from reportlab.lib.styles import getSampleStyleSheet
        from reportlab.lib.pagesizes import A4

        styles = getSampleStyleSheet()

        page_width, page_height = A4
        rmargin = lmargin = tmargin = 72
        bmargin = 18
        doc = RL.SimpleDocTemplate(report_out, pagesize=A4,
                                   rightMargin=rmargin, leftMargin=lmargin,
                                   topMargin=tmargin, bottomMargin=bmargin)

        elements = [RL.Paragraph("Project: " + self.project, styles["Title"])]

        for h in sorted(self.sections):
            header = h.lstrip("0123456789 ")
            e = self.sections[h]
            elements.append(RL.Paragraph(header, styles["Heading1"]))
            for t in e:
                if t[0] == "p":
                    elements.append(RL.Paragraph(t[1], styles["Normal"]))
                if t[0] == "png/pdf":
                    if t[3] == "page":
                        img_width = page_width - rmargin - lmargin
                        img_height = img_width*t[4]
                    PdfImage = self.PdfImage()
                    if PdfImage:
                        img = PdfImage(t[2], width=img_width, height=img_height)
                    else:
                        log.debug("pdfrw not found, using png images in output report.")
                        t[1].seek(0)
                        img = RL.Image(t[1], width=img_width, height=img_height)
                    elements.append(img)

        doc.build(elements)

    def PdfImage(self):
        """PDF flowable for report lab

        Wrapped in a try statement since it can fail silently if pdfrw is
        not installed, and wrapped inside a def statement to shield the logger
        from the pdfrw logger, which overwrites a logger config.

        """
        try:
            from pdfrw import PdfReader
            from pdfrw.buildxobj import pagexobj
            from pdfrw.toreportlab import makerl

            from reportlab.platypus import Flowable
            from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_RIGHT
        except ImportError:
            return False

        class PdfImageInner(Flowable):
            """PdfImage wraps the first page from a PDF file as a Flowable
            which can be included into a ReportLab Platypus document.
            Based on the vectorpdf extension in rst2pdf code.google.com/p/rst2pdf

            Hijacked from: http://stackoverflow.com/a/13870512

            """

            def __init__(self, filename_or_object, width=None, height=None, kind='direct'):
                # If using StringIO buffer, set pointer to begining
                if hasattr(filename_or_object, 'read'):
                    filename_or_object.seek(0)

                page = PdfReader(filename_or_object, decompress=False).pages[0]
                self.xobj = pagexobj(page)
                self.imageWidth = width
                self.imageHeight = height
                x1, y1, x2, y2 = self.xobj.BBox

                self._w, self._h = x2 - x1, y2 - y1
                if not self.imageWidth:
                    self.imageWidth = self._w
                if not self.imageHeight:
                    self.imageHeight = self._h
                self.__ratio = float(self.imageWidth)/self.imageHeight
                if kind in ['direct','absolute'] or width==None or height==None:
                    self.drawWidth = width or self.imageWidth
                    self.drawHeight = height or self.imageHeight
                elif kind in ['bound','proportional']:
                    factor = min(float(width)/self._w,float(height)/self._h)
                    self.drawWidth = self._w*factor
                    self.drawHeight = self._h*factor

            def wrap(self, aW, aH):
                return self.drawWidth, self.drawHeight

            def drawOn(self, canv, x, y, _sW=0):
                if _sW > 0 and hasattr(self, 'hAlign'):
                    a = self.hAlign
                    if a in ('CENTER', 'CENTRE', TA_CENTER):
                        x += 0.5*_sW
                    elif a in ('RIGHT', TA_RIGHT):
                        x += _sW
                    elif a not in ('LEFT', TA_LEFT):
                        raise ValueError("Bad hAlign value " + str(a))

                xobj = self.xobj
                xobj_name = makerl(canv._doc, xobj)

                xscale = self.drawWidth/self._w
                yscale = self.drawHeight/self._h

                x -= xobj.BBox[0] * xscale
                y -= xobj.BBox[1] * yscale

                canv.saveState()
                canv.translate(x, y)
                canv.scale(xscale, yscale)
                canv.doForm(xobj_name)
                canv.restoreState()

        return PdfImageInner


if __name__ == "__main__":
    import doctest
    doctest.testmod()
