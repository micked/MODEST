MODEST general
==============

Input
-----
The program needs 3 input files or 4 if using barcodes:

 * A file with gene adjustments
 * A genbank file with the whole genome
 * A genome config file with additional information about the organism.
 * Optional: A file with barcode information

.. _gene-operations-file:

Gene operations file
--------------------
The gene operations file is a space or tab separated file.  Where every line
corresponds to an operation and look as follows::

    Gene-Name    Operation   Options

**Example**: ::

    thiM start_codon_optimal barcodes=ID1+ID2;ID2+ID3
    yaaA custom_mutation     mut=TATCAACGCC[GCTC=]GCTTTCATGACT,barcodes=ID2       
    thiC translational_KO    KO_frame=10,barcodes=ID1
    thiE RBS_library         barcodes=ID1
    thiF RBS_library         max_mutations=5,barcodes=ID3+ID4
    thiS RBS_library         n=15,passes=3,barcodes=ID1

It is important to notice that there are **NO** spaces in the options
parameter. **Any option after a space will be omitted.**

Some operations will suppert the ``genome`` "gene" as input, such as ``mutation``,
``deletion`` and ``insertion``.

Multiple barcodes are supported, they should be separated by a semicolon (\ ``;``\
) or a plus (\ ``+``\ ). Semicolon separated ID’s give multiple oligos with
different barcodes, whereas plus separated ID’s result in multiple barcodes on
the same oligo.

I.e. thiM in the list above will result in two identical copies of the created
oligo, with different barcodes. One with barcodes specified as ID1 and ID2 in
the barcode file. And one with ID2 and ID3.

A list of all the operations and corrosponding documentation can be found
here: :doc:`/operations`


Genome and genome configuration files
-------------------------------------
MODEST can read genome files in GenBank format
(`specification <http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html>`_).
Additional information is however required, such as replication
origin/termination.

Genome configuration files (usually ending in .gbcfg), are YAML format text
files. MODEST includes tools to generate genome configuration files with
minimal input. Read more about :doc:`/genome_config`.


Barcoding
---------
MODEST supports optional barcoding of oligos. Barcodes are specified as barcode
ids using the ``barcodes=<ID>`` option.

Read more about :doc:`/barcoding_libraries`.