Contents:

.. toctree::
   :maxdepth: 2

User manual
===========

Dependencies
------------
 * Biopython 1.6.0 or higher.
 * ViennaRNA 2.0 or higher.


Program input
-------------
The program needs 4 input files:
 * A file with gene adjustments
 * A genbank file with the whole genome
 * A file with barcode information
 * A config file with information about the organism.


Gene adjustment file
~~~~~~~~~~~~~~~~~~~~
The gene adjustment file is a space or tab separated file.  Where every line corresponds to an operation and look as follows::

    Gene-Name    Operation   Barcode-ID   Options

**Example**: ::

    thiM start_codon_optimal ID1+ID2,ID2+ID3
    yaaA custom_mutation     ID2              mut=TATCAACGCC[GCTC=]GCTTTCATGACT
    thiC translational_KO    ID1              KO_frame=10
    thiE RBS_library         ID1
    thiF RBS_library         ID3+ID4          max_mutations=5
    thiS RBS_library         ID1              n=15,passes=3

It is important to notice that there are **NO** spaces in the options
parameter. **Any option after a space will be omitted.**

Multiple barcodes are supported, they should be separated by a comma (\ ``,``\
) or a plus (\ ``+``\ ). Comma separated ID’s give multiple oligos with
different barcodes, whereas plus separated ID’s result in multiple barcodes on
the same oligo.

I.e. thiM in the list above will result in two identical copies of the created
oligo, with different barcodes. One with barcodes specified as ID1 and ID2 in
the barcode file. And one with ID2 and ID3.


List of operations
~~~~~~~~~~~~~~~~~~

 - Translational modifications

   - ``RBS_library``
   - ``translational_knockout``
   - ``start_codon_optimal``
 - Custom mutations

   - ``dna_mutation``
   - ``residue_mutation``

RBS_library
^^^^^^^^^^^
.. autosimple:: mage_tool.operations.translation.RBSLibrary

translational_knockout
^^^^^^^^^^^^^^^^^^^^^^
.. autosimple:: mage_tool.operations.translation.TranslationalKnockout

start_codon_optimal
^^^^^^^^^^^^^^^^^^^
.. autosimple:: mage_tool.operations.translation.StartCodonOptimal

dna_mutation
^^^^^^^^^^^^
.. autosimple:: mage_tool.operations.manual.DNAMutation

residue_mutation
^^^^^^^^^^^^^^^^
.. autosimple:: mage_tool.operations.manual.ResidueMutation


Barcoding file
~~~~~~~~~~~~~~

The barcoding file consists of a list of primers followed by a list of barcode
ID’s combining a forward and reverse primer.  They must follow the headers
``>PRIMERS`` and ``>LIBRARY`` respectively. The ID’s are used in the gene
adjustments file to add the desired barcodes. ID’s should only contain one
forward and one reverse primer.  ::

    >PRIMERS
    F1  Sequence
    F2  Sequence
    R1  Sequence
    R2  Sequence
    ..  ..

    >LIBRARY
    ID1 F1  R1
    ID2 F2  R2
    ID3 F1  R1
    ..  ..   ..


Running MODEST.py
-----------------
MODEST.py take three required arguments: ``adjustments``, ``barcodes`` and ``config``.
``adjustments`` is is the adjustment list, and ``barcodes`` is
the barcoding library mentioned above. ``config`` is the genome configuration file.
Since a genome and a genbank file is linked, MODEST.py will automatically attempt
to find a ``.gb`` file named after the config file.

The full search path for locating a genome file is:

 * ``<cfg>.(gb|genbank)``
 * ``<cwd>/<locus>[.gb|.genbank]``
 * ``<cfg_path>/<locus>[.gb|.genbank]``

Where:

 * ``<cfg>`` is the basename (anything but the extension) of the config file. (I.e. given config file ``path/e_coli.config``, ``<cfg>`` is ``path/e_coli``)
 * ``(gb|genbank)`` means MODEST.py will look for the two file extensions ``.gb`` or ``.genbank``.
 * ``<cwd>`` is current working directory
 * ``<locus>`` is the genbank ID described in the locus field in the config file.
 * ``<cfg_path>`` is the directory where the config file resides.
 * ``[.gb|.genbank]`` means either no extension, ``.gb``, or ``.genbank`` extension.

The ``-p`` or ``--project`` flag gives the project a name. This name is used in
the output id's as well as output filenames.

To output a PDF report, use the ``--PDF`` flag, optionally with an argument which
will be the pdf file that is written.

Similarly, to output MASC primers, use the ``--MASC`` flag with an optional output file.

How output is logged can be configered with the ``--log`` parameter. Default is
to output log to a file called <project>.log. Specifying a filename will log output
to that file. Specifying ``-`` will disable logging and ``stdout`` will print log
events to screen.

A typical procject folder will look like this::

  Project
  ├─ data
  │  ├─ organism.config
  │  ├─ organism.gb
  │  ├─ adjustments.txt
  │  ├─ barcodes.txt
  │  └─ operon_file.opr
  └─ <output files>

As ``Project`` as your current working directory, MODEST.py would be invoked like so::

  $ MODEST.py data/adjustments.txt data/organism.config data/barcodes.txt -p project_name