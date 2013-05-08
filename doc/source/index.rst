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

**Example**:
::
    thiM start_codon_optimal ID1+ID2,ID2+ID3
    yaaA custom_mutation     ID2              mut=TATCAACGCC[GCTC=]GCTTTCATGACT
    thiC translational_KO    ID1              KO_frame=10
    thiE RBS_library         ID1
    thiF RBS_library         ID3+ID4          max_mutations=5
    thiS RBS_library         ID1              n=15,passes=3

It is important to notice that there are **NO** spaces in the options parameter. **Any option after a space will be omitted.**

Multiple barcodes are supported, they should be separated by a comma (\ ``,``\ ) or a plus (\ ``+``\ ). Comma separated ID’s give multiple oligos with different barcodes, whereas plus separated ID’s result in multiple barcodes on the same oligo.

I.e. thiM in the list above will result in two identical copies of the created oligo, with different barcodes. One with barcodes specified as ID1 and ID2 in the barcode file. And one with ID2 and ID3.


List of operations
~~~~~~~~~~~~~~~~~~

Translational modifications
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: mage_tool.interface.RBS_library
.. autofunction:: mage_tool.interface.translational_KO
.. autofunction:: mage_tool.interface.start_codon_optimal

Custom mutations
^^^^^^^^^^^^^^^^

.. autofunction:: mage_tool.interface.gene_mutation
.. autofunction:: mage_tool.interface.residue_mutation


Barcoding file
~~~~~~~~~~~~~~
The barcoding file consists of a list of primers followed by a list of barcode ID’s combining a forward and reverse primer.  They must follow the headers ``>PRIMERS`` and ``>LIBRARY`` respectively. The ID’s are used in the gene adjustments file to add the desired barcodes. ID’s should only contain one forward and one reverse primer.
::

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


Running the program
-------------------
TDB.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
