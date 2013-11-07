MODEST commandline tool
=======================

.. contents::
   :local:
   :depth: 1

Installation
------------

Installing ``mage_tool`` has the following dependencies:

 * Biopython 1.6.0 or higher.
 * ViennaRNA 2.0 or higher.


ViennaRNA
~~~~~~~~~

ViennaRNA is a required component of MODEST.py..


Running MODEST.py
-----------------
MODEST.py take two required arguments: ``adjustments`` and ``config``.
``adjustments`` is is the adjustment list, and ``config`` is the genome configuration file.
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

A typical project folder will look like this::

 Project
 |-- data
 |   |-- organism.gbcfg
 |   |-- organism.gb
 |   |-- adjustments.txt
 |   |-- barcodes.txt
 |   ·-- operon_file.opr
 ·-- <output files>

As ``Project`` as your current working directory, MODEST.py would be invoked like so::

  $ MODEST.py data/adjustments.txt data/organism.config data/barcodes.txt -p project_name

Configuring mage_tool
~~~~~~~~~~~~~~~~~~~~~
``mage_tool`` can be configured via rc files. rc files are either placed in
``~/.conf/mage_tool.rc`` or ``./mage_tool.rc``. rc files are YAML files which
overrides the default program parameters::

 "oligo_length": 90
 "processes": None
 "operations": []
