MODEST web interface
====================

The MODEST frontpage can be split in several sections:

The first text field is used to specify the name of the project. This name will be
used for the files generated as well as for naming the oligos designed with
MODEST.

Oligos created with MODEST is optimal for use with MAGE and recombineering methods
where the methyl-directed mismatch repair (MMR) system is unfunctional e.g. by
using a MutS mutant strain or a system where the MMR system can be temporarily
knocked out. Please refer to the following articles for more information about
strains, effeciency and the MMR system:

    1. Wang, H. H., Isaacs, F. J., Carr, P. A., Sun, Z. Z., Xu, G., Forest, C. R.
       and Church, G. M. (2009) *Programming cells by multiplex genome engineering
       and accelerated evolution.* Nature, 460, 894–8.
    2. Nyerges, A., Csörgo, B., Nagy, I., Latinovics, D., Szamecz, B., Pósfai, G.
       and Pál, C. (2014) *Conditional DNA repair mutants enable highly precise
       genome engineering.* Nucleic Acids Res., 10.1093/nar/gku105.
    3. Sawitzke, J. A., Costantino, N., Li, X.-T., Thomason, L. C., Bubunenko, M.,
       Court, C. and Court, D. L. (2011) *Probing cellular processes with
       oligo-mediated recombination and using the knowledge gained to optimize
       recombineering.* J. Mol. Biol., 407, 45–59.

Selecting organism
------------------

Step 1 is selecting the desired organism. By default the organism E. coli str.
K12 substr. MG1655 is selected. Other organisms can be chosen either by
uploading a genbank file or by specifying a genbank identifier. When using other
organisms than the default, you are required to specify replication origin and
termination regions. These are necessary in order for MODEST to target the
lagging strand of the replication fork when designing oligos. This can also be
specified by uploading a genome configuration file. MODEST will load a list of
genes from the genbank file. These are used to specify mutations in Step 2.

Specifying desired changes
--------------------------

Step 2 can be used to select genes and operations. When an operation
has been selected, the appropriate options are shown. These can then be
specified. Clicking the "Add line" button will then append this line in the
text field.

The main input text field contains the input file for MODEST. This is where operations
are written. It is also possible to write the operations manually or copy from a
text document. This input field can be downloaded as a text file.

More about the different: :doc:`/operations`.

More about the: :ref:`gene-operations-file`.

It is possible to download the input field as a text file for later use by
clicking the "Download" button. It is also possible to upload a text file to use
as input. *NOTICE*: Uploading a text file will override the operations
specified in the input field!

Advanced options such as oligo length can be set be expanding the advanced options in step 3.
Here it is possible to upload or paste a barcoding library as well as specify the length of the designed oligos
(standard is 90 bp).

Extra information
-----------------

The side bar contains documentation about the chosen operation. Here examples of
the operations are given and what the different options are used for
can be found.

Gene information is also displayed in the sidebar. The coding DNA
sequence can be found, as well as information about the genome location
and predicted expression level. This level can be used to determine a target
for RBS libraries.

After the inputs have been given. The job can be submitted by clicking the
"Submit" button.

After submission you will be taken to the job page. This page contains
information about the job. While the job is running, you can see the logfile.
Here you can see how many of the operations have been completed. The
page will update itself every 30 seconds, as long as the job is running.