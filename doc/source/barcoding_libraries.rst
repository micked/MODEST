Barcoding
=========

Barcoding file
--------------

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



