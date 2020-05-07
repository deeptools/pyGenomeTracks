bed
==========

Description
-----------

A track for all bed-like files (with first column chromosome, second start and third end). If the other columns fit the requirement as defined in `UCSC <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_ additional fields can be used.
For example, the 5th and 9th column can be used to change the color of intervals, the 6th column indicate the strand...
By default, intervals without strand are displayed as rectangle and for intervals with strand an arrow is added at the extremity (not included in the interval). In case of BED12 format, the introns are displayed into another color. But other styles are available.

Parameters
----------

.. include:: auto/bed_deduced_from_code.txt

