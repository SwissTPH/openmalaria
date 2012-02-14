Tab-deliminated fields format
==================

Author: Diggory Hardy
Copyright: 2010 Swiss Tropical and Public Health Institute (Swiss TPH)
Licence:
Date: Jan 2010

This is a formal description of a data format for future data output by openmalaria. It is intended:

*   to be similar to existing formats used
*   to be easily read into a spread-sheet program and analysed
*   to convey 2 or more dimensional data content (i.e. per measure, per time-point and optionally
    per measure-specific extra dimensions)
*   to be space efficient
*   to be easy to generate


Lexical format
---------------------

The file is separated into lines, by \n, \r or \r\n. Blank lines are ignored. (This yields a set of
 lines.)

For each line yielded by the above, it is separated into fields by tab (\t) characters. These form
deliminators; fields may be empty. Each line must have the same number of fields.

These line-by-line fields form a table (each line being one row).


Table header
--------------------

The first row of the table is a header. Each item is a textual description of the contents of the
column below it. Names may consist of alphanumeric characters, spaces and underscores:
"a-zA-Z_ ". This name in each column's header describes the "measure" dimension.

To allow extra dimensions, the "measure" name of a column header may be appended by one or more
extra names within square brackets, describing the coordinate within the extra dimensions of the
column. That is, the header cell contains contents of the form:

    NAME ( '[' NAME ']' )*

where each occurance of NAME matches the regular expression `[a-zA-Z_ ]+` (thus the position
within each dimension is an alpha-numeric identifier).


Time dimension
-----------------------

The first column within the table should be the time-point, and have name "time" (unit system is not
specified). Each value within this column should be unique, and greater than the previous value, to
allow it to be used as a row key instead of the row number.


Table values
------------------

All values within the table excluding the header line (but including the "time" column) should be
floating point numbers.


Example
------------

    time avg_temp temp[0][0] temp[0][1] temp[1][0] temp[1][1]
    0
