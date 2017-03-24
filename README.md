bltools: Unix style tools for biological sequences
======================================================

I've often wanted to be able to use "grep" to search a FASTA file, but
the standard Unix text-manipulation tools are line-based rather than
sequence-record-based. The goal of this software is to create a set of
tools that perform the same functions as head, tail, grep, sort, etc.,
but are aware of the structure of sequence files.

blgrep: Grep for biological sequences
--------------------------------------

blgrep is a C++ program based on the Seqan library that performs regex
searches on biological sequences. It is meant to act very much like
grep, except that it works on sequence records instead of lines in a
file.

blhead and bltail
-------------------

Very much like head and tail, but for sequences. Note that there has to
be a space between `-n' and the value, unlike the real tail and head
commands.

blwc
-----

Counts the number of records in a file (by default), or the length
of each record.
