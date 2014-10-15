perfectread
===========

Filter out erroneous reads from a collection of reads with substitution errors, keeping only perfect (error-free) reads. This project is a basic implementation of the algorithm in [].

Installation
------------

The project require the Jellyfish software (version 2.1.4) to be already installed. The details of Jellyfish can be found at

Update the variable JELLYFISH in ``src/Makefile`` to point to the top level directory of Jellyfish 2.1.4. Compile using

    $ make

Run
---
Update the variable JELLYFISH in ``src/perfectread`` to point to the top level directory of Jellyfish 2.1.4.

    $ src/perfectread -d -t 8 -k 24 <fastq file>
