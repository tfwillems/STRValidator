STRValidator
============

STRValidator is a simple web application designed to facilitate the comparison of STR call sets.
Given two sets of STR calls which should be identical, it facilitates their comparison
by generating a bubble plot of length_in_call_set_a vs length_in_call_set_b. Each bubble is linked to the resulting calls 
such that the user can click on a bubble and view the alignments of the underlying reads. This tool should be useful 
when debugging errors in STR genotyping algorithms. For instance, some useful comparisons may involve:
    1. Pairs of father-son Y-chromosome calls
    2. Autosomal calls generated using capillary electrophoresis and whole-genome sequencing data for the same samples

To view the command line syntax required to run the program, simply type
   ./webserver.py -h

