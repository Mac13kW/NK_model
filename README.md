NK_model
========

New set of NK model exercises, created on July 2, 2018
Files:
* 1_landscape_creation.py
* 2_local_search.py
* 3_decentralized.py

The three files form a set of exercises which are aimed at introducing
the user to the NK model. The first file creates NK landscapes and reports
simple summary statistics. The other two files are different versions of
search on those landscapes and report overall fitness level produced by
that search. The code has been adjusted to work on Python 2.7 and 3.6 without
additional changes.



NK_model_basic (ver 3.1) (old version)

This is an example of an NK model. The simulation uses an NK Model
to represent a task environment of varying degree of complexity (K).
It produces i NK landscapes, together with all combinations, contributions
of each combination of decision variables and information whether a
given location is a local peak (all combinations, which differ by only
one decision from the focal one have lower overall fitness value)

You can change the N and K of the model, but be careful with the custom
interaction matrices, as you need to adjust the N and K manually.

I set N = 4 for simplicity.
K must be an integer from 0 to N-1

The model was rewritten on February 18th 2013 and later on
April 13th 2013 and December 11th, 2014

- The model has custom modules for interaction matrices
- different hamming distance settings

Updated on February 19th 2013 (ver 2.0)
Updated on April 13th 2013 (ver 3.0)
Updated on December 11th, 2014 (ver 3.1)

Options:
- Press CNTRL+C to interrupt


Installation
------------------------------------
The script has been tested with Python 2.7 (Anaconda Python Distribution - Windows 8.1 64-bit)

Required modules:
- numpy
- itertools
- time
