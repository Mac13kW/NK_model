NK_model
========

NK model (Version 3.1)

This is a simple example of an NK model. The simulation uses an NK Model
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
The script has been written and tested using Python 2.7 (Anaconda Python Distribution)
