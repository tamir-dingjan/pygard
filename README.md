# pygard
A Python implementation of the Graded Autocatalysis Replication Domain (GARD) algorithm.

## Background
The GARD algorithm simulates the association of amphiphiles in solution, exploring the
composition of micelles over generations of growth and splitting. The algorithm was 
developed by Prof. Doron Lancet at the Weizmann Institute of Science; please see the following
publications:


## Input Parameters
This version of the algorithm implements a number of basic parameters:
- N: starting number of amphiphiles in a micelle. Micelle splitting occurs at 2N.
- Ng: total number of distinct types of amphiphile in the system.
- K: intrinsic affinity matrix of each amphiphile type toward each other amphiphile type.
- B: enhancement matrix for each amphiphile type toward each other amphiphile type.


## Algorithm
After providing the above input parameters, the algorithm begins by instantiating 
a micelle of size N with a random composition. Each step of the algorithm calculates
the weighted affinity of the micelle toward each of the amphiphile types in the system,
and subsequently adds to the micelle a single amphiphile of the most favoured amphiphile type.

Once the micelle reaches a population of 2N, it will split into two daughter micelles of size N.
The composition of the daughter micelles is a random sample of the composition of the parent:
that is, each amphiphile in the parent is assigned to one of the daughters.

The algorithm can be run for as many steps as you have memory to contain the ever-growing micelles.

## Analysis
One of the most interesting properties to analyse from the behaviour of the GARD model is
the generation of 'composomes', generationally-consistent regions in amphiphile composition space.
That is, the property that certain compositions of micelles reinforce their own composition as they
grow and divide. To enable ease of analysis for this property, one of the features to be added
will be a facile way to keep track of micellar compositions over algorithm steps (and hopefully
also pretty ways to graph it).