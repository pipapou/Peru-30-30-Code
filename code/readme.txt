This python program aims to propose maps extending the network of Peruvian protected and conserved areas from 17% (current) to 30%, 
by taking into account biodiversity, carbon, water, human impact, ecoregion representativeness (during step 1), structural connectivity (during step 2). 

It is a 2-step program that integrates optimization methods: Integer Linear Programming (ILP) during step 1 and Constraint Programming (CP) during step 2.

The "data" and "results" folders contain the program's input data and the generated outputs, respectively.

There are 5 python files:

- main_3030.py: the file to execute to run the entire program

- preprocessing: formatting of data used as input

- reference_values: this script generates the best and worst possible values ​​for each variable to optimize, for the reference point method used in step 1. 
We individually optimize each model variable under the same parameters and constraints as step 1.

- a_ilp: step 1, using ILP

- b_cp: step 2, using CP

