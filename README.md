# Overview
 The uploaded files are to simulate and demonstrate the capability of an auxiliary method called the prediction-correction algorithm. The algorithm helps to avoid a local minimum solution of a multi-objective optimization problem and to search for one of its near-global optimum points. The key idea of the auxiliary algorithm is to investigate the effect of the sign of the derivative of each sub-objective function and to change the sign if necessary by considering the difference between former and current sub-objective functions.


# Descriptions
 There are 6 Matlab files to run the simulation.

### Test_function_2D_CG.m
The Matlab file for solving 4 multi-objective problems composed of two design variables with the conjugate gradient algorithm.

### Test_function_2D_MMA.m
The Matlab file for solving 4 multi-objective problems composed of two design variables with the MMA algorithm.

### Test_function_2D_MMAPC.m
The Matlab file for solving 4 multi-objective problems composed of two design variables with the MMA and prediction-correction algorithms.

### Test_function_5D_CG.m
The Matlab file for solving 4 multi-objective problems composed of five design variables with the conjugate gradient algorithm.

### Test_function_5D_MMA.m
The Matlab file for solving 4 multi-objective problems composed of five design variables with the MMA algorithm.

### Test_function_5D_MMAPC.m
The Matlab file for solving 4 multi-objective problems composed of five design variables with the MMA and prediction-correction algorithms.


# Miscellaneous
There are multiple txt files about optimization outputs as well as minor Matlab files for the subroutines to run the simulation.
