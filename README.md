# Lab 13: Final Project
Zane Rothe

## Overview
This program uses an iterative approach to determine the optimal stacking sequence for a composite material layup based on the loads applied. This program is limited to non-woven, unidirectional lamina. Composite properties are determined using the Halpin-Tsai method based on costituent properties. A stress state is determined by the user-supplied loading conditions. A change of basis is performed for each lamina in the stacking sequence accoridng to the current angle interation. The Tsai-Hill failure criteria is employed for finding the safety factor.

## Features
This program allows the user to populate a csv file containing material properties, laminate dimensions, applied loads, and calculation constraints. the user can simply issue the program with the number of layers as the first argument and the file path as the second argument. For example: `$ ./composite 4 data.csv` will find the four-layer stacking sequence that provides the largest factor of safety for the parameters given in data.csv. A progress indicator is updated in the terminal for every 1% of program progress. Provided below is a flowchart detailing the computations.

## Code Structure
The main program first reads in the provided csv file. Intermediate calculations are made and then data is populated within a struct. The heart of the program lies within a resursive for-loop function. The number of nested for loops is equal to the number of layers in the composite. The for loops iterate through all of the possible layer angle combinations, calculating the safety factor for each configuration. The best safety factor and associated stacking sequence is displayed. To speed up the computations, multiprocessing is employed, using the maximum number of computer cores that satisfy resolution requirements. Mutex locks are used to ensure data is not lost across processes.

## Performance
A similar program was developed in MATLAB side-by-side with this C program to verify correct operation. The purpose of creating the program in C was to decrease the processing time needed to acheive the same result. This was accomplished as shown below. The same input parameters were given to both programs. the MATLAB program took 19.38 seconds whereas this C program took only 0.36 seconds.

## Techniques Used
- File IO
- Multiprocessing
- Shared Memory
- Mutex Locks (not used in any other labs)

