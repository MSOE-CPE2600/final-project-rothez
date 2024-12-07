# Lab 13: Final Project
Zane Rothe

## Overview
This program uses an iterative approach to determine the optimal stacking sequence for a composite material layup based on the loads applied. This program is limited to non-woven, unidirectional lamina. Composite properties are determined using the Halpin-Tsai method based on costituent properties. A stress state is determined by the user-supplied loading conditions. A change of basis is performed for each lamina in the stacking sequence according to the current angle iteration. The Tsai-Hill failure criteria is employed for finding the safety factor. The diagram below shows a rough outline of the computations performed within the program.

![Screenshot 2024-12-07 140204](https://github.com/user-attachments/assets/84081ab9-c2c3-4696-b9a2-5b8c14c7a04f)

## About Composites
The term "composite material" is a broad term that encompasses a wide variety of materials. A composite material is defined as having two distinct, bonded materials, one of which acts as a matrix, and the other which acts as a reinforcement. Common matrix materials are epoxies and thermoplastics. Common reinforcements are carbon fiber, glass fiber, and Kevlar. Reinforcement fibers are often arranged in sheets, called lamina. A lamina can either be woven or unidirectional. This program will be limited to unidirectional lamina, in which all fibers in a lammina are parallel. Multiple laminae can be stacked to create a laminate. By varying the angle between the fiber directions in the laminae, different laminate properties can be acheived. It is often desirable for composites engineers to prescribe these angles such that certain loads can safely be applied to the lamina. This program will do exacly that. By providing properties and a loading scenario, the program will provide the stacking sequence that results in the highest factor fo safety for the given load. The figure below shows an example of interpreting the stacking sequence.

![image](https://github.com/user-attachments/assets/a65dbd49-7613-4395-a171-6599f1754b82)


## Features
This program allows the user to populate a csv file containing material properties, laminate dimensions, applied loads, and calculation constraints. the user can simply issue the program with the number of layers as the first argument and the file path as the second argument. For example: `$ ./composite 4 data.csv` will find the four-layer stacking sequence that provides the largest factor of safety for the parameters given in data.csv. A progress indicator is updated in the terminal for every 1% of program progress. 

## Code Structure
The main program first reads in the provided csv file. Intermediate calculations are made and then data is populated within a struct. The heart of the program lies within a recursive for-loop function. The number of nested for loops is equal to the number of layers in the composite. The for loops iterate through all of the possible layer angle combinations, calculating the safety factor for each configuration. The best safety factor and associated stacking sequence is displayed. To speed up the computations, multiprocessing is employed, using the maximum number of computer cores that satisfy resolution requirements. Mutex locks are used to ensure data is not lost across processes.

## Performance
A similar program was developed in MATLAB side-by-side with this C program to verify correct operation. The purpose of creating the program in C was to decrease the processing time needed to acheive the same result. This was accomplished as shown below. The same input parameters were given to both programs. the MATLAB program took 19.38 seconds whereas this C program took only 0.36 seconds.

![Screenshot 2024-12-05 133721](https://github.com/user-attachments/assets/a4db8b6e-6cbe-46dd-addf-697cbfebf4f4)
![Screenshot 2024-12-05 133848](https://github.com/user-attachments/assets/8e8a58d8-e9fb-4c91-b873-757229a140b7)


## Techniques Used
- File IO
- Multiprocessing
- Shared Memory
- Mutex Locks (not used in any other labs)

