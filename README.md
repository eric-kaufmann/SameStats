# SameStats

As in the Paper [*Same Stats, Different Graphs: Generating Datasets with Varied Appearance and Identical Statistics through Simulated Annealing*](https://www.autodesk.com/research/publications/same-stats-different-graphs) by Justin Matejka and George Fitzmaurice introduced, there are datasets with the same statistical properties but dissimilar graphs. The current approach is still very slow because it is written in Python. In this project we wan't to optimize and parallelize the code using C++ and OpenMP.

# Requirements:
1. Python3
2. OpenMP library
3. g++ Compiler

# Building project
1. create directory e.g.: mkdir build
2. copy scatter_plot.py into build: cp scatter_plot.py build
3. go to directory: cd build 
4. execute: cmake -DCMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=g++ ..
        or: cmake -DCMAKE_BUILD_TYPE=Debug -D CMAKE_CXX_COMPILER=g++ ..
        or: cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -D CMAKE_CXX_COMPILER=g++ ..
5. execute: cmake --build .
6. run executable: ./test

# Output
This Version creates png files which display the plot of our pointcloud before the transformation and after. 
Plot before transformation: input_<shape>_scatter_plot.png
                     after: output_<shape>_scatter_plot.png

# Speed up
Compared to the python version of the above mentioned autodesk paper, this version is much faster.
When running with 10^6 iterations, the python version needs more than 10 minutes whereas this version is terminating after a few seconds on the same machine. 

