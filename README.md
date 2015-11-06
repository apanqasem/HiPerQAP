

Introduction
------------

HiperQAP is a high-performance CUDA library for QAP solvers. Currently
implements three heuristic search algorithm: Tabu Search, 2opt and Simulated
Annealing  


Build Instructions
------------------

To install, in the top-level directory type

   make

add the bin directory to your PATH if you want to perfom experiments



Directories
-----------

src/
	source code for search algorithms, I/O routines etc. 

utils/
	scripts to run QAP instances, check for correctness, accuracy.

datasets/
	Taillard and Lipa datasets from http://anjos.mgi.polymtl.ca/qaplib 

bin/
	location of scripts after build  
