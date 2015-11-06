

Introduction
------------

HiperQAP is a high-performance CUDA library for QAP solvers. Currently
implements three heuristic search algorithm: Tabu Search, 2opt and Simulated
Annealing  


Build Instructions
------------------

To install, in the top-level directory type

   make

set QAPHOME to be the install directory

   export QAPHOME=path_to_install_dir

add ${QAPHOME}/bin to your PATH 

   export PATH=$PATH:${QAPHOME}/bin

add above to your .bashrc, so that you don't have to type the above at every
login



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
