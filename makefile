all: src utils 

.PHONY: src
.PHONY: utils

src:
	cd src && make

utils:	
	cd utils && make
	/bin/cp utils/costcheck bin/costcheck

	/bin/ln -fs  ../utils/qap.sh bin/qap.sh 
	/bin/ln -fs  ../utils/verify_sol.sh  bin/verify_sol.sh 
	/bin/ln -fs  ../utils/calc_gap.sh  bin/calc_gap.sh 
	/bin/ln -fs  ../utils/exp_driver.sh bin/exp_driver.sh




