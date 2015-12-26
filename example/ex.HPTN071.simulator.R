##--------------------------------------------------------------------------------------------------------
##	create command line string to simulate epidemic scenarios under the 
##	regional HPTN071/PopART model
##	code by Mike Pickles and Anne Cori, version 150119
##--------------------------------------------------------------------------------------------------------
cat(cmd.HPTN071.simulator('.', seed=NA, opt.acute='high', opt.intervention='fast'))
