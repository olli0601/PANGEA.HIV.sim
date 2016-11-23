CWD= $$(pwd)

all:	
	echo "current directory is "$(CWD)	
	@echo "*********************************************** "
	@echo "********* BUILDING seq-gen ******************** "
	@echo "*********************************************** "	
	make seq-gen
	@echo "*********************************************** "
	@echo "********* BUILDING HPTN071/HighAcute ********** "
	@echo "*********************************************** "	
	make hptn-high
	@echo "*********************************************** "
	@echo "********* BUILDING HPTN071/LowAcute *********** "
	@echo "*********************************************** "	
	make hptn-low	
	
		
seq-gen:
	echo $(CWD)	
	tar -zxvf ../inst/ext/seq-gen_v161121.tar.gz
	make -C Seq-Gen/source 
	mv -f $(CWD)/Seq-Gen/source/seq-gen ../inst/ext/seq-gen
	rm -rf Seq-Gen	

hptn-high:
	echo $(CWD)	
	tar -zxvf ../inst/ext/SCENARIO_HIGHACUTE_150129.tar.gz
	make -C SCENARIO_HIGHACUTE 
	mv -f $(CWD)/SCENARIO_HIGHACUTE/popart-simul.exe ../inst/ext/popart-highacute
	cp -fR $(CWD)/SCENARIO_HIGHACUTE/PangeaParams ../inst/ext/PangeaParamsHighAcute
	rm -rf SCENARIO_HIGHACUTE
	
hptn-low:
	echo $(CWD)	
	tar -zxvf ../inst/ext/SCENARIO_LOWACUTE_150129.tar.gz
	make -C SCENARIO_LOWACUTE 
	mv -f $(CWD)/SCENARIO_LOWACUTE/popart-simul.exe ../inst/ext/popart-lowacute
	cp -fR $(CWD)/SCENARIO_LOWACUTE/PangeaParams ../inst/ext/PangeaParamsLowAcute
	rm -rf SCENARIO_LOWACUTE	