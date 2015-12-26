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
	tar -zxvf ../inst/ext/Seq-Gen.v1.3.3.tgz
	make -C Seq-Gen.v1.3.3/source 
	mv -f $(CWD)/Seq-Gen.v1.3.3/source/seq-gen ../inst/ext/seq-gen
	rm -rf Seq-Gen.v1.3.3	

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