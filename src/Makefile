CWD= $$(pwd)

all:	
	echo "current directory is "$(CWD)	
	@echo "*********************************************** "
	@echo "********* BUILDING HPTN071 sim model ********** "
	@echo "*********************************************** "	
	make hptn	
	@echo "*********************************************** "
	@echo "********* BUILDING seq-gen ******************** "
	@echo "*********************************************** "	
	make seq-gen	

hptn:		
	echo $(CWD)		
	curl "https://raw.githubusercontent.com/mrehp2/POPART-IBM/PangeaSim/PANGEA_SIM.tar.gz" -o PANGEA_SIM.tar.gz
	tar -zxvf PANGEA_SIM.tar.gz
	make -C PANGEA_SIM/SCENARIO_HIGHACUTE
	mv -f $(CWD)/PANGEA_SIM/SCENARIO_HIGHACUTE/popart-simul.exe ../inst/ext/popart-highacute
	cp -fR $(CWD)/PANGEA_SIM/SCENARIO_HIGHACUTE/PangeaParams ../inst/ext/PangeaParamsHighAcute	
	make -C PANGEA_SIM/SCENARIO_LOWACUTE 
	mv -f $(CWD)/PANGEA_SIM/SCENARIO_LOWACUTE/popart-simul.exe ../inst/ext/popart-lowacute
	cp -fR $(CWD)/PANGEA_SIM/SCENARIO_LOWACUTE/PangeaParams ../inst/ext/PangeaParamsLowAcute
	rm -rf PANGEA_SIM
		
seq-gen:
	echo $(CWD)	
	tar -zxvf ../inst/ext/seq-gen_v161121.tar.gz
	make -C Seq-Gen/source 
	mv -f $(CWD)/Seq-Gen/source/seq-gen ../inst/ext/seq-gen
	rm -rf Seq-Gen	
