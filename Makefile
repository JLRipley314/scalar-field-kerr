all:
	$(MAKE) --directory=Source/

test:
	$(MAKE) --directory=Test/

clean_Bin:
	$(MAKE) clean --directory=Bin/

clean_Output:
	$(MAKE) clean --directory=Output/

clean_Source:
	$(MAKE) clean --directory=Source/

clean:
	$(MAKE) clean_Bin
	$(MAKE) clean_Output
	$(MAKE) clean_Source
