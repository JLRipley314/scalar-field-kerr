all:
	$(MAKE) --directory=Source/

test:
	$(MAKE) --directory=Test/

clean_bin:
	$(MAKE) clean --directory=Bin/

clean_out:
	$(MAKE) clean --directory=Output/
