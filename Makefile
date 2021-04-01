all:
	$(MAKE) --directory=Source/

test:
	$(MAKE) --directory=Test/

clean:
	$(MAKE) clean --directory=Bin/
	$(MAKE) clean --directory=Output/
