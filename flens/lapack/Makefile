all :
	$(MAKE) -C interface

dev :
	$(MAKE) -C interface dev

test/ :
	tar xfz test.tgz

dev-check : test/
	$(MAKE) -C test/MATGEN
	$(MAKE) -C test cleanup
	$(MAKE) -C test zsg.out
	$(MAKE) -C test dtest.out
	$(MAKE) -C test ztest.out
	$(MAKE) -C test dnep.out
	$(MAKE) -C test znep.out

clean : test/
	$(MAKE) -C interface clean

distclean :
	$(MAKE) -C interface clean
	$(RM) -rf test/
