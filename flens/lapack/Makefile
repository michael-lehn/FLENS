all :
	$(MAKE) -C interface

test/ :
	tar xfz test.tgz

check : test/
	$(MAKE) -C test/MATGEN
	$(MAKE) -C test cleanup
	$(MAKE) -C test dtest.out

clean : test/
	$(MAKE) -C interface clean

distclean :
	$(RM) -rf test/
