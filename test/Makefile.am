#############################################################################
SUBDIRS 		=

#############################################################################

EXTRA_DIST                =			\
	boulet.dat				\
	boulet.ll				\
	bouleti.dat				\
	bouleti.ll				\
	challenges				\
	challenges/pipFile_0			\
	challenges/pipFile_1			\
	challenges/system_sysmo_for_pipMP	\
	challenges/vivien32.dat			\
	crescat.dat				\
	crescat.ll				\
	dirk.dat				\
	dirk.ll					\
	discr.dat				\
	discr.ll				\
	equus.dat				\
	equus.ll				\
	invert.dat				\
	invert.ll				\
	linear.dat				\
	linear.ll				\
	lineri.dat				\
	lineri.ll				\
	loz.dat					\
	loz.ll					\
	Makefile				\
	Makefile.am				\
	max.dat					\
	max.ll					\
	maxb.dat				\
	maxb.ll					\
	pairi.dat				\
	pairi.ll				\
	petit.dat				\
	petit.ll				\
	rairo.dat				\
	rairo.ll				\
	rairoi.dat				\
	rairoi.ll				\
	test.dat				\
	test.ll					\
	test10i.dat				\
	test10i.ll				\
	test11.dat				\
	test11.ll				\
	test11i.dat				\
	test11i.ll				\
	test12.dat				\
	test12.ll				\
	test12i.dat				\
	test12i.ll				\
	test2.dat				\
	test2.ll				\
	test2i.dat				\
	test2i.ll				\
	test3.dat				\
	test3.ll				\
	test3i.dat				\
	test3i.ll				\
	test4.dat				\
	test4.ll				\
	test4i.dat				\
	test4i.ll				\
	test5.dat				\
	test5.ll				\
	test5i.dat				\
	test5i.ll				\
	test6.dat				\
	test6.ll				\
	test6i.dat				\
	test6i.ll				\
	test7.dat				\
	test7.ll				\
	test7i.dat				\
	test7i.ll				\
	test8i.dat				\
	test8i.ll				\
	test9i.dat				\
	test9i.ll

MAINTAINERCLEANFILES 	=			\
	Makefile.in

#############################################################################

PIPTEST =       \
	crescat \
	discr   \
	equus   \
	invert  \
	linear  \
	lineri  \
	loz     \
	max     \
	maxb    \
	pairi   \
	petit   \
	rairo   \
	rairoi  \
	test    \
	test2   \
	test2i  \
	test3   \
	test3i  \
	test4   \
	test4i  \
	test5   \
	test5i  \
	test6   \
	test6i  \
	test7   \
	test7i  \
	test8i  \
	test9i  \
	test10i \
	test11  \
	test11i \
	test12  \
	test12i

BITS = 64

generate:
	@echo "             /*-----------------------------------------------*"
	@echo "              *                 Generate files                *"
	@echo "              *-----------------------------------------------*/"
	for x in $(PIPTEST) ; do \
		echo "Generate $$x.ll" ; \
		$(top_builddir)/source/pip$(BITS) $$x.dat $$x.ll ; \
	done

check:
	@echo "             /*-----------------------------------------------*"
	@echo "              *                  Testing PIP                  *"
	@echo "              *-----------------------------------------------*/"
	failedtest=0; \
	for x in $(PIPTEST) ; do \
		echo "Verify file $$x" ; \
		$(top_builddir)/source/pip$(BITS) $$x.dat > xyz ;\
		diff -w xyz $$x.ll ; \
		result=$$?; \
		if [ "$$result" -eq "1" ]; then \
			echo "Error: $$x is not the same"; \
			failedtest=`expr $$failedtest + 1`; \
		else \
			echo "$$x passed"; \
	        fi; \
	done ; \
	rm xyz ; \
	if [ $$failedtest != 0 ]; then \
		echo "$$failedtest tests failed"; \
	else \
	        echo "PIP works correctly :-)"; \
        fi

#############################################################################
