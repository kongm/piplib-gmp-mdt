#############################################################################
SUBDIRS 		= doc source include test example

#############################################################################
ACLOCAL_AMFLAGS		= -I config

m4datadir		= $(datadir)/aclocal

AUX_DIST                =			\
	$(ac_aux_dir)/config.guess		\
	$(ac_aux_dir)/config.sub		\
	$(ac_aux_dir)/install-sh		\
	$(ac_aux_dir)/ltmain.sh			\
	$(ac_aux_dir)/missing			\
	$(ac_aux_dir)/depcomp


#EXTRA_DIST		= 			\
#	bootstrap.sh


MAINTAINERCLEANFILES 	=			\
	Makefile.in				\
	aclocal.m4				\
	configure				\
	source/stamp-h.in			\
	$(AUX_DIST)

dist-hook:
	(cd $(distdir) && mkdir -p $(ac_aux_dir))
	for file in $(AUX_DIST); do \
	  cp $$file $(distdir)/$$file; \
	done
#############################################################################
