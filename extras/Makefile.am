SUBDIRS =
DIST_SUBDIRS = corrlen gap localsus parallel simple top transmag xyz

check-extras :
	for d in $(DIST_SUBDIRS); do test -d $$d && make -C $$d check; done
