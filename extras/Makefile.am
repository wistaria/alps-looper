SUBDIRS =
DIST_SUBDIRS = corrlen gap localsus qwl top transmag xyz

check-extras :
	for d in $(DIST_SUBDIRS); do test -d $$d && make -C $$d check; done
