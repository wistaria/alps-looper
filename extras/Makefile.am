SUBDIRS =
if BUILD_EXTRAS
SUBDIRS += gap localsus transmag
endif
DIST_SUBDIRS = gap localsus transmag

check-extras :
	for d in $(DIST_SUBDIRS); do test -d $$d && make -C $$d check; done
