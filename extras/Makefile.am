SUBDIRS =
if BUILD_EXTRAS
SUBDIRS += gap localsus top transmag
endif
DIST_SUBDIRS = gap localsus top transmag

check-extras :
	for d in $(DIST_SUBDIRS); do test -d $$d && make -C $$d check; done
