SUBDIRS =
if BUILD_EXTRAS
SUBDIRS += corrlen gap localsus top transmag xyz
endif
DIST_SUBDIRS = corrlen gap localsus top transmag xyz

check-extras :
	for d in $(DIST_SUBDIRS); do test -d $$d && make -C $$d check; done
