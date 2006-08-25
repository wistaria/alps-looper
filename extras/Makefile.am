SUBDIRS =
if BUILD_EXTRAS
SUBDIRS += gap localsus transmag
endif
DIST_SUBDIRS = gap localsus transmag

dist_noinst_HEADERS = \
	looper/gap_measurement.h \
	looper/localsus_measurement.h \
	looper/transmag_measurement.h

check-extras :
	for d in $(DIST_SUBDIRS); do \
	  if test -d $$d; then \
	    make -C $$d check; \
	  fi; \
	done
