dist_noinst_SCRIPTS = run-test.sh
dist_noinst_DATA = plot-energy plot-m2 plot-sheat plot-sus

.PHONY : run

run: chain dterm xyz

#

dist_noinst_DATA += chain-1 chain-2 chain-3 chain-4 chain-5
.PHONY : chain
chain:
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp chain-1 plot-energy plot-sheat plot-sus plot-ssus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp chain-2 plot-energy plot-sheat plot-sus plot-ssus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp chain-3 plot-energy plot-sheat plot-sus plot-ssus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp chain-4 plot-energy plot-sheat plot-sus plot-ssus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp chain-5 plot-energy plot-sheat plot-sus plot-ssus plot-m2

#

dist_noinst_DATA += dterm-1 dterm-2 dterm-3

.PHONY : dterm

dterm:
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp dterm-1 plot-energy plot-sheat plot-sus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp dterm-2 plot-energy plot-sheat plot-sus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp dterm-3 plot-energy plot-sheat plot-sus plot-m2

#

dist_noinst_DATA += xyz-1 xyz-2 xyz-3

.PHONY : xyz

xyz:
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/xyz/xyz_pp xyz-1 plot-energy plot-sheat plot-sus plot-ssus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/xyz/xyz_pp xyz-2 plot-energy plot-sheat plot-sus plot-ssus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/xyz/xyz_pp xyz-3 plot-energy plot-sheat plot-sus plot-ssus plot-m2

#

clean-local :
	rm -rf *.run
	rm -f *.dat *.db *.text
