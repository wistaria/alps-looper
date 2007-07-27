dist_noinst_SCRIPTS = run-test.sh
dist_noinst_DATA = plot-energy plot-sheat plot-m plot-m2 plot-sus plot-sm plot-sm2 plot-ssus

.PHONY : run

run: chain dterm site triangle xyz

#

dist_noinst_DATA += chain-1 chain-2 chain-3 chain-4 chain-5
.PHONY : chain
chain:
	(cd $(top_builddir) && $(MAKE) loop_pp)
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp chain-1 plot-energy plot-sheat plot-m2 plot-m4 plot-sus plot-sm2 plot-sm4 plot-ssus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp chain-2 plot-energy plot-sheat plot-m2 plot-m4 plot-sus plot-sm2 plot-sm4 plot-ssus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp chain-3 plot-energy plot-sheat plot-m2 plot-m4 plot-sus plot-sm2 plot-sm4 plot-ssus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp chain-4 plot-energy plot-sheat plot-m2 plot-m4 plot-sus plot-sm2 plot-sm4 plot-ssus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp chain-5 plot-energy plot-sheat plot-m2 plot-m4 plot-sus plot-sm2 plot-sm4 plot-ssus

#

dist_noinst_DATA += dterm-1 dterm-2 dterm-3

.PHONY : dterm

dterm:
	(cd $(top_builddir) && $(MAKE) loop_pp)
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp dterm-1 plot-energy plot-sheat plot-m2 plot-m4 plot-sus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp dterm-2 plot-energy plot-sheat plot-m2 plot-m4 plot-sus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp dterm-3 plot-energy plot-sheat plot-m2 plot-m4 plot-sus

#

dist_noinst_DATA += site-1 site-2 site-3 site-4 site-5
.PHONY : site
site:
	(cd $(top_builddir) && $(MAKE) loop_pp)
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp site-1 plot-energy plot-sheat plot-m plot-sus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp site-2 plot-energy plot-sheat plot-m plot-sus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp site-3 plot-energy plot-sheat plot-m plot-sus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp site-4 plot-energy plot-sheat plot-m plot-sus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp site-5 plot-energy plot-sheat plot-m plot-sus

#

dist_noinst_DATA += triangle-1 triangle-2 triangle-3 triangle-4 triangle-5
.PHONY : triangle
triangle:
	(cd $(top_builddir) && $(MAKE) loop_pp)
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp triangle-1 plot-energy plot-sheat plot-sus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp triangle-2 plot-energy plot-sheat plot-sus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp triangle-3 plot-energy plot-sheat plot-sus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_pp triangle-4 plot-energy plot-sheat plot-sus plot-m2

#

dist_noinst_DATA += xyz-1 xyz-2 xyz-3

.PHONY : xyz

xyz:
	(cd $(top_builddir)/extras/xyz && $(MAKE) xyz_pp)
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/xyz/xyz_pp xyz-1 plot-energy plot-sheat plot-sus plot-ssus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/xyz/xyz_pp xyz-2 plot-energy plot-sheat plot-sus plot-ssus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/xyz/xyz_pp xyz-3 plot-energy plot-sheat plot-sus plot-ssus plot-m2

#

clean-local :
	rm -rf *.run
	rm -f *.dat *.db *.text
