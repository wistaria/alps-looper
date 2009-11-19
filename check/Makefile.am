dist_noinst_SCRIPTS = run-test.sh
dist_noinst_DATA = plot-energy plot-sheat plot-m plot-m2 plot-sus plot-sm plot-sm2 plot-ssus plot-tm plot-sign

.PHONY : run

run: chain dterm kagome site transmag triangle xyz

#

dist_noinst_DATA += chain-1 chain-2 chain-3 chain-4 chain-5 chain-6
.PHONY : chain
chain:
	(cd $(top_builddir) && $(MAKE) loop_mpi)
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi chain-1 plot-energy plot-sheat plot-m2 plot-m4 plot-sus plot-sm2 plot-sm4 plot-ssus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi chain-2 plot-energy plot-sheat plot-m2 plot-m4 plot-sus plot-sm2 plot-sm4 plot-ssus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi chain-3 plot-energy plot-sheat plot-m2 plot-m4 plot-sus plot-sm2 plot-sm4 plot-ssus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi chain-4 plot-energy plot-sheat plot-m2 plot-m4 plot-sus plot-sm2 plot-sm4 plot-ssus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi chain-5 plot-energy plot-sheat plot-m2 plot-m4 plot-sus plot-sm2 plot-sm4 plot-ssus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi chain-6 plot-energy plot-sheat plot-m2 plot-m4 plot-sus plot-sm2 plot-sm4 plot-ssus

#

dist_noinst_DATA += dterm-1 dterm-2 dterm-3

.PHONY : dterm

dterm:
	(cd $(top_builddir) && $(MAKE) loop_mpi)
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi dterm-1 plot-energy plot-sheat plot-m2 plot-m4 plot-sus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi dterm-2 plot-energy plot-sheat plot-m2 plot-m4 plot-sus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi dterm-3 plot-energy plot-sheat plot-m2 plot-m4 plot-sus

#

dist_noinst_DATA += kagome-1 kagome-2 kagome-3 kagome-4
.PHONY : kagome
kagome:
	(cd $(top_builddir) && $(MAKE) loop_mpi)
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi kagome-1 plot-energy plot-sheat plot-sus plot-m2 plot-sign plot-zeromeron plot-signmeron
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi kagome-2 plot-energy plot-sheat plot-sus plot-m2 plot-sign plot-zeromeron plot-signmeron
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi kagome-3 plot-energy plot-sheat plot-sus plot-m2 plot-sign plot-zeromeron plot-signmeron
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi kagome-4 plot-energy plot-sheat plot-sus plot-m2 plot-sign plot-zeromeron plot-signmeron

#

dist_noinst_DATA += site-1 site-2 site-3 site-4 site-5
.PHONY : site
site:
	(cd $(top_builddir) && $(MAKE) loop_mpi)
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi site-1 plot-energy plot-sheat plot-m plot-sus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi site-2 plot-energy plot-sheat plot-m plot-sus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi site-3 plot-energy plot-sheat plot-m plot-sus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi site-4 plot-energy plot-sheat plot-m plot-sus
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi site-5 plot-energy plot-sheat plot-m plot-sus

#

dist_noinst_DATA += transmag-1 transmag-2 transmag-3 transmag-4
.PHONY : transmag
transmag:
	(cd $(top_builddir)/extras/transmag && $(MAKE) transmag_mpi)
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/transmag/transmag_mpi transmag-1 plot-energy plot-sheat plot-tm
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/transmag/transmag_mpi transmag-2 plot-energy plot-sheat plot-m
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/transmag/transmag_mpi transmag-3 plot-energy plot-sheat plot-tm
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/transmag/transmag_mpi transmag-4 plot-energy plot-sheat plot-sm

#

dist_noinst_DATA += triangle-1 triangle-2 triangle-3 triangle-4 triangle-5 triangle-6
.PHONY : triangle
triangle:
	(cd $(top_builddir) && $(MAKE) loop_mpi)
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi triangle-1 plot-energy plot-sheat plot-sus plot-m2 plot-sign plot-zeromeron plot-signmeron
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi triangle-2 plot-energy plot-sheat plot-sus plot-m2 plot-sign plot-zeromeron plot-signmeron
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi triangle-3 plot-energy plot-sheat plot-sus plot-m2 plot-sign plot-zeromeron plot-signmeron
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi triangle-4 plot-energy plot-sheat plot-sus plot-m2 plot-sign plot-zeromeron plot-signmeron
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi triangle-5 plot-energy plot-sheat plot-sus plot-m2 plot-sign
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/loop_mpi triangle-6 plot-energy plot-sheat plot-sus plot-m2 plot-sign

#

dist_noinst_DATA += xyz-1 xyz-2 xyz-3 xyz-4 xyz-5 xyz-6 xyz-7 xyz-8 xyz-9 xyz-10

.PHONY : xyz

xyz:
	(cd $(top_builddir)/extras/xyz && $(MAKE) xyz_mpi)
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/xyz/xyz_mpi xyz-1 plot-energy plot-sheat plot-sus plot-ssus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/xyz/xyz_mpi xyz-2 plot-energy plot-sheat plot-sus plot-ssus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/xyz/xyz_mpi xyz-3 plot-energy plot-sheat plot-sus plot-ssus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/xyz/xyz_mpi xyz-4 plot-energy plot-sheat plot-sus plot-ssus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/xyz/xyz_mpi xyz-5 plot-energy plot-sheat plot-sus plot-ssus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/xyz/xyz_mpi xyz-6 plot-energy plot-sheat plot-sus plot-ssus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/xyz/xyz_mpi xyz-7 plot-energy plot-sheat plot-sus plot-ssus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/xyz/xyz_mpi xyz-8 plot-energy plot-sheat plot-sus plot-ssus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/xyz/xyz_mpi xyz-9 plot-energy plot-sheat plot-sus plot-ssus plot-m2
	$(SHELL) $(srcdir)/run-test.sh $(top_builddir)/extras/xyz/xyz_mpi xyz-10 plot-energy plot-sheat plot-sus plot-ssus plot-m2

#

clean-local :
	rm -rf *.run
	rm -f *.dat *.db *.text
