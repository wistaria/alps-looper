dist_noinst_SCRIPTS = evaluate-chain.sh evaluate-site.sh evaluate-triangle.sh
dist_noinst_DATA = chain site triangle

run: run-site run-chain run-triangle

run-site: site-d.db site-p.db site-s.db
site-d.db:
	$(SHELL) $(srcdir)/run.sh site-d
site-p.db:
	$(SHELL) $(srcdir)/run.sh site-p
site-s.db:
	$(SHELL) $(srcdir)/run.sh site-s

clean-local :
	rm -rf *.run
	rm -f *.xml *.run* *.out *.dat
