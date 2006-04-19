dist_noinst_SCRIPTS = evaluate-chain.sh evaluate-site.sh evaluate-triangle.sh
dist_noinst_DATA = chain site triangle

run: run-site run-chain run-triangle

../loop:
	(cd .. && make loop)
../loop_evaluate:
	(cd .. && make loop_evaluate)

run-site: run-site-p run-site-e run-site-g ../loop_evaluate
	$(SHELL) $(srcdir)/evaluate-site.sh @bindir@

site-d.db: ../loop ../loop_evaluate
	$(SHELL) $(srcdir)/run.sh site-d
site-p.db: ../loop ../loop_evaluate
	$(SHELL) $(srcdir)/run.sh site-p
site-s.db: ../loop ../loop_evaluate
	$(SHELL) $(srcdir)/run.sh site-s

clean-local :
	rm -rf *.run
	rm -f *.xml *.run* *.out *.dat
