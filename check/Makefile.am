run: run-site run-chain run-triangle

../loop:
	(cd .. && make loop)
../loop_evaluate:
	(cd .. && make loop_evaluate)
../diag:
	(cd .. && make diag)

run-site: run-site-p run-site-e run-site-g ../loop_evaluate
	$(SHELL) $(srcdir)/evaluate-site.sh @bindir@

run-site-p: ../loop
	if test ! -f site-p.out.xml; then \
	  echo 'REPRESENTATION = "path integral"' > site-p; \
	  cat $(srcdir)/site >> site-p; \
	  $(DESTDIR)$(bindir)/parameter2xml site-p; \
	  cp -f site-p.in.xml site-p.out.xml; \
	fi
	../loop --Tmin 1 site-p.out.xml
run-site-e: ../loop
	if test ! -f site-e.out.xml; then \
	  echo 'REPRESENTATION = "SSE"' > site-e; \
	  cat $(srcdir)/site >> site-e; \
	  $(DESTDIR)$(bindir)/parameter2xml site-e; \
	  cp -f site-e.in.xml site-e.out.xml; \
	fi
	../loop --Tmin 1 site-e.out.xml
run-site-g: ../diag
	if test ! -f site-g.out; then \
	  ../diag < $(srcdir)/site > site-g.out; \
	fi

run-chain: run-chain-p run-chain-e run-chain-g ../loop_evaluate
	$(SHELL) $(srcdir)/evaluate-chain.sh @bindir@

run-chain-p: ../loop
	if test ! -f chain-p.out.xml; then \
	  echo 'REPRESENTATION = "path integral"' > chain-p; \
	  cat $(srcdir)/chain >> chain-p; \
	  $(DESTDIR)$(bindir)/parameter2xml chain-p; \
	  cp -f chain-p.in.xml chain-p.out.xml; \
	fi
	../loop --Tmin 1 chain-p.out.xml
run-chain-e: ../loop
	if test ! -f chain-e.out.xml; then \
	  echo 'REPRESENTATION = "SSE"' > chain-e; \
	  cat $(srcdir)/chain >> chain-e; \
	  $(DESTDIR)$(bindir)/parameter2xml chain-e; \
	  cp -f chain-e.in.xml chain-e.out.xml; \
	fi
	../loop --Tmin 1 chain-e.out.xml
run-chain-g: ../diag
	if test ! -f chain-g.out; then \
	  ../diag < $(srcdir)/chain > chain-g.out; \
	fi

run-triangle: run-triangle-p run-triangle-e run-triangle-g ../loop_evaluate
	$(SHELL) $(srcdir)/evaluate-triangle.sh @bindir@

run-triangle-p: ../loop
	if test ! -f triangle-p.out.xml; then \
	  echo 'REPRESENTATION = "path integral"' > triangle-p; \
	  cat $(srcdir)/triangle >> triangle-p; \
	  $(DESTDIR)$(bindir)/parameter2xml triangle-p; \
	  cp -f triangle-p.in.xml triangle-p.out.xml; \
	fi
	../loop --Tmin 1 triangle-p.out.xml
run-triangle-e: ../loop
	if test ! -f triangle-e.out.xml; then \
	  echo 'REPRESENTATION = "SSE"' > triangle-e; \
	  cat $(srcdir)/triangle >> triangle-e; \
	  $(DESTDIR)$(bindir)/parameter2xml triangle-e; \
	  cp -f triangle-e.in.xml triangle-e.out.xml; \
	fi
	../loop --Tmin 1 triangle-e.out.xml
run-triangle-g: ../loop
	if test ! -f triangle-g.out; then \
	  ../diag < $(srcdir)/triangle > triangle-g.out; \
	fi

clean-local :
	rm -f *.xml *.run* *.out *.dat
