run: run-chain run-triangle

run-chain: run-chain-p run-chain-e run-chain-g ../loop_evaluate
	$(SHELL) $(srcdir)/evaluate-chain.sh @bindir@

run-chain-p: ../loop
	if test ! -f chain-p.out.xml; then \
	  $(DESTDIR)$(bindir)/parameter2xml $(srcdir)/chain-p chain-p; \
	  cp -f chain-p.in.xml chain-p.out.xml; \
	fi
	../loop --Tmin 1 chain-p.out.xml
run-chain-e: ../loop
	if test ! -f chain-e.out.xml; then \
	  $(DESTDIR)$(bindir)/parameter2xml $(srcdir)/chain-e chain-e; \
	  cp -f chain-e.in.xml chain-e.out.xml; \
	fi
	../loop --Tmin 1 chain-e.out.xml
run-chain-g: ../ed/diag
	if test ! -f chain-g.out; then \
	  ../ed/diag < $(srcdir)/chain-p > chain-g.out; \
	fi

run-triangle: run-triangle-p run-triangle-e run-triangle-g ../loop_evaluate
	$(SHELL) $(srcdir)/evaluate-triangle.sh @bindir@

run-triangle-p: ../loop
	if test ! -f triangle-p.out.xml; then \
	  $(DESTDIR)$(bindir)/parameter2xml $(srcdir)/triangle-p triangle-p; \
	  cp -f triangle-p.in.xml triangle-p.out.xml; \
	fi
	../loop --Tmin 1 triangle-p.out.xml
run-triangle-e: ../loop
	if test ! -f triangle-e.out.xml; then \
	  $(DESTDIR)$(bindir)/parameter2xml $(srcdir)/triangle-e triangle-e; \
	  cp -f triangle-e.in.xml triangle-e.out.xml; \
	fi
	../loop --Tmin 1 triangle-e.out.xml
run-triangle-g: ../loop
	if test ! -f triangle-g.out; then \
	  ../ed/diag < $(srcdir)/triangle-p > triangle-g.out; \
	fi

clean-local :
	rm -f *.xml *.run* *.out *.dat
