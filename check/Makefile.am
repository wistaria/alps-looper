run: run-chain run-triangle

run-chain: run-chain-p run-chain-e run-chain-g
	$(SHELL) $(srcdir)/evaluate-chain.sh @bindir@

run-chain-p:
	if test ! -f chain-p.out.xml; then \
	  $(DESTDIR)$(bindir)/parameter2xml chain-p; \
	  cp -f chain-p.in.xml chain-p.out.xml; \
	fi
	../loop --Tmin 2 chain-p.out.xml
run-chain-e:
	if test ! -f chain-e.out.xml; then \
	  $(DESTDIR)$(bindir)/parameter2xml chain-e; \
	  cp -f chain-e.in.xml chain-e.out.xml; \
	fi
	../loop --Tmin 2 chain-e.out.xml
run-chain-g:
	if test ! -f chain-g.out; then \
	  ../ed/diag < chain-p > chain-g.out; \
	fi

run-triangle: run-triangle-p run-triangle-e run-triangle-g
	$(SHELL) $(srcdir)/evaluate-triangle.sh @bindir@

run-triangle-p:
	if test ! -f triangle-p.out.xml; then \
	  $(DESTDIR)$(bindir)/parameter2xml triangle-p; \
	  cp -f triangle-p.in.xml triangle-p.out.xml; \
	fi
	../loop --Tmin 2 triangle-p.out.xml
run-triangle-e:
	if test ! -f triangle-e.out.xml; then \
	  $(DESTDIR)$(bindir)/parameter2xml triangle-e; \
	  cp -f triangle-e.in.xml triangle-e.out.xml; \
	fi
	../loop --Tmin 2 triangle-e.out.xml
run-triangle-g:
	if test ! -f triangle-g.out; then \
	  ../ed/diag < triangle-p > triangle-g.out; \
	fi

clean-local :
	rm -f *.xml *.run* *.out *.dat
