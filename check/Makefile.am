run-chain-p:
	$(DESTDIR)$(bindir)/parameter2xml chain-p
	../loop -Tmin 2 chain-p.in.xml
run-chain-e:
	$(DESTDIR)$(bindir)/parameter2xml chain-e
	../loop -Tmin 2 chain-e.in.xml
run-chain-g:
	../ed/diag < chain-p > chain-g.out

run-triangle-p:
	$(DESTDIR)$(bindir)/parameter2xml triangle-p
	../loop -Tmin 2 triangle-p.in.xml
run-triangle-e:
	$(DESTDIR)$(bindir)/parameter2xml triangle-e
	../loop -Tmin 2 triangle-e.in.xml
run-triangle-g:
	../ed/diag < triangle-p > triangle-g.out

clean-local :
	rm -f *.xml *.run* plot*.dat *.out
