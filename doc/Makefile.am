install_sh_HTML = $(install_sh_DATA)
pkgdocdir = @prefix@/doc/@PACKAGE_NAME@
dist_noinst_DATA = $(HTML)

HTML = index.html

install: install-html
install-html:
	@list='$(HTML)'; for p in $$list; do \
	  if test -d $(pkgdocdir); then :; else \
	    echo "$(mkinstalldirs) $(pkgdocdir)"; \
	    $(mkinstalldirs) $(pkgdocdir); \
	  fi; \
	  sd=`dirname $$p`; \
	  if test "$$sd" != "."; then \
	    if test -d "$(pkgdocdir)/$$sd"; then :; else \
	      echo "$(mkinstalldirs) $(pkgdocdir)/$$sd"; \
	      $(mkinstalldirs) $(pkgdocdir)/$$sd; \
	    fi; \
	  fi; \
	  if test -f "$$p"; then d=; else d="$(srcdir)/"; fi; \
	  echo "$(install_sh_HTML) $$d$$p $(pkgdocdir)/$$p"; \
	  $(install_sh_HTML) $$d$$p $(pkgdocdir)/$$p; \
	done
uninstall: uninstall-html
uninstall-html:
	@list='$(HTML)'; for p in $$list; do \
	  echo "rm -f $(pkgdocdir)/$$p"; \
	  rm -f $(pkgdocdir)/$$p; \
	done
