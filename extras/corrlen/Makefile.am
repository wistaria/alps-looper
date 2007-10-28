#
# compiler options
#

DEFAULT_INCLUDES = -I$(srcdir) -I$(top_srcdir) -DLOOP_CONFIG_HEADER=\"corrlen_config.h\"

#
# programs
#

noinst_PROGRAMS = corrlen corrlen_evaluate
if HAVE_MPI
noinst_PROGRAMS += corrlen_mpi
endif
if HAVE_PARAPACK
noinst_PROGRAMS += corrlen_ps corrlen_pevaluate
if HAVE_MPI
noinst_PROGRAMS += corrlen_pp
endif
endif

SRC = $(top_srcdir)/loop_factory.C $(top_srcdir)/loop_model.C $(top_srcdir)/path_integral.C

  dist_corrlen_SOURCES = corrlen_config.h
nodist_corrlen_SOURCES = $(top_srcdir)/loop.C $(SRC)
nodist_corrlen_mpi_SOURCES = $(top_srcdir)/loop.C $(SRC)
       corrlen_mpi_LDFLAGS = @LDFLAGS_MPI@
       corrlen_mpi_LDADD = @LIBS_MPI@
nodist_corrlen_evaluate_SOURCES = $(top_srcdir)/loop_evaluate.C $(SRC)

nodist_corrlen_ps_SOURCES = $(top_srcdir)/loop_pp.C $(SRC)
       corrlen_ps_LDFLAGS = @LDFLAGS_PS@
       corrlen_ps_LDADD = @LIBS_PS@
nodist_corrlen_pp_SOURCES = $(top_srcdir)/loop_pp.C $(SRC)
       corrlen_pp_LDFLAGS = @LDFLAGS_PP@
       corrlen_pp_LDADD = @LIBS_PP@
nodist_corrlen_pevaluate_SOURCES = $(top_srcdir)/loop_pevaluate.C $(SRC)
       corrlen_pevaluate_LDFLAGS = @LDFLAGS_PP@
       corrlen_pevaluate_LDADD = @LIBS_PP@


#
# tests
#

dist_noinst_DATA = corrlen_ns.ip corrlen_ns.op

check_PROGRAMS = corrlen_ns

nodist_corrlen_ns_SOURCES =  $(top_srcdir)/loop_ns.C $(SRC)

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests
test tests : check
