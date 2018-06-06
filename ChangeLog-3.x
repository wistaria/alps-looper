2009-06-23  Synge Todo <wistaria@comp-phys.org>
	
	* Fixed improved estimators for |Magnetization|, |Magnetization Density|,
	  |Staggered Magnetization|, |Staggered Magnetization Density| in
	  looper/susceptibility.h

2009-06-19  Synge Todo <wistaria@comp-phys.org>

	* Reduced memory usage of estimators
	
2009-06-17  Synge Todo <wistaria@comp-phys.org>

	* [pro] fixed correlation length estimator
	
2009-03-30  Synge Todo <wistaria@comp-phys.org>
	
	* version 3.2b9 for ALPS Applications 1.3

2009-03-30  Synge Todo <wistaria@comp-phys.org>

	* Fixed estimator initialization
	
2009-03-24  Synge Todo <wistaria@comp-phys.org>

	* version 3.2b7 for ALPS Applications 1.3

2009-03-24  Synge Todo <wistaria@comp-phys.org>
	
	* Removed old alps scheduler support

2008-11-08  Synge Todo <wistaria@comp-phys.org>

	* Adapted for parapack scheduler inclusion into ALPS library
	  (To-Do: remove dependency on conventional ALPS scheduler)

	* Removed looper/integer_range.h (use alps/parapack/integer_range.h instead)

2008-10-08  Synge Todo <wistaria@comp-phys.org>

	* Update for alps-parapack new factory class

2008-01-12  Synge Todo <wistaria@comp-phys.org>
	
	* [pro] Added stand-alone reference implementations (serial and
	parallel versions) in 'standalone' directory.
	
2008-01-08  Synge Todo <wistaria@comp-phys.org>

	* version 3.2b5 for ALPS Applications 1.3

2007-10-28  Synge Todo <wistaria@comp-phys.org>
	
	* [pro] Added correlation length estimator.
	
2007-10-09  Synge Todo <wistaria@comp-phys.org>

	* Introduced SCATTERING_RATIO parameter to xyz_bond_weight_helper
	class
	
2007-10-03  Synge Todo <wistaria@comp-phys.org>

	* Fixed unimproved estimator of correlation functions for S>1/2 in
	looper/correlation.h
	
2007-07-31  Synge Todo <wistaria@comp-phys.org>

	* Added support for XYZ interaction.

2007-07-30  Synge Todo <wistaria@comp-phys.org>

	* Fixed (improved) negative sign estimator.

2007-07-01  Synge Todo <wistaria@comp-phys.org>

	* Fixed a bug in calculating weights for longitudinal external field.

2007-05-30  Synge Todo <wistaria@comp-phys.org>

	* Added linear_distance_helper in looper/lattice.h for calculating
	the 'real' linear distance between the vertices

2007-03-23  Synge Todo <wistaria@comp-phys.org>

	* Added pre_evaluator.
	
2007-03-12  Synge Todo <wistaria@comp-phys.org>
	
	* [pro] Added string order parameter.
	
2007-02-24  Synge Todo <wistaria@comp-phys.org>

	* version 3.2b3

	- prerequisites: boost-1.31.1, alps-1.3b3, lapack

	- Implemented all the feasures planed for alps-application release 1.3

	- Support for ALPS/parapack
	
	- ToDo: fixed for negative signs, external field for SSE, QWL,
	more tests, documentation
	
2007-02-24  Synge Todo <wistaria@comp-phys.org>
	
	* Added custom measurements (average, local, correlations,
	structure factor)
	
2006-12-15  Synge Todo <wistaria@comp-phys.org>

	* Added support for ALPS/parapack.
	
2006-01-26  Synge Todo <wistaria@comp-phys.org>
	
	* Added extras subdirectory
	
2006-01-23  Synge Todo <wistaria@comp-phys.org>

	* Fixed evaluation of non-linear quantities.
	
2006-01-20  Synge Todo <wistaria@comp-phys.org>

        * version 3.2b1

	- prerequisites: boost-1.31.1, alps snapshot (2006/01/20), lapack

	- new features: new internal tree-like data structure for
	clusters, two-pass cluster update, support for longitudinal field,
	transverse field, single-ion anisotropy (D-term), variable length SSE
	
	- ToDo: fixed for negative signs, external field for SSE, QWL, add
	measurement of stiffness, more tests, documentation
	
2006-01-06  Synge Todo <wistaria@comp-phys.org>

	* First working version of new alps-looper (candidate of version 3.2)

	* Old version is moved into v3.1 subdirectory
	
2005-02-23  Synge Todo <wistaria@comp-phys.org>

	* Added support for disordered lattices
	
2004-12-10  Synge Todo <wistaria@comp-phys.org>

	* Added support for MKL 7.2
	
2004-11-25  Synge Todo <wistaria@comp-phys.org>

	* Adopted boost::program_options
	
2004-11-07  Synge Todo <wistaria@comp-phys.org>

	* Moved all the operator definitions into site basis (for ALPS 1.3).
	
2004-10-08  Synge Todo <wistaria@comp-phys.org>

        * version 3.1.2

	* Rearranged source files to reduce memory size for build.

	* Added master header file 'looper.h'.

	* Removed virtual_graph class template and virtual_graph.h.
	
2004-10-06  Synge Todo <wistaria@comp-phys.org>

        * version 3.1.1

	* Fixed improved estimators for models with negative signs.

	* Fixed (or improved) default labeling probability for models on
	non-bipartite lattices.
	
	* Added diagonal energy measurement.
	
2004-09-23  Synge Todo <wistaria@comp-phys.org>

        * version 3.1

	- prerequisites: boost-1.31, alps-1.2, lapack

        - TODO: optimization, more measurements, bug fix for improved
	estimator with negative sign problem

2004-09-21  Synge Todo <wistaria@comp-phys.org>

        * version 3.1b3

	* Force "ergodic" loop solutions also for non-bipartite lattices.

2004-09-17  Synge Todo <wistaria@comp-phys.org>

	* Removed unused parameter from generate_loops.
	
2004-09-14  Synge Todo <wistaria@comp-phys.org>

        * version 3.1b2

	* Fixed for non-bipartite lattices.

	* Added (hidden) simulation parameters STRICT_MCS and FIXED_SEED
	  for debugging.

	* Fixed sign calculation.

2004-09-13  Synge Todo <wistaria@comp-phys.org>

        * version 3.1b1

	* Fixed energy offset for SSE with "ergodic" loop solutions.

	* Added negative sign support.
	
2004-08-25  Synge Todo <wistaria@comp-phys.org>

	* Fixed xxz.h for systems with different spin sizes.

	* Added support for extenal field in site_parameter.

	* Added an extra parameter to the default weight classes for
	supporting "erdogic" loop solutions.
	
2004-08-20  Synge Todo <wistaria@comp-phys.org>

	* Changed solve_llsp so that the SVD is used instead of QL
	factorization.  Added test programs gesvd and llsp.  Removed
	the test program gels.
	
2004-08-19  Synge Todo <wistaria@comp-phys.org>

	* Added support for MKL 7.0
	
2004-07-22  Synge Todo  <wistaria@comp-phys.org>

        * version 3.0.1 released.

2004-07-09  Synge Todo  <wistaria@comp-phys.org>

	* Fixed random_choice class template.
	
2004-06-25  Synge Todo  <wistaria@comp-phys.org>

	* Added alternating_tensor (Levi-Civita or permutation symbol)
	function in looper/util.h.
	
2004-06-11  Synge Todo  <wistaria@comp-phys.org>

        * version 3.0 released.

2004-06-08  Synge Todo  <wistaria@comp-phys.org>
	
	* Turned off maintainer mode for GNU autotools.
	
2004-06-03  Synge Todo  <wistaria@comp-phys.org>

	* Added new measurement staggered_generalized_susceptibility_imp.
	Renamed generalized_susceptibility_imp as
	uniform_generalized_susceptibility_imp.

2004-05-17  Synge Todo  <wistaria@comp-phys.org>
	
	* Updated automake version to 1.8.5.

2004-05-14  Synge Todo  <wistaria@comp-phys.org>

	* Changed excecutable names (qmc* -> loop*).
	
2004-04-21  Synge Todo  <wistaria@comp-phys.org>

        * Fixed a bug in static_correlation in looper/path_integral.h

2004-04-20  Synge Todo  <wistaria@comp-phys.org>

        * Optimized measurement of staggered susceptibility in SSE.

2004-04-14  Synge Todo  <wistaria@comp-phys.org>

        * Added wrapper and test program for linear least squares solver
        (gels) in LAPACK.

        * Replaced internal matrix representation in xxz_matrix by
        boost::multiarray<>.  Improved fitting algorithm in fit2xxz by
        using the linear least squares solver.
        
2004-04-07  Synge Todo  <wistaria@comp-phys.org>

        * version 3.0b5 released.
        - prerequisites: boost_1_31_0, alps-20040407
        - TODO: optimization, improved estimators, more measurements,
          license issue.  Bug fixes for some mesurements...

2004-04-06  Synge Todo  <wistaria@comp-phys.org>
        
        * Added exact diagonalization checking program.
        
        * Moved percolation programs into an independent package.
        
2004-04-03  Synge Todo  <wistaria@comp-phys.org>

        * Introduced site_parameter<T> class template.
        
2004-04-01  Synge Todo  <wistaria@comp-phys.org>

        * Add new test generate_01, in which icc 8.0 with -xW fails, and
        changed looper/path_integral.h as a workaround for this bug.
        
        * Changed looper/path_integral.h as a workaround for a bug in icc
        8.0 with -O1 or higher.
        
2004-03-23  Synge Todo  <wistaria@comp-phys.org>
        
        * Changed looper/graph.h so that the generate_graph function uses
          ALPS/Lattice Library.

2004-03-06  Synge Todo  <wistaria@comp-phys.org>

        * New configure script and Makefiles using GNU automake tool.

2004-03-04  Synge Todo  <wistaria@comp-phys.org>

        * Adopted parameter and measurement names to the new 
          ALPS/QMC standard.  (See applications/qmc/sse/conventions.txt)
        
2004-02-29  Matthias Troyer <troyer@comp-phys.org>

        * Added support for default parameters in Hamiltonian
        
2004-02-19  Synge Todo  <wistaria@comp-phys.org>

        * Supported installing headers etc.
        
2004-01-16  Synge Todo  <wistaria@comp-phys.org>
        
        * Fixed worker<QMC_WORKER>::work_done() function.

        * Changed measurement type in qmc_impl.h and percolation_impl.h.
        
2003-11-20  Synge Todo  <wistaria@comp-phys.org>

        * version 3.0b4 released.
        - prerequisites: boost-ss-20031116, boost-ss-20031116.patch,
          alps-20031116
        - TODO: optimization, more measurements, license issue.  Bug fixes
          for some mesurements (specific heat in path integral?)...
        
2003-11-20  Synge Todo  <wistaria@comp-phys.org>

        * Fixed susceptibility calculation by means of improved estimator
          in SSE representation.
        
2003-11-16  Synge Todo  <wistaria@comp-phys.org>

        * Fixed specific heat calculation in SSE representation.
        
2003-11-14  Synge Todo  <wistaria@comp-phys.org>
        
        * version 3.0b3 released.
        - prerequisites: boost-ss-20031109, boost-ss-20031109.patch,
          alps-20031114
        - TODO: optimization, improved estimators, more measurements,
          license issue.  Bug fixes for some mesurements...

2003-11-13  Synge Todo  <wistaria@comp-phys.org>
        
        * Added check directory, which contains sample parameter for
          qmc (qmc_mpi) for checking results.
        
2003-11-12  Synge Todo  <wistaria@comp-phys.org>
        
        * Updated copyright and file headers and added program `copyright'.

        * Added exact diagonalization code (from dpQLM) for checking QMC
          results.
        
2003-11-07  Synge Todo  <wistaria@comp-phys.org>
        
        * Added improved estimators for uniform/staggered magnetization^2,
          etc.
        
2003-11-06  Synge Todo  <wistaria@comp-phys.org>
        
        * version 3.0b2 released.
        - prerequisites: boost-ss-20031102, boost-ss-20031102.patch,
          alps-20031105
        - TODO: optimization, improved estimators, more measurements,
          license issue.

2003-11-06  Synge Todo  <wistaria@comp-phys.org>
        
        * Moved measurements into looper/measuerment.h

        * Added is_bipartite member in [path_integral,sse]::parameter_type.

        * Removed check_and_resize from path_integral.h.

        * Added specific heat measurement.

        * Fixed staggered susceptibility measurement.
        
2003-11-05  Synge Todo  <wistaria@comp-phys.org>
        
        * Modified Makefile.in so that MPI version of programs (qmc_mpi,
          percolation_mpi) will also be built.

        * Added measurement for staggered susceptibility (unimproved
          estimator).  (SSE output seems still incorrect.)
        
2003-11-05  Synge Todo  <wistaria@comp-phys.org>
        
        * version 3.0b1 released.
        - prerequisites: boost-ss-20031102, boost-ss-20031102.patch,
          alps-20031105
