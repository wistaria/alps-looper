# ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems

Copyright (C) 1997-2015 by Synge Todo <wistaria@comp-phys.org>

This software is published under the ALPS Application License; you can use,
redistribute and/or modify this software under the terms of the license,
either version 1 or (at your option) any later version.

You should have received a copy of the ALPS Application License along with
the ALPS Library; see the file LICENSE. If not, the license is also
available from http://alps.comp-phys.org/.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
DEALINGS IN THE SOFTWARE.

## Overview

The ALPS/looper package is one of the applications of the ALPS Project.  The ALPS/looper is the successor of the Looper Library version 2, which is the Fortran 90 library for quantum Monte Carlo method.

The ALPS/looper application provides multi-cluster quantum Monte Carlo algorithms for generic spin systems in the path-integral or SSE representations.  it supports generic spin models with arbitrary spin size *S* on arbitrary lattices. The supported interactions are

* XXZ two-spin interaction (Ising, Heisenberg, XY)
* XYZ two-spin interaction (requires recompilation)
* single-ion anisotropy (D-term) for *S* > 1/2
* longitudinal magnetic field
* transverse matnetic field

Since it supports so-called *freezing graphs*, models with easy-axis
(Ising-like) anisotropy can also be simulated very efficiently.

From version 4, it includes non-trivial hybrid parallelization by using MPI and OpenMP, which has been tested on several state-of-the-art supercomputers (The K computer, SGI Altix, etc).

A stripped version of the ALPS/looper is distributed as a part of ALPS.  The full version of ALPS/looper library (including test suites, etc) is available from [http://wistaria.github.io/alps-looper/](http://wistaria.github.io/alps-looper/).

## History

* Version 1
   * 1997: written in FORTRAN 77, only works for S=1/2, bipartite lattice, Heisenberg and easy plane anisotropy
* Version 2.0
   * 1998/10/09 version 2.0: Fortran 90, large S support
* Version 2.1
   * 1998/10/14 version 2.1.0
   * 1998/10/19 version 2.1.4
* Version 2.2
   * 1998/10/23 version 2.2.0:  Parallelization (MPI) support
   * 2001/02/01 version 2.2.9
* Version 2.3
   * 2000/10/31 version 2.3: Freezing (easy-axis anisotropy) support
* Version 2.9
   * 2001/08/16 version 2.9.1: C++ version, non-bipartite lattice support
   * 2001/10/26 version 2.9.6
* Version 3.0
   * 2004/06/11 version 3.0: Based on ALPS framework
   * 2004/09/23 version 3.1: Negative sign support
* Version 4
   * March 1, 2012: Alpha version for ALPS/looper 4 (4.0a1)
   * March 5, 2012: First release of ALPS/looper 4 (4.0.0)

## Installation

Note: This section applies only for the full version of ALPS/looper.

### Prerequisites

* [Boost C++ Library](http://www.boost.org/) -- version 1.48.0 or higher
* [ALPS](http://alps-comp-phys.org) -- version 2.0 or higher
* [LAPACK Library](http://www.netlib.org/lapack/)
  
### Download
  
[https://github.com/wistaria/alps-looper/releases](https://github.com/wistaria/alps-looper/releases)

### Compilation

```sh
$ cmake -DCMAKE_INSTALL_PREFIX=$ALPS_HOME $LOOPER_SRC_DIR
$ make
```

where `$ALPS_HOME` denotes the top directory in which the ALPS package is installed, `$LOOPER_SRC_DIR` is the directory where the source code of the ALPS/looper is located.

CMake accepts the following configuration options:

| Name | Default value | Descriptin |
|:-----|:--------------|:-----------|
| LOOPER\_BUILD\_TESTS | ON | Build looper tests |
| LOOPER\_BUILD\_EXTRAS | OFF | Build in extras subdirectory |
| LOOPER\_ENABLE\_DIAG | ON | Enable exact diagonalization method |
| LOOPER\_ENABLE\_ISING | ON | Enable Swendsen-Wang algorithm for classical Ising model |
| LOOPER\_ENABLE\_PI | ON | Enable loop algorithm in path integral reprentation |
| LOOPER\_ENABLE\_SSE | ON | Enable loop algorithm in high-temperature expansion representation |
| LOOPER\_ENABLE\_OPENMP | ON | OpenMP (hybrid parallelization) support |
| LOOPER\_ENABLE\_TIMER | OFF | Perform time measurement of each section |
| LOOPER\_TIMER\_TRACE | OFF | Turn on debug trace output of timer |
| LOOPER\_TIMER\_DETAILED | OFF | Turn on detailed timer report |

Note: OpenMP (and hybrid parallelization) support will be enabled only if the ALPS library has been configured with `-DALPS_ENABLE_OPENMP=ON` and
`-DALPS_ENABLE_OPENMP_WORKER=ON`.

### Running tests

```sh
$ ctest
```
  
## Running a simulation

### Command line options

| Option | Default value | Descriptin |
|:-------|:--------------|:-----------|
| -h [ --help ] | | produce help message |
| -l [ --license ] | | print license conditions |
| --auto-evaluate | evaluate observables upon halting [default] |
| --check-parameter | | perform parameter checking [default| |
| --check-interval arg | 100 | time between internal status check [unit = millisec] |
| --checkpoint-interval arg | 3600 | time between checkpointing [unit = sec] |
| --evaluate | | only perform evaluation |
| --mpi | | run in parallel using MPI |
| --Nmin | | [obsolete] |
| --Nmax | | [obsolete] |
| --no-evaluate | | prevent evaluating observables upon halting |
| --report-interval arg | 600 | time between progress report of clones [unit = sec] |
| -T [ --time-limit ] arg | unlimited | time limit for the simulation [unit = sec] |
| --Tmin | | [obsolete] |
| --Tmax | | [obsolete] |
| -p [ --threads-per-clone ] arg | 1 | number of threads for each clone |
| -r [ --total-threads ] arg | # of MPI processes | total number of threads [integer or 'auto'] |
| --input-file arg | | input master XML files |

### Parameters

| Name | Default Value | Description |
|:-----|:--------------|:------------|
| NUM\_CLONES | 1 number of clones for each task |
| LATTICE\_LIBRARY | ``$PREFIX/lib/xml/lattices.xml`` | path to a file containing lattice descriptions |
| LATTICE | *none* | name of the lattice |
| MODEL\_LIBRARY | ``$PREFIX/lib/xml/models.xml`` | path to a file containing model descriptions |
| MODEL | *none* | name of the model |
| ALGORITHM | *none* | type of algorithm ("loop", "loop; path integral", "loop; sse", "ising", or "diagonalization") |
| | |
| T | *none* | temperature |
| T\_START|_\# | *undefined* | [optional] see section "Temperature Annealing" |
| T\_DURATION\_\# | *undefined* | [optional] see section "Temperature Annealing" |
| | |
| SWEEPS | 65536 | Monte Carlo steps after thermalization |
| THERMALIZATION | SWEEPS/8 | Monte Carlo steps for thermalization |
| | |
| USE\_SITE\_INDICES\_AS\_TYPES | false | if true, all sites will have distinct site types, which are identical to site indices (starting from 0) |
| USE\_BOND\_INDICES\_AS\_TYPES | false | if true, all bonds will have distinct bond types, which are identical to bond indices (starting from 0) |
| | |
| PARTITION | *none* | [parallel] list of number of processes (separated by *colons*) used in each stage of cluster unification (e.g. 2:3:3:2).  The product of all the factors should be equal to the total number of processes for each clone |
| DUPLEX | true | [parallel] if true, bidirectional communication is enabled |
| | |
| MEASURE[Correlations] | false | if true, correlation function will be calculated |
| MEASURE[Green Function] | false | if true, Green's function will be calculated |
| INITIAL\_SITE | *undefined* | initical site from which correlation function and Green's function are measured.  If not defined, correlation between all possible site pairs will be calculated |
| MEASURE[Local Susceptibility] | false | if true, local suscepbility and local magnetization will be calculated |
| MEASURE[Structure Factor] | false | if true, structure factor for all possible k-values will be calculated |
| | |
| DISABLE\_IMPROVED\_ESTIMATOR | false | [optional] use normal (i.e. unimproved) estimator for measurements (will be set to true automatically in the presence of longitudinal magnetic field) |
| FORCE\_SCATTER | 0 | [optional] minimum probability for forward scattering graph (will be set to 0.1 automatically for classically frustrated models in order to ensure the ergodicity |
| | |
| LOOPER\_DEBUG[MODEL OUTPUT] | *undefined* | if defined, all the coupling constants will be printed out to *standard output* (std::cout) before calculation.  If the value is ``cerr``, the output will be made for *standard error* (std::cerr).

Note: when PARTITION is not specified, the number of MPI processes per clone should be factorable by 2, 3, 5, 7, 11, and 13.

In addition, the lattice/model descriptions can require further parameters (e.g. L or W) as specified in the lattice (model) description file.

## Temperature Annealing

Some models (typically models with competing interactions) have very long equilibration time.  In such cases, one might want to lower the temperature slowly during the thermalization steps.  Such temperature annealing process can be specified by using the parameters T\_START\_\# and T\_DURATION\_\# (\# = 0, 1, 2,...), where T\_START\_\# and \T\_DURATION\_\# specify the starting temperature and the duration of the \# th block, respectively.

![anneal](anneal.jpg)

For example the protocol shown in the above figure can be realized by specifing the following parameter set.

```
T_START_0 = 2.0;
T_DURATION_0 = 100;
T_START_1 = 1.5;
T_DURATION_1 = 400;
THERMALIZATION = 1000;
T = 1.2;
```

Note that the sum of T\_DURATION\_\#'s must be less than or equal to THERMALIZATION.

## Measurements

The following observables are measured by the loop application:

| Name | Description |
|:-----|:------------|
| Energy | total energy |
| Energy Density | energy per spin |
| Energy^2 | square of total energy |
| Specific Heat | specific heat |
| | |
| Stiffness | stiffness constant |
| | |
| Magnetization | total uniform magnetization |
| Magnetization^2 | square of total uniform magnetization |
| Magnetization^4 | 4th power of total uniform magnetization |
| Binder Ratio of Magnetization | Binder ratio of uniform magnetization |
| Magnetization Density | uniform magnetization per spin |
| Magnetization Density^2 | square of uniform magnetization per spin |
| Magnetization Density^4 | 4th power of uniform magnetization per spin |
| Susceptibility | uniform susceptibility |
| | |
| Staggered Magnetization | total staggered magnetization [1] |
| Staggered Magnetization^2 | square of total staggered magnetization [1] |
| Staggered Magnetization^4 | 4th power of total staggered magnetization [1] |
| Staggered Magnetization Density | staggered magnetization per spin [1] |
| Staggered Magnetization Density^2 | square of staggered magnetization per spin [1] |
| Staggered Magnetization Density^4 | 4th power of staggered magnetization per spin [1] |
| Staggered Susceptibility | staggered susceptibility [1] |
| | |
| Generalized Magnetization | total generalized magnetization [2] |
| Generalized Magnetization^2 | square of total generalized magnetization [2] |
| Generalized Magnetization^4 | 4th power of total generalized magnetization [2] |
| Generalized Binder Ratio of Magnetization | Binder ratio of generalized magnetization [2] |
| Generalized Magnetization Density | generalized magnetization per spin [2] |
| Generalized Magnetization Density^2 | square of generalized magnetization per spin [2] |
| Generalized Magnetization Density^4 | 4th power of generalized magnetization per spin [2] |
| Generalized Susceptibility | generalized susceptibility [2] |
| | |
| Generalized Staggered Magnetization | total generalized staggered magnetization [1,2] |
| Generalized Staggered Magnetization^2 | square of total generalized staggered magnetization [1,2] |
| Generalized Staggered Magnetization^4 | 4th power of total generalized staggered magnetization [1,2] |
| Generalized Staggered Magnetization Density | generalized staggered magnetization per spin [1,2] |
| Generalized Staggered Magnetization Density^2 | square of generalized staggered magnetization per spin [1,2] |
| Generalized Staggered Magnetization Density^4 | 4th power of generalized staggered magnetization per spin [1,2] |
| Generalized Staggered Susceptibility | generalized staggered susceptibility [1,2] |
| | |
| Spin Correlations | correlation functions [3] |
| Staggered Spin Correlations | staggered correlation functions [1,3] |
| Generalized Spin Correlations | generalized correlation functions [2,3] |
| Generalized Staggered Spin Correlations | generalized staggered correlation function [1,2,3] |
| | |
| Green's Function | Green's function [2, 4] |
| | |
| Structure Factor | static structure factor [5] |
| | |
| Local Magnetization | local magnetization [6] |
| Local Susceptibility | local susceptibility (i.e., response of uniform magnetization against local magnetic field or that of local magnetization against uniform magnetic field) [6] |
| Staggered Local Susceptibility | statggered local susceptibility (i.e., response of staggered magnetization against local magnetic field or that of local magnetization against staggered magnetic field) [1,6] |
| Local Field Susceptibility |local-field susceptibility (i.e., response of local magnetization against local magnetic field) [6] |

[1] measured only for bipartite lattices

[2] measured only by improved estimator

[3] set MEASURE[Correlations] to true

[4] set MEASURE[Green Function] to true

[5] set MEASURE[Structure Factor] to true
[6] set MEASURE[Local Susceptiblity] to true

A few remarks concerning the generalized susceptibilities: The general idea is that the uniform generalized susceptibility is the susceptibility for the magnetic order supposed to be stronger.  The staggered generalized susceptibility is the staggered susceptibility of this specific order.  To be more precise, we have the following correspondance:

| Model | Generalized Susceptibility | Generalized Staggered Susceptibility |
|:------|:---------------------------|:-----------------------------------------|
| Ferromagnetic model with Ising anisotropy (in the *Z* axis) | uniform *Z-Z* susceptibility | staggered *Z-Z* susceptibility |
Ferromagnetic model (Heisenberg point) | uniform *Z-Z* (or *X-X*) susceptibility | staggered *Z-Z* (or *X-X*) susceptibility |
| Ferromagnetic model with XY anisotropy (in the *XY* plane) | uniform *X-X* susceptibility | staggered *X-X* susceptibility |
| Antiferromagnetic model with XY anisotropy (in the *XY* plane) | staggered *X-X* susceptibility             uniform *X-X* susceptibility |
| Antiferromagnetic model (Heisenberg point) | staggered *Z-Z* (or *X-X*) susceptibility | uniform *Z-Z* (or *X-X*) susceptibility |
| Antiferromagnetic model with Ising anisotropy (in the *Z* axis) | staggered *Z-Z* susceptibility | uniform *Z-Z* susceptibility |

Same considerations apply for the uniform and staggered magnetization to the square.

### Tutorials

ALPS/looper tutorials can be found in [ALPS Tutorials Page](http://alps.comp-phys.org/mediawiki/index.php/ALPS_2_Tutorials:MC-02_Susceptibilities).

## License

The [license](http://alps.comp-phys.org/software/applications/LICENSE.txt) allows the use of the applications for non-commercial scientific use provided that the use of the ALPS/looper Library and the ALPS Libraries is acknowledged, and the papers listed below are referenced in any scientific publication.  For detail please see the [ALPS Applications Licence](http://alps.comp-phys.org/software/applications/LICENSE.txt).

* reference

  - http://alps.comp-phys.org/ 
  - http://wistaria.comp-phys.org/alps-looper/

* publications

  - S. Todo and K. Kato, *Cluster Algorithms for General-S Quantum Spin Systems*, Phys. Rev. Lett., **87**, 047203 (2001).

  - A.F. Albuquerque, F. Alet, P. Corboz, P. Dayal, A. Feiguin, L. Gamper, E. Gull, S Gürtler, A. Honecker, R. Igarashi, M. Körner, A. Kozhevnikov, A. Läuchli, S.R. Manmana, M. Matsumoto, I.P. McCulloch, F. Michel, R.M. Noack, G. Pawlowski, L. Pollet, T. Pruschke, U. Schollwöck, S. Todo, S. Trebst, M. Troyer, P. Werner, S. Wessel (ALPS collaboration), *The ALPS project release 1.3: open source software for strongly correlated systems*, [J. Mag. Mag. Mat. 310, 1187 (2007)](http://dx.doi.org/10.1016/j.jmmm.2006.10.304).

  - B. Bauer, L. D. Carr, A. Feiguin, J. Freire, S. Fuchs, L. Gamper, J. Gukelberger, E. Gull, S. Guertler, A. Hehn, R. Igarashi, S.V. Isakov, D. Koop, P.N. Ma, P. Mates, H. Matsuo, O. Parcollet, G. Pawlowski, J.D. Picon, L. Pollet, E. Santos, V.W. Scarola, U. Schollwöck, C. Silva, B. Surer, S. Todo, S. Trebst, M. Troyer, M.L. Wall, P. Werner, S. Wessel, *The ALPS project release 2.0: Open source software for strongly correlated systems*, [Journal of Statistical Mechanics: Theory and Experiment, P05001 (2011)](http://iopscience.iop.org/1742-5468/2011/05/P05001).

## Developers

* Version 1: Synge Todo
* Version 2: Synge Todo, Kiyoshi Kato, and Shigeru Koikegami
* Version 3: Synge Todo
* Version 4: Synge Todo, Haruhiko Matuo, and Hideyuki Shitara

## Questions and request for support

can be addressed to Synge Todo <wistaria@comp-phys.org>.

## Acknowledgement

I wish to thank M. Troyer and F. Alet for many useful comments and suggestions.  The development of ALPS/looper version 4 was supported by Grand Challenges in Next-Generation Integrated Nanoscience, Next-Generation Supercomputer Project, MEXT, Japan.

## References

* H.G. Evertz, *The loop algorithm*, Adv. in Physics, **52**, 1 (2003).

* S. Todo and K. Kato, *Cluster Algorithms for General-S Quantum Spin Systems*, Phys. Rev. Lett., **87**, 047203 (2001).

* S. Todo, *Parallel Quantum Monte Carlo Simulation of S=3 Antiferromagnetic Heisenberg Chain*, Computer Simulation Studies in Condensed-Matter Physics XV, ed D.P. Landau, S.P. Lewis, and H.-B. Schuettler (Springer-Verlag, Berlin, 2003) pp. 89-94.

* F. Alet, P. Dayal, A. Grzesik, A. Honecker, M. Körner, A. Läuchli, S.R. Manmana, I.P. McCulloch, F. Michel, R.M. Noack, G. Schmid, U. Schollwöck, F. Stöckli, S. Todo, S. Trebst, M. Troyer, P. Werner, S. Wessel (ALPS collaboration), *The ALPS project: open source software for strongly correlated systems*,
  [J. Phys. Soc. Jpn. Suppl. 74, 30 (2005)]
  (http://jpsj.ipap.jp/link?JPSJS/74S/30).

* A.F. Albuquerque, F. Alet, P. Corboz, P. Dayal, A. Feiguin, L. Gamper, E. Gull, S Gürtler, A. Honecker, R. Igarashi, M. Körner, A. Kozhevnikov, A. Läuchli, S.R. Manmana, M. Matsumoto, I.P. McCulloch, F. Michel, R.M. Noack, G. Pawlowski, L. Pollet, T. Pruschke, U. Schollwöck, S. Todo, S. Trebst, M. Troyer, P. Werner, S. Wessel (ALPS collaboration), *The ALPS project release 1.3: open source software for strongly correlated systems*, [J. Mag. Mag. Mat. 310, 1187 (2007)](http://dx.doi.org/10.1016/j.jmmm.2006.10.304).

* B. Bauer, L. D. Carr, A. Feiguin, J. Freire, S. Fuchs, L. Gamper, J. Gukelberger, E. Gull, S. Guertler, A. Hehn, R. Igarashi, S.V. Isakov, D. Koop, P.N. Ma, P. Mates, H. Matsuo, O. Parcollet, G. Pawlowski, J.D. Picon, L. Pollet, E. Santos, V.W. Scarola, U. Schollwöck, C. Silva, B. Surer, S. Todo, S. Trebst, M. Troyer, M.L. Wall, P. Werner, S. Wessel, *The ALPS project release 2.0: Open source software for strongly correlated systems*, [Journal of Statistical Mechanics: Theory and Experiment, P05001 (2011)](http://iopscience.iop.org/1742-5468/2011/05/P05001).

-------------------------------------

Copyright (c) 1997-2015 by Synge Todo <wistaria@comp-phys.org>