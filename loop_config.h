/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
*
* This software is published under the ALPS Application License; you
* can use, redistribute it and/or modify it under the terms of the
* license, either version 1 or (at your option) any later version.
* 
* You should have received a copy of the ALPS Application License
* along with this software; see the file LICENSE. If not, the license
* is also available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifdef LOOP_CONFIG_HEADER

#include LOOP_CONFIG_HEADER

#else

#ifndef LOOP_CONFIG_H
#define LOOP_CONFIG_H

#include <looper/lattice.h>
#include <looper/location.h>
#include <looper/graph.h>
#include <looper/measurement.h>

struct loop_config
{
  // lattice structure
  typedef looper::graph_type lattice_graph_t;

  // imaginary time
  typedef double time_t;

  // model, weights, and local_graph
  typedef looper::local_graph<looper::location> loop_graph_t;

  // measurements
  // typedef looper::composite_estimator<
  //   looper::susceptibility_estimator,
  //   looper::stiffness_estimator<3> >
  //   estimator_t;
  typedef looper::susceptibility_estimator estimator_t;
};

#endif // LOOP_CONFIG_H

#endif // LOOP_CONFIG_HEADER
