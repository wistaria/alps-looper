/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2006 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef TRANSMAG_CONFIG_H
#define TRANSMAG_CONFIG_H

#include <alps/lattice.h>
#include <looper/location.h>
#include <looper/graph.h>

#include <looper/transmag_measurement.h>

struct loop_config
{
  // lattice structure
  typedef alps::coordinate_graph_type lattice_graph_t;
  typedef looper::virtual_lattice<lattice_graph_t> virtual_lattice_t;

  // imaginary time
  typedef double time_t;

  // model, weights, and local_graph
  typedef looper::local_graph<looper::location> loop_graph_t;

  // measurements
  typedef looper::composite_estimator<
          looper::composite_estimator<
            looper::susceptibility_estimator<virtual_lattice_t, time_t>,
            looper::stiffness_estimator<virtual_lattice_t, time_t, 3> >,
            transverse_magnetization_estimator<virtual_lattice_t, time_t>
  > estimator_t;
};

#endif // TRANSMAG_CONFIG_H
