//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2021, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file box_gravity.hpp
//  \brief Class for defining gravitational acceleration for a cluster+bcg+smbh
#ifndef BOX_GRAVITY_HPP_
#define BOX_GRAVITY_HPP_

// Parthenon headers
#include <parameter_input.hpp>

// AthenaPK headers
#include "../../units.hpp"

namespace box {

/************************************************************
 *  Cluster Gravity Class, for computing gravitational acceleration
 *    Lightweight object for inlined computation within kernels
 ************************************************************/
class ClusterGravity {

  // Radius underwhich to truncate
  parthenon::Real smoothing_z_;

 public:
  ClusterGravity(parthenon::ParameterInput *pin) {
    Units units(pin);

    smoothing_z_ = pin->GetOrAddReal("problem/cluster", "g_smoothing_height", 0.0);
  }

  // Inline functions to compute gravitational acceleration
  KOKKOS_INLINE_FUNCTION parthenon::Real g_from_z(const parthenon::Real z_in) const
      __attribute__((always_inline)) {

    const parthenon::Real z = std::abs(z_in);

    parthenon::Real g_z = 0;

    return g_z;
  }  
};

} // namespace box

#endif // BOX_GRAVITY_HPP_
