//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2021, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file entropy_profiles_cartesian.hpp
//  \brief Classes defining initial entropy profile
#ifndef BOX_ENTROPY_PROFILES_HPP_
#define BOX_ENTROPY_PROFILES_HPP_

// Parthenon headers
#include <parameter_input.hpp>

// AthenaPK headers
#include "../../units.hpp"

namespace box {

class ACCEPTEntropyProfile {
 private:
  // Entropy Profile
  parthenon::Real K_0_, K_100_, z_K_, alpha_K_;

 public:
  ACCEPTEntropyProfile(parthenon::ParameterInput *pin) {
    Units units(pin);

    K_0_ = pin->GetOrAddReal("problem/cluster", "K_0",
                             20 * units.kev() * units.cm() * units.cm());
    K_100_ = pin->GetOrAddReal("problem/cluster", "K_100",
                               120 * units.kev() * units.cm() * units.cm());
    z_K_ = pin->GetOrAddReal("problem/cluster", "z_K", 100 * units.kpc());
    alpha_K_ = pin->GetOrAddReal("problem/cluster", "alpha_K", 1.75);
  }

  // Get entropy from radius, using broken power law profile for entropy
  parthenon::Real K_from_z(const parthenon::Real z) const {
    const parthenon::Real K = K_0_ + K_100_ * pow(z / z_K_, alpha_K_);
    return K;
  }
};

} // namespace box

#endif // CLUSTER_ENTROPY_PROFILES_HPP_
