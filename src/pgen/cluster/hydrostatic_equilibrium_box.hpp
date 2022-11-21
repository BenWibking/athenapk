//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2021, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hydrostatic_equilbirum_box.hpp
//  \brief Class for initializing a sphere in hydrostatic equiblibrium

#ifndef BOX_HYDROSTATIC_EQUILIBRIUM_BOX_HPP_
#define BOX_HYDROSTATIC_EQUILIBRIUM_BOX_HPP_

// Parthenon headers
#include <mesh/domain.hpp>
#include <parameter_input.hpp>

// AthenaPK headers
#include "../../units.hpp"

namespace box {

/************************************************************
 *  Hydrostatic Equilbrium Sphere Class,
 *  for initializing a sphere in hydrostatic equiblibrium
 *
 *
 *  GravitationField:
 *    Graviational field class with member function g_from_r(parthenon::Real r)
 *  EntropyProfile:
 *    Entropy profile class with member function g_from_r(parthenon::Real r)
 ************************************************************/
template <typename GravitationalField, typename EntropyProfile>
class HydrostaticEquilibriumBox {
 private:
  // Graviational field and entropy profile
  const GravitationalField gravitational_field_;
  const EntropyProfile entropy_profile_;

  // Physical constants
  parthenon::Real atomic_mass_unit_, k_boltzmann_;

  // Density to fix baryons at a radius (change to temperature?)
  parthenon::Real z_fix_, rho_fix_;

  // Molecular weights
  parthenon::Real mu_, mu_e_;

  // R mesh sampling parameters
  parthenon::Real z_sampling_, max_dz_;

  /************************************************************
   * Functions to build the cluster model
   *
   * Using lambda functions to make picking up model parameters seamlessly
   * Read A_from_B as "A as a function of B"
   ************************************************************/

  // Get pressure from density and entropy, using ideal gas law and definition
  // of entropy
  parthenon::Real P_from_rho_K(const parthenon::Real rho, const parthenon::Real K) const {
    const parthenon::Real P =
        K * pow(mu_ / mu_e_, 2. / 3.) * pow(rho / (mu_ * atomic_mass_unit_), 5. / 3.);
    return P;
  }

  // Get density from pressure and entropy, using ideal gas law and definition
  // of entropy
  parthenon::Real rho_from_P_K(const parthenon::Real P, const parthenon::Real K) const {
    const parthenon::Real rho =
        pow(P / K, 3. / 5.) * mu_ * atomic_mass_unit_ / pow(mu_ / mu_e_, 2. / 5);
    return rho;
  }

  parthenon::Real rho_from_z(const parthenon::Real z) const {
    const parthenon::Real rho = rho_fix_ * std::tanh(z_fix_ / z);
    return rho;
  }

  // Get total number density from density
  parthenon::Real n_from_rho(const parthenon::Real rho) const {
    const parthenon::Real n = rho / (mu_ * atomic_mass_unit_);
    return n;
  }

  // Get electron number density from density
  parthenon::Real ne_from_rho(const parthenon::Real rho) const {
    const parthenon::Real ne = mu_ / mu_e_ * n_from_rho(rho);
    return ne;
  }

  // Get the temperature from density and pressure
  parthenon::Real T_from_rho_P(const parthenon::Real rho, const parthenon::Real P) const {
    const parthenon::Real T = P / (n_from_rho(rho) * k_boltzmann_);
    return T;
  }

  // Get the dP/dz from radius and pressure, which we will     //
  /************************************************************
   *  dP_dz_from_z_P_functor
   *
   *  Functor class giving dP/dz (z,P)
   *  which is used to compute the HSE profile
   ************************************************************/
  class dP_dz_from_z_P_functor {
    const HydrostaticEquilibriumBox<GravitationalField, EntropyProfile> &box_;

   public:
    dP_dz_from_z_P_functor(
        const HydrostaticEquilibriumBox<GravitationalField, EntropyProfile> &box)
        : box_(box) {}
    parthenon::Real operator()(const parthenon::Real z, const parthenon::Real P) const {

      const parthenon::Real g = box_.gravitational_field_.g_from_z(z);
      const parthenon::Real rho = box_.rho_from_z(z);
      const parthenon::Real dP_dz = -rho * g;
      return dP_dz;
    }
  };

  // Takes one rk4 step from t0 to t1, taking y0 and returning y1, using y'(t) = f(t,y)
  template <typename Function>
  parthenon::Real step_rk4(const parthenon::Real t0, const parthenon::Real t1,
                           const parthenon::Real y0, Function f) const {
    const parthenon::Real h = t1 - t0;
    const parthenon::Real k1 = f(t0, y0);
    const parthenon::Real k2 = f(t0 + h / 2., y0 + h / 2. * k1);
    const parthenon::Real k3 = f(t0 + h / 2., y0 + h / 2. * k2);
    const parthenon::Real k4 = f(t0 + h, y0 + h * k3);
    const parthenon::Real y1 = y0 + h / 6. * (k1 + 2 * k2 + 2 * k3 + k4);
    return y1;
  }

  static constexpr parthenon::Real kRTol = 1e-15;

 public:
  HydrostaticEquilibriumBox(parthenon::ParameterInput *pin,
                               GravitationalField gravitational_field,
                               EntropyProfile entropy_profile);

  template <typename View1D>
  class PRhoProfile {
   private:
    const View1D z_;
    const View1D P_;
    const HydrostaticEquilibriumBox &box_;

    const int n_z_;
    const parthenon::Real z_start_, z_end_;

   public:
    PRhoProfile(const View1D z, const View1D P,
                const HydrostaticEquilibriumBox &box)
        : z_(z), P_(P), box_(box), n_z_(z_.extent(0)), z_start_(z_(0)),
          z_end_(z_(n_z_ - 1)) {}

    parthenon::Real P_from_z(const parthenon::Real z) const;
    parthenon::Real rho_from_z(const parthenon::Real z) const;
    std::ostream &write_to_ostream(std::ostream &os) const;
  };

  template <typename View1D, typename Coords>
  PRhoProfile<View1D>
  generate_P_rho_profile(parthenon::IndexRange ib, parthenon::IndexRange jb,
                         parthenon::IndexRange kb, Coords coords) const;

  template <typename View1D>
  PRhoProfile<View1D> generate_P_rho_profile(const parthenon::Real z_start,
                                             const parthenon::Real z_end,
                                             const unsigned int n_z) const;
};

} // namespace box

#endif // BOX_HYDROSTATIC_EQUILIBRIUM_BOX_HPP_
