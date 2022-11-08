//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2021, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hydrostatic_equilbirum_sphere.cpp
//  \brief Creates pressure profile in hydrostatic equilbrium
//
// Setups up a pressure profile in hydrostatic equilbrium given an entropy
// profile and gravitational field
//========================================================================================

// C++ headers
#include <fstream>

// Parthenon headers
#include <coordinates/uniform_cartesian.hpp>
#include <globals.hpp>
#include <mesh/domain.hpp>
#include <parameter_input.hpp>

// AthenaPK headers
#include "../../units.hpp"

// Cluster headers
#include "box_gravity.hpp"
#include "entropy_profiles_cartesian.hpp"
#include "hydrostatic_equilibrium_box.hpp"

namespace box {
using namespace parthenon;

/************************************************************
 * HydrostaticEquilibriumBox constructor
 ************************************************************/
template <typename GravitationalField, typename EntropyProfile>
HydrostaticEquilibriumBox<GravitationalField, EntropyProfile>::
    HydrostaticEquilibriumBox(ParameterInput *pin,
                                 GravitationalField gravitational_field,
                                 EntropyProfile entropy_profile)
    : gravitational_field_(gravitational_field), entropy_profile_(entropy_profile) {
  Units units(pin);

  atomic_mass_unit_ = units.atomic_mass_unit();
  k_boltzmann_ = units.k_boltzmann();

  const Real He_mass_fraction = pin->GetReal("hydro", "He_mass_fraction");
  const Real H_mass_fraction = 1.0 - He_mass_fraction;

  mu_ = 1 / (He_mass_fraction * 3. / 4. + (1 - He_mass_fraction) * 2);
  mu_e_ = 1 / (He_mass_fraction * 2. / 4. + (1 - He_mass_fraction));

  z_fix_ =
      pin->GetOrAddReal("problem/cluster", "z_fix", 1953.9724519818478 * units.kpc());
  rho_fix_ = pin->GetOrAddReal("problem/cluster", "rho_fix",
                               8.607065015897638e-30 * units.g() / pow(units.kpc(), 3));
  const Real gam = pin->GetReal("hydro", "gamma");
  const Real gm1 = (gam - 1.0);

  z_sampling_ = pin->GetOrAddReal("problem/cluster", "z_sampling", 4.0);
  max_dz_ = pin->GetOrAddReal("problem/cluster", "max_dz", 1e-3);

  // Test out the HSE sphere if requested
  const bool test_he_sphere =
      pin->GetOrAddBoolean("problem/cluster", "test_he_sphere", false);
  if (test_he_sphere) {
    const Real test_he_sphere_z_start = pin->GetOrAddReal(
        "problem/cluster", "test_he_sphere_z_start_kpc", 1e-3 * units.kpc());
    const Real test_he_sphere_z_end = pin->GetOrAddReal(
        "problem/cluster", "test_he_sphere_z_end_kpc", 4000 * units.kpc());
    const int test_he_sphere_n_z =
        pin->GetOrAddInteger("problem/cluster", "test_he_sphere_n_z", 4000);
    if (Globals::my_rank == 0) {
      typedef Kokkos::View<Real *, Kokkos::LayoutRight, HostMemSpace> View1D;

      auto P_rho_profile = generate_P_rho_profile<View1D>(
          test_he_sphere_z_start, test_he_sphere_z_end, test_he_sphere_n_z);

      std::ofstream test_he_file;
      test_he_file.open("test_he_sphere.dat");
      P_rho_profile.write_to_ostream(test_he_file);
      test_he_file.close();
    }
  }
}

/************************************************************
 * PRhoProfile::P_from_z
 ************************************************************/
template <typename EntropyProfile, typename GravitationalField>
template <typename View1D>
Real HydrostaticEquilibriumBox<EntropyProfile, GravitationalField>::PRhoProfile<
    View1D>::P_from_z(const Real z) const {

  // Determine indices in R bounding r
  const int i_r =
      static_cast<int>(floor((n_z_ - 1) / (z_end_ - z_start_) * (z - z_start_)));

  if (z < z_(i_r) - kRTol || z > z_(i_r + 1) + kRTol) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [HydrostaticEquilibriumBox::PRhoProfile]"
        << std::endl
        << "R(i_r) to R_(i_r+1) does not contain r" << std::endl
        << "R(i_r) R_r R(i_r+1):" << z_(i_r) << " " << z << " " << z_(i_r + 1)
        << std::endl;
    PARTHENON_FAIL(msg);
  }

  // Linearly interpolate Pressure from P
  const Real P_z = (P_(i_r) * (z_(i_r + 1) - z) + P_(i_r + 1) * (z - z_(i_r))) /
                   (z_(i_r + 1) - z_(i_r));

  return P_z;
}

/************************************************************
 * PRhoProfile::rho_from_r
 ************************************************************/
template <typename EntropyProfile, typename GravitationalField>
template <typename View1D>
Real HydrostaticEquilibriumBox<EntropyProfile, GravitationalField>::PRhoProfile<
    View1D>::rho_from_z(const Real z) const {

  // Get pressure first
  const Real P_z = P_from_z(z);
  // Compute entropy and pressure here
  const Real K_z = box_.entropy_profile_.K_from_z(z);
  const Real rho_z = box_.rho_from_P_K(P_z, K_z);
  return rho_z;
}

/************************************************************
 * PRhoProfile::write_to_ostream
 ************************************************************/
template <typename EntropyProfile, typename GravitationalField>
template <typename View1D>
std::ostream &
HydrostaticEquilibriumBox<EntropyProfile, GravitationalField>::PRhoProfile<
    View1D>::write_to_ostream(std::ostream &os) const {

  const dP_dz_from_z_P_functor dP_dz_func(box_);

  for (int i = 0; i < z_.extent(0); i++) {
    const Real z = z_(i);
    const Real P = P_(i);
    const Real K = box_.entropy_profile_.K_from_z(z);
    const Real rho = box_.rho_from_P_K(P, K);
    const Real n = box_.n_from_rho(rho);
    const Real ne = box_.ne_from_rho(rho);
    const Real T = box_.T_from_rho_P(rho, P);
    const Real g = box_.gravitational_field_.g_from_z(z);
    const Real dP_dz = dP_dz_func(z, P);

    os << z << " " << P << " " << K << " " << rho << " " << n << " " << ne << " " << T
       << " " << g << " " << dP_dz << std::endl;
  }
  return os;
}

/************************************************************
 * HydrostaticEquilibriumBox::generate_P_rho_profile(x,y,z)
 ************************************************************/
template <typename EntropyProfile, typename GravitationalField>
template <typename View1D, typename Coords>
typename HydrostaticEquilibriumBox<EntropyProfile,
                                      GravitationalField>::template PRhoProfile<View1D>
HydrostaticEquilibriumBox<EntropyProfile, GravitationalField>::generate_P_rho_profile(
    IndexRange ib, IndexRange jb, IndexRange kb, Coords coords) const {

  /************************************************************
   * Define R mesh to integrate pressure along
   *
   * R mesh should adapt with requirements of MeshBlock
   ************************************************************/

  // Determine spacing of grid (WARNING assumes equispaced grid in x,y,z)
  PARTHENON_REQUIRE(coords.dx1v(0) == coords.dx1v(1), "No equidistant grid in x1dir");
  PARTHENON_REQUIRE(coords.dx2v(0) == coords.dx2v(1), "No equidistant grid in x2dir");
  PARTHENON_REQUIRE(coords.dx3v(0) == coords.dx3v(1), "No equidistant grid in x3dir");
  PARTHENON_REQUIRE(coords.dx1v(0) == coords.dx2v(1),
                    "No equidistant grid between x1 and x2 dir");
  PARTHENON_REQUIRE(coords.dx2v(0) == coords.dx3v(1),
                    "No equidistant grid between x2 and x3 dir");
  const Real dz = std::min(coords.dx1v(0) / z_sampling_, max_dz_);

  // Loop through mesh for minimum and maximum radius
  // Make sure to include R_fix_
  Real z_start = z_fix_;
  Real z_end = z_fix_;
  for (int k = kb.s; k <= kb.e; k++) {
    for (int j = jb.s; j <= jb.e; j++) {
      for (int i = ib.s; i <= ib.e; i++) {
        const Real z = std::abs(coords.x3v(k));
        z_start = std::min(z, z_start);
        z_end = std::max(z, z_end);
      }
    }
  }

  // Add some room for R_start and R_end
  z_start = std::max(0.0, z_start - z_sampling_ * dz);
  z_end += z_sampling_ * dz;

  // Compute number of cells needed
  const unsigned int n_z = static_cast<unsigned int>(ceil((z_end - z_start) / dz));
  // Make R_end  consistent
  z_end = z_start + dz * (n_z - 1);

  return generate_P_rho_profile<View1D>(z_start, z_end, n_z);
}

/************************************************************
 * HydrostaticEquilibriumBox::generate_P_rho_profile(Ri,Re,nR)
 ************************************************************/
template <typename EntropyProfile, typename GravitationalField>
template <typename View1D>
typename HydrostaticEquilibriumBox<EntropyProfile,
                                      GravitationalField>::template PRhoProfile<View1D>
HydrostaticEquilibriumBox<EntropyProfile, GravitationalField>::generate_P_rho_profile(
    const Real z_start, const Real z_end, const unsigned int n_z) const {

  // Array of radii along which to compute the profile
  View1D z("z", n_z);
  const Real dz = (z_end - z_start) / (n_z - 1.0);

  // Use a linear R - possibly adapt if using a mesh with logrithmic r
  for (int i = 0; i < n_z; i++) {
    z(i) = z_start + i * dz;
  }

  /************************************************************
   * Integrate Pressure inward and outward from virial radius
   ************************************************************/
  // Create array for pressure
  View1D P("P", n_z);

  const Real K_fix = entropy_profile_.K_from_z(z_fix_);
  const Real P_fix = P_from_rho_K(rho_fix_, K_fix);

  // Integrate P inward from R_fix_
  Real zi = z_fix_; // Start zi at z_fix_ first
  Real Pi = P_fix;  // Start with pressure at z_fix_

  // Find the index in R right before R_fix_
  int i_fix = static_cast<int>(floor((n_z - 1) / (z_end - z_start) * (z_fix_ - z_start)));
  if (z_fix_ < z(i_fix) - kRTol || z_fix_ > z(i_fix + 1) + kRTol) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function "
           "[HydrostaticEquilibriumBox::generate_P_rho_profile]"
        << std::endl
        << "R(i_fix) to R_(i_fix+1) does not contain R_fix_" << std::endl
        << "R(i_fix) R_fix_ R(i_fix+1):" << z(i_fix) << " " << z_fix_ << " "
        << z(i_fix + 1) << std::endl;
    PARTHENON_FAIL(msg);
  }

  dP_dz_from_z_P_functor dP_dz_from_z_P(*this);

  // Make is the i right before R_fix_
  for (int i = i_fix + 1; i > 0; i--) { // Move is up one, to account for initial R_fix_
    P(i - 1) = step_rk4(zi, z(i - 1), Pi, dP_dz_from_z_P);
    zi = z(i - 1);
    Pi = P(i - 1);
  }

  // Integrate P outward from R_fix_
  zi = z_fix_; // Start Ri at R_fix_ first
  Pi = P_fix;  // Start with pressure at R_fix_

  // Make is the i right after R_fix_
  for (int i = i_fix; i < n_z - 1;
       i++) { // Move is back one, to account for initial R_fix_
    P(i + 1) = step_rk4(zi, z(i + 1), Pi, dP_dz_from_z_P);
    zi = z(i + 1);
    Pi = P(i + 1);
  }

  return PRhoProfile<View1D>(z, P, *this);
}

// Instantiate HydrostaticEquilibriumBox
template class HydrostaticEquilibriumBox<ClusterGravity, box::ACCEPTEntropyProfile>;

// Instantiate PRhoProfile
template class HydrostaticEquilibriumBox<ClusterGravity, box::ACCEPTEntropyProfile>::
    PRhoProfile<Kokkos::View<parthenon::Real *, parthenon::LayoutWrapper,
                             parthenon::DevMemSpace>>;

#if (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
template class HydrostaticEquilibriumBox<ClusterGravity, box::ACCEPTEntropyProfile>::
    PRhoProfile<Kokkos::View<parthenon::Real *, LayoutWrapper, HostMemSpace>>;
#endif

// Instantiate generate_P_rho_profile
template HydrostaticEquilibriumBox<ClusterGravity, box::ACCEPTEntropyProfile>::PRhoProfile<
    Kokkos::View<parthenon::Real *, parthenon::LayoutWrapper, parthenon::DevMemSpace>>
    HydrostaticEquilibriumBox<ClusterGravity, box::ACCEPTEntropyProfile>::
        generate_P_rho_profile<Kokkos::View<parthenon::Real *, parthenon::LayoutWrapper,
                                            parthenon::DevMemSpace>,
                               parthenon::UniformCartesian>(
            parthenon::IndexRange, parthenon::IndexRange, parthenon::IndexRange,
            parthenon::UniformCartesian) const;

template HydrostaticEquilibriumBox<ClusterGravity, box::ACCEPTEntropyProfile>::PRhoProfile<
    Kokkos::View<parthenon::Real *, parthenon::LayoutWrapper, parthenon::DevMemSpace>>
HydrostaticEquilibriumBox<ClusterGravity, box::ACCEPTEntropyProfile>::
    generate_P_rho_profile<Kokkos::View<parthenon::Real *, parthenon::LayoutWrapper,
                                        parthenon::DevMemSpace>>(
        const parthenon::Real, const parthenon::Real, const unsigned int) const;

#if (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
template HydrostaticEquilibriumBox<ClusterGravity, box::ACCEPTEntropyProfile>::PRhoProfile<
    Kokkos::View<parthenon::Real *, LayoutWrapper, HostMemSpace>>
    HydrostaticEquilibriumBox<ClusterGravity, box::ACCEPTEntropyProfile>::
        generate_P_rho_profile<
            Kokkos::View<parthenon::Real *, LayoutWrapper, HostMemSpace>,
            parthenon::UniformCartesian>(parthenon::IndexRange, parthenon::IndexRange,
                                         parthenon::IndexRange,
                                         parthenon::UniformCartesian) const;

template HydrostaticEquilibriumBox<ClusterGravity, box::ACCEPTEntropyProfile>::PRhoProfile<
    Kokkos::View<parthenon::Real *, LayoutWrapper, HostMemSpace>>
HydrostaticEquilibriumBox<ClusterGravity, ACCEPTEntropyProfile>::
    generate_P_rho_profile<Kokkos::View<parthenon::Real *, LayoutWrapper, HostMemSpace>>(
        const parthenon::Real, const parthenon::Real, const unsigned int) const;
#endif

} // namespace box
