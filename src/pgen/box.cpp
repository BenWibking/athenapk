//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2021, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file box.cpp
//  \brief Idealized galaxy cluster problem generator
//
// Setups up an idealized galaxy cluster box in hydrostatic equilbrium
// with a nearly-constant vertical gravitational acceleration.
// Includes tabular cooling and "magic" heating.
//========================================================================================

// C headers

// C++ headers
#include <algorithm> // min, max
#include <cmath>     // sqrt()
#include <cstdio>    // fopen(), fprintf(), freopen()
#include <iostream>  // endl
#include <limits>
#include <sstream>   // stringstream
#include <stdexcept> // runtime_error
#include <string>    // c_str()

// Parthenon headers
#include "mesh/mesh.hpp"
#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>

// AthenaPK headers
#include "../hydro/hydro.hpp"
#include "../hydro/srcterms/gravitational_field.hpp"
#include "../main.hpp"

// Cluster headers
#include "cluster/box_gravity.hpp"
#include "cluster/entropy_profiles_cartesian.hpp"
#include "cluster/hydrostatic_equilibrium_box.hpp"

namespace box {
using namespace parthenon::driver::prelude;
using namespace parthenon::package::prelude;

void ClusterSrcTerm(MeshData<Real> *md, const parthenon::SimTime, const Real beta_dt) {
  auto hydro_pkg = md->GetBlockData(0)->GetBlockPointer()->packages.Get("Hydro");

  const bool &gravity_srcterm = hydro_pkg->Param<bool>("gravity_srcterm");

  if (gravity_srcterm) {
    const ClusterGravity &cluster_gravity =
        hydro_pkg->Param<ClusterGravity>("cluster_gravity");

    cluster::GravitationalFieldSrcTerm(md, beta_dt, cluster_gravity);
  }
}

//========================================================================================
//! \fn void InitUserMeshData(Mesh *mesh, ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void ProblemGenerator(MeshBlock *pmb, parthenon::ParameterInput *pin) {
  auto hydro_pkg = pmb->packages.Get("Hydro");
  if (pmb->lid == 0) {
    /************************************************************
     * Read Uniform Gas
     ************************************************************/

    const bool init_uniform_gas =
        pin->GetOrAddBoolean("problem/cluster", "init_uniform_gas", false);
    hydro_pkg->AddParam<>("init_uniform_gas", init_uniform_gas);

    if (init_uniform_gas) {
      const Real uniform_gas_rho = pin->GetReal("problem/cluster", "uniform_gas_rho");
      const Real uniform_gas_ux = pin->GetReal("problem/cluster", "uniform_gas_ux");
      const Real uniform_gas_uy = pin->GetReal("problem/cluster", "uniform_gas_uy");
      const Real uniform_gas_uz = pin->GetReal("problem/cluster", "uniform_gas_uz");
      const Real uniform_gas_pres = pin->GetReal("problem/cluster", "uniform_gas_pres");

      hydro_pkg->AddParam<>("uniform_gas_rho", uniform_gas_rho);
      hydro_pkg->AddParam<>("uniform_gas_ux", uniform_gas_ux);
      hydro_pkg->AddParam<>("uniform_gas_uy", uniform_gas_uy);
      hydro_pkg->AddParam<>("uniform_gas_uz", uniform_gas_uz);
      hydro_pkg->AddParam<>("uniform_gas_pres", uniform_gas_pres);
    }

    /************************************************************
     * Read Cluster Gravity Parameters
     ************************************************************/

    // Build cluster_gravity object
    ClusterGravity cluster_gravity(pin);
    hydro_pkg->AddParam<>("cluster_gravity", cluster_gravity);

    // Include gravity as a source term during evolution
    const bool gravity_srcterm = pin->GetBoolean("problem/cluster", "gravity_srcterm");
    hydro_pkg->AddParam<>("gravity_srcterm", gravity_srcterm);

    /************************************************************
     * Read Initial Entropy Profile
     ************************************************************/

    // Build entropy_profile object
    box::ACCEPTEntropyProfile entropy_profile(pin);

    /************************************************************
     * Build Hydrostatic Equilibrium Sphere
     ************************************************************/

    HydrostaticEquilibriumBox hse_box(pin, cluster_gravity, entropy_profile);
    hydro_pkg->AddParam<>("hydrostatic_equilibrium_box", hse_box);
  }

  IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);

  // initialize conserved variables
  auto &rc = pmb->meshblock_data.Get();
  auto &u_dev = rc->Get("cons").data;
  auto &coords = pmb->coords;

  // Initialize the conserved variables
  auto u = u_dev.GetHostMirrorAndCopy();

  // Get Adiabatic Index
  const Real gam = pin->GetReal("hydro", "gamma");
  const Real gm1 = (gam - 1.0);

  const auto &init_uniform_gas = hydro_pkg->Param<bool>("init_uniform_gas");
  if (init_uniform_gas) {
    const Real rho = hydro_pkg->Param<Real>("uniform_gas_rho");
    const Real ux = hydro_pkg->Param<Real>("uniform_gas_ux");
    const Real uy = hydro_pkg->Param<Real>("uniform_gas_uy");
    const Real uz = hydro_pkg->Param<Real>("uniform_gas_uz");
    const Real pres = hydro_pkg->Param<Real>("uniform_gas_pres");

    const Real Mx = rho * ux;
    const Real My = rho * uy;
    const Real Mz = rho * uz;
    const Real E = rho * (0.5 * (ux * ux + uy * uy + uz * uz) + pres / (gm1 * rho));

    for (int k = kb.s; k <= kb.e; k++) {
      for (int j = jb.s; j <= jb.e; j++) {
        for (int i = ib.s; i <= ib.e; i++) {

          u(IDN, k, j, i) = rho;
          u(IM1, k, j, i) = Mx;
          u(IM2, k, j, i) = My;
          u(IM3, k, j, i) = Mz;
          u(IEN, k, j, i) = E;
        }
      }
    }

  } else {
    /************************************************************
     * Initialize a HydrostaticEquilibriumBox
     ************************************************************/
    const auto &he_box =
        hydro_pkg
            ->Param<HydrostaticEquilibriumBox<ClusterGravity, box::ACCEPTEntropyProfile>>(
                "hydrostatic_equilibrium_box");

    const auto P_rho_profile = he_box.generate_P_rho_profile<
        Kokkos::View<parthenon::Real *, parthenon::LayoutWrapper,
                     parthenon::HostMemSpace>,
        parthenon::UniformCartesian>(ib, jb, kb, coords);

    // initialize conserved variables
    for (int k = kb.s; k <= kb.e; k++) {
      for (int j = jb.s; j <= jb.e; j++) {
        for (int i = ib.s; i <= ib.e; i++) {

          // Calculate radius
          const Real z = std::abs(coords.x3v(k));

          // Get pressure and density from generated profile
          const Real P_z = P_rho_profile.P_from_z(z);
          const Real rho_z = P_rho_profile.rho_from_z(z);

          // Fill conserved states, 0 initial velocity
          u(IDN, k, j, i) = rho_z;
          u(IM1, k, j, i) = 0.0;
          u(IM2, k, j, i) = 0.0;
          u(IM3, k, j, i) = 0.0;
          u(IEN, k, j, i) = P_z / gm1;
        }
      }
    }
  }

  // copy initialized cons to device
  u_dev.DeepCopy(u);
}

} // namespace box
