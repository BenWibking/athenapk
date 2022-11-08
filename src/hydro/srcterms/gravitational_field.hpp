//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2021, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gravitational_field.hpp
//  \brief Defines GravitationalFieldSrcTerm
// GravitationalFieldSrcTerm is templated function to apply an arbitrary
// gravitational field as a source term
//========================================================================================
#ifndef HYDRO_SRCTERMS_GRAVITATIONAL_FIELD_HPP_
#define HYDRO_SRCTERMS_GRAVITATIONAL_FIELD_HPP_

// Parthenon headers
#include <interface/mesh_data.hpp>
#include <interface/variable_pack.hpp>
#include <mesh/domain.hpp>
#include <mesh/meshblock_pack.hpp>

// AthenaPK headers
#include "../../main.hpp"

namespace cluster {

template <typename GravitationalField>
void GravitationalFieldSrcTerm(parthenon::MeshData<parthenon::Real> *md,
                               const parthenon::Real beta_dt,
                               GravitationalField gravitationalField) {
  using parthenon::IndexDomain;
  using parthenon::IndexRange;
  using parthenon::Real;

  // Grab some necessary variables
  const auto &prim_pack = md->PackVariables(std::vector<std::string>{"prim"});
  const auto &cons_pack = md->PackVariables(std::vector<std::string>{"cons"});
  IndexRange ib = md->GetBlockData(0)->GetBoundsI(IndexDomain::interior);
  IndexRange jb = md->GetBlockData(0)->GetBoundsJ(IndexDomain::interior);
  IndexRange kb = md->GetBlockData(0)->GetBoundsK(IndexDomain::interior);

  parthenon::par_for(
      DEFAULT_LOOP_PATTERN, "GravitationalFieldSrcTerm", parthenon::DevExecSpace(), 0,
      cons_pack.GetDim(5) - 1, kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
      KOKKOS_LAMBDA(const int &b, const int &k, const int &j, const int &i) {
        auto &cons = cons_pack(b);
        auto &prim = prim_pack(b);
        const auto &coords = cons_pack.GetCoords(b);

        const Real x = std::abs(coords.x3v(i));
        const Real y = std::abs(coords.x3v(j));
        const Real z = std::abs(coords.x3v(k));

        const Real g_x = 0; // gravitationalField.g_from_x(x);
        const Real g_y = 0; // gravitationalField.g_from_y(y);
        const Real g_z = gravitationalField.g_from_z(z);

        // Apply g as a source term
        const Real dens = prim(IDN, k, j, i);
        Real p1 = cons(IM1, k, j, i);
        Real p2 = cons(IM2, k, j, i);
        Real p3 = cons(IM3, k, j, i);
        const Real KE0 = 0.5 * p1*p1 + p2*p2 + p3*p3 / dens;

        // compute new momenta
        p1 += beta_dt * dens * g_x;
        p2 += beta_dt * dens * g_y;
        p3 += beta_dt * dens * g_z;
        const Real KE1 = 0.5 * p1*p1 + p2*p2 + p3*p3 / dens;
        const Real dKE = KE1 - KE0;

        cons(IM1, k, j, i) = p1;
        cons(IM2, k, j, i) = p2;
        cons(IM3, k, j, i) = p3;
        cons(IEN, k, j, i) += dKE;
      });
}

} // namespace cluster

#endif // HYDRO_SRCTERMS_GRAVITATIONAL_FIELD_HPP_
