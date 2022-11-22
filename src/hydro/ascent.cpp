//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2020-2022, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================

#include <parthenon/parthenon.hpp>

#include "ascent.hpp"

void RenderAscent(parthenon::MeshData<parthenon::Real> *md) {
  // call ascent to render the data at this timestep
  CreateBlueprintFromMeshBlock(md);
}

void CreateBlueprintFromMeshBlock(parthenon::MeshData<parthenon::Real> *md) {
  Kokkos::Profiling::pushRegion("CreateBlueprintFromMeshBlock");
  
  auto pm = md->GetParentPointer();
  parthenon::IndexRange ib = md->GetBoundsI(parthenon::IndexDomain::interior);
  parthenon::IndexRange jb = md->GetBoundsJ(parthenon::IndexDomain::interior);
  parthenon::IndexRange kb = md->GetBoundsK(parthenon::IndexDomain::interior);

  std::vector<std::string> vars({"U"});
  auto &v = md->PackVariables(vars);
  const int ndim = pm->ndim;
}
