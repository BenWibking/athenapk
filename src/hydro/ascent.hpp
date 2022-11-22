#ifndef ASCENT_HPP_
#define ASCENT_HPP_
//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2020-2022, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================

#include "interface/mesh_data.hpp"

void RenderAscent(parthenon::MeshData<parthenon::Real> *md);

void CreateBlueprintFromMeshBlock(parthenon::MeshData<parthenon::Real> *md);


#endif // ASCENT_HPP_