//========================================================================================
// AthenaPK - a performance portable block structured AMR astrophysical MHD code.
// Copyright (c) 2020-2022, Athena-Parthenon Collaboration. All rights reserved.
// Licensed under the BSD 3-Clause License (the "LICENSE").
//========================================================================================

#include <parthenon/parthenon.hpp>

#include "ascent.hpp"

void RenderAscent(parthenon::MeshData<parthenon::Real> *md) {
    // call ascent to render the data at this timestep
    CreateBlueprintForAscent(md);
    
}
