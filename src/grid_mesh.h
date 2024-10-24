//
// Created by Charles Du on 1/15/24.
//
#pragma once

#include <vector>
#include <string>
#include <cassert>
#include "3rd/mtet/mtet.h"
#include <nlohmann/json.hpp>

namespace grid_mesh {

    enum GridStyle {
        TET5, // 5 tetrahedrons per grid cell
        TET6  // 6 tetrahedrons per grid cell
    };

    mtet::MTetMesh generate_tet_mesh(const std::array<size_t, 3> &resolution,
                                     const std::array<double, 3> &bbox_min,
                                     const std::array<double, 3> &bbox_max,
                                     GridStyle style = TET5);

    mtet::MTetMesh generate_from_kuhn_mesh(const std::array<size_t, 3> &resolution,
                                     const std::array<double, 3> &bbox_min,
                                     const std::array<double, 3> &bbox_max,
                                     GridStyle style = TET6);
    // load tet mesh from json file
    mtet::MTetMesh load_tet_mesh(const std::string &filename) ;


} // namespace grid_mesh


