#pragma once

#include "MshSpec.h"

#include <ostream>

namespace mshio {

void save_mesh_format(std::ostream& out, const MshSpec& spec);

}
