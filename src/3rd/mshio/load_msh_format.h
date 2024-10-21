#pragma once
#include "MshSpec.h"
#include <iostream>

namespace mshio {

void load_mesh_format(std::istream& in, MshSpec& spec);

}
