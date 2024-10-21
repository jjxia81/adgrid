#pragma once

#include "MshSpec.h"

#include <iostream>

namespace mshio {

void save_nodes(std::ostream& out, const MshSpec& spec);

}
