#pragma once

#include "MshSpec.h"
#include <iostream>

namespace mshio {

void save_entities(std::ostream& out, const MshSpec& spec);

} // namespace mshio
