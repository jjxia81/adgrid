//#define Check_Flip_Tets
//#include <mtet/mtet.h>
//#include <mtet/io.h>
//#include <ankerl/unordered_dense.h>
#include <span>
#include <queue>
#include <optional>
#include <CLI/CLI.hpp>

#include "io.h"
#include "timer.h"
#include "csg.h"
#include "grid_mesh.h"
#include "grid_refine.h"
#include "3rd/implicit_functions/implicit_functions.h"

//using namespace mtet;

int test_ad(int argc, const char *argv[]);
