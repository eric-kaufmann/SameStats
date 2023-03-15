#include "utils.h"
#include "line_shapes.h"
#include "sim_annealing.h"
#include "catch.hpp"


// TEST_CASE("x", "[correctness]"){

//         REQUIRE(sim_annealing("x") <= 0.05);
// }


TEST_CASE("star", "[correctness]"){
        REQUIRE(sim_annealing("star", 1, 1e6) <= 0.05);
}



