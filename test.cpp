#include "utils.h"
#include "line_shapes.h"
#include "sim_annealing.h"
#include "catch.hpp"


TEST_CASE("x", "[correctness]"){

        REQUIRE(sim_annealing("x", 8, 1e6) <= 0.05); // shape: "x", 8 threads, 1e6 iterations
}
//error in statistics should be less than 0.05

TEST_CASE("star", "[correctness]"){
        REQUIRE(sim_annealing("star", 8, 1e6) <= 0.05); // shape: "star", 8 threads, 1e6 iterations
}
//error in statistics should be less than 0.05



