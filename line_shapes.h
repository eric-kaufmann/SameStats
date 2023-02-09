#include <vector>
#include <string>
#include <math.h>
#include <iostream>

std::vector<std::vector<double>> get_points_by_lines(std::vector<std::vector<std::vector<double>>> lines, int num_samples);
std::vector<std::vector<double>> get_points_for_shape(std::string shape, int num_samples);
std::vector<std::vector<double>> get_datasaurus_data();