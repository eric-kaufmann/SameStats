#ifndef LINE_SHAPES_H
#define LINE_SHAPES_H

#include <vector>
#include <iostream>
#include <random>

std::vector<std::vector<double>> get_points_by_lines(std::vector<std::vector<std::vector<double>>> lines, int num_samples);
std::vector<std::vector<double>> get_points_for_shape(std::string shape, int num_samples);
std::vector<std::vector<double>> get_datasaurus_data();
std::vector<std::vector<double>> generate_point_cloud(int num_samples);
double easeInOutQuad(double n);

#endif