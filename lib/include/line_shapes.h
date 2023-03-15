#ifndef LINE_SHAPES_H
#define LINE_SHAPES_H

#include <vector>
#include <iostream>
#include <random>

std::vector<std::vector<double>> get_points_by_lines(std::vector<std::vector<std::vector<double>>> lines, int num_samples);
//function that creates a point-wise shape from lines

std::vector<std::vector<double>> get_points_for_shape(std::string shape, int num_samples);
//function that creates pointcloud which corresponds to a specified shape

std::vector<std::vector<double>> get_datasaurus_data();
//function that creates a 'dinosaur' pointcloud like that from the autodesk paper (see README)

std::vector<std::vector<double>> generate_point_cloud(int num_samples);
//function that creates a random pointcloud from a uniform distribution

double easeInOutQuad(double n);
//smoothing function needed for temperature adjustment of simulated annealing (adopted from autodesk paper)

#endif