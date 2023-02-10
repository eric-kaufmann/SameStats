#ifndef SIM_ANN_H
#define SIM_ANN_H


#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <omp.h>
#include <chrono>


std::vector<std::vector<double>> get_points_by_lines(std::vector<std::vector<std::vector<double>>> lines, int num_samples);
std::vector<std::vector<double>> get_points_for_shape(std::string shape, int num_samples);
std::vector<std::vector<double>> get_datasaurus_data();
std::vector<std::vector<double>> generate_point_cloud(int num_samples);
void print_matrix(std::vector<std::vector<double>> data); 
std::vector<std::vector<double>> transpose_data(std::vector<std::vector<double>> data);
void data_to_csv(std::vector<std::vector<double>>* data, std::string filename);
void generate_scatter_plot();
void calc_stats(std::vector<std::vector<double>>* data, std::vector<double>* stats);
double minDist(const std::vector<double> v1, const std::vector<std::vector <double>> shape, int shape_size);
bool check_stats(std::vector<double>* stats1, std::vector<double>* stats2, double error);








#endif