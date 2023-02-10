#ifndef UTILS_H
#define UTILS_H


#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <omp.h>
#include <chrono>

void print_matrix(std::vector<std::vector<double>> data); 
void data_to_csv(std::vector<std::vector<double>> data, std::string filename);
void generate_scatter_plot(std::string filename);
void calc_stats(std::vector<std::vector<double>>* data, std::vector<double>* stats);
double minDist(const std::vector<double> v1, const std::vector<std::vector <double>> shape, int shape_size);
bool check_stats(std::vector<double>* stats1, std::vector<double>* stats2, double error);
std::vector<std::vector<double>> transpose_data(std::vector<std::vector<double>> data);








#endif