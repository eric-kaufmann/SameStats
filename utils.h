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
void calc_partitionedStats(std::vector<std::vector<double>> data, std::vector<double>* stats, int num_points, double global_Xmean, double global_Ymean, double global_Xstd, double global_Ystd); //num_points: globale anzahl datenpunkte
double minDist(const std::vector<double> v1, const std::vector<std::vector <double>> shape, int shape_size);
bool check_stats(std::vector<double>* stats1, std::vector<double>* stats2, double error);
bool check_Partitionedstats(std::vector<double>* stats1, std::vector<double>* stats2, double error);
std::vector<std::vector<double>> transpose_data(std::vector<std::vector<double>> data);
std::vector<std::vector<std::vector<double>>> partitionData(int num_threads, int padding, std::vector<std::vector<double>> pointCloud);
std::vector<std::vector<double>> refactorData(std::vector<std::vector<std::vector<double>>> new_working_data, std::vector<std::vector<double>> working_data, int num_threads, int padding);
void printVector(std::vector<std::vector<std::vector<double>>> working_data);
int count_nonzeros(std::vector<std::vector<double>> working_data, int padding);
void compareData(std::vector<std::vector<double>> vec1, std::vector<std::vector<double>>vec2);












#endif