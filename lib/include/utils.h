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



void data_to_csv(std::vector<std::vector<double>> data, std::string filename);
//writes data (2d vector) into a csv file in the same directory

void generate_scatter_plot(std::string filename);
//generates a scatter plot using python with matplotlib from csv data 

void calc_stats(std::vector<std::vector<double>>* data, std::vector<double>* stats);
//calculates the statistcs of data: mu_x, mu_y, sigma_x, sigma_y, corr_xy

void calc_partitionedStats(std::vector<std::vector<double>> data, std::vector<double>* stats, int num_points, double global_Xmean, double global_Ymean, double global_Xstd, double global_Ystd); //num_points: globale anzahl datenpunkte
//calculates the contribution of a subgroup to the global statistics of the initial data set. 

double minDist(const std::vector<double> v1, const std::vector<std::vector <double>> shape, int shape_size);
//calculates the smallest distance of a point to all points of a pointcloud

bool check_stats(std::vector<double>* stats1, std::vector<double>* stats2, double error);
//checks if the difference of two statistics is within a specified error bound 

bool check_Partitionedstats(std::vector<double>* stats1, std::vector<double>* stats2, double error);
//checks if the difference of two statistic contributions of a subgroup is within a specified error bound

std::vector<std::vector<double>> transpose_data(std::vector<std::vector<double>> data);
//transposes a 2d vector from {{x1,y1},{x2,y2},{x3.y3},...,{xn,yn}} to {{x1,x2,x3,...,xn},{y1,y2,y3,...,yn}}

std::vector<std::vector<std::vector<double>>> partitionData(int num_threads, int padding, std::vector<std::vector<double>> pointCloud);
//creates thread-wise groups of data which are seperated by a specified padding size. By that, false sharing should be prevented.
//num_threads = z, pointcloud_size = n, padding_size = 2
//{{x1,x2,x3,...,xn},{y1,y2,y3,...,yn}} -------> {{{0,0,x11,x12,x13,...,x1n/z},{y11,y12,y13,...,y1n/z,0,0}}, {{0,0,x21,x22,x23,...,x2n/z},{y21,y22,y23,...,y2n/z,0,0}}},...,{{0,0,xz1,xz2,xn3,...,xzn/z},{yz1,yz2,yz3,...,yzn/z,0,0}}}}

std::vector<std::vector<double>> refactorData(std::vector<std::vector<std::vector<double>>> new_working_data, std::vector<std::vector<double>> working_data, int num_threads, int padding);
//merge the thread-wise partitioned data into one data set of the initial form

int count_nonzeros(std::vector<std::vector<double>> working_data, int padding);
//counts how many elements of the data are already considered

double get_MaxError(std::vector<double> init_stats, std::vector<double> final_stats);
//calculates the maximum error between two given statistics


#endif