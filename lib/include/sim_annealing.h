#ifndef SIM_ANNEALING_H
#define SIM_ANNEALING_H

#include <string>

double sim_annealing(std::string shape, int threads, int serial_steps); 
//main function which shifts a pointcloud into a specified shape while maintainig the same statistics of the initial pointcloud
#endif