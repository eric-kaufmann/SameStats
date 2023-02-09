#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <string>
#include <random>
#include "line_shapes.h"


void data_to_csv(std::vector<std::vector<double>>* data, std::string filename)
{
    std::ofstream outfile;
    outfile.open(filename);

    for(const auto point : *data){
        outfile << point[0] << "," << point[1] << "\n";
    } 

    outfile.close();
}

void generate_scatter_plot(){
    system("python3 scatter_plot.py");
}


void calc_stats(std::vector<std::vector<double>>* data, std::vector<double>* stats){

    double x_mean = 0;
    double y_mean = 0;
    double x_std = 0;
    double y_std = 0;
    double pearson = 0;

    int num_points = (data->at(0)).size();

    for(int i=0; i<num_points; ++i){
        x_mean += (data->at(0))[i];
        y_mean += (data->at(1))[i];
    } 
    x_mean = x_mean / num_points;
    y_mean = y_mean / num_points;

    for(int i=0; i<num_points; ++i){
        x_std += pow((data->at(0))[i] - x_mean, 2);
        y_std += pow((data->at(1))[i] - x_mean, 2);
    } 
    x_std = sqrt(x_std / num_points);
    y_std = sqrt(y_std / num_points);


    for(int i=0; i<num_points; ++i){
        pearson += ((data->at(0))[i] - x_mean)*((data->at(1))[i] - y_mean) / (num_points * x_std * y_std);
    } 

    stats->at(0) = x_mean;
    stats->at(1) = y_mean;
    stats->at(2) = x_std;
    stats->at(3) = y_std;
    stats->at(4) = pearson;

    std::cout << "x_mean: " << stats->at(0) << std::endl;
    std::cout << "y_mean: " << stats->at(1) << std::endl;
    std::cout << "x_std: " << stats->at(2) << std::endl;
    std::cout << "y_std: " << stats->at(3) << std::endl;
    std::cout << "pearson: " << stats->at(4) << std::endl;
}

//shape daten in 2d vector: 1. vektor mit x koordinaten, 2. vektor mit y koordinaten 
//shape_size: anzahl punkte
double minDist(const std::vector<double> v1, const std::vector<std::vector <double>> shape, int shape_size){ 
    double result=0.0, x_diff, y_diff;
    double x_vec, y_vec;
    double buff[shape_size];
#pragma omp parallel for 
    for(int i = 0; i < shape_size; i++){
        x_diff = std::pow((v1[0] - shape[0][i]),2); //quadratischer abstand
        y_diff = std::pow((v1[1] - shape[1][i]),2); //quadratischer abstand
        buff[i] = std::sqrt(x_diff+y_diff); //norm
    }
#pragma omp parallel for reduction(min:result)
    for(int i = 0; i < shape_size; i++){
        // result = std::min(result, buff[i]);
        if(buff[i]<result)
            result = buff[i];   
    }
    return result;    
}


bool check_stats(std::vector<double>* stats1, std::vector<double>* stats2, double error){
    for(int i = 0; i < stats1->size(); ++i){
        if( (stats1->at(i) - stats2->at(i)) > error ){
            return false;
        }
    }
    return true;
}

int main(int argc, char const *argv[]){

    bool save_data = true; // generates csv files and plots data

    int num_steps = 10; // number of iteration steps
    double error = 10e-2; // maximum statistical diffence
    double max_shift = 10e-3; // std of rand-norm for data point shift
    int num_samples = 100000; // #datapoints if a random point cloud is generated
    
    // set data format by transposing data
    // true -> {{x1,y1}, {x2,y2}, ...} :: false -> {{x1, x2, ...}, {y1, y2, ...}}
    bool transpose = true;
    
    double init_temperature = 1.0; // initial temperature (gets smaller every iteration)
    double min_temperature = 0.0; // minimal temperature (always temperature > min_temperature)
   
    // All possible init shapes 
    // {"datasaurus", "random"};
    std::string init_shape = "datasaurus";
    std::vector<std::vector<double>> working_data;
    if (init_shape == "datasaurus"){
        working_data = get_datasaurus_data();
    } 
    else if (init_shape == "random"){
        working_data = generate_point_cloud(num_samples);
    }

    // All possible target shapes 
    // {"x", "h_lines", "v_lines", "wide_lines", "high_lines", "slant_up", "slant_down", "center", "star", "down_parab"};
    std::string target_shape = "x";
    std::vector<std::vector<double>> target_data = get_points_for_shape(target_shape, working_data.size());

    // Transpose data
    if(transpose){
        working_data = transpose_data(working_data);
        target_data = transpose_data(target_data);
    }

    // // target points
    // std::vector<std::vector<double>> target_data = {
    //     {1.0, 2.0, -1.0, 0.0}, // x values
    //     {1.0, 2.0, -1.0, 0.0} // y values
    // };

    // Initialize the data using a 2D vector (x,y)


    // std::vector<std::vector<double>> working_data = {
    //     {3.0, 2.0, 0.0, 1.0}, // x values
    //     {2.0, 0.0, -3.0, -2.0} // y values
    // };


    // Calculate initial stats
    std::vector<double> init_stats(5); 
    std::vector<double> working_stats(5);
    std::vector<double> target_stats(5);

    calc_stats(&working_data, &init_stats); // set initial stats
    calc_stats(&target_data, &target_stats); // set target stats


    // create radom number generator
    std::default_random_engine generator;

    // init distributions
    std::uniform_int_distribution<int> rand_point_dist(0,working_data[0].size()); // selects a rondom point idx
    std::normal_distribution<double> rand_shift(0.0,max_shift); // generates shift for points
    std::uniform_real_distribution<> rand_break(0.0, 1.0); // gives a chance to break out of loop

    double temperature = init_temperature;

    // Main loop over num_steps steps
    for(int step=0; step<num_steps; ++step){
        // get random point
        int rand_point_idx = rand_point_dist(generator);
        std::cout << "random point " << working_data[0][rand_point_idx] << "  " << working_data[1][rand_point_idx] << std::endl;
        double rpx = working_data[0][rand_point_idx]; // x value of random point
        double rpy = working_data[1][rand_point_idx]; // y value of random point
        
        // TODO: calculate initial min dist from (rpx, rpy) to target_data
        double init_min_dist = 0.3; 

        // Get new point with smaller distance
        double new_min_dist = INFINITY;
        double rpx_shift; // x value of random point with shift
        double rpy_shift; // y value of random point with shift
        while(new_min_dist > init_min_dist){
            // shift point
            rpx_shift = rpx + rand_shift(generator);
            rpy_shift = rpy + rand_shift(generator);
            std::cout << "random shift " << rpx_shift << "  " << rpy_shift << std::endl;
            std::cout <<  std::endl;
            
            // TODO: calculate new min dist from (rpx_shift, rpy_shift) to target_data
            new_min_dist = 0.001;
            
            // simmulated annealing: break out of loop with a specific chance
            if(rand_break(generator)<temperature){
                std::cout << "BREAK" << std::endl;
                break;
            }
        }

        // embed new point into vector
        working_data[0][rand_point_idx] = rpx_shift;
        working_data[1][rand_point_idx] = rpy_shift;

        // calc new stats
        calc_stats(&working_data, &working_stats);

        // check if statistics are ca. the same to our inital stats
        if(!check_stats(&working_stats, &init_stats, error)){
            // not the same -> return points to old value
            working_data[0][rand_point_idx] = rpx;
            working_data[1][rand_point_idx] = rpy;
        }

        // Adjust temperature for simulated annealing
        if(init_temperature/(step+1) > min_temperature){
            temperature = init_temperature/(step+1);
        }
    }

    if(save_data){
        data_to_csv(&working_data, "generated_data.csv");
        data_to_csv(&working_data, "target_data.csv");
        generate_scatter_plot();
    }

    return 0;
}
