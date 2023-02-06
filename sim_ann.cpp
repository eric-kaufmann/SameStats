#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <string>
#include <random>


void matrixToCSV(std::vector<std::vector<double>>* data, std::string filename)
{
    std::ofstream outfile;
    outfile.open(filename);

    for(const auto point : *data){
        outfile << point[0] << "," << point[1] << "\n";
    } 

    outfile.close();
}

void generateScatterPlots(){
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

bool check_stats(std::vector<double>* stats1, std::vector<double>* stats2, double error){
    for(int i = 0; i < stats1->size(); ++i){
        if( (stats1->at(i) - stats2->at(i)) > error ){
            return false;
        }
    }
    return true;
}

/***
 * Generate a point cloud based on given statistics
*/
void generate_point_cloud(std::vector<double>* stats){
    
}

int main(int argc, char const *argv[]){

    int num_steps = 10;
    double error = 10e-2;
    double max_shift = 10e-3;
    double temperature = 10e-4;

    // target points
    std::vector<std::vector<double>> target_data = {
        {1.0, 2.0, -1.0, 0.0}, // x values
        {1.0, 2.0, -1.0, 0.0} // y values
    };

    // calc stats of target points
    std::vector<double> target_stats(5);

    // Initialize the data using a 2D vector (x,y)
    std::vector<std::vector<double>> working_data = {
        {3.0, 2.0, 0.0, 1.0}, // x values
        {2.0, 0.0, -3.0, -2.0} // y values
    };

    // // data for checking
    // std::vector<std::vector<double>> test_data = {
    //     {5.0, -2.0, 1.5, 1.0}, // x values
    //     {1.5, 3.0, -2.0, -1.5} // y values
    // };


    std::vector<double> init_stats(5); // initial stats
    std::vector<double> working_stats(5); // stats of p

    calc_stats(&working_data, &init_stats);
    calc_stats(&target_data, &target_stats);

    // std::cout << "stats equal: " << check_stats(&init_stats, &test_stats, 0.5) << std::endl;

    std::default_random_engine generator;
    std::uniform_int_distribution<int> rand_point_dist(0,working_data[0].size());

    // do num_steps iterations
    for(int step=0; step<num_steps; ++step){
        // get random point
        int rand_point_idx = rand_point_dist(generator);
        std::cout << "random point " << working_data[0][rand_point_idx] << "  " << working_data[1][rand_point_idx] << std::endl;
        double rpx = working_data[0][rand_point_idx]; // x value of random point
        double rpy = working_data[1][rand_point_idx]; // y value of random point
        
        // TODO: calculate initial min dist
        double init_min_dist = 0.3;

        // init normal distributions
        std::normal_distribution<double> rand_shift(0.0,max_shift); // generates shift for points
        std::uniform_real_distribution<> rand_break(0.0, 1.0); // gives a chance to break out of loop

        double new_min_dist = INFINITY;
        double rpx_shift; // x value of random point with shift
        double rpy_shift; // y value of random point with shift
        while(new_min_dist > init_min_dist){
            // shift point
            rpx_shift = rpx + rand_shift(generator);
            rpy_shift = rpy + rand_shift(generator);
            std::cout << "random shift " << rpx_shift << "  " << rpy_shift << std::endl;
            std::cout <<  std::endl;
            
            // TODO: calculate new min dist from (rpx_shift, rpy_shift) to other points
            new_min_dist = 0.001;
            
            // break even if new_min_dist is still bigger
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

        // check if statistics are to an error the same
        if(!check_stats(&working_stats, &target_stats, error)){
            // not equal -> return points to old value
            working_data[0][rand_point_idx] = rpx;
            working_data[1][rand_point_idx] = rpy;
        }
    }

    return 0;
}
