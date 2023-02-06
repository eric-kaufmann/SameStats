#include <iostream>
#include <vector>
#include <math.h>


void calc_stats(std::vector<std::vector<double>>* data, std::vector<double>* stats){

    double x_mean = 0;
    double y_mean = 0;
    double x_std = 0;
    double y_std = 0;
    double pearson = 0;

    for(const auto point : *data){
        x_mean += point[0];
        y_mean += point[1];
    }  
    x_mean = x_mean / data->size();
    y_mean = y_mean / data->size();

    for(const auto point : *data){
        x_std += pow(point[0] - x_mean, 2);
        y_std += pow(point[1] - y_mean, 2);
    }  
    x_std = sqrt(x_std / data->size());
    y_std = sqrt(y_std / data->size());

    for(const auto point : *data){
        pearson += (point[0] - x_mean)*(point[1] - y_mean) / (data->size() * x_std * y_std);
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

int main(int argc, char const *argv[])
{
    // Initialize the data using a 2D vector (x,y)
    std::vector<std::vector<double>> init_data = {
        {1000.0, 2.0},
        {2.0, 0.0},
        {0.0, -2.0},
        {-2.0, 0.0}
    };

    // Initialize the target data using a 2D vector (x,y)
    std::vector<std::vector<double>> target_data = {
        {-2.0, -2.0},
        {-2.0, 2.0},
        {2.0, 2.0},
        {2.0, -2.0}
    };

    // Initialize the test data using a 2D vector (x,y)
    std::vector<std::vector<double>> test_data = {
        {1.0, 2.0},
        {3.0, 4.0},
        {5.0, -2.0},
        {-2.0, 7.0}
    };

    std::vector<double> init_stats(5);

    calc_stats(&init_data, &init_stats);


    return 0;
}
