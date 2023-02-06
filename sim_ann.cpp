#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <string>


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

int main(int argc, char const *argv[])
{
    // Initialize the data using a 2D vector (x,y)
    std::vector<std::vector<double>> init_data = {
        {3.0, 2.0, 0.0, 1.0}, // x values
        {2.0, 0.0, -3.0, -2.0} // y values
    };

    std::vector<std::vector<double>> test_data = {
        {5.0, -2.0, 1.5, 1.0}, // x values
        {1.5, 3.0, -2.0, -1.5} // y values
    };

    std::vector<double> init_stats(5);
    std::vector<double> test_stats(5);

    calc_stats(&init_data, &init_stats);
    calc_stats(&test_data, &test_stats);

    std::cout << "stats equal: " << check_stats(&init_stats, &test_stats, 0.5) << std::endl;

    //matrixToCSV(&init_stats, "test");

    return 0;
}
