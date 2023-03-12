#include "utils.h"

void print_matrix(std::vector<std::vector<double>> data){ // gesammelte x,y koordinaten falls transpose = false
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[i].size(); j++)
            std::cout << data[i][j] << " ";
        std::cout << std::endl;
    }
}

std::vector<std::vector<double>> transpose_data(std::vector<std::vector<double>> data){
    std::vector<double> x_vals;
    std::vector<double> y_vals;

    for (auto p : data){
        x_vals.push_back(p[0]);
        y_vals.push_back(p[1]);
    }
    std::vector<std::vector<double>> data_T;
    data_T.push_back(x_vals);
    data_T.push_back(y_vals);

    return data_T;
}


void data_to_csv(std::vector<std::vector<double>> data, std::string filename) //assuming data in format {{x1, x2, ...}, {y1, y2, ...}}
{
    std::ofstream outfile;
    outfile.open(filename, std::ios::trunc);
    //std::cout << "data[0].size() " << data[0].size() << " data[1].size() " << data[1].size() << std::endl;
    for(int i = 0; i < data[0].size(); i++){
        outfile << data[0][i] << "," << data[1][i] << "\n";
        //std::cout << "write " << data[0][i] << "   " << data[1][i] << std::endl;
    }
    outfile.close();
}

void generate_scatter_plot(std::string filename){
    std::string sysout = "python3 scatter_plot.py " + filename;
    system(sysout.data());
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
        x_std += std::pow((data->at(0))[i] - x_mean, 2);
        y_std += std::pow((data->at(1))[i] - y_mean, 2);
    } 
    x_std = std::sqrt(x_std / (num_points-1));
    y_std = std::sqrt(y_std / (num_points-1));


    for(int i=0; i<num_points; ++i){
        pearson += ((data->at(0))[i] - x_mean)*((data->at(1))[i] - y_mean);
    } 
    pearson = pearson / (num_points * x_std * y_std);

    stats->at(0) = x_mean;
    stats->at(1) = y_mean;
    stats->at(2) = x_std;
    stats->at(3) = y_std;
    stats->at(4) = pearson;

//     std::cout << "x_mean: " << stats->at(0) << std::endl;
//     std::cout << "y_mean: " << stats->at(1) << std::endl;
//     std::cout << "x_std: " << stats->at(2) << std::endl;
//     std::cout << "y_std: " << stats->at(3) << std::endl;
//     std::cout << "pearson: " << stats->at(4) << std::endl;
}

//shape daten in 2d vector: 1. vektor mit x koordinaten, 2. vektor mit y koordinaten 
//shape_size: anzahl punkte
double minDist(const std::vector<double> v1, const std::vector<std::vector <double>> shape, int shape_size){ 
    double result=INFINITY;
    double x_diff, y_diff;
    double buff[shape_size];
// #pragma omp parallel for shared(buff)
    for(int i = 0; i < shape_size; i++){
        x_diff = std::pow((v1[0] - shape[0][i]),2); //quadratic difference
        y_diff = std::pow((v1[1] - shape[1][i]),2); //quadratic difference
        buff[i] = std::sqrt(x_diff+y_diff); //norm
    }
// #pragma omp parallel for reduction(min:result)
    for(int i = 0; i < shape_size; i++){
        // result = std::min(result, buff[i]);
        if(buff[i]<result)
            result = buff[i];   
    }
    return result;    
}


bool check_stats(std::vector<double>* stats1, std::vector<double>* stats2, double error){
    for(int i = 0; i < stats1->size(); ++i){
        if(fabs( (stats1->at(i) - stats2->at(i))) > error ){
            return false;
        }
    }
    return true;
}
