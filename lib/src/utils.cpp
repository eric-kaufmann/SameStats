#include "utils.h"


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
    for(int i = 0; i < (int)data[0].size(); i++){
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

}

void calc_partitionedStats(std::vector<std::vector<double>> data, std::vector<double>* stats, int num_points, double global_Xmean, double global_Ymean, double global_Xstd, double global_Ystd){ //num_points: globale anzahl datenpunkte

    double x_mean = 0;
    double y_mean = 0;
    double x_std = 0;
    double y_std = 0;
    double pearson = 0;
    for(auto &it : data[0]){ //x values
        x_mean += it;
    }
    x_mean = pow(x_mean/num_points,2); //pow because we compare with sqrt(error)

    for(auto &it : data[1]){ //y values
        y_mean += it;
    }
    y_mean = pow(y_mean/num_points,2); //pow because we compare with sqrt(error)

    for(auto &it : data[0]){ //x values
        if(it != 0.0) //filter zeros
            x_std += pow(it - global_Xmean, 2); 
    }
    x_std = x_std/num_points; //no sqrt() because we compare with sqrt(error)  

    for(auto &it : data[1]){ //y values
        if(it != 0.0) //filter zeros
            y_std += pow(it - global_Ymean, 2); 
    }
    y_std = y_std/num_points; //no sqrt() because we compare with sqrt(error)

    int local_points = data[0].size();
    if(data[0].size() != data[1].size())
        std::cout << "error: invalid format (in function calc_partitionedStats)" << std::endl;
    for(int i = 0; i < local_points; i++){
        if(data[0][i] != 0.0 && data[1][i] != 0.0) //filter zeros
            pearson += (data[0][i] - global_Xmean)*(data[1][i] - global_Ymean);
    }
    pearson = pearson / (num_points * global_Xstd * global_Ystd); 

    stats->at(0) = x_mean;
    stats->at(1) = y_mean;
    stats->at(2) = x_std;
    stats->at(3) = y_std;
    stats->at(4) = pearson;

}


//shape daten in 2d vector: 1. vektor mit x koordinaten, 2. vektor mit y koordinaten 
//shape_size: anzahl punkte
double minDist(const std::vector<double> v1, const std::vector<std::vector<double>> shape, int shape_size){ 
    double result=INFINITY;
    double x_diff, y_diff;
    double buff[shape_size];

    for(int i = 0; i < shape_size; i++){
        x_diff = std::pow((v1[0] - shape[0][i]),2); //quadratic difference
        y_diff = std::pow((v1[1] - shape[1][i]),2); //quadratic difference
        buff[i] = std::sqrt(x_diff+y_diff); //norm
    }

    for(int i = 0; i < shape_size; i++){

        if(buff[i]<result)
            result = buff[i];   
    }
    return result;    
}


bool check_stats(std::vector<double>* stats1, std::vector<double>* stats2, double error){
    for(int i = 0; i < (int)stats1->size(); ++i){
        if(fabs( (stats1->at(i) - stats2->at(i)) ) > error ){
            return false;
        }
    }
    return true;
}

bool check_Partitionedstats(std::vector<double>* stats1, std::vector<double>* stats2, double error){
    for(int i = 0; i < (int)stats1->size()-1; ++i){
        if(sqrt(fabs( stats1->at(i) - stats2->at(i)) ) > error ){
            return false;
        }
    }
    if(fabs(stats1->at(4) - stats2->at(4)) > pow(error,2)){ //we compare with sqrt(error) --> to obtain error we need to take error^2
        return false;
    }
    return true;
}

//partition data thread-wise, insert padding into data to prevent false sharing
std::vector<std::vector<std::vector<double>>> partitionData(int num_threads, int padding, std::vector<std::vector<double>> pointCloud){ //pointCloud is assumed to be transposed
    int size_x = pointCloud[0].size();
    int size_y = pointCloud[1].size();
    if(size_x!=size_y){
        std::cout << "error: pointcloud invalid format" << std::endl;
    }
    int num_points = size_x/num_threads;
    if(size_x%num_threads != 0){
        num_points++; //wenn es einen rest gibt, soll jeder vektor um 1 vergrößert werden -> genug platz für alle punkte
    }
    std::vector<std::vector<std::vector<double>>> paddingData(num_threads, std::vector<std::vector<double>>(2, std::vector<double>(num_points + padding, 0))); //3d: each thread gets a 2d vector, initialized with zeros 
    //TODO originale daten in gepaddeten vektor laden
    int j = 0;
    int elements = 0;
    for(auto &element : paddingData){ //loop over all 2d vectors
        for(int i = padding; i < num_points+padding; i++){ //fill new vector with data from old vector
            element[0][i] = pointCloud[0][i-padding+j*num_points]; // x values: {0,0,0,0,x11,x12,...,x1n}
            element[1][i-padding] = pointCloud[1][i-padding+j*num_points]; // y values: {y11,y12,...,y1n,0,0,0,0}
            elements++;
            if(elements == size_x) //müssen abbrechen wenn alle elemente aufgeteilt wurden --> geht nicht immer auf, daher manuelles ausbrechen aus loop
                break;
        }
        j++;
    }
    return paddingData;
}



int count_nonzeros(std::vector<std::vector<double>> working_data, int padding){ //debugging
    int counts = 0;
    for(auto &element : working_data[0]){
        if(element != 0){
            counts++;
        }
    }
    return (counts-2*padding);
}

std::vector<std::vector<double>> refactorData(std::vector<std::vector<std::vector<double>>> new_working_data, std::vector<std::vector<double>> working_data, int num_threads, int padding){
    int j = 0;
    int size = working_data[0].size();
    int toggle = 0;
    int vec_size = size/num_threads;
    int gap = 0;
    if(size%num_threads != 0){
        vec_size++;
        gap = vec_size - size%vec_size; 
    }

    for(auto &element : new_working_data){
        if(j==(num_threads-1)){
            toggle = gap;
        }
        for(int i = padding; i < (int)element[0].size()-toggle; i++){
            working_data[0][i-padding+j*vec_size] = element[0][i]; //x value
            working_data[1][i-padding+j*vec_size] = element[1][i-padding]; //y value
        }
        j++;
    }
    return working_data;    
}


double get_MaxError(std::vector<double> init_stats, std::vector<double> final_stats){
    double max = 0.0;
    double diff;
    for(int i = 0; i < (int)init_stats.size(); i++){
        diff = fabs(init_stats[i]-final_stats[i]);
        if(diff > max)
            max = diff;
    }
    return max;
}