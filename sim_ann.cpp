#include "sim_ann.h"

/***
 * calculates distributed points on defined lines
*/
std::vector<std::vector<double>> get_points_by_lines(std::vector<std::vector<std::vector<double>>> lines, int num_samples){
    
    std::vector<std::vector<double>> new_data;

    int num_lines = lines.size();
    int samples_per_line = num_samples / num_lines;

    for(auto l : lines){
        auto p1 = l[0]; // anfangspunkt linie
        auto p2 = l[1]; //endpunkt linie

        // lam * (x1,y1) + (1-lam) * (x2,y2)
        for(double lam = 0; lam <= 1.0; lam = lam + 1.0/(samples_per_line-1)){ //generiert äquidistante punkte entlang der linie
            std::vector<double> new_point = {
                lam * p1[0] + (1.0-lam) * p2[0], //p1[0] ist eine koordinate des anfangspunktes
                lam * p1[1] + (1.0-lam) * p2[1]
            };
            new_data.push_back(new_point);
        }
    }

    return new_data;
}

std::vector<std::vector<double>> get_points_for_shape(std::string shape, int num_samples){
    
    // Define lines
    std::vector<std::vector<std::vector<double>>> lines;

    if(shape == "x"){
        lines = {
            {{20.0,0.0},{100.0,100.0}}, // one line starting at (0,0) and ending at (100,100)
            {{100.0,0.0},{20.0,100.0}}
        };
    }
    else if (shape == "h_lines"){
        lines = {
            {{0.0,10.0},{100.0,10.0}}, 
            {{0.0,30.0},{100.0,10.0}}, //y koordiante anpassen? 
            {{0.0,50.0},{100.0,10.0}}, //y koordiante anpassen? 
            {{0.0,70.0},{100.0,10.0}}, //y koordiante anpassen? 
            {{0.0,90.0},{100.0,10.0}}, //y koordiante anpassen? 
        };
    }
    else if (shape == "v_lines"){
        lines = {
            {{10.0,0.0},{10.0,100.0}},
            {{30.0,0.0},{30.0,100.0}},
            {{50.0,0.0},{50.0,100.0}},
            {{70.0,0.0},{70.0,100.0}},
            {{90.0,0.0},{90.0,100.0}},
        };
    }
    else if (shape == "wide_lines"){
        lines = {
            {{10.0,0.0},{10.0,100.0}},
            {{90.0,0.0},{90.0,100.0}},
        };
    }
    else if (shape == "high_lines"){
        lines = {
            {{0.0,10.0},{100.0,10.0}},
            {{0.0,90.0},{100.0,90.0}},
        };
    }
    else if (shape == "slant_up"){
        lines = {
            {{0.0,0.0},{100.0,100.0}},
            {{0.0,30.0},{70.0,100.0}},
            {{30.0,0.0},{100.0,70.0}},
            {{50.0,0.0},{100.0,50.0}},
            {{0.0,50.0},{50.0,100.0}},
        };
    }
    else if (shape == "slant_down"){
        lines = {
            {{0.0,100.0},{100.0,0.0}},
            {{0.0,70.0},{70.0,0.0}},
            {{30.0,100.0},{100.0,30.0}},
            {{0.0,50.0},{50.0,0.0}},
            {{50.0,100.0},{100.0,50.0}},
        };
    }
    else if (shape == "center"){
        lines = {
            {{54.26,47.83},{54.26,47.83}},
        };
    }
    else if (shape == "star"){
        lines = {
            {{28.0, 60}, {52.0, 60}},
            {{52.0, 60}, {60.0, 90}},
            {{60.0, 90}, {68.0, 60}},
            {{68.0, 60}, {92.0, 60}},
            {{92.0, 60}, {72.0, 40}},
            {{72.0, 40}, {80.0, 10}},
            {{80.0, 10}, {60.0, 30}},
            {{60.0, 30}, {40.0, 10}},
            {{40.0, 10}, {48.0, 40}},
            {{48.0, 40}, {28.0, 60}}
        };
    }
    else if (shape == "down_parab"){
        lines = {
            {{0, -66.25}, {3, -48.0625}},
            {{3, -48.0625}, {6, -31.0}},
            {{6, -31.0}, {9, -15.0625}},
            {{9, -15.0625}, {12, -0.25}},
            {{12, -0.25}, {15, 13.4375}},
            {{15, 13.4375}, {18, 26.0}},
            {{18, 26.0}, {21, 37.4375}},
            {{21, 37.4375}, {24, 47.75}},
            {{24, 47.75}, {27, 56.9375}},
            {{27, 56.9375}, {30, 65.0}},
            {{30, 65.0}, {33, 71.9375}},
            {{33, 71.9375}, {36, 77.75}},
            {{36, 77.75}, {39, 82.4375}},
            {{39, 82.4375}, {42, 86.0}},
            {{42, 86.0}, {45, 88.4375}},
            {{45, 88.4375}, {48, 89.75}},
            {{48, 89.75}, {51, 89.9375}},
            {{51, 89.9375}, {54, 89.0}},
            {{54, 89.0}, {57, 86.9375}},
            {{57, 86.9375}, {60, 83.75}},
            {{60, 83.75}, {63, 79.4375}},
            {{63, 79.4375}, {66, 74.0}},
            {{66, 74.0}, {69, 67.4375}},
            {{69, 67.4375}, {72, 59.75}},
            {{72, 59.75}, {75, 50.9375}},
            {{75, 50.9375}, {78, 41.0}},
            {{78, 41.0}, {81, 29.9375}},
            {{81, 29.9375}, {84, 17.75}},
            {{84, 17.75}, {87, 4.4375}},
            {{87, 4.4375}, {90, -10.0}},
            {{90, -10.0}, {93, -25.5625}},
            {{93, -25.5625}, {96, -42.25}},
            {{96, -42.25}, {99, -60.0625}}
        };
    }

    // generate target data
    std::vector<std::vector<double>> target_data = get_points_by_lines(lines, num_samples);

    // for(auto p : target_data){
    //     std::cout << p[0] << "  " << p[1] << std::endl;
    // }
    return target_data;
}

std::vector<std::vector<double>> get_datasaurus_data(){
    std::vector<std::vector<double>> data = {
        {51.5385, 96.0256},
        {46.1538, 94.4872},
        {42.8205, 91.4103},
        {40.7692, 88.3333},
        {38.7179, 84.8718},
        {35.641 , 79.8718},
        {33.0769, 77.5641},
        {28.9744, 74.4872},
        {26.1538, 71.4103},
        {23.0769, 66.4103},
        {22.3077, 61.7949},
        {22.3077, 57.1795},
        {23.3333, 52.9487},
        {25.8974, 51.0256},
        {29.4872, 51.0256},
        {32.8205, 51.0256},
        {35.3846, 51.4103},
        {40.2564, 51.4103},
        {44.1026, 52.9487},
        {46.6667, 54.1026},
        {50.    , 55.2564},
        {53.0769, 55.641 },
        {56.6667, 56.0256},
        {59.2308, 57.9487},
        {61.2821, 62.1795},
        {61.5385, 66.4103},
        {61.7949, 69.1026},
        {57.4359, 55.2564},
        {54.8718, 49.8718},
        {52.5641, 46.0256},
        {48.2051, 38.3333},
        {49.4872, 42.1795},
        {51.0256, 44.1026},
        {45.3846, 36.4103},
        {42.8205, 32.5641},
        {38.7179, 31.4103},
        {35.1282, 30.2564},
        {32.5641, 32.1795},
        {30.    , 36.7949},
        {33.5897, 41.4103},
        {36.6667, 45.641 },
        {38.2051, 49.1026},
        {29.7436, 36.0256},
        {29.7436, 32.1795},
        {30.    , 29.1026},
        {32.0513, 26.7949},
        {35.8974, 25.2564},
        {41.0256, 25.2564},
        {44.1026, 25.641 },
        {47.1795, 28.718 },
        {49.4872, 31.4103},
        {51.5385, 34.8718},
        {53.5897, 37.5641},
        {55.1282, 40.641 },
        {56.6667, 42.1795},
        {59.2308, 44.4872},
        {62.3077, 46.0256},
        {64.8718, 46.7949},
        {67.9487, 47.9487},
        {70.5128, 53.718 },
        {71.5385, 60.641 },
        {71.5385, 64.4872},
        {69.4872, 69.4872},
        {46.9231, 79.8718},
        {48.2051, 84.1026},
        {50.    , 85.2564},
        {53.0769, 85.2564},
        {55.3846, 86.0256},
        {56.6667, 86.0256},
        {56.1538, 82.9487},
        {53.8462, 80.641 },
        {51.2821, 78.718 },
        {50.    , 78.718 },
        {47.9487, 77.5641},
        {29.7436, 59.8718},
        {29.7436, 62.1795},
        {31.2821, 62.5641},
        {57.9487, 99.4872},
        {61.7949, 99.1026},
        {64.8718, 97.5641},
        {68.4615, 94.1026},
        {70.7692, 91.0256},
        {72.0513, 86.4103},
        {73.8462, 83.3333},
        {75.1282, 79.1026},
        {76.6667, 75.2564},
        {77.6923, 71.4103},
        {79.7436, 66.7949},
        {81.7949, 60.2564},
        {83.3333, 55.2564},
        {85.1282, 51.4103},
        {86.4103, 47.5641},
        {87.9487, 46.0256},
        {89.4872, 42.5641},
        {93.3333, 39.8718},
        {95.3846, 36.7949},
        {98.2051, 33.718 },
        {56.6667, 40.641 },
        {59.2308, 38.3333},
        {60.7692, 33.718 },
        {63.0769, 29.1026},
        {64.1026, 25.2564},
        {64.359 , 24.1026},
        {74.359 , 22.9487},
        {71.2821, 22.9487},
        {67.9487, 22.1795},
        {65.8974, 20.2564},
        {63.0769, 19.1026},
        {61.2821, 19.1026},
        {58.7179, 18.3333},
        {55.1282, 18.3333},
        {52.3077, 18.3333},
        {49.7436, 17.5641},
        {47.4359, 16.0256},
        {44.8718, 13.718 },
        {48.7179, 14.8718},
        {51.2821, 14.8718},
        {54.1026, 14.8718},
        {56.1538, 14.1026},
        {52.0513, 12.5641},
        {48.7179, 11.0256},
        {47.1795,  9.8718},
        {46.1538,  6.0256},
        {50.5128,  9.4872},
        {53.8462, 10.2564},
        {57.4359, 10.2564},
        {60.    , 10.641 },
        {64.1026, 10.641 },
        {66.9231, 10.641 },
        {71.2821, 10.641 },
        {74.359 , 10.641 },
        {78.2051, 10.641 },
        {67.9487,  8.718 },
        {68.4615,  5.2564},
        {68.2051,  2.9487},
        {37.6923, 25.7692},
        {39.4872, 25.3846},
        {91.2821, 41.5385},
        {50.    , 95.7692},
        {47.9487, 95.    },
        {44.1026, 92.6923}
    };
    return data;
}

/***
 * Generate a point cloud based on given statistics
*/
std::vector<std::vector<double>> generate_point_cloud(int num_samples){
    std::vector<std::vector<double>> data;

    std::default_random_engine generator;
    std::uniform_real_distribution<> rand_value(0.0, 100.0);
    for (int i = 0; i<num_samples; ++i){
        std::vector<double> new_point = {
            rand_value(generator),
            rand_value(generator),
        };
        data.push_back(new_point);
    }
    return data; //Zufällige gleichverteilte Daten x aus [0,100] und y aus [0,100] 
}

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
        if( (stats1->at(i) - stats2->at(i)) > error ){
            return false;
        }
    }
    return true;
}
