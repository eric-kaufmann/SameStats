#include "utils.h"
#include "line_shapes.h"


int main(int argc, char const *argv[]){

    bool save_data = true; // generates csv files and plots data

    int num_steps = 1e6; // number of iteration steps
    double error = 1; // maximum statistical diffence (1e-2)
    double max_shift = 0.1; // std of rand-norm for data point shift (0.1 adopted from autodesk paper) 
    int num_samples = 100; // #datapoints if a random point cloud is generated
    double max_temp = 0.4; //param for sim annealing (0.4 adopted from paper)
    double min_temp = 0.01; // param for sim annealing (0.0 or 0.01 in paper)
    // set data format by transposing data
    // true -> {{x1, x2, ...}, {y1, y2, ...}} :: false -> {{x1,y1}, {x2,y2}, ...} 
    bool transpose = true;
    
    // double init_temperature = 1.0; // initial temperature (gets smaller every iteration)
    // double min_temperature = 0.0; // minimal temperature (always temperature > min_temperature)
   
    // All possible init shapes 
    std::vector<std::string> shape_name = {"datasaurus", "random"};
    std::string init_shape = shape_name[0];
    std::vector<std::vector<double>> working_data;
    if (init_shape == "datasaurus"){
        working_data = get_datasaurus_data();
    } 
    else if (init_shape == "random"){
        working_data = generate_point_cloud(num_samples);
    }

    std::cout << "size working data: " << working_data.size() << std::endl;

    // All possible target shapes -> target_shape[i]
    // 0:"x"   1: "h_lines"    2: "v_lines"    3: "wide_lines" 4: "high_lines" 5: "slant_up" 6: "slant_down" 7: "center" 8: "star" 9: "down_parab"     
    std::vector<std::string> target_shape = {"x", "h_lines", "v_lines", "wide_lines", "high_lines", "slant_up", "slant_down", "center", "star", "down_parab"};
    std::vector<std::vector<double>> target_data = get_points_for_shape(target_shape[0], working_data.size());

    std::cout << "size target data: " << target_data.size() << std::endl;
    int num_points;

    // Transpose data
    if(transpose){
        working_data = transpose_data(working_data);
        target_data = transpose_data(target_data);
        num_points = working_data[0].size(); //number of points in data set
    }

    // Save Input data
    if(save_data){
        std::string fn = "input";
        data_to_csv(working_data, fn+"_generated_data.csv");
        data_to_csv(target_data, fn+"_target_data.csv");
        generate_scatter_plot(fn);
    }
    int padding = 4;
    int threads=2; 
#pragma omp parallel num_threads(threads)
{
        threads = omp_get_num_threads();
}
    std::cout << "num threads: " << threads << std::endl;

    std::vector<std::vector<std::vector<double>>> new_working_data;
    new_working_data = partitionData(threads, padding, working_data);
    printVector(new_working_data);
    auto start = std::chrono::high_resolution_clock::now();
    int check_portion = 0; //debug
    int rest = working_data[0].size()%threads;
    std::vector<double> global_stats(5); //x_mean, y_mean, x_std, y_std, pearson
    calc_stats(&working_data, &global_stats);

    // int j = 0;
// #pragma omp parallel for num_threads(threads)
    for(auto &data_portion : new_working_data){
        // bool rest_check;
        int toggle = 0;
        print_matrix(data_portion);
        // Calculate initial stats
        std::vector<double> init_stats(5); 
        std::vector<double> working_stats(5);
        // std::vector<double> target_stats(5);
        int vec_size = num_points/threads;
        int gap = 0;
        if(num_points%threads != 0){
            vec_size++;
            gap = vec_size - num_points%vec_size; 
        }
        if(count_nonzeros(data_portion, padding) < vec_size)
            toggle = 1;
        calc_partitionedStats(data_portion, &init_stats, num_points, global_stats[0], global_stats[1], global_stats[2], global_stats[3]); // set initial stats
        // calc_stats(&target_data, &target_stats); // set target stats

        // create radom number generator
        std::default_random_engine generator;

        // init distributions
        std::uniform_int_distribution<int> rand_point_dist(padding, vec_size+padding-1-toggle*gap); // selects a random point idx
        std::normal_distribution<double> rand_shift(0.0,1.0); // generates shift for points (standard normal distribution -> adopted from autodesk paper)
        std::uniform_real_distribution<> rand_break(0.0, 1.0); // gives a chance to break out of loop (continous uniform dist in [0.0, 2.0) adopted from autodesk paper)
        // double temperature = init_temperature;
        double temperature;
        
        bool check = true;
        bool loop =  true;
        double tol = 1.0;
        int counter = 0;
        double old_max = INFINITY;
        // Main loop over num_steps steps
        std::vector<double> max_point = {0.0,0.0};      

        for(int step=0; loop; step++){
            // omp_get_thread_num() == 1 && 
            if(step%1000 == 0)
                std::cout << "iter: " << check_portion << " step: " << step <<"  maxnorm: " << old_max << " point (x,y): " << max_point[0] << " , " << max_point[1] << std::endl;
            temperature = (max_temp - min_temp) * easeInOutQuad(((num_steps-step)/num_steps)) + min_temp;
            // get random point
            int rand_point_idx = rand_point_dist(generator);
            double rpx = data_portion[0][rand_point_idx]; // x value of random point
            double rpy = data_portion[1][rand_point_idx-padding]; // y value of random point
            //funktioniert der zugriff so auch wenn transpose = false ist? -> Nein
            std::vector<double> new_point = {rpx, rpy};

            double init_min_dist = minDist(new_point, target_data, target_data[0].size()); 

            // Get new point with smaller distance
            double new_min_dist = INFINITY;
            double rpx_shift = 0.0; // x value of random point with shift
            double rpy_shift = 0.0; // y value of random point with shift
            while(new_min_dist >= init_min_dist){
                // shift point
                rpx_shift = rpx + max_shift*rand_shift(generator); //scaled standard normal distribution --> std = max_shift, mu = 0
                rpy_shift = rpy + max_shift*rand_shift(generator);
                std::vector<double> mod_point = {rpx_shift, rpy_shift};

                new_min_dist = minDist(mod_point, target_data, target_data[0].size());
                
                // simmulated annealing: break out of loop with a specific chance
                if(rand_break(generator)<temperature){
                    break;
                }
            }

            // embed new point into vector
            data_portion[0][rand_point_idx] = rpx_shift;
            data_portion[1][rand_point_idx-padding] = rpy_shift;

            // calc new stats
            calc_partitionedStats(data_portion, &working_stats, num_points, global_stats[0], global_stats[1], global_stats[2], global_stats[3]); // set initial stats
            check = true;
            // check if statistics are ca. the same to our inital stats
            if(!check_stats(&working_stats, &init_stats, error)){ //kann sein dass error verändert werden muss --> vergleichen größere werte
                // not the same -> return points to old value
                check = false;
                data_portion[0][rand_point_idx] = rpx;
                data_portion[1][rand_point_idx-padding] = rpy;
            }
            if(check){
                counter++;
                //std::cout << "hits: " << counter << std::endl;
            }     
            if (step % 100 ==0){
                double max = 0.0;
                for(int i = padding; i < vec_size+padding-toggle*gap; i++){ //maximalen abstand eines punktes zum shape finden
                    std::vector<double> point = {data_portion[0][i], data_portion[1][i-padding]};
                    if(minDist(point, target_data, target_data[0].size()) > max){
                        max = minDist(point, target_data, target_data[0].size());
                        max_point = point;
                    }
                }
                if(max < old_max)
                    old_max = max;
                if(old_max < tol)
                    loop = false;
            }
            if(step > 200000)
                loop = false;
            // // Adjust temperature for simulated annealing
            // if(init_temperature/(step+1) > min_temperature){ 
            //     temperature = init_temperature/(step+1);
            // }
        }
        check_portion++;
    }

    working_data = refactorData(new_working_data, working_data, threads, padding);
    //time measurement of data transformation process
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "elapsed time: " << duration.count()*1e-9 << "s" << std::endl;

    std::vector<double> end_stats(5);
    calc_stats(&working_data, &end_stats);

    std::cout << "initial stats: " << std::endl;
    for(auto &stats : global_stats){
        std::cout << stats << std::endl;
    }
    std::cout << "end stats: " << std::endl;
    for(auto &stats : end_stats){
        std::cout << stats << std::endl;
    }

    if(save_data){
        std::string fn = "output";
        std::cout << "save working data.." << std::endl;
        data_to_csv(working_data, fn+"_generated_data.csv");
        std::cout << "save target data.." << std::endl;
        data_to_csv(target_data, fn+"_target_data.csv");
        generate_scatter_plot(fn);
    }
    return 0;
}
