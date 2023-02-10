#include "sim_ann.h"


int main(int argc, char const *argv[]){

    bool save_data = true; // generates csv files and plots data

    int num_steps = 1e2; // number of iteration steps
    double error = 1e-2; // maximum statistical diffence
    double max_shift = 5; // std of rand-norm for data point shift 
    //warum 0.01?
    int num_samples = 100000; // #datapoints if a random point cloud is generated
    
    // set data format by transposing data
    // true -> {{x1, x2, ...}, {y1, y2, ...}} :: false -> {{x1,y1}, {x2,y2}, ...} 
    bool transpose = true;
    
    double init_temperature = 1.0; // initial temperature (gets smaller every iteration)
    double min_temperature = 0.0; // minimal temperature (always temperature > min_temperature)
   
    // All possible init shapes 
    // {"datasaurus", "random"};
    std::string init_shape = "datasaurus";
    // std::string init_shape = "random";
    std::vector<std::vector<double>> working_data;
    if (init_shape == "datasaurus"){
        working_data = get_datasaurus_data();
    } 
    else if (init_shape == "random"){
        working_data = generate_point_cloud(num_samples);
    }

    // All possible target shapes 
    // {"x", "h_lines", "v_lines", "wide_lines", "high_lines", "slant_up", "slant_down", "center", "star", "down_parab"};
    std::vector<std::string> target_shape = {"x", "h_lines", "v_lines", "wide_lines", "high_lines", "slant_up", "slant_down", "center", "star", "down_parab"};
    std::vector<std::vector<double>> target_data = get_points_for_shape(target_shape[0], working_data.size());
    // std::cout << "working data: " << std::endl;
    // print_matrix(working_data);
    // std::cout << std::endl << "target data: " << std::endl;
    // print_matrix(target_data);
    // Transpose data
    if(transpose){
        working_data = transpose_data(working_data);
        target_data = transpose_data(target_data);
    }

    // Calculate initial stats
    std::vector<double> init_stats(5); 
    std::vector<double> working_stats(5);
    std::vector<double> target_stats(5);

    calc_stats(&working_data, &init_stats); // set initial stats
    calc_stats(&target_data, &target_stats); // set target stats


    // create radom number generator
    std::default_random_engine generator;

    // init distributions
    std::uniform_int_distribution<int> rand_point_dist(0,working_data[0].size()); // selects a random point idx
    std::normal_distribution<double> rand_shift(0.0,max_shift); // generates shift for points
    std::uniform_real_distribution<> rand_break(0.0, 1.0); // gives a chance to break out of loop

    double temperature = init_temperature;
    auto start = std::chrono::high_resolution_clock::now();
    bool check = true;
    int counter = 0;
    // Main loop over num_steps steps
    for(int step=0; step<num_steps; ++step){
        std::cout << "step: " << step << std::endl;
        // get random point
        int rand_point_idx = rand_point_dist(generator);
        std::cout << "random point " << working_data[0][rand_point_idx] << "  " << working_data[1][rand_point_idx] << std::endl;
        double rpx = working_data[0][rand_point_idx]; // x value of random point
        double rpy = working_data[1][rand_point_idx]; // y value of random point
        //funktioniert der zugriff so auch wenn transpose = false ist?
        std::vector<double> new_point = {rpx, rpy};
        // TODO: calculate initial min dist from (rpx, rpy) to target_data
        double init_min_dist = minDist(new_point, target_data, target_data[0].size()); 
        std::cout << "dist: " << init_min_dist << std::endl; //debugg
        // std::cout << "size: " << target_data[0].size() << std::endl; //debugg


        // Get new point with smaller distance
        double new_min_dist = INFINITY;
        double rpx_shift; // x value of random point with shift
        double rpy_shift; // y value of random point with shift
        while(new_min_dist >= init_min_dist){
            // shift point
            rpx_shift = rpx + rand_shift(generator);
            rpy_shift = rpy + rand_shift(generator);
            // std::cout << "random shift x: " << rpx_shift << "  y: " << rpy_shift << std::endl;
            // std::cout <<  std::endl;
            std::vector<double> mod_point = {rpx_shift, rpy_shift};
            // TODO: calculate new min dist from (rpx_shift, rpy_shift) to target_data
            new_min_dist = minDist(mod_point, target_data, target_data[0].size());
            
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
        check = true;
        // check if statistics are ca. the same to our inital stats
        if(!check_stats(&working_stats, &init_stats, error)){
            // not the same -> return points to old value
            check = false;
            working_data[0][rand_point_idx] = rpx;
            working_data[1][rand_point_idx] = rpy;
        }
        if(check){
            counter++;
            std::cout << "hits: " << counter << std::endl;
        }
            

        // Adjust temperature for simulated annealing
        if(init_temperature/(step+1) > min_temperature){ 
            temperature = init_temperature/(step+1);
            std::cout << "temp: " << temperature << std::endl;
        }
    }
    //time measurement of data transformation process
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "elapsed time: " << duration.count()*1e-9 << "s" << std::endl;
    std::cout << "working data: " << std::endl;
    print_matrix(working_data);
    std::cout << std::endl << "target data: " << std::endl;
    print_matrix(target_data);
    if(save_data){
        data_to_csv(&working_data, "generated_data.csv");
        data_to_csv(&working_data, "target_data.csv");
        generate_scatter_plot();
    }
    return 0;
}
