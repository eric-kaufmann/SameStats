#include "utils.h"
#include "line_shapes.h"


int main(int argc, char const *argv[]){

    bool save_data = true; // generates csv files and plots data

    int num_steps = 1e6; // number of iteration steps
    double error = 1e-2; // maximum statistical diffence
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
    // {"datasaurus", "random"};
    std::string init_shape = "random";
    init_shape = "datasaurus";
    std::vector<std::vector<double>> working_data;
    if (init_shape == "datasaurus"){
        working_data = get_datasaurus_data();
    } 
    else if (init_shape == "random"){
        working_data = generate_point_cloud(num_samples);
    }

    // All possible target shapes -> target_shape[i]
    // 0:"x"   1: "h_lines"    2: "v_lines"    3: "wide_lines" 4: "high_lines" 5: "slant_up" 6: "slant_down" 7: "center" 8: "star" 9: "down_parab"     
    std::vector<std::string> target_shape = {"x", "h_lines", "v_lines", "wide_lines", "high_lines", "slant_up", "slant_down", "center", "star", "down_parab"};
    std::vector<std::vector<double>> target_data = get_points_for_shape(target_shape[8], working_data.size());

    // Transpose data
    if(transpose){
        working_data = transpose_data(working_data);
        target_data = transpose_data(target_data);
    }

    // Save Input data
    if(save_data){
        std::string fn = "input";
        data_to_csv(working_data, fn+"_generated_data.csv");
        data_to_csv(target_data, fn+"_target_data.csv");
        generate_scatter_plot(fn);
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
    std::normal_distribution<double> rand_shift(0.0,1.0); // generates shift for points (standard normal distribution -> adopted from autodesk paper)
    std::uniform_real_distribution<> rand_break(0.0, 1.0); // gives a chance to break out of loop (continous uniform dist in [0.0, 2.0) adopted from autodesk paper)

    // double temperature = init_temperature;
    double temperature;
    auto start = std::chrono::high_resolution_clock::now();
    bool check = true;
    int counter = 0;
    // Main loop over num_steps steps
    for(int step=0; step<num_steps; ++step){
        if (step % 1000==0){
            std::cout << "step: " << step << std::endl;
        }
        temperature = (max_temp - min_temp) * easeInOutQuad(((num_steps-step)/num_steps)) + min_temp;
        // get random point
        int rand_point_idx = rand_point_dist(generator);
        double rpx = working_data[0][rand_point_idx]; // x value of random point
        double rpy = working_data[1][rand_point_idx]; // y value of random point
        //funktioniert der zugriff so auch wenn transpose = false ist? -> Nein
        std::vector<double> new_point = {rpx, rpy};

        double init_min_dist = minDist(new_point, target_data, target_data[0].size()); 

        // Get new point with smaller distance
        double new_min_dist = INFINITY;
        double rpx_shift; // x value of random point with shift
        double rpy_shift; // y value of random point with shift
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
            //std::cout << "hits: " << counter << std::endl;
        }           

        // // Adjust temperature for simulated annealing
        // if(init_temperature/(step+1) > min_temperature){ 
        //     temperature = init_temperature/(step+1);
        // }
    }
    //time measurement of data transformation process
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "elapsed time: " << duration.count()*1e-9 << "s" << std::endl;
    if(save_data){
        std::string fn = "output";
        data_to_csv(working_data, fn+"_generated_data.csv");
        data_to_csv(target_data, fn+"_target_data.csv");
        generate_scatter_plot(fn);
    }
    return 0;
}
