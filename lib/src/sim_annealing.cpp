#include "utils.h"
#include "line_shapes.h"

double sim_annealing(std::string shape, int threads, int serial_steps){  // all possible shapes: {"x", "h_lines", "v_lines", "wide_lines", "high_lines", "slant_up", "slant_down", "center", "star", "down_parab"};


    int max_threads;
    // int threads=8; //seems to be the best trade off between speed and accuracy. --> more threads hinder shift to a specific shape due to lower degree of freedom for each point in a smaller group. 
#pragma omp parallel
{
    max_threads = omp_get_max_threads();
}
    if(threads > max_threads)
        threads = max_threads;
    std::cout << "max threads: " << max_threads << " involved threads: " << threads << std::endl;
    bool save_data = true; // generates csv files and plots data
    int padding = 4;

    // int serial_steps = 1e6;
    int num_steps = serial_steps/threads; // number of iteration steps --> to shift one point the same times in parallelized version, we can divide serial_steps by number of threads.
    double error = sqrt(1e-2); // maximum statistical diffence (1e-2)     we take sqrt() because statistics in local groups is squared
    double max_shift = 0.1; // std of rand-norm for data point shift (0.1 adopted from autodesk paper) 
    int num_samples = 1000; // #datapoints if a random point cloud is generated
    double max_temp = 0.4; //param for sim annealing (0.4 adopted from paper)
    double min_temp = 0.01; // param for sim annealing (0.0 or 0.01 in paper)

    // set data format by transposing data
    bool transpose = true;    // true -> {{x1, x2, ...}, {y1, y2, ...}} :: false -> {{x1,y1}, {x2,y2}, ...} 

    
    // All possible init shapes 
    std::vector<std::string> shape_name = {"datasaurus", "random"};
    std::string init_shape = shape_name[1];
    std::vector<std::vector<double>> working_data;
    if (init_shape == "datasaurus"){
        working_data = get_datasaurus_data();
    } 
    else if (init_shape == "random"){
        working_data = generate_point_cloud(num_samples);
    }

    std::cout << "size working data: " << working_data.size() << std::endl;

    std::vector<std::vector<double>> target_data = get_points_for_shape(shape, working_data.size());

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
        std::string fn = "input_"+shape;
        data_to_csv(working_data, fn+"_generated_data.csv");
        data_to_csv(target_data, fn+"_target_data.csv");
        generate_scatter_plot(fn);
    }
    std::vector<double> global_stats(5); //x_mean, y_mean, x_std, y_std, corr_xy
    calc_stats(&working_data, &global_stats);


    std::vector<std::vector<std::vector<double>>> new_working_data; 
    new_working_data = partitionData(threads, padding, working_data); //create thread-wise data with padding in between to prevent false sharing
    // printVector(new_working_data);
    auto start = std::chrono::high_resolution_clock::now(); 
    int check_portion = 0; //debug
    

    // int j = 0;
#pragma omp parallel for num_threads(threads) //do sim_annealing independently in each thread
    for(auto &data_portion : new_working_data){
        // bool rest_check;

        int toggle = 0;
        // print_matrix(data_portion);
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
        double temperature;
        
        bool check = true;
        double tol = 0.5;
        int counter = 0;
        double old_max = INFINITY;
        // Main loop over num_steps steps
        std::vector<double> max_point = {0.0,0.0};      

        for(int step=0; step < num_steps; step++){
            
            temperature = (max_temp - min_temp) * easeInOutQuad(((num_steps-step)/num_steps)) + min_temp;

            // get random point
            int rand_point_idx = rand_point_dist(generator);
            double rpx = data_portion[0][rand_point_idx]; // x value of random point
            double rpy = data_portion[1][rand_point_idx-padding]; // y value of random point

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

                // simmulated annealing: break out of loop with a specific chance
                if(rand_break(generator)<temperature){
                    break;
                }

                std::vector<double> mod_point = {rpx_shift, rpy_shift};
                new_min_dist = minDist(mod_point, target_data, target_data[0].size());
                
                
            }

            // embed new point into vector
            data_portion[0][rand_point_idx] = rpx_shift;
            data_portion[1][rand_point_idx-padding] = rpy_shift;

            // calc new stats
            calc_partitionedStats(data_portion, &working_stats, num_points, global_stats[0], global_stats[1], global_stats[2], global_stats[3]); // set initial stats
            check = true;
            // check if statistics are ca. the same to our inital stats
            if(!check_Partitionedstats(&working_stats, &init_stats, error)){ //kann sein dass error verändert werden muss --> vergleichen größere werte
                // not the same -> return points to old value
                check = false;
                data_portion[0][rand_point_idx] = rpx;
                data_portion[1][rand_point_idx-padding] = rpy;
            }
           
        }

    }
    
    working_data = refactorData(new_working_data, working_data, threads, padding); //merge data from each thread to one data set
    //time measurement of data transformation process
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

    std::vector<double> end_stats(5);
    calc_stats(&working_data, &end_stats);


    std::vector<std::string> stat_names = {"mu_x: ", "mu_y: ", "sigma_x: ", "sigma_y: ", "corr_xy: "};
    std::cout << std::endl << "initial stats: " << std::endl;
    for(int i = 0; i < (int)global_stats.size(); i++){
        std::cout << stat_names[i] << global_stats[i] << std::endl;
    }

    std::cout << std::endl << "end stats: " << std::endl;
    for(int i = 0; i < (int)end_stats.size(); i++){
        std::cout << stat_names[i] << end_stats[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "elapsed time sim annealing: " << duration.count()*1e-9 << "s" << std::endl;

    // print_matrix(result_data);
    // compareData(working_data, result_data);
    double maxError = get_MaxError(global_stats, end_stats);
    if(save_data){
        std::string fn = "output_"+shape;
        std::cout << "save working data.." << std::endl;
        data_to_csv(working_data, fn+"_generated_data.csv");
        std::cout << "save target data.." << std::endl;
        data_to_csv(target_data, fn+"_target_data.csv");
        generate_scatter_plot(fn);
    }


    return maxError;
}