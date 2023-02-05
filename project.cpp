#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <numeric>
#include <random>

double calc_p2p(std::vector<double> p1, std::vector<double> p2) {
  //calculte the euclidean distance between two points
  double sum = 0.0;
  for (int i = 0; i < p1.size(); ++i) {
    sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
  }
  return sqrt(sum);
}

std::map<std::string, double> calc_stats(const std::vector<std::vector<double>> &data) {
  // Calculate the mean and standard deviation of the x and y coordinates

  double xmean = std::accumulate(data.begin(), data.end(), 0.0, [](double a, const auto &point) { return a + point[0]; }) / data.size();
  double ymean = std::accumulate(data.begin(), data.end(), 0.0, [](double a, const auto &point) { return a + point[1]; }) / data.size();
  double xstd = std::sqrt(std::accumulate(data.begin(), data.end(), 0.0, [xmean](double a, const auto &point) { return a + std::pow(point[0] - xmean, 2); }) / (data.size() - 1));
  double ystd = std::sqrt(std::accumulate(data.begin(), data.end(), 0.0, [ymean](double a, const auto &point) { return a + std::pow(point[1] - ymean, 2); }) / (data.size() - 1));

  // Calculate the Pearson correlation coefficient of the data set
  double pearson = std::accumulate(data.begin(), data.end(), 0.0, [xmean, &ymean, xstd](double a, const auto &point) { return a + (point[0] - xmean) * (point[1] - ymean); }) / (data.size() - 1);
  pearson /= (xstd * ystd);

  return {
    {"xmean", xmean},
    {"ymean", ymean},
    {"xstd", xstd},
    {"ystd", ystd},
    {"pearson", pearson},
  };
}

bool check_stats(std::map<std::string, double> init_stats, std::map<std::string, double> target_stats, int d){
  bool result = true;
  double val1;
  for(const auto& [key1, value1] : init_stats){
  // if the key is not in target_stats, skip it
    if (target_stats.count(key1) == 0) continue;
      // get the corresponding value from target_stats
      double value2 = target_stats.at(key1);
      val1 = value1;
      val1 = std::trunc(val1*pow(10,d))/pow(10,d); //take statistic values with two decimal places 
      value2 = std::trunc(value2*pow(10,d))/pow(10,d);
      // std::cout << key1 << ": val1 " << val1 << " value2 " << value2 << std::endl;
      // std::cout << key1 << ": " << value1 << " value2 " << target_stats.at(key1) << std::endl;
      if(val1!=value2){ //compare statistic values 
        result = false;
        break;
      }
        
  }  
  return result;
}

void print_stats(std::map<std::string, double> stats){
  for(auto it = stats.cbegin(); it != stats.cend(); ++it)
  {
    std::cout << it->first << ": " << it->second << std::endl;
  }  
}

void print_data(std::vector<std::vector<double>> data){
  for (const auto &point : data) {
    std::cout << "x: " << point[0] << ", y: " << point[1] << std::endl;
  }  
}

int main() {

  // Initialize the data using a 2D vector (x,y)
  std::vector<std::vector<double>> init_data = {
    {0.0, 2.0},
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


  // Print the x and y coordinates of each point
  print_data(init_data);

  std::cout << std::endl << std::endl;
  
  // Print the x and y coordinates of each point
  print_data(target_data);

  std::cout << std::endl << std::endl;

  //calculate statistics of init_data and print them
  std::map<std::string, double> init_stats = calc_stats(init_data);  
  print_stats(init_stats);

  std::cout << std::endl << std::endl;
  
  //calculate statistics of target_data and print them
  std::map<std::string, double> target_stats = calc_stats(target_data);  
  print_stats(target_stats);
  
  int d = 2;
  //check if the statistics are approx the same
  if(check_stats(init_stats, target_stats, d)){
    std::cout << std::endl << "gleiche stats!" << std::endl;
  }else{
    std::cout << std::endl << "nicht gleiche stats!" << std::endl;
  }




  int steps = 1e6;

  double init_data_distance = INFINITY;
  double value;
  // calculate the smallest aggregated distance between all points of data to target
  for (const auto &point1 : init_data) {
    for (const auto &point2 : target_data) {
      value = calc_p2p(point1, point2); 
      if(value < init_data_distance)
        init_data_distance = value;
    }  
  }
  std::cout << std::endl << "init_data_distance: " << init_data_distance << std::endl;

  int hit_counter = 0; // gets number of changes of the data points
  int max_rd = init_data.size();
  std::cout << "num data points: " << max_rd << std::endl;

  init_stats = calc_stats(init_data);
  
  for(int i = 0; i < steps; i++){
    std::cout << i << std::endl;
    std::random_device rd;
    std::uniform_int_distribution<int> dist(0, max_rd-1);
    int r_num = dist(rd); // get random int between 0 and num_data_points-1
    std::vector<double> x = init_data.at(r_num);
    // calculate the smallest distance between the point to all target data points
    double curr_min_dist = INFINITY, value;
    for (const auto &point2 : target_data) {
      value = calc_p2p(x, point2); 
      if(value < curr_min_dist)
        curr_min_dist = value;
    }
    //add a random number between -0.1, 0.1 to the respective point and calculate smallest distance to target data. Loop until smallest distance gets smaller than initiate smallest distance. 
    double new_min_dist = curr_min_dist+1.0;
    std::uniform_real_distribution<float> dist2(-0.1,0.1);
    std::vector<double> new_x = {0.0, 0.0};
    while (new_min_dist>curr_min_dist && curr_min_dist != 0){
      std::random_device rd2;
      for (int i = 0; i < new_x.size(); i++) {
        new_x[i] = x[i] + dist2(rd2);
      }
      double buff = INFINITY;
      for (const auto &point2 : target_data) {
        value = calc_p2p(new_x, point2); 
        if(value < buff){
          buff = value;
          new_min_dist = buff;
        }
      } 
    }

    // replace respective data point in test data with new_x
    test_data = init_data;
    test_data.at(r_num) = new_x;

    // check if the statistics are approx correct
    if(check_stats(init_stats, calc_stats(test_data), d)){
      // new data point gets integrated in init_data
      init_data.at(r_num) = new_x;
      hit_counter++;
      // std::cout << i << std::endl;
    }
  }

  std::cout << "hit counter: " << hit_counter << std::endl;
  std::cout << "target data:" << std::endl;
  print_data(target_data);
  std::cout << std::endl << "found data set: "<< std::endl;
  print_data(init_data);

  std::cout << std::endl << "initial statistics:" << std::endl;
  print_stats(init_stats);
  std::cout << std::endl << "final statistics: " << std::endl;
  print_stats(calc_stats(test_data));

  return 0;
}