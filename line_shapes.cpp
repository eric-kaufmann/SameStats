#include <vector>
#include <string>
#include <math.h>
#include <iostream>


/***
 * calculates distributed points on defined lines
*/
std::vector<std::vector<double>> get_points_by_lines(std::vector<std::vector<std::vector<double>>> lines, int num_samples){
    
    std::vector<std::vector<double>> new_data;

    int num_lines = lines.size();
    int samples_per_line = num_samples / num_lines;

    for(auto l : lines){
        auto p1 = l[0];
        auto p2 = l[1];

        // lam * (x1,y1) + (1-lam) * (x2,y2)
        for(double lam = 0; lam <= 1.0; lam = lam + 1.0/(samples_per_line-1)){
            std::vector<double> new_point = {
                lam * p1[0] + (1.0-lam) * p2[0],
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
            {{0.0,30.0},{100.0,10.0}},
            {{0.0,50.0},{100.0,10.0}},
            {{0.0,70.0},{100.0,10.0}},
            {{0.0,90.0},{100.0,10.0}},
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

int main(){
    std::vector<std::string> shapes = {"x", "h_lines", "v_lines", "wide_lines", "high_lines", "slant_up", "slant_down", "center", "star", "down_parab"};
    for(auto s : shapes){
        std::vector<std::vector<double>> temp = get_points_for_shape("h_lines", 100000);
        std::cout << "worked for " << s << " num datapoints: " << temp.size() <<std::endl;
    }
}