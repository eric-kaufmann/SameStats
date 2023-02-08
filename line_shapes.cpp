// def get_points_for_shape(line_shape):
//     lines = []
//     if line_shape == 'x':
//         l1 = [[20, 0], [100, 100]]
//         l2 = [[20, 100], [100, 0]]
//         lines = [l1, l2]
//     elif line_shape == "h_lines":
//         lines = [[[0, y], [100, y]] for y in [10, 30, 50, 70, 90]]
//     elif line_shape == 'v_lines':
//         lines = [[[x, 0], [x, 100]] for x in [10, 30, 50, 70, 90]]
//     elif line_shape == 'wide_lines':
//         l1 = [[10, 0], [10, 100]]
//         l2 = [[90, 0], [90, 100]]
//         lines = [l1, l2]
//     elif line_shape == 'high_lines':
//         l1 = [[0, 10], [100, 10]]
//         l2 = [[0, 90], [100, 90]]
//         lines = [l1, l2]
//     elif line_shape == 'slant_up':
//         l1 = [[0, 0], [100, 100]]
//         l2 = [[0, 30], [70, 100]]
//         l3 = [[30, 0], [100, 70]]
//         l4 = [[50, 0], [100, 50]]
//         l5 = [[0, 50], [50, 100]]
//         lines = [l1, l2, l3, l4, l5]
//     elif line_shape == 'slant_down':
//         l1 = [[0, 100], [100, 0]]
//         l2 = [[0, 70], [70, 0]]
//         l3 = [[30, 100], [100, 30]]
//         l4 = [[0, 50], [50, 0]]
//         l5 = [[50, 100], [100, 50]]
//         lines = [l1, l2, l3, l4, l5]
//     elif line_shape == 'center':
//         cx = 54.26
//         cy = 47.83
//         l1 = [[cx, cy], [cx, cy]]
//         lines = [l1]
//     elif line_shape == 'star':
//         star_pts = [10,40,40,40,50,10,60,40,90,40,65,60,75,90,50,70,25,90,35,60]
//         pts = [star_pts[i:i+2] for i in range(0, len(star_pts), 2)]
//         pts = [[p[0]*0.8 + 20, 100 - p[1]] for p in pts]
//         pts.append(pts[0])
//         lines = [pts[i:i+2] for i in range(0, len(pts)-1, 1)]
//     elif line_shape == 'down_parab':
//         curve = [[x, -((x-50)/4)**2 + 90] for x in np.arange(0, 100, 3)]
//         lines = [curve[i:i+2] for i in range(0, len(curve)-1, 1)]

//     return lines


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
        for(double lam = 0; lam < 1.0; lam = lam + 1.0/samples_per_line){
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
            {{0.0,0.0},{100.0,100.0}}, // one line starting at (0,0) and ending at (100,100)
            {{100.0,0.0},{0.0,100.0}}
        };
    }
    else if (shape == "h_lines"){
        //lines = [[[0, y], [100, y]] for y in [10, 30, 50, 70, 90]]
        lines = {
            {{0.0,10.0},{100.0,10.0}},
            {{0.0,30.0},{100.0,10.0}},
            {{0.0,50.0},{100.0,10.0}},
            {{0.0,70.0},{100.0,10.0}},
            {{0.0,90.0},{100.0,10.0}},
        };
    }

    // generate target data
    std::vector<std::vector<double>> target_data = get_points_by_lines(lines, num_samples);

    for(auto p : target_data){
        std::cout << p[0] << "  " << p[1] << std::endl;
    }
    return target_data;
}

int main(){
    std::vector<std::vector<double>> temp = get_points_for_shape("h_lines", 50);
}