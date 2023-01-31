#include <fstream>
#include <string>

void matrixToCSV(double m[SIZE][DIMENSIONS], std::string filename)
{
    std::ofstream outfile;
    outfile.open(filename);

    for (int i = 0; i < SIZE; i++)
    {
        for (int j = 0; j < DIMENSIONS; j++)
        {
            outfile << m[i][j] << ",";
        }
        outfile << "\n";
    }

    outfile.close();
}

void generateScatterPlots(){
    system("python3 scatter_plot.py");
    return 0;
}