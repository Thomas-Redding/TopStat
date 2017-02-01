#include <fstream>
#include <iostream>
#include <unistd.h>
#include <vector>
#include <string>

#include "Graph.hpp"
#include "Matrix.hpp"

/*
 * utility function to read the contents of a text file
 */
std::string read_file(std::string filepath) {
    std::string line;
    std::string contents = "";

    std::ifstream file(filepath);
    if (file.is_open()) {
        while(getline(file, line)) {
            contents += line;
            contents += "\n";
        }
        file.close();
    }
    else {
        throw std::invalid_argument("File \"" + filepath + "\" not found.");
    }
    return contents;
}

int main(int argc, const char * argv[]) {
    // make sure we have the correct number of arguments
	if (argc != 3) {
		std::cout << "Error: need to pass in two arguments (input-file and epsilon)" << std::endl;
		return 0;
	}

    // process arguments
	std::string file_path = argv[1];
    double epsilon = atof(argv[2]);
    double epislon_squared = epsilon * epsilon;
    std::string file_contents = read_file(file_path);
    Matrix<double> points = string_to_matrix(file_contents);

    // construct graph
    int num_points = points.get_width();
    int num_dim = points.get_height();
    double dx;
    double dist_squred;
    IntGraph my_graph(num_points);
    for (uint i = 0; i < num_points; ++i) {
        for (uint j = 0; j < i; ++j) {
            dist_squred = 0;
            for (uint k = 0; k < num_dim; ++k) {
                dx = points.get(i, k) - points.get(j, k);
                dist_squred += dx*dx;
            }
            if (dist_squred < epislon_squared) {
                my_graph.add_edge(i, j);
            }
        }
    }

    // print graph
    std::cout << my_graph;
    return 0;
}

