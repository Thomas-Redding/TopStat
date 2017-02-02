#include <fstream>
#include <iostream>
#include <unistd.h>
#include <vector>
#include <string>

#include "Graph.hpp"
#include "Matrix.hpp"
#include "Algorithms.hpp"

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
    // let n = number of points; d = number of dimensions; e = edges in graph
	std::string file_path = argv[1];
    double epsilon = atof(argv[2]);

    // read file
    // O(n * d)
    std::string file_contents = read_file(file_path);

    // convert points as text into a matrix
    // O(n * d)
    Matrix<double>* points = string_to_matrix(file_contents);

    // convert points as a matrix into a graph
    // O(n^2 * d)
    // can probably speed up to O(e * d)
    IntGraph *my_graph = matrix_to_graph(points, epsilon);

    // print graph
    std::cout << *my_graph;
    return 0;
}
