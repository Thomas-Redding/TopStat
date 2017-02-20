#include <fstream>
#include <iostream>
#include <unistd.h>
#include <vector>
#include <string>

#include "IntGraph.hpp"
#include "Array.hpp"
#include "Matrix.hpp"
#include "Algorithms.hpp"

/*
 * utility function to read the contents of a text file
 */
std::string *read_file(std::string filepath) {
    std::string line;
    std::string *contents = new std::string("");

    std::ifstream file(filepath);
    if (file.is_open()) {
        while(getline(file, line)) {
            *contents += line;
            *contents += "\n";
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
	if (argc != 4) {
		std::cout << "Error: need to pass in three arguments (input-file, epsilon, and betti #s)" << std::endl;
		return 0;
	}

    // process arguments
    // let n = number of points
    // let d = number of dimensions
    // let e = edges in graph of average point
    // let b = number of betti #s
	std::string file_path = argv[1];
    double epsilon = atof(argv[2]);
    double max_betti = atoi(argv[3]);

    if (epsilon <= 0) {
        std::cout << "Error: epsilon must be positive" << std::endl;
        return 0;
    }

    if (max_betti <= 1) {
        std::cout << "Error: max_betti must be greater than 1" << std::endl;
        return 0;
    }

    // read file
    // O(n * d)
    std::string *file_contents = read_file(file_path);

    // convert points as text into a matrix
    // O(n * d)
    Matrix<double>* points = string_to_matrix(file_contents);

    // convert points as a matrix into a graph
    // O(n^2 * d)
    // can probably speed up to O(n * e * d)
    IntGraph *my_graph = matrix_to_graph(points, epsilon);
    std::cout << *my_graph;

    std::cout << "START";

    // TODO: find boundry matrices
    // O(n * d * e^b)
    // creates O(n * e^b) cells
    Array<Array<std::set<uint>>*> k_cells = Array<Array<std::set<uint>>*>(max_betti+1);

    std::cout << "MID1";

    // 0-cells
    k_cells[0] = new Array<std::set<uint>>(my_graph->size());
    for (int i = 0; i < k_cells[0]->size(); ++i)
        k_cells[0]->at(i) = std::set<uint>();

    std::cout << "MID2";
    
    k_cells[1] = new Array<std::set<uint>>(my_graph->edge_count());

    std::cout << "END";

    for (int i = 2; i <= max_betti; ++i) {
        // get_k_cells(const IntGraph *graph, uint dim, Array<std::set<uint>>* lower_cells);
    }

    // TODO: row-reduce matrices [slowest step]
    // O(n^3 * e^b)
    return 0;
}
