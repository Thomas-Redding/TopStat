#include <fstream>
#include <iostream>
#include <unistd.h>
#include <vector>
#include <string>

#include "Graph.hpp"
#include "Matrix.hpp"

std::string read_file(std::string path) {
    std::string line;
    std::string contents = "";

    std::ifstream file(path);
    if (file.is_open()) {
        while(getline(file, line)) {
            contents += line;
            contents += "\n";
        }
        file.close();
    }
    else {
        throw std::invalid_argument("File \"" + path + "\" not found.");
    }
    return contents;
}

Matrix<double> string_to_matrix(std::string str) {
    uint from = 0;
    std::vector<std::vector<double>> answer;
    answer.push_back(std::vector<double>());
    for (int i = 0; i < str.length(); ++i) {
        if (str[i] == ' ') {
            std::string float_str = str.substr(from, i - from);
            answer.back().push_back(stof(float_str));
            from = i + 1;
        }
        else if (str[i] == '\n') {
            std::string float_str = str.substr(from, i - from);
            answer.back().push_back(stof(float_str));
            from = i + 1;
            answer.push_back(std::vector<double>());
        }
    }
    answer.pop_back();
    return Matrix<double>(answer);
}

int main(int argc, const char * argv[]) {
	auto timeStart = std::chrono::high_resolution_clock::now();
	if (argc != 3) {
		std::cout << "Error: need to pass in two arguments (input-file and epsilon)" << std::endl;
		return 0;
	}

	std::string file_path = argv[1];
    double epsilon = atof(argv[2]);
    double epislon_squared = epsilon * epsilon;
    std::string file_contents = read_file(file_path);
    Matrix<double> points = string_to_matrix(file_contents);

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
    std::cout << my_graph;
    return 0;
}

