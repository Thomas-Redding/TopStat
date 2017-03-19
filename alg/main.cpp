
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>

#include "Matrix.hpp"
#include "util.cpp"
#include "SimplexInfo.hpp"


Matrix* create_distance_matrix(std::string *str, std::string format) {
    std::vector<std::vector<float>>* table = split(str);
    if (format == "points") {
        Matrix* rtn = new Matrix(table->size(), table->size());
        for (int i = 0; i < table->size(); ++i) {
            for (int j = 0; j < table->size(); ++j)
                rtn->set(i, j, dist(table->at(i), table->at(j)));
        }
        return rtn;
    }
    else if (format == "distances") {
        // distance
        Matrix* rtn = new Matrix(table->size(), table->size());
        for (int i = 0; i < table->size(); ++i) {
            for (int j = 0; j < table->size(); ++j) {
                rtn->set(i, j, table->at(i)[j]);
            }
        }

        // verify correctness
        for (int i = 0; i < rtn->get_width(); ++i) {
            for (int j = 0; j < i; ++j) {
                if (rtn->get(i, j) != rtn->get(j, i))
                    terminate("distance matrix not symmetric");
            }
            if (rtn->get(i,i) != 0)
                terminate("distance matrix has a non-zero diagonal entry");
        }
        return rtn;
    }
    else {
        return nullptr;
    }
}

Simplex points_to_edge(Simplex P, Simplex Q, unsigned int num_points) {
    Simplex from = std::min(P, Q);
    Simplex to = std::max(P, Q);
    return num_points*from - from*(from+1.0)/2.0 + (to - from - 1) + num_points;
}

/*
 * ./bars input-file.txt output-file.txt
 */
int main(int argc, const char * argv[]) {
    if (argc < 3) terminate("too few arguments given (2 or 3 expected)");
    if (argc > 4) terminate("too many arguments given (2 or 3 expected)");

    std::string input_path = argv[1];
    std::string output_path = argv[2];
    std::string format = "points";
    if (argc == 4) format = argv[3];

    if (format != "points" && format != "distances")
        terminate("expected 'points', 'distances', or nothing as third argument");

    std::cout << "loading file..." << std::endl;
    std::string *input = read_file(input_path);

    std::cout << "creating distance matrix..." << std::endl;
    Matrix *dist = create_distance_matrix(input, format);
    std::cout << *dist;

    std::cout << "creating list of simplicies..." << std::endl;
    unsigned int num_points = dist->get_width();
    unsigned int num_simplices = 0;
    num_simplices += choose(num_points, 1);            // points
    num_simplices += choose(num_points, 2);            // edges
    num_simplices += choose(num_points, 3);            // triangles
    SimplexInfo *simplices = new SimplexInfo[num_simplices];

    std::cout << "adding points to list of simplicies..." << std::endl;
    int counter = 0;
    for (int i = 0; i < num_points; ++i) {
        simplices[counter] = SimplexInfo();
        ++counter;
    }

    std::cout << "adding edges to list of simplicies..." << std::endl;
    for (int i = 0; i < num_points; ++i) {
        for (int j = i+1; j < num_points; ++j) {
            simplices[counter] = SimplexInfo(dist->get(i, j), i, j);
            ++counter;
        }
    }

    std::cout << "adding triangles to list of simplicies..." << std::endl;
    for (int i = 0; i < num_points; ++i) {
        for (int j = i+1; j < num_points; ++j) {
            for (int k = j+1; k < num_points; ++k) {
                float ep = std::max(std::max(dist->get(i, j), dist->get(i, k)), dist->get(j, k));
                Simplex edge1 = points_to_edge(i, j, num_points);
                Simplex edge2 = points_to_edge(i, k, num_points);
                Simplex edge3 = points_to_edge(j, k, num_points);
                simplices[counter] = SimplexInfo(ep, edge1, edge2, edge3);
                ++counter;
            }
        }
    }

    std::cout << "sorting simplicies in order of creation..." << std::endl;
    Simplex *simplices_by_epsilon = new Simplex[num_simplices];
    for (int i = 0; i < num_simplices; ++i) simplices_by_epsilon[i] = i;
    sort_simplices_by_epsilon(simplices_by_epsilon, simplices, num_simplices);

    std::cout << "coloring points..." << std::endl;
    std::map<Simplex, unsigned int> point_to_color;
    std::map<unsigned int, std::vector<Simplex>> color_to_points;
    for (int i = 0; i < dist->get_width(); ++i) {
        point_to_color[i] = i;
        color_to_points[i] = std::vector<Simplex>();
        color_to_points[i].push_back(i);
    }

    std::cout << "starting Betti-0 bars..." << std::endl;
    std::set<Simplex> creations;              // (index, dimension)
    std::map<Simplex, Simplex> killers;       // (index, dimension) -> (index, dimension)
    for (int i = 0; i < dist->get_width(); ++i) {
        creations.insert(i);
    }

    std::cout << "expanding epsilons..." << std::endl;
    for (int it = 0; it < num_simplices; ++it) {
        Simplex new_simplex = simplices_by_epsilon[it];
        SimplexInfo new_simplex_info = simplices[new_simplex];
        forfeit
        // todo
    }

    std::cout << "writing out bars..." << std::endl;
    std::string output = "stuff";
    write_file(output_path, output);

    // for (int i = 0; i < num_simplices; ++i) {
    //     std::cout << i << ": ";
    //     print_simplex(simplices_by_epsilon[i], simplices);
    // }

    std::cout << "cleaning up memory..." << std::endl;
    delete input;
    delete dist;

    std::cout << "Done." << std::endl;
    return 0;
}
