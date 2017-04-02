
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>
#include <limits>

#include "Matrix.hpp"
#include "Basis.hpp"
#include "util.cpp"
#include "SimplexInfo.hpp"
#include "SimplexPriorityQueue.hpp"

/*
 * Note, the 2-dimensional bars are definitely not good,
 * because we need to consider only 3-simplices, which we currently don't.
 * I don't guarantee the lower-dimensional bars work, but they might :)
 */

class Bar {
public:
    Bar(unsigned int d, float s, float e) {
        dimension = d;
        start = s;
        end = e;
    }
    unsigned int dimension;
    float start;
    float end;
};


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
    return num_points*from - from*(from+1.0)/2.0 + (to-from-1) + num_points;
}


// https://i.stack.imgur.com/M6L8y.png
float cubic_solver(float a, float b, float c, float d) {
    float term_1 = b/(-3*a);
    float term_2 = 1/(-3*a) * cbrt(0.5*(2*b*b*b - 9*a*b*c + 27*a*a*d + sqrt(pow(2*b*b*b-9*a*b*c+27*a*a*d, 2) - 4*pow(b*b-3*a*c, 3))));
    float term_3 = 1/(-3*a) * cbrt(0.5*(2*b*b*b - 9*a*b*c + 27*a*a*d - sqrt(pow(2*b*b*b-9*a*b*c+27*a*a*d, 2) - 4*pow(b*b-3*a*c, 3))));
    return term_1 + term_2 + term_3;
}


Simplex points_to_triangle(Simplex P, Simplex Q, Simplex R, unsigned int num_points, unsigned int num_edges) {
    unsigned int n = num_points;
    Simplex p0 = std::min(P, std::min(Q, R));
    Simplex p2 = std::max(P, std::max(Q, R));
    Simplex p1 = 0;
    if (P != p0 && P != p2) p1 = P;
    else if (Q != p0 && Q != p2) p1 = Q;
    else if (R != p0 && R != p2) p1 = R;

    unsigned int counter = 0;
    for (int i = 0; i < num_points; ++i) {
        for (int j = i+1; j < num_points; ++j) {
            for (int k = j+1; k < num_points; ++k) {
                if (i == p0 && j == p1 && k == p2)
                    return counter + num_edges + num_points;
                ++counter;
            }
        }
    }

    return 0;
}


/*
 * ./bars input-file.txt output-file.txt
 */
int main(int argc, const char * argv[]) {
    // general error checking
    std::cout << "making sure you didn't mess up..." << std::endl;
    if (argc < 2) terminate("too few arguments given (2 or 3 expected)");
    if (argc > 3) terminate("too many arguments given (2 or 3 expected)");
    std::string input_path = argv[1];
    std::string format = "points";
    if (argc == 3) format = argv[2];
    if (format != "points" && format != "distances")
        terminate("expected 'points', 'distances', or nothing as third argument");


    // read file with input data
    std::cout << "counting your data really fast..." << std::endl;
    std::string *input = read_file(input_path);


    // compute distance matrix
    std::cout << "using my ruler really fast..." << std::endl;
    Matrix *dist = create_distance_matrix(input, format);


    // create empty array of simplicies
    std::cout << "asking your computer for memory..." << std::endl;
    unsigned int num_points = dist->get_width();
    unsigned int num_edges = choose(num_points, 2);
    unsigned int num_triangles = choose(num_points, 3);
    unsigned int num_tetrahedron = choose(num_points, 4);
    unsigned int num_simplices = 0;
    num_simplices += num_points;
    num_simplices += num_edges;
    num_simplices += num_triangles;
    num_simplices += num_tetrahedron;
    SimplexInfo *simplices = new SimplexInfo[num_simplices];

    std::cout << num_points << ":" << num_edges << ":" << num_triangles << ":" << num_tetrahedron << std::endl;


    // add points to array of simplices
    std::cout << "looking at all your seeds..." << std::endl;
    int counter = 0;
    for (int i = 0; i < num_points; ++i) {
        simplices[counter] = SimplexInfo();
        ++counter;
    }


    // add edges to array of simplices
    std::cout << "plowing land..." << std::endl;
    for (int i = 0; i < num_points; ++i) {
        for (int j = i+1; j < num_points; ++j) {
            simplices[counter] = SimplexInfo(dist->get(i, j), i, j);
            ++counter;
        }
    }


    // add triangles to array of simplices
    std::cout << "buying more land..." << std::endl;
    for (int i = 0; i < num_points; ++i) {
        for (int j = i+1; j < num_points; ++j) {
            for (int k = j+1; k < num_points; ++k) {
                float ep = std::max(std::max(dist->get(i, j), dist->get(i, k)), dist->get(j, k));
                Simplex edge1 = points_to_edge(i, j, num_points);
                Simplex edge2 = points_to_edge(i, k, num_points);
                Simplex edge3 = points_to_edge(j, k, num_points);
                simplices[counter] = SimplexInfo(ep, i, j, k, edge1, edge2, edge3);
                ++counter;
            }
        }
    }


    // add tetrahedrons to array of simplices
    std::cout << "buying 3d land..." << std::endl;
    for (int i = 0; i < num_points; ++i) {
        for (int j = i+1; j < num_points; ++j) {
            for (int k = j+1; k < num_points; ++k) {
                for (int l = k+1; l < num_points; ++l) {
                    float ep = std::max(
                        std::max(dist->get(i, j), std::max(dist->get(i, k), dist->get(i, l))),
                        std::max(dist->get(j, k), std::max(dist->get(j, l), dist->get(k, l)))
                    );
                    Simplex face1 = points_to_triangle(i, j, k, num_points, num_edges);
                    Simplex face2 = points_to_triangle(i, j, l, num_points, num_edges);
                    Simplex face3 = points_to_triangle(i, k, l, num_points, num_edges);
                    Simplex face4 = points_to_triangle(j, k, l, num_points, num_edges);
                    simplices[counter] = SimplexInfo(ep, i, j, k, l, face1, face2, face3, face4);
                    ++counter;
                }
            }
        }
    }

    // print simplices for debuggin
    // for (int i = 0; i < num_simplices; ++i)
    //     std::cout << simplices[i];

    // sort simplices by creation date (epsilon)
    std::cout << "sorting shapes by distance to retirement..." << std::endl;
    Simplex *simplices_by_epsilon = new Simplex[num_simplices];
    for (int i = 0; i < num_simplices; ++i) simplices_by_epsilon[i] = i;
    sort_simplices_by_epsilon(simplices_by_epsilon, simplices, num_simplices);


    // we "color" simplices to detect whether new simplicies are positive or negative
    std::cout << "coloring inside the lines..." << std::endl;
    Basis boundary_operator2(num_edges);
    int *point_to_color = new int[num_points]; // point_to_color[point_index]
    std::map<unsigned int, std::vector<Simplex>> color_to_points; // color_to_points[color] = [point_indices]
    for (int i = 0; i < num_points; ++i) {
        point_to_color[i] = i;
        color_to_points[i] = std::vector<Simplex>();
        color_to_points[i].push_back(i);
    }


    // create the bars associated with the 0-simplices (points)
    std::cout << "planting seeds..." << std::endl;
    std::set<Simplex> creations;              // (index, dimension)
    std::map<Simplex, Simplex> killers;       // (index, dimension) -> (index, dimension)
    for (int i = 0; i < dist->get_width(); ++i) {
        creations.insert(i);
    }


    // continue adding simplices as per
    // Topological Persistence and Simplification
    // - Herbert Edelsbrunner, David Letscher, and Afra Zomorodian
    std::cout << "feeding and watering our plants..." << std::endl;
    for (int it = num_points; it < num_simplices; ++it) {
        Simplex new_simplex = simplices_by_epsilon[it];
        SimplexInfo new_simplex_info = simplices[new_simplex];
        unsigned int dim = new_simplex_info.dim();
        std::vector<Simplex> boundary = new_simplex_info.boundary();
        std::vector<Simplex> points = new_simplex_info.points();

        bool is_positive;
        if (dim == 1) {
            bool all_same_color = true;
            unsigned int col = point_to_color[points[0]];
            for (int i = 1; i < points.size(); ++i) {
                if (point_to_color[points[i]] < col)
                    col = point_to_color[points[i]];
            }

            for (int i = 0; i < points.size(); ++i) {
                if (point_to_color[points[i]] != col) {
                    all_same_color = false;
                    break;
                }
            }

            is_positive = all_same_color;
            if (!is_positive) {
                // color the points the same color
                for (int i = 0; i < points.size(); ++i) {
                    if (point_to_color[points[i]] != col) {
                        std::vector<Simplex> simplices_we_need_to_color = color_to_points[points[i]];
                        for (int j = 0; j < simplices_we_need_to_color.size(); ++j) {
                            point_to_color[simplices_we_need_to_color[j]] = col;
                            color_to_points[col].push_back(simplices_we_need_to_color[j]);
                        }
                    }
                }
            }
        }
        else if (dim == 2) {
            float* vec = new float[num_edges]();
            unsigned int min_point = *std::min_element(points.begin(), points.end());
            unsigned int max_point = *std::max_element(points.begin(), points.end());

            std::vector<Simplex> bound0 = simplices[boundary[0]].boundary();
            std::vector<Simplex> bound1 = simplices[boundary[1]].boundary();
            std::vector<Simplex> bound2 = simplices[boundary[2]].boundary();

            vec[boundary[0] - num_points] = 1;
            vec[boundary[1] - num_points] = 1;
            vec[boundary[2] - num_points] = 1;

            if (std::find(bound0.begin(), bound0.end(), min_point) != bound0.end() && std::find(bound0.begin(), bound0.end(), max_point) != bound0.end())
                vec[boundary[0] - num_points] = -1;
            else if (std::find(bound1.begin(), bound1.end(), min_point) != bound1.end() && std::find(bound1.begin(), bound1.end(), max_point) != bound1.end())
                vec[boundary[1] - num_points] = -1;
            else if (std::find(bound2.begin(), bound2.end(), min_point) != bound2.end() && std::find(bound2.begin(), bound2.end(), max_point) != bound2.end())
                vec[boundary[2] - num_points] = -1;

            bool is_independent = boundary_operator2.add(vec);
            is_positive = !is_independent;
        }
        else {
            //
        }

        if (is_positive) {
            // positive - belongs to a cycle
            creations.insert(new_simplex);
        }
        else {
            // negative - connects components
            // color newly connected components together
            
            // find "victim" (simplex whose bar I end)
            SimplexPriorityQueue spq(simplices, num_simplices);
            std::set<Simplex> considerd;
            for (int i = 0; i < boundary.size(); ++i) {
                spq.add(boundary[i]);
                considerd.insert(boundary[i]);
            }

            while (true) {
                if (spq.size() == 0) break;
                Simplex youngest = spq.pop();
                if (creations.find(youngest) != creations.end()) {
                    // "youngest" is still alive, so kill it
                    creations.erase(creations.find(youngest));
                    killers[youngest] = new_simplex;
                    break;
                }
                else {
                    // "youngest" is already dead; add it's killer's boundary
                    if (killers.find(youngest) == killers.end())
                        continue;
                    std::vector<Simplex> bound = simplices[killers[youngest]].boundary();
                    for (int i = 0; i < bound.size(); ++i) {
                        if (considerd.find(bound[i]) == considerd.end()) {
                            spq.add(bound[i]);
                            considerd.insert(bound[i]);
                        }
                    }
                }
            }
        }
    }


    // eliminate bars of length 0
    std::cout << "pruning pitiful puny bars..." << std::endl;
    for (auto it = killers.begin(); it != killers.end();) {
        SimplexInfo start = simplices[(*it).first];
        SimplexInfo end = simplices[(*it).second];
        if (start.epsilon() == end.epsilon()) {
            auto old_it = it;
            ++it;
            killers.erase(old_it);
        }
        else {
            ++it;
        }
    }


    // combine infinite and finite bars
    std::cout << "combining gods and humans..." << std::endl;
    std::map<unsigned int, std::vector<Bar>> bars;
    bars[0] = std::vector<Bar>();
    bars[1] = std::vector<Bar>();
    bars[2] = std::vector<Bar>();
    bars[3] = std::vector<Bar>();
    float inf = std::numeric_limits<float>::max();
    for (auto it = creations.begin(); it != creations.end(); ++it) {
        SimplexInfo simp = simplices[*it];
        bars[simp.dim()].push_back(Bar(simp.dim(), simp.epsilon(), inf));
    }
    for (auto it = killers.begin(); it != killers.end(); ++it) {
        SimplexInfo start = simplices[(*it).first];
        SimplexInfo end = simplices[(*it).second];
        bars[start.dim()].push_back(Bar(start.dim(), start.epsilon(), end.epsilon()));
    }

    // print bars to terminal
    std::cout << "showing you the best bars..." << std::endl;
    for (auto it = bars.begin(); it != bars.end(); ++it) {
        unsigned int dim = (*it).first;
        std::cout << "dim: " << dim << std::endl;
        std::vector<Bar> bar_vec = (*it).second;
        for (unsigned int i = 0; i < bar_vec.size(); ++i) {
            if (bar_vec[i].end == inf)
                std::cout << "(" << bar_vec[i].start << ", inf)" << std::endl;
            else
                std::cout << "(" << bar_vec[i].start << ", " << bar_vec[i].end << ")" << std::endl;
        }
    }


    std::cout << "patching up memory leaks..." << std::endl;
    delete input;
    delete dist;
    delete[] simplices;


    std::cout << "You're welcome." << std::endl;
    return 0;
}
