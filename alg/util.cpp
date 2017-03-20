#include <map>
#include <set>
#include <fstream>
#include <iostream>
#include "SimplexInfo.hpp"

template <class V>
std::ostream& operator << (std::ostream& os, std::set<V>& set) {
    os << "{";
    int counter = 0;
    for (auto it = set.begin(); it != set.end(); ++it) {
        if (counter != 0) os << ", ";
        os << *it;
        ++counter;
    }
    os << "}";
    return os;  
}

template <class K, class V>
std::ostream& operator << (std::ostream& os, std::map<K,V>& map) {
    os << "{" << std::endl;
    for (auto it = map.begin(); it != map.end(); ++it) {
        os << (*it).first;
        os << ": ";
        os << (*it).second;
        os << std::endl;
    }
    os << "}" << std::endl;
    return os;  
}

template <class V>
std::ostream& operator << (std::ostream& os, std::vector<V>& vec) {
    os << "[";
    for (int i = 0; i < vec.size(); ++i) {
        if (i != 0) os << ", ";
        os << vec[i];
    }
    os << "]" << std::endl;
    return os;  
}

void print_simplex(Simplex s, SimplexInfo* info) {
    if (info[s].dim() == 0) std::cout << "{" << s << "}";
    else std::cout << info[s];
}

int choose(int n, int k) {
    int rtn = 1;
    for (int i = n-k+1; i <= n; ++i) rtn *= i;
    for (int i = 1; i <= k; ++i) rtn /= i;
    return rtn;
}


/*
 * utility function to read the contents of a file as a string
 * @param file_path - path to file to read
 * @return contents of given fiel
 */
std::string* read_file(std::string &file_path) {
    std::string line;
    std::string *contents = new std::string("");

    std::ifstream file(file_path);
    if (file.is_open()) {
        while(getline(file, line)) {
            *contents += line;
            *contents += "\n";
        }
        file.close();
    }
    else {
        throw std::invalid_argument("File \"" + file_path + "\" not found.");
    }
    return contents;
}

/*
 * utility function to write to the contents of a file as a string
 * if no file of the name exists, one will be created
 * if a file of the name exists, it will be overwritten
 * @param file_path - path to file to write to
 * @param contents - new contents of file
 */
void write_file(std::string &file_path, std::string &contents) {
    std::ofstream file;
    file.open(file_path);
    file << contents;
    file.close();
}

/*
 * this function prints out the given message and then terminates the program
 * @param message - message to be printed
 */
void terminate(std::string message) {
    std::cout << message << "\n";
    exit(1);
}

std::vector<std::vector<float>>* split(std::string *str) {
    std::vector<std::vector<float>> *rtn = new std::vector<std::vector<float>>();
    rtn->push_back(std::vector<float>());
    int last_split = -1;
    for (int i = 0; i < str->length(); ++i) {
        if (str->at(i) == ',' || str->at(i) == '\n') {
            rtn->back().push_back(stof(str->substr(last_split+1, i-last_split-1)));
            last_split = i;
        }
        if (str->at(i) == '\n' && i+1 != str->length()) rtn->push_back(std::vector<float>());
    }
    return rtn;
}

float dist(std::vector<float> &x, std::vector<float> &y) {
    if (x.size() != y.size()) terminate("file format incorrect");

    float rtn = 0;
    for (int i = 0; i < x.size(); ++i) {
        rtn += (x[i] - y[i]) * (x[i] - y[i]);
    }
    return sqrt(rtn);
}
