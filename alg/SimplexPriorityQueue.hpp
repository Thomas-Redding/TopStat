
#ifndef SIMPLEXPRIORITYQUEUE_HPP
#define SIMPLEXPRIORITYQUEUE_HPP

#include <vector>

#include "SimplexInfo.hpp"

class SimplexPriorityQueue {
public:
    SimplexPriorityQueue(SimplexInfo* i, unsigned int l);
    void add(Simplex val);
    Simplex pop();
    Simplex peek();
    unsigned int size();
    void print();
private:
    SimplexInfo *info = nullptr;
    unsigned int len = 0;
    std::vector<Simplex> arr;
};

SimplexPriorityQueue::SimplexPriorityQueue(SimplexInfo* i, unsigned int l) {
    info = i;
    len = l;
}

unsigned int SimplexPriorityQueue::size() {
    return arr.size();
}

void SimplexPriorityQueue::print() {
    std::cout << "[";
    for (int i = 0; i < arr.size(); ++i) {
        std::cout << arr[i] << ">";
        print_simplex(arr[i], info);
    }
    std::cout << "]" << std::endl;
}

void SimplexPriorityQueue::add(Simplex val) {
    // don't add duplicates
    for (int i = 0; i < arr.size(); ++i) {
        if (arr[i] == val)
            return;
    }

    // add to queue
    int i;
    for (i = 0; i < arr.size(); ++i) {
        if (info[val].epsilon() > info[arr[i]].epsilon()) {
            arr.insert(arr.begin()+i, val);
            break;
        }
    }
    if (i == arr.size()) arr.push_back(val);
}

Simplex SimplexPriorityQueue::pop() {
    Simplex rtn = arr[0];
    arr.erase(arr.begin(), arr.begin() + 1);
    return rtn;
}

Simplex SimplexPriorityQueue::peek() {
    return arr.front();
}

#endif
