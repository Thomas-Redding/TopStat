
#ifndef SIMPLEXINFO_HPP
#define SIMPLEXINFO_HPP

#include <vector>

typedef unsigned int Simplex;

class SimplexInfo {
public:
    SimplexInfo() {
        // point
        ep = 0;
    }
    SimplexInfo(float e, Simplex from, Simplex to) {
        // edge
        ep = e;
        _boundary.push_back(from);
        _boundary.push_back(to);
    }
    SimplexInfo(float e, Simplex edge1, Simplex edge2, Simplex edge3) {
        // triangle
        ep = e;
        _boundary.push_back(edge1);
        _boundary.push_back(edge2);
        _boundary.push_back(edge3);
    }

    float epsilon() const {
        return ep;
    }

    unsigned int dim() const {
        if (_boundary.size() == 0) return 0;
        else return _boundary.size() - 1;
    }

    std::vector<Simplex> boundary() const {
        return _boundary;
    }
private:
    float ep;
    std::vector<Simplex> _boundary;
};

bool operator < (const SimplexInfo &x, const SimplexInfo &y) {
    if (x.epsilon() < y.epsilon()) return true;
    if (x.epsilon() > y.epsilon()) return false;
    if (x.dim() < y.dim()) return true;
    if (x.dim() > y.dim()) return false;
    return false; // equal
}

bool operator > (const SimplexInfo &x, const SimplexInfo &y) {
    return y < x;
}

void sort_simplices_by_epsilon(Simplex *indices, SimplexInfo *info, int len) {
    for (int i = 0; i < len - 1; ++i) {
        if (info[indices[i]] > info[indices[i+1]]) {
            Simplex sto = indices[i];
            indices[i] = indices[i+1];
            indices[i+1] = sto;
        }
    }
}


/*
 * allow simplicies to printed out nicely using std::cout
 * {
 *   ep: 0.32
 *   bound: 0, 1
 * }
 * 
 * 
 */
std::ostream& operator << (std::ostream& os, const SimplexInfo& simp) {
    os << "{\n  epsilon: " << simp.epsilon() << "\n  boundry: ";
    std::vector<unsigned int> bound = simp.boundary();
    for (int i = 0; i < bound.size(); ++i) {
        if (i != 0) os << ", ";
        os << bound[i];
    }
    os << "\n}\n";
    return os;  
}

#endif
