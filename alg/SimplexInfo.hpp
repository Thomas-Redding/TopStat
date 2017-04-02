
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
        _points.push_back(from);
        _points.push_back(to);
        _boundary.push_back(from);
        _boundary.push_back(to);
    }
    SimplexInfo(float e, Simplex p1, Simplex p2, Simplex p3, Simplex edge1, Simplex edge2, Simplex edge3) {
        // triangle
        ep = e;
        _points.push_back(p1);
        _points.push_back(p2);
        _points.push_back(p3);
        _boundary.push_back(edge1);
        _boundary.push_back(edge2);
        _boundary.push_back(edge3);
    }

    SimplexInfo(float e, Simplex p1, Simplex p2, Simplex p3, Simplex p4, Simplex tri1, Simplex tri2, Simplex tri3, Simplex tri4) {
        ep = e;
        _points.push_back(p1);
        _points.push_back(p2);
        _points.push_back(p3);
        _points.push_back(p4);
        _boundary.push_back(tri1);
        _boundary.push_back(tri2);
        _boundary.push_back(tri3);
        _boundary.push_back(tri4);
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

    std::vector<Simplex> points() const {
        return _points;
    }

private:
    float ep;
    std::vector<Simplex> _boundary;
    std::vector<Simplex> _points;
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
    for (int j = 0; j < len; ++j) {
        for (int i = 0; i < len - 1; ++i) {
            if (info[indices[i]] > info[indices[i+1]]) {
                Simplex sto = indices[i];
                indices[i] = indices[i+1];
                indices[i+1] = sto;
            }
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
    std::vector<Simplex> bound = simp.boundary();
    for (int i = 0; i < bound.size(); ++i) {
        if (i != 0) os << ", ";
        os << bound[i];
    }
    os << "\n  points: ";
    std::vector<Simplex> pts = simp.points();
    for (int i = 0; i < pts.size(); ++i) {
        if (i != 0) os << ", ";
        os << pts[i];
    }
    os << "\n}\n";
    return os;  
}

#endif
