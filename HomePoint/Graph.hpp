#include <map>
#include <set>

/*
 * A graph with undirected, unweighted edges, and vertices labeled by integers.
 */
class IntGraph {
public:
    IntGraph(uint num) {
        vertex_num = num;
    }
    void add_edge(uint from, uint to) {
        if (edges.find(from) == edges.end())
            edges[from] = std::set<uint> ();
        if (edges.find(to) == edges.end())
            edges[to] = std::set<uint> ();
        edges[from].insert(to);
        edges[to].insert(from);
    }
    std::set<uint> get_neighbors(uint vertex) {
        if (vertex >= vertex_num)
            return std::set<uint>();
        if (edges.find(vertex) == edges.end())
            return std::set<uint>();
        return edges[vertex];
    }
    uint size() {
        return vertex_num;
    }
private:
    uint vertex_num;
    std::map<uint, std::set<uint>> edges;
};

/*
 * Allows IntGraph to printed out nicely using std::cout
 * For instance, the graph 0 -- 1 -- 2 -- 3 would be displayed:
 * 0: 1
 * 1: 0, 2
 * 2: 1, 3
 * 3: 2
 * or, more generally,
 * vertex: list-of-neighbors
 */
std::ostream& operator << (std::ostream& os, IntGraph& graph) {
    for (uint i = 0; i < graph.size(); ++i) {
        std::set<uint> neighbors = graph.get_neighbors(i);
        os << i << ": ";
        bool is_first = true;
        for (std::set<uint>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
            if (is_first)
                is_first = false;
            else
                os << ", ";
            os << *it;
        }
        os << "\n";
    }
    return os;  
}