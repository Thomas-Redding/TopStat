
#ifndef BASIS_HPP
#define BASIS_HPP

#include <vector>
#include <algorithm>

class Basis {
public:
    /*
     * Create the zero matrix with the given width and height
     */
    Basis(unsigned int d) {
        dim = d;
    }

    ~Basis() {
        for (unsigned int i = 0; i < vectors.size(); ++i) {
            delete[] vectors[i];
        }
    }

    /*
     */
    bool add(float *vec) {
        unsigned int original_basis_size = vectors.size();
        vectors.push_back(vec);
        reduce_basis();
        eliminate_zero_columns();
        return vectors.size() != original_basis_size;
    }

    std::string to_string() {
        std::string rtn = "{";
        for (unsigned int col = 0; col < vectors.size(); ++col) {
            rtn += '[';
            for (unsigned int row = 0; row < dim; ++row) {
                if (row != 0) rtn += ", ";
                rtn += std::to_string(vectors[col][row]);
            }
            rtn += ']';
            if (col != vectors.size() - 1)
                rtn += '\n';
        }
        rtn += "}\n";
        return rtn;
    }

private:
    unsigned int dim;

    // vectors[i][j] = the jth component of the ith vector
    // vectors[i][j] = column i, row j
    std::vector<float*> vectors;

    void swap_columns(unsigned int col1, unsigned int col2) {
        iter_swap(vectors.begin() + col1, vectors.begin() + col2);
    }

    void scale_column(unsigned int col, float k) {
        for (int i = 0; i < dim; ++i)
            vectors[col][i] *= k;
    }

    void add_multiple_of_column(unsigned int col_from, unsigned int col_to, float k) {
        for (int i = 0; i < dim; ++i)
            vectors[col_to][i] += k * vectors[col_from][i];
    }

    void reduce_basis() {
        for (unsigned int row = 0; row < dim; ++row) {
            // find a non-zero entry in this row
            if (row >= vectors.size())
                break;
            unsigned int col;
            for (col = row; col < vectors.size(); ++col) {
                if (vectors[col][row] != 0)
                    break;
            }
            if (col == vectors.size()) {
                // all entries are zero in this row, so skip
                continue;
            }
            // we found a nonzero entry
            // swap columns so that the row^th column has the nonzero entry
            if (col != row)
                swap_columns(row, col);

            // scale column to make entry 1
            scale_column(row, 1.0 / vectors[row][row]);
            for (col = 0; col < vectors.size(); ++col) {
                if (col != row)
                    add_multiple_of_column(row, col, -1 * vectors[col][row]);
            }
        }
    }

    void eliminate_zero_columns() {
        for (unsigned int col = 0; col < vectors.size(); ++col) {
            unsigned int row;
            for (row = 0; row < dim; ++row) {
                if (vectors[col][row] != 0)
                    break;
            }
            if (row == dim) {
                vectors.erase(vectors.begin() + col);
                --col;
            }
        }
    }

};

/*
 * allow matrices to printed out nicely using std::cout
 * 1  2  3
 * 4  5  6
 * 7  8  9
 */
std::ostream& operator << (std::ostream& os, Basis& basis) {
    os << basis.to_string();
    return os;  
}

#endif
