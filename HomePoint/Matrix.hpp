#ifndef MATRIX_HPP
#define MATRIX_HPP

template <class T>
class Matrix {
public:
    /*
     * Create the zero matrix with the given width and height
     */
    Matrix(uint w, uint h) {
        width = w;
        height = h;
        arr = new T*[w];
        for (uint i = 0; i < w; ++i) {
            arr[i] = new T[h];
            for (uint j = 0; j < h; ++j) {
                arr[i][j] = T();
            }
        }
    }

    /*
     * Create a matrix from a vector of vectors
     */
    Matrix (std::vector<std::vector<double>> input) {
        width = input.size();
        height = input[0].size();
        arr = new T*[width];
        for (uint i = 0; i < width; ++i) {
            arr[i] = new T[height];
            if (input[i].size() != height)
                throw std::invalid_argument("Matrix constructor passed 2d vector with non-equal subvectors.");
            for (uint j = 0; j < height; ++j)
                arr[i][j] = input[i][j];
        }
    }
    ~Matrix() {
        delete arr;
    }
    void set(uint x, uint y, T value) {
        arr[x][y] = value;
    }
    T get(uint x, uint y) {
        return arr[x][y];
    }
    uint get_width() {
        return width;
    }
    uint get_height() {
        return height;
    }
private:
    uint width = 0;
    uint height = 0;
    T** arr;
};

/*
 * allow matrices to printed out nicely using std::cout
 * 1  2  3
 * 4  5  6
 * 7  8  9
 */
template <class T>
std::ostream& operator << (std::ostream& os, Matrix<T>& mat) {
    for (int i = 0; i < mat.get_width(); ++i) {
        for (int j = 0; j < mat.get_height(); ++j)
            os << mat.get(i, j) << "\t";
        os << "\n";
    }
    return os;  
}

#endif
