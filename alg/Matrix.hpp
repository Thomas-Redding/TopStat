
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>

class Matrix {
public:
    /*
     * Create the zero matrix with the given width and height
     */
    Matrix(int w, int h) {
        width = w;
        height = h;
        arr = new float*[w];
        for (int i = 0; i < w; ++i) {
            arr[i] = new float[h];
            for (int j = 0; j < h; ++j) {
                arr[i][j] = float();
            }
        }
    }

    /*
     * Create a matrix from a vector of vectors
     */
    Matrix (std::vector<std::vector<float>> *input) {
        width = input->size();
        height = input->at(0).size();
        arr = new float*[width];
        for (int i = 0; i < width; ++i) {
            arr[i] = new float[height];
            if (input->at(i).size() != height)
                throw std::invalid_argument("Matrix constructor passed 2d vector with non-equal subvectors.");
            for (int j = 0; j < height; ++j)
                arr[i][j] = input->at(i)[j];
        }
    }
    ~Matrix() {
        delete arr;
    }
    void set(int x, int y, float value) {
        arr[x][y] = value;
    }
    float get(int x, int y) const {
        return arr[x][y];
    }
    int get_width() const {
        return width;
    }
    int get_height() const {
        return height;
    }
private:
    int width = 0;
    int height = 0;
    float** arr;
};

/*
 * allow matrices to printed out nicely using std::cout
 * 1  2  3
 * 4  5  6
 * 7  8  9
 */
std::ostream& operator << (std::ostream& os, Matrix& mat) {
    for (int i = 0; i < mat.get_width(); ++i) {
        for (int j = 0; j < mat.get_height(); ++j)
            os << mat.get(i, j) << "\t";
        os << "\n";
    }
    return os;  
}

#endif
