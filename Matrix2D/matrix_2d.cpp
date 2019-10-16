/**
* @file matrix_2d.cpp
* Contains implementations of the Matrix2D class.
*/
#include "stdafx.h"
#include "matrix_2d.h"

using namespace std;

namespace m2d {
	void InputMatrix(ifstream &ifs, Matrix2D &mat) {
		double temp;
		mat.m.lock();
		for (size_t l = 0; l < mat.getSizeX(); l++) {
			for (size_t c = 0; c < mat.getSizeY(); c++) {
				ifs >> temp;
				mat.setAt(l, c, temp);
			}
		}
		mat.m.unlock();
	}

	// Constructor and destructor
	Matrix2D::Matrix2D(size_t size_x, size_t size_y) :
		size_x(size_x), size_y(size_y) {
		elem = new double*[size_x];
		for (size_t i = 0; i < size_x; i++) {
			elem[i] = new double[size_y] {0};
		}
	}
	Matrix2D::~Matrix2D() {
		for (size_t i = 0; i < size_x; i++) delete[] elem[i];
		delete[] elem;
	}

	// Getter and setter
	double Matrix2D::getAt(size_t pos_x, size_t pos_y) const {
		if (pos_x >= size_x || pos_y >= size_y) throw out_of_range("Indices exceeded Matrix2D range.");
		return elem[pos_x][pos_y];
	}
	void Matrix2D::setAt(size_t pos_x, size_t pos_y, double val) {
		m.lock();
		try {
			if (pos_x >= size_x || pos_y >= size_y) throw out_of_range("Indices exceeded Matrix2D range.");
			elem[pos_x][pos_y] = val;
			m.unlock();
			return;
		}
		catch (out_of_range &e) {
			cerr << e.what();
			m.unlock();
			return;
		}
	}

	// Display
	void Matrix2D::print() const {
		for (size_t x = 0; x < size_x; x++) {
			for (size_t y = 0; y < size_y; y++) {
				cout << this->getAt(x, y) << ' ';
			}
			cout << endl;
		}
	}

	// Concatenators
	Matrix2D &ConcatenateHorizontally(const Matrix2D &left, const Matrix2D &right) {
		if (left.getSizeX() != right.getSizeX())
			throw range_error("Cannot horizontally concatenate two matrices with different row counts.");
		Matrix2D result(left.getSizeX(), left.getSizeY() + right.getSizeY()); // TODO: warn if overflow
		for (size_t y = 0; y < left.getSizeY(); y++) { // put left matrix in result
			for (size_t x = 0; x < left.getSizeX(); x++) {
				result.setAt(x, y, left.getAt(x, y));
			}
		}
		for (size_t y = left.getSizeY(); y < left.getSizeY() + right.getSizeY(); y++) { // put right matrix in result
			for (size_t x = 0; x < left.getSizeX(); x++) {
				result.setAt(x, y, right.getAt(x, y));
			}
		}
		return result;
	}
	Matrix2D &ConcatenateVertically(const Matrix2D &top, const Matrix2D &bottom) {
		if (top.getSizeY() != bottom.getSizeY())
			throw range_error("Cannot vertically concatenate two matrices with different column counts.");
		Matrix2D result(top.getSizeX() + bottom.getSizeX(), top.getSizeY()); // TODO: warn if overflow
		for (size_t x = 0; x < top.getSizeX(); x++) { // put top matrix in result
			for (size_t y = 0; y < top.getSizeY(); y++) {
				result.setAt(x, y, top.getAt(x, y));
			}
		}
		for (size_t x = top.getSizeX(); x < top.getSizeX() + bottom.getSizeX(); x++) { // put top matrix in result
			for (size_t y = 0; y < top.getSizeY(); y++) {
				result.setAt(x, y, top.getAt(x, y));
			}
		}
		return result;
	}

	// Submatrix extractor
	Matrix2D &Matrix2D::subMatrix(size_t pos_x, size_t pos_y, size_t sub_size_x, size_t sub_size_y) {
		if (pos_x + sub_size_x - 1 > size_x || pos_y + sub_size_y - 1 > size_y)
			throw range_error("Submatrix is out of bounds.");
		Matrix2D result(sub_size_x, sub_size_y);
		for (size_t x = 0; x < sub_size_x; x++) {
			for (size_t y = 0; y < sub_size_y; y++) {
				result.setAt(x, y, this->getAt(pos_x + x, pos_y + y));
			}
			return result;
		}
	}

	// Transposer
	void Matrix2D::transpose() noexcept {
		double** transposed = new double* [size_y];
		for (size_t i = 0; i < size_y; i++) {
			transposed[i] = new double[size_x];
			for (size_t j = 0; j < size_x; j++) {
				transposed[i][j] = getAt(i, j);
			}
		}
		for (size_t i = 0; i < size_x; i++) delete[] elem[i];
		delete[] elem;
		elem = transposed;
		size_t temp = size_x;
		size_x = size_y;
		size_y = temp;
		return;
	}

	// Bool checks
	bool Matrix2D::isUpperTriangular() const {
		if (!isSquare()) return false;
		for (size_t x = 1; x < size_x; x++) {
			for (size_t y = 0; y < x; y++) {
				if (getAt(x, y) != 0) return false;
			}
		}
		return true;
	}
	bool Matrix2D::isLowerTriangular() const {
		if (!isSquare()) return false;
		for (size_t x = 0; x < size_x; x++) {
			for (size_t y = x + 1; y < size_y; y++) {
				if (getAt(x, y) != 0) return false;
			}
		}
		return true;
	}
	bool Matrix2D::isDiagonal() const {
		if (!isSquare()) return false;
		for (size_t x = 0; x < size_x; x++) {
			for (size_t y = 0; y < size_y; y++) {
				if (x != y && getAt(x, y) != 0) return false;
			}
		}
		return true;
	}

}




